%%% Paul Adkisson
%%% 9.2.21
%%% Main Body of Wang (2002) Biophysical Attractor Model + early stopping
%%% when firing rate hits 15Hz
clear;
tic;
load("bam_constants.mat")
load("brain.mat")
load("conductances.mat")
mkdir("data")
for I_dc_idx = 1:size(I_dcs, 2)
    I_dc = I_dcs(:, I_dc_idx)';
    fprintf("I_dc: %0.0fpA \n", I_dc(1)*1e12)
    output_dcpath = sprintf("data/I_dc=%0.0fpA", I_dc(1)*1e12);
    mkdir(output_dcpath)
    for c = coherences
        fprintf("Coherence: %0.1f%% \n", c*100)
        input_coherentpath = sprintf("spikes/c=%0.3f", c);
        output_coherentpath = strcat(output_dcpath, sprintf("/c=%0.3f", c));
        mkdir(output_coherentpath)
        parfor trial = start_trial:end_trial
        %for trial = start_trial:end_trial
            t_earlystop = t_span;
            fprintf("trial: %0.0f \n", trial)
            input_trialpath = strcat(input_coherentpath, sprintf("/trial%0.0f/input.mat", trial));
            output_trialpath = strcat(output_coherentpath, sprintf("/trial%0.0f.mat", trial));

            %external input
            input_spikes = load(input_trialpath, 'spikes');
            spikes = input_spikes.spikes;

            %LIF variables
            RP_ind = zeros(1, N);
            Vm = EL*ones(length(t), N);
            s_ampa = zeros(length(t), N);
            s_ampa_ext = zeros(length(t), N);
            s_nmda = zeros(length(t), N);
            x_nmda = zeros(length(t), N);
            s_gaba = zeros(length(t), N);

            %Record Spikes
            rec_spikes = zeros(length(t), N);

            for i = 1:length(t)-1
                % Synaptic current
                if i <= delay_ind
                    I_ch = zeros(1, N);
                else
                    V_ch = Vm(i, :);
                    s_ampa_ch = s_ampa(i-delay_ind, :);
                    s_ampa_ext_ch = s_ampa_ext(i, :);
                    g_ampa_ext = s_ampa_ext_ch.*G_ampa_ext(pop_type);
                    s_nmda_ch = s_nmda(i-delay_ind, :);
                    s_gaba_ch = s_gaba(i-delay_ind, :);
                    I_temp = synapse_current(V_ch, s_ampa_ch, g_ampa_ext, s_nmda_ch, s_gaba_ch, ...
                        adja, AMPA, NMDA, GABA);
                    I_ch = sum(I_temp, 1);
                end

                %External input
                spikes_ch = spikes(i, :);
                if any(spikes_ch)
                    s_ampa_ext(i, spikes_ch) = s_ampa_ext(i, spikes_ch) + 1;
                end

                %Spiking Behavior
                is_spike = (Vm(i, :)>= Vs);
                if any(is_spike)
                    Vm(i, is_spike) = 0; %draw spike
                    rec_spikes(i, is_spike) = 1; % record spike
                    Vm(i+1, is_spike) = Vr;
                    RP_ind(is_spike) = refract_ind(pop_type(is_spike));
                    s_ampa(i, is_spike) = s_ampa(i, is_spike) + 1;
                    x_nmda(i, is_spike) = x_nmda(i, is_spike) + 1;
                    s_gaba(i, is_spike) = s_gaba(i, is_spike) + 1;
                end

                %Synapse Update
                s_ampa(i+1, :) = s_ampa(i, :) - dt*s_ampa(i, :) / tau_AMPA;
                s_ampa_ext(i+1, :) = s_ampa_ext(i, :) - dt*s_ampa_ext(i, :)/tau_AMPA;
                s_nmda(i+1, :) = s_nmda(i, :) + dt*(alpha*x_nmda(i, :).*(1-s_nmda(i, :)) ...
                    - s_nmda(i, :)/tau_NMDA_2);
                x_nmda(i+1, :) = x_nmda(i, :) - dt*x_nmda(i, :)/tau_NMDA_1;
                s_gaba(i+1, :) = s_gaba(i, :) - dt*s_gaba(i, :)/tau_GABA;

                %Voltage update (for non-refractory neurons)
                non_rp = (RP_ind==0);
                in_rp = (~non_rp)&(~is_spike);
                dvdt = (-gL(pop_type).*(Vm(i, :)-EL) + I_ch + I_dc(pop_type))./C(pop_type);
                Vm(i+1, non_rp) = Vm(i, non_rp) + dvdt(non_rp)*dt;
                Vm(i+1, in_rp) = Vr; %hold refractory neurons at Vr
                RP_ind(in_rp) = RP_ind(in_rp) - 1; %decrement remaining refractory time
                
                %Real-time firing rate
                if i > win_index
                    subspikes = rec_spikes(i-win_index:i, 1:2*num_group);
                    frs = zeros(2, 1);
                    frs(1) = sum(subspikes(:, 1:num_group), 'all')/(win_size*num_group);
                    frs(2) = sum(subspikes(:, num_group+1:2*num_group), 'all')/(win_size*num_group);
                    if any(frs>earlystop_thresh) && t_earlystop==t_span 
                        t_earlystop = t(i) + 2*avg_win_size;
                    end
                    if t(i) > t_earlystop
                        fprintf("Stopping Early at t=%0.1fs", t(i))
                        break
                    end
                end
            end
            pop_frs = spikes2popfrs(rec_spikes, dt, p, f, N_E);
            parsave(output_trialpath, pop_frs);
            %save(output_trialpath, 'pop_frs', "rec_rtfrs")
        end
    end
end
toc

function I = synapse_current(V, s_ampa, g_ampa_ext, s_nmda, s_gaba, adja, AMPA, NMDA, GABA)
    Mg = 1; %mM
    E_AMPA = 0e-3; %mV
    E_NMDA = 0e-3; %mV
    E_GABA = -70e-3; %mV
    
    %First, find the total synaptic conductance (summed over all synapses)
    %for each neuron (for AMPA, NMDA, and GABA)
    g_AMPA = s_ampa*(adja.*AMPA) + g_ampa_ext;
    g_NMDA = s_nmda*(adja.*NMDA);
    g_GABA = s_gaba*(adja.*GABA);
    
    %Then, use post-synaptic membrane potential to calculate total current
    I = zeros(3, length(V));
    I(1, :) = g_AMPA.*(E_AMPA-V);
    I(2, :) = g_NMDA.*(E_NMDA-V) ./ (1 + Mg*exp(-0.062*V*1000)/3.57);
    I(3, :) = g_GABA.*(E_GABA-V);
end

function frs = spikes2popfrs(spikes, dt, p, f, N_E)
    win_size = 5e-3;
    avg_win_size = 50e-3;
    num_group = floor(f*N_E);
    w = ones(floor(win_size/dt), 1);
    w = w ./ length(w);
    neuron_frs = filter(w, 1, spikes) ./ dt;
    w = ones(floor(avg_win_size/dt), 1);
    w = w ./ length(w);
    neuron_frs = filter(w, 1, neuron_frs);
    frs = zeros(size(spikes, 1), p+2);
    for i = 1:p
        group_idx = ((i-1)*num_group+1):i*num_group;
        frs(:, i) = mean(neuron_frs(:, group_idx), 2);
    end
    frs(:, end-1) = mean(neuron_frs(:, f*p*N_E+1:N_E), 2);
    frs(:, end) = mean(neuron_frs(:, N_E+1:end), 2);
end

function parsave(savepath, pop_frs)
    save(savepath, 'pop_frs')
end