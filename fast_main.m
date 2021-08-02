clear;
tic
load("tdcs_constants.mat")
%% setups for different experiment
% fixed k
fixed_k = 3;
if ismember(fixed_k, 0:4)
    loop_k = fixed_k;
else
    loop_k = 0:4;
end

% fixed E/I balance
% index of the E/I balance wanted in the list above
% anything else put the E/I balance loop through the list
fixed_ei = 1; 
if ismember(fixed_ei, 1:length(ei_balance))
    loop_ei = fixed_ei;
else
    loop_ei = 1 : length(ei_balance);
end

% this part determines whether or not there is DC in the experiment, and
% to which part DC is applied, options are +/-1) whole cortical population; 
% +/-2) excitatory neurons in cortical population; +/-3) inhibitory neurons in 
% cortical population; 4) no DC
dc_type = 4;
dc_amp = 2e-12;
switch dc_type
    case 1
        dc = dc_amp*repmat((population_type==1), length(ei_balance), 1);
    case -1
        dc = -dc_amp*repmat((population_type==1), length(ei_balance), 1);
    case 2
        dc = dc_amp*repmat((population_type==1), length(ei_balance), 1).*(neuron_type)';
    case -2
        dc = -dc_amp*repmat((population_type==1), length(ei_balance), 1).*(neuron_type)';
    case 3
        dc = dc_amp*repmat((population_type==1), length(ei_balance), 1).*(~neuron_type)';
    case -3
        dc = -dc_amp*repmat((population_type==1), length(ei_balance), 1).*(~neuron_type)';
    case 4
        dc = zeros(length(ei_balance), num);
end

% the subjects
% you can choose a certain subject, or the code will loop through all of
% them
fixed_brains = 1; 
if ismember(fixed_brains, 1:num_sub)
    loop_brains = fixed_brains;
    % Clear brains you don't need
    for i = 1:num_sub
        if i ~= fixed_brains
            clear(sprintf("brain%0.0f", i))
        end
    end
else
    loop_brains = 1 : num_sub;
end
%% main body
for brain = loop_brains
    load(sprintf("brains/%0.0f.mat", brain), "adja");
    for ei = loop_ei
        eir = ei_balance(ei);
        load(sprintf("conductance_maps/ei%0.1f.mat", eir), "AMPA", "NMDA", "GABA");
        for k = loop_k
            save_path = ['data/brain=', num2str(brain), ...
                ' ei=', num2str(ei_balance(ei)), ' k=', num2str(k), ...
                ' dc_type=', num2str(dc_type)];
            mkdir(save_path);
            parfor trial = 1 : num_trials
                % background input
                bgSP = pSpike(200, length(t), 2*num_cor, dt, I_sti);
                bgSI = pSpike(0, length(t), num_inter, dt, I_sti);
                bgSpike = [bgSP; bgSI];
                % task-related input
                task_freq = mapping(step_freq, std_freq, k, rand_on);
                p1Spike = pSpike(task_freq(1), length(t), num_cor, dt, I_sti);
                p2Spike = pSpike(task_freq(2), length(t), num_cor, dt, I_sti);
                spike = [p1Spike; p2Spike; zeros(num_inter, length(t))];
                % loop through time to calculate membrane potential
                RP_ind = zeros(1, num);
                Vm = -70*ones(length(t), num)*1e-3;
                t_channel = inf*ones(length(t), num); % time after the last spike
                for i = 1 : length(t)-1
                    if i <= delay_ind
                        % the following three are variables considering the
                        % delay
                        t_ch = inf*ones(1, num);
                        V_ch = zeros(1, num);
                        I_ch = zeros(1, num);
                    else
                        t_ch = t_channel(i-delay_ind, :);
                        V_ch = Vm(i-delay_ind, :);
                        I_temp = synapse_current(V_ch, t_ch);
                        I_ch = I_temp(1, :)*(adja.*AMPA)+ ...
                            I_temp(2, :)*(adja.*NMDA)+ ...
                            I_temp(3, :)*(adja.*GABA);
                    end
                    % calculate external input
                    if (t(i)>t_task)
                        I_ext = bgSpike(:, i)'+spike(:, i)'+dc(ei, :);
                    else
                        I_ext = bgSpike(:, i)'+dc(ei, :);
                    end
                    % if there is a spike, reset the membrane potential and
                    % put the neuron in refractory period, also update
                    % t_channel
                    is_spike = (Vm(i, :)>=VS);
                    Vm(i+1, is_spike) = Vr;
                    RP_ind(is_spike) = RP(sign(population_type(is_spike)-3)+2);
                    t_channel(i, is_spike) = 0;
                    t_channel(i+1, :) = t_channel(i, :)+dt;
                    % choose the neurons that are not in refractory period,
                    % calculate dvdt and update membrane potential
                    non_rp = (RP_ind==0);
                    in_rp = (RP_ind~=0)&(~is_spike);
                    dvdt = (-gL(sign(population_type-3)+2).*(Vm(i, :)-EL)+ ...
                                gL(sign(population_type-3)+2)*dT.*exp((Vm(i, :)-VT)/dT) ...
                                +I_ch+I_ext)./C(sign(population_type-3)+2);
                    Vm(i+1, non_rp) = Vm(i, non_rp)+dvdt(non_rp)*dt;
                    % for those in refractory period, reduce the time they
                    % will still be in it
                    RP_ind(in_rp) = RP_ind(in_rp)-1;
                end
                % save membrane potential change over time
                parsave([save_path, '/trial', num2str(trial), '.mat'], Vm);
            end
        end
    end
end
toc
%% 
function I = synapse_current(V, t)
    tau_AMPA = 2e-3;
    tau_NMDA_1 = 2e-3;
    tau_NMDA_2 = 100e-3;
    tau_GABA = 5e-3;
    Mg = 1e-3;
    E_AMPA = -70e-3;
    E_NMDA = 0;
    E_GABA = 0;
    I = zeros(length(V), 3);
    I(:, 1) = exp(-t/tau_AMPA).*(V-E_AMPA);
    g_NMDA = tau_NMDA_2/(tau_NMDA_2-tau_NMDA_1) ...
        *(exp(-t/tau_NMDA_1)-exp(-t/tau_NMDA_2));
    NMDA = g_NMDA.*(V-E_NMDA)./(1+Mg*exp(-0.062*V)/3.57);
    I(:, 2) = NMDA;
    I(:, 3) = exp(-t/tau_GABA).*(V-E_GABA);
    I = I';
end
%%
function [sTrain] = pSpike(fr, bins, trials, dt, I_sti)
    s = double(rand(trials, bins) < fr*dt);
    sTrain = zeros(trials, bins);
    for i = 1 : trials
        temp = conv(s(i, :), I_sti);
        sTrain(i, :) = temp(1:length(temp)-length(I_sti)+1);
    end
end
%%
function [freq] = mapping(step, std, k, rand_on)
    freq = zeros(1, 2);
    freq(1) = step*(4+k)+rand_on*std*(rand-0.5);
    freq(2) = step*(4-k)+rand_on*std*(rand-0.5);
    freq(freq<0) = 0;
end
%%
function [] = parsave(dir, Vm)
    save(dir, 'Vm');
end