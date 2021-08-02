%%% Paul Adkisson
%%% 7.28.21
%%% Purpose: Generate Poisson Spike Trains for various trials of a given
%%% frequency.
function GenerateSpikes(fr, t, num_neurons, trials, G_AMPA, G_NMDA)
    mkdir("spikes")
    dt = t(2) - t(1);
    for trial = 1:trials
        spikes = rand(num_neurons, length(t)) < (dt*fr);
        tau_AMPA = 2e-3;
        tau_NMDA_1 = 2e-3;
        tau_NMDA_2 = 100e-3;
        g_AMPA = exp(-t/tau_AMPA);
        g_NMDA = tau_NMDA_2/(tau_NMDA_2-tau_NMDA_1) ...
        *(exp(-t/tau_NMDA_2)-exp(-t/tau_NMDA_1));
        g__AMPA = zeros(num_neurons, length(t));
        g__NMDA = zeros(num_neurons, length(t));
        for j = 1:num_neurons
            temp_ampa = conv(g_AMPA, spikes(j, :));
            temp_nmda = conv(g_NMDA, spikes(j, :));
            g__AMPA(j, :) = temp_ampa(1:length(t));
            g__NMDA(j, :) = temp_nmda(1:length(t));
        end
        g_AMPA = g__AMPA*G_AMPA;
        g_NMDA = g__NMDA*G_NMDA;
        mkdir(sprintf("spikes/trial%0.0f", trial))
        save(sprintf("spikes/trial%0.0f/%0.0fHz_N=%0.0f", [trial, fr, num_neurons]), "spikes", "g_AMPA", "g_NMDA")
    end
end