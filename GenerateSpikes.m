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
        temp_ampa = cnv(g_AMPA', spikes');
        temp_nmda = cnv(g_NMDA', spikes');
        g_AMPA = temp_ampa(1:length(t), :)*G_AMPA;
        g_NMDA = temp_nmda(1:length(t), :)*G_NMDA;
        mkdir(sprintf("spikes/trial%0.0f", trial))
        save(sprintf("spikes/trial%0.0f/%0.0fHz_N=%0.0f", [trial, fr, num_neurons]), "spikes", "g_AMPA", "g_NMDA")
        disp("Ran Generate Spikes")
    end
end