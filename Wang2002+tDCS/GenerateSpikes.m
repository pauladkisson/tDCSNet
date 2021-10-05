%%% Paul Adkisson
%%% 9.2.21
%%% Purpose: Generate Poisson Spike Trains for various trials of a given
%%% frequency.
function GenerateSpikes(fr_bg, max_fr_task, coherences, f, N_E, N_I, t_task, ...
    t, start_trial, end_trial)
    mkdir("spikes")
    dt = t(2) - t(1);
    N = N_E + N_I;
    num_group = floor(f*N_E);
    g1_idx = 1:num_group;
    g2_idx = num_group+1:num_group*2;
    time_idx = t>=t_task;
    for c = coherences
        mkdir(sprintf("spikes/c=%0.3f%", c))
        fr_task1 = max_fr_task/2*(1 + c);
        fr_task2 = max_fr_task/2*(1 - c);
        fr_tot1 = fr_bg + fr_task1;
        fr_tot2 = fr_bg + fr_task2;
        for trial = start_trial:end_trial
            spikes = rand(length(t), N) < (dt*fr_bg);
            spikes(time_idx, g1_idx) = rand(sum(time_idx), num_group) < (dt*fr_tot1);
            spikes(time_idx, g2_idx) = rand(sum(time_idx), num_group) < (dt*fr_tot2);
            basepath = sprintf("spikes/c=%0.3f/trial%0.0f", [c, trial]);
            mkdir(basepath)
            save(strcat(basepath, "/input.mat"), "spikes", "-v7.3")
        end
    end
end