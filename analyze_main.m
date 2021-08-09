load("tdcs_constants.mat")
do_save = true;
do_load = false;

% Analyze different ei/dc types
%dc_types = [1, 4];
dc_types = [4];
%eis = 0:0.2:0.8;
eis = [0];
num_trials = 1;
num_groups = 2;
dt_avgs = zeros(length(eis), length(dc_types), 2);
dt_stds = zeros(length(eis), length(dc_types), 2);
decisions = zeros(length(eis), length(dc_types), num_trials); %1 for C1, 2 for C2, 3 for no decision
for i = 1:length(dc_types)
    dc_type = dc_types(i);
    fprintf("dc_type=%0.0f \n", dc_type)
    for j = 1:length(eis)
        ei = eis(j);
        fprintf("ei=%0.1f \n", ei)
        setting_dts = zeros(num_trials, num_groups);
        wait = waitbar(1/num_trials, sprintf("trial 1 / %0.0f", num_trials));
        for trial = 1:num_trials
            wait = waitbar(trial/num_trials, wait, compose("trial %0.0f / %0.0f", [trial, num_trials]));
            if ei==0
                trialname = sprintf("data/brain=1 ei=%0.0f k=4 dc_type=%0.0f/trial%0.0f.mat", [ei, dc_type, trial]);
            else
                trialname = sprintf("data/brain=1 ei=%0.1f k=4 dc_type=%0.0f/trial%0.0f.mat", [ei, dc_type, trial]);
            end
            if do_load
                frname = split(trialname, ".mat");
                frname = append(frname(1), "fr.mat");
                load(frname)
            else
                tic;
                [fr1_avg, fr2_avg, ~] = get_firing_rates(trialname, false, true);
                toc
            end
            fr_avgs = [fr1_avg, fr2_avg];
            setting_dts(trial, :) = get_dt_exp(t, fr_avgs);
            if all(setting_dts(trial, :) == Inf)
                decisions(j, i, trial) = Inf;
            elseif setting_dts(trial, 2) == Inf
                decisions(j, i, trial) = 1;
            else
                decisions(j, i, trial) = 2;
            end
        end
        setting_decisions = decisions(j, i, :);
        if do_save
            savepath = split(trialname, "trial");
            savepath = append(savepath(1), "dt_data.mat");
            save(savepath, "setting_decisions", "setting_dts")
        end
        for group = 1:num_groups
            dt_avgs(j, i, group) = mean(setting_dts(setting_decisions==group, group));
            dt_stds(j, i, group) = std(setting_dts(setting_decisions==group, group));
        end
    end
end

if do_save
    %save("dt_varying_ei.mat", "dt_avgs", "dt_stds", "decisions")
end

function trial_dts = get_dt_thresh(t, fr_avgs)
    fr_thresh = 10;
    trial_dts = zeros(size(fr_avgs, 2), 1);
    for group = 1:size(fr_avgs, 2)
        t_over_thresh = t(fr_avgs(:, group)>fr_thresh);
        if isempty(t_over_thresh)
            trial_dts(group) = Inf;
        else
            trial_dts(group) = t_over_thresh(1);
        end
    end
end

function trial_dts = get_dt_exp(t, fr_avgs)
    trial_dts = zeros(size(fr_avgs, 2), 1);
    max_fr1 = max(fr_avgs(:, 1));
    max_fr2 = max(fr_avgs(:, 2));
    if max_fr1 >= max_fr2
        trial_dts(1) = t(find(fr_avgs(:, 1)>=max_fr1*(1-1/exp(1)), 1));
        trial_dts(2) = Inf;
    else
        trial_dts(2) = t(find(fr_avgs(:, 2)>=max_fr2*(1-1/exp(1)), 1));
        trial_dts(1) = Inf;
    end
end