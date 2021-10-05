%%% Paul Adkisson
%%% 9.6.21
%%% Purpose: Calculate decision time and accuracy from population firing
%%% rates
load("bam_constants.mat")
start_trial = 73;
end_trial = 84;
num_trials = end_trial - start_trial + 1;
I_dcs_pyr = [-4, 0, 4]*1e-12; %pA
I_dcs_int = I_dcs_pyr / (-2);
I_dcs = [I_dcs_pyr; I_dcs_int];
batch_size = 12;
num_batches = end_trial / batch_size;
num_rand_iter = 0;
for I_dc = I_dcs
    fprintf("I_dc: %0.0fpA \n", I_dc(1)*1e12)
    output_dcpath = sprintf("data/I_dc=%0.0fpA", I_dc(1)*1e12);
    decisions = zeros(num_trials, length(coherences)); %0 for no decision, 1 for C1, 2 for C2
    decision_times = zeros(num_trials, length(coherences));
    avg_dts = zeros(length(coherences), 1);
    std_dts = zeros(length(coherences), 1);
    avg_correct_dts = zeros(length(coherences), 1);
    std_correct_dts = zeros(length(coherences), 1);
    avg_incorrect_dts = zeros(length(coherences), 1);
    std_incorrect_dts = zeros(length(coherences), 1);
    avg_acc = zeros(length(coherences), 1);
    percent_nodec = zeros(length(coherences), 1);
    percent_earlydec = zeros(length(coherences), 1);
    tot_frs = zeros(length(coherences), num_trials, length(t), 4);
    avg_correct_frs = zeros(length(coherences), 1, length(t), 4);
    std_correct_frs = zeros(length(coherences), 1, length(t), 4);
    avg_incorrect_frs = zeros(length(coherences), 1, length(t), 4);
    std_incorrect_frs = zeros(length(coherences), 1, length(t), 4); 
    avg_nodec_frs = zeros(length(coherences), 1, length(t), 4);
    std_nodec_frs = zeros(length(coherences), 1, length(t), 4);
    avg_earlydec_frs = zeros(length(coherences), 1, length(t), 4);
    std_earlydec_frs = zeros(length(coherences), 1, length(t), 4);
    avg_correct_peaks = zeros(length(coherences), 4);
    std_correct_peaks = zeros(length(coherences), 4);
    avg_incorrect_peaks = zeros(length(coherences), 4);
    std_incorrect_peaks = zeros(length(coherences), 4);
    avg_nodec_peaks = zeros(length(coherences), 4);
    std_nodec_peaks = zeros(length(coherences), 4);
    dec_thresholds_lb = zeros(length(coherences), num_trials);
    dec_thresholds_ub = zeros(length(coherences), num_trials);
    for i = 1:length(coherences)
        c = coherences(i);
        fprintf("Coherence: %0.1f%% \n", c*100)
        output_coherentpath = strcat(output_dcpath, sprintf("/c=%0.3f", c));
        for trial = start_trial:end_trial
            output_trialpath = strcat(output_coherentpath, sprintf("/trial%0.0f.mat", trial));
            load(output_trialpath, "pop_frs")
            relative_trial = trial - start_trial + 1;
            [decision, decision_idx] = get_decision(pop_frs);
            decisions(relative_trial, i) = decision;
            if decision_idx == -1
                decision_times(relative_trial, i) = Inf;
            else
                decision_times(relative_trial, i) = t(decision_idx) - t_task;
            end
            tot_frs(i, relative_trial, :, :) = pop_frs;
        end
        %omit trials with no decision & those that decide before onset of task-related input
        coherent_times = decision_times(decisions(:, i)~=0 & decision_times(:, i)>0, i);
        avg_dts(i) = mean(coherent_times);
        std_dts(i) = std(coherent_times);
        %omit trials with no decision & those that decide before onset of task-related input
        coherent_decisions = decisions(decisions(:, i)~=0 & decision_times(:, i)>0, i);
        avg_acc(i) = sum(coherent_decisions == 1) / length(coherent_decisions);
        percent_nodec(i) = sum(decisions(:, i)==0, 'all') / end_trial;
        percent_earlydec(i) = sum(decision_times(:, i)<=0, 'all') / end_trial;
        
        %Breakdown DTs by outcome
        avg_correct_dts(i) = mean(coherent_times(coherent_decisions==1));
        std_correct_dts(i) = std(coherent_times(coherent_decisions==1));
        avg_incorrect_dts(i) = mean(coherent_times(coherent_decisions==2));
        std_incorrect_dts(i) = std(coherent_times(coherent_decisions==2));
        
        %Breakdown by outcome (correct, incorrect, no-decision, early-decision)
        correct_frs = tot_frs(i, decisions(:, i)==1, :, :);
        incorrect_frs = tot_frs(i, decisions(:, i)==2, :, :);
        nodec_frs = tot_frs(i, decisions(:, i)==0, :, :);
        earlydec_frs = tot_frs(i, decision_times(:, i)<=0, :, :); 
        avg_correct_frs(i, :, :, :) = mean(correct_frs, 2);
        std_correct_frs(i, :, :, :) = std(correct_frs, 0, 2);
        avg_incorrect_frs(i, :, :, :) = mean(incorrect_frs, 2);
        std_incorrect_frs(i, :, :, :) = std(incorrect_frs, 0, 2);
        avg_incorrect_frs(i, :, :, :) = mean(incorrect_frs, 2);
        std_incorrect_frs(i, :, :, :) = std(incorrect_frs, 0, 2);
        avg_nodec_frs(i, :, :, :) = mean(nodec_frs, 2);
        std_nodec_frs(i, :, :, :) = std(nodec_frs, 0, 2);
        avg_earlydec_frs(i, :, :, :) = mean(earlydec_frs, 2);
        std_nodec_frs(i, :, :, :) = std(earlydec_frs, 0, 2);
        
        %Record peak frs split by outcome
        peak_frs = max(tot_frs(i, :, :, :), [], 3);
        correct_peaks = peak_frs(:, decisions(:, i)==1, :, :);
        incorrect_peaks = peak_frs(:, decisions(:, i)==2, :, :);
        nodec_peaks = peak_frs(:, decisions(:, i)==0, :, :);
        if ~isempty(correct_peaks)
            avg_correct_peaks(i, :) = mean(correct_peaks, 2);
            std_correct_peaks(i, :) = std(correct_peaks, 0, 2);
        end
        if ~isempty(incorrect_peaks)
            avg_incorrect_peaks(i, :) = mean(incorrect_peaks, 2);
            std_incorrect_peaks(i, :) = std(incorrect_peaks, 0, 2);
        end
        if ~isempty(nodec_peaks)
            avg_nodec_peaks(i, :) = mean(nodec_peaks, 2);
            std_nodec_peaks(i, :) = std(nodec_peaks, 0, 2);
        end
        
        %Decision Thresholds
        dec_thresholds_lb(i, decisions(:, i)==1) = max(tot_frs(i, decisions(:, i)==1, :, 2), [], 3);
        dec_thresholds_lb(i, decisions(:, i)==2) = max(tot_frs(i, decisions(:, i)==2, :, 1), [], 3);
        dec_thresholds_ub(i, decisions(:, i)==1) = max(tot_frs(i, decisions(:, i)==1, :, 1), [], 3);
        dec_thresholds_ub(i, decisions(:, i)==2) = max(tot_frs(i, decisions(:, i)==2, :, 2), [], 3);
        dec_thresholds_ub(i, decisions(:, i)==0) = max(dec_thresholds_ub(i, :), [], 'all'); %ensure that no-dec trials are not affecting threshold upper bound
    end
    
    %Batch trials
    batch_acc = zeros(length(coherences), num_rand_iter*num_batches);
    batch_nodec = zeros(length(coherences), num_rand_iter*num_batches);
    batch_earlydec = zeros(length(coherences), num_rand_iter*num_batches);
    for j = 1:num_rand_iter
        redo_iter = false;
        dowhile = true;
        while redo_iter || dowhile
            redo_iter = false;
            shuffled_trials = randperm(end_trial);
            for b = 1:num_batches
                batch_trials = shuffled_trials((b-1)*batch_size+1:b*batch_size);
                for i = 1:length(coherences)
                    c = coherences(i);
                    batch_decs = decisions(batch_trials, i);
                    coherent_decisions = batch_decs(batch_decs~=0 & decision_times(batch_trials, i)>0);
                    if isempty(coherent_decisions)
                        disp("No coherent decisions")
                        redo_iter = true;
                        break
                    end
                    batch_acc(i, b+(j-1)*num_rand_iter) = sum(coherent_decisions == 1) / length(coherent_decisions);
                    batch_nodec(i, b+(j-1)*num_rand_iter) = sum(decisions(batch_trials, i)==0, 'all') / batch_size;
                    batch_earlydec(i, b+(j-1)*num_rand_iter) = sum(decision_times(batch_trials, i)<=0, 'all') / batch_size;
                end
            end
            dowhile = false;
        end
    end
    batch_avg_acc = mean(batch_acc, 2);
    batch_std_acc = std(batch_acc, 0, 2);
    batch_avg_nodec = mean(batch_nodec, 2);
    batch_std_nodec = std(batch_nodec, 0, 2);
    batch_avg_earlydec = mean(batch_earlydec, 2);
    batch_std_earlydec = std(batch_earlydec, 0, 2);
    
    %Pre-stimulus bias
    prestim_bias = tot_frs(:, :, abs(t-t_task)<dt/2, 1) - tot_frs(:, :, abs(t-t_task)<dt/2, 2);
    
    %Decision Thresholds
    dec_thresh = zeros(2, 1);
    dec_thresh_idx = zeros(2, 2);
    [dec_thresh(1), idx] = max(dec_thresholds_lb, [], 'all', 'linear');
    [c_idx, trial_idx] = ind2sub(size(dec_thresholds_lb), idx);
    dec_thresh_idx(1, :) = [c_idx, trial_idx];
    [dec_thresh(2), idx] = min(dec_thresholds_ub, [], 'all', 'linear');
    [c_idx, trial_idx] = ind2sub(size(dec_thresholds_ub), idx);
    dec_thresh_idx(2, :) = [c_idx, trial_idx];
    
    %Fit to weibull FN
    %coeffs = lsqcurvefit(@weibull, [1e-9, 0], coherences, avg_acc');
    coeffs = [1e-9, 0];
    avg_frs = mean(tot_frs, 2);
    std_frs = std(tot_frs, 0, 2);
    decisionpath = strcat(output_dcpath, "/decisions.mat");
    save(decisionpath, "decisions", "decision_times", ...
        "avg_dts", "std_dts", "avg_acc", "percent_nodec", "coeffs", ...
        "avg_frs", "std_frs", "avg_correct_frs", "std_correct_frs", ...
        "avg_incorrect_frs", "std_incorrect_frs", ...
        "avg_nodec_frs", "std_nodec_frs", ...
        "percent_earlydec", "avg_earlydec_frs", "std_earlydec_frs", ...
        "avg_correct_dts", "std_correct_dts", "avg_incorrect_dts", "std_incorrect_dts", ...
        "batch_avg_acc", "batch_std_acc", "batch_avg_nodec", "batch_std_nodec", ...
        "batch_avg_earlydec", "batch_std_earlydec", "avg_correct_peaks", "std_correct_peaks", ...
        "avg_incorrect_peaks", "std_incorrect_peaks", ...
        "avg_nodec_peaks", "std_nodec_peaks", "prestim_bias", ...
        "dec_thresh", "dec_thresh_idx", "dec_thresholds_lb", "dec_thresholds_ub")
end

function [decision, decision_idx] = get_decision(pop_frs)
    decision_thresh = 15; %Hz
    decision1_idx = find(pop_frs(:, 1)>=decision_thresh, 1);
    decision2_idx = find(pop_frs(:, 2)>=decision_thresh, 1);
    if isempty(decision1_idx) && isempty(decision2_idx) %no decision
        decision = 0;
        decision_idx = -1;
    elseif isempty(decision1_idx)
        decision = 2;
        decision_idx = decision2_idx;
    elseif isempty(decision2_idx) || decision1_idx <= decision2_idx
        decision = 1;
        decision_idx = decision1_idx;
    else
        decision = 2;
        decision_idx = decision2_idx;
    end
end

function [decision, decision_idx] = get_diff_decision(pop_frs)
    %decision_thresh = 15; %Hz
    decision_thresh = 7.5;
    decision_idx = find(abs(pop_frs(:, 1)-pop_frs(:, 2))>=decision_thresh, 1);
    if isempty(decision_idx) %no decision
        decision = 0;
        decision_idx = -1;
    elseif pop_frs(decision_idx, 1) > pop_frs(decision_idx, 2)
        decision = 1;
    elseif pop_frs(decision_idx, 1) < pop_frs(decision_idx, 2)
        decision = 2;
    else
        assert(1==2, "Error!")
    end
end