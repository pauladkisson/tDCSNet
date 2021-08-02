function [fr1_avg, fr2_avg, fr3_avg] = get_firing_rates(trial_fullpath, do_plot, do_save)
    % Constants
    load("tdcs_constants.mat")
    
    % load trial
    load(trial_fullpath)
    
    % Compute firing rate for each neuron (50ms window)
    win_size = 0.050;
    frs = zeros(length(t), num);
    for i = 1 : num
        [~, locs] = findpeaks(Vm(:, i), 'MinPeakHeight', -25e-3);
        peak_times = t(locs);
        spike_count = zeros(length(t), 1);
        
        for j = 1:length(peak_times)
            peak_time = peak_times(j);
            window = t>=(peak_time-win_size/2) & t<= (peak_time+win_size/2);
            spike_count(window) = spike_count(window) + 1;
        end
        frs(:, i) = spike_count / win_size; 
        frs(1:floor(win_size/dt), i) = sum(peak_times <= win_size) / win_size;
        frs(end-floor(win_size/dt):end, i) = sum(peak_times >= (t(end) - win_size)) / win_size;
        %{
        % Shuming's old version
        if isempty(locs)
            locs = [locs; length(t)];
        end
        fr_ins = 1./(diff([0; locs])*dt);
        fr = interp1([0; locs*dt], [0; fr_ins], t, 'linear', 0);
        frs(:, i) = fr;
        %}
    end
    fr1_avg = mean(frs(:, population_type==1), 2); %Average instaneous firing rate of group 1
    fr1_std = std(frs(:, population_type==1), [], 2); 
    fr2_avg = mean(frs(:, population_type==2), 2); %Average instaneous firing rate of group 2
    fr2_std = std(frs(:, population_type==2), [], 2);
    fr3_avg = mean(frs(:, population_type==3), 2); %Average instaneous firing rate of group 3
    fr3_std = std(frs(:, population_type==3), [], 2);
    
    %{
    % Using gaussian window filter
    w = gausswin(floor(win_size/dt));
    w = w ./ (w'*w);
    %}
    % Simple Moving Average filter
    win_size = 0.050;
    w = ones(floor(win_size/dt), 1);
    w = w ./ length(w);
    fr1_avg = filter(w, 1, fr1_avg);
    fr2_avg = filter(w, 1, fr2_avg);
    fr3_avg = filter(w, 1, fr3_avg);
    
    if do_plot
        figure;
        hold on
        errorcloud(t, fr1_avg, fr1_std, [0.8, 0, 0])
        errorcloud(t, fr2_avg, fr2_std, [0, 0.8, 0])
        errorcloud(t, fr3_avg, fr3_std, [0, 0, 0.8])
        children = get(gca, "Children");
        reordered_ch = [children(1:2:end-1), children(2:2:end)];
        set(gca, "Children", reordered_ch)
        hold off
        legend(["group 1 SD", "group 2 SD", "group 3 SD", "group 1 mean", "group 2 mean", "group 3 mean"])
    end
    
    if do_save
        savepath = split(trial_fullpath, ".mat");
        savepath = append(savepath(1), "fr_srf5.mat");
        save(savepath, "fr1_avg", "fr1_std", "fr2_avg", "fr2_std", "fr3_avg", "fr3_std")
    end
end

function errorcloud(x, y, err, color)
    y_high = y + err;
    y_low = y - err;
    [c_max, idx] = max(color);
    fill_color = c_max*ones(1, 3)/2;
    fill_color(idx) = fill_color(idx)*1.5;
    fill([x; flipud(x)], [y_high; flipud(y_low)], fill_color, "linestyle", "none");
    plot(x, y, "color", color)
end