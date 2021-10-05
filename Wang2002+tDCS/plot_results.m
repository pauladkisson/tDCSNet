clear;
load("bam_constants.mat")
start_trial = 49;
end_trial = 60;
num_trials = end_trial - start_trial + 1;
I_dcs_pyr = [-4, 0, 4]*1e-12; %pA
I_dcs_int = I_dcs_pyr / (-2);
I_dcs = [I_dcs_pyr; I_dcs_int];
dc_colors = ["g", "k", "r"];
%{
% Membrane Potential
for I_dc = I_dcs
    for c = coherences
        load(sprintf("data/I_dc=%0.0fpA/c=%0.3f/trial1.mat", [I_dc(1)*1e12, c]))
        figure;
        hold on
        plot(t, Vm(:, 1)*1000)
        plot(t, Vm(:, num_group+1)*1000)
        plot(t, Vm(:, num_group*2+1)*1000)
        plot(t, Vm(:, end)*1000)
        hold off
        legend([compose("Selective Population %0.0f", 1:p), "Non-Selective", "Inhibitory"])
        title(sprintf("I-dc: %0.0fpA, Coherence: %0.1f%%", [I_dc(1)*1e12, c*100]))
        xlabel("Time (s)")
        ylabel("Membrane Potential (mV)")
    end
end
%}


%Example FRs
for I_dc = I_dcs
    for c = coherences
        if c ~= 0.064
            continue
        end
        load(sprintf("data/I_dc=%0.0fpA/c=%0.3f/trial54.mat", [I_dc(1)*1e12, c]))
        figure;
        hold on
        for i = 1:p+2
            plot(t, pop_frs(:, i))
        end
        hold off
        legend([compose("Selective Population %0.0f", 1:p), "Non-Selective", "Inhibitory"])
        title(sprintf("I-dc: %0.0fpA, Coherence: %0.1f%%", [I_dc(1)*1e12, c*100]))
    end
end
%}

%{
%Example FRs by outcome vs population type
for I_dc = I_dcs
    load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12))
    for c_idx = 1:length(coherences)
        c = coherences(c_idx);
        figure;
        for i = 1:p+2
            ax = subplot(p+2, 1, i);
            hold on
            for trial = start_trial:end_trial
                load(sprintf("data/I_dc=%0.0fpA/c=%0.3f/trial%0.0f.mat", [I_dc(1)*1e12, c, trial]));
                if decisions(trial, c_idx) == 0
                    plot(t, pop_frs(:, i), "k-")
                elseif decisions(trial, c_idx)==1
                    plot(t, pop_frs(:, i), "g-")
                else
                    plot(t, pop_frs(:, i), "r-")
                end
            end
            hold off
            title(sprintf("I-dc: %0.0fpA, Coherence: %0.1f%%, Pop %0.0f", [I_dc(1)*1e12, c*100, i]))
        end
    end
end
%}

%{
%Example FRs by outcome vs. DC
for i = 1:1
    for c_idx = 1:length(coherences)
        c = coherences(c_idx);
        figure;
        for I_dc_idx = 1:size(I_dcs, 2)
            I_dc = I_dcs(:, I_dc_idx);
            load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12), 'decisions')
            ax = subplot(3, 1, I_dc_idx);
            hold on
            for trial = start_trial:end_trial
                load(sprintf("data/I_dc=%0.0fpA/c=%0.3f/trial%0.0f.mat", [I_dc(1)*1e12, c, trial]), 'pop_frs');
                if decisions(trial, c_idx) == 0
                    plot(t, pop_frs(:, i), "k-")
                elseif decisions(trial, c_idx)==1
                    plot(t, pop_frs(:, i), "g-")
                else
                    plot(t, pop_frs(:, i), "r-")
                end
            end
            hold off
            title(sprintf("I-dc: %0.0fpA, Coherence: %0.1f%%, Pop %0.0f", [I_dc(1)*1e12, c*100, i]))
        end
    end
end
%}

%{
%Example FRs by trial vs. DC 
for i = 1:p+2
    for c_idx = 1:length(coherences)
        c = coherences(c_idx);
        figure;
        for trial = start_trial:end_trial
            ax = subplot(3, 4, trial-start_trial+1);
            title(sprintf("Trial %0.0f", trial-start_trial+1))
            if trial - start_trial == 0
                title({sprintf("Coherence: %0.1f%%, Pop %0.0f", [c*100, i]), ...
                    sprintf("Trial %0.0f", trial-start_trial+1)})
            end
            hold on
            for I_dc_idx = 1:size(I_dcs, 2)
                I_dc = I_dcs(:, I_dc_idx);
                load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12), "decisions") 
                load(sprintf("data/I_dc=%0.0fpA/c=%0.3f/trial%0.0f.mat", [I_dc(1)*1e12, c, trial]), 'pop_frs');
                if decisions(trial, c_idx) == 0
                    plot(t, pop_frs(:, i), dc_colors(I_dc_idx))
                elseif decisions(trial, c_idx)==1
                    plot(t, pop_frs(:, i), dc_colors(I_dc_idx))
                else
                    plot(t, pop_frs(:, i), dc_colors(I_dc_idx))
                end
            end
            hold off
        end
    end
end
%}

%{
% Avg FRs
default_colors = get(gca, 'colororder');
for I_dc = I_dcs
    load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12))
    for c_idx = 1:length(coherences)
        c = coherences(c_idx);
        figure;
        hold on
        for i = 1:p+2
            shadedErrorBar(t, avg_frs(c_idx, :, :, i), std_frs(c_idx, :, :, i), ...
                'lineprops', {'color', default_colors(i, :)})
        end
        hold off
        legend([compose("Selective Population %0.0f", 1:p), "Non-Selective", "Inhibitory"])
        title(sprintf("I-dc: %0.0fpA, Coherence: %0.1f%%", [I_dc(1)*1e12, c*100]))
    end
end
%}

%{
% Avg FRs by outcome (correct, incorrect, no-decision, early-decision)
figure;
default_colors = get(gca, 'colororder');
plot([], []) % dummy plot for default_colors
plot_coherences = coherences(1:2:5);
for I_dc = I_dcs
    load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12))
    figure;
    for c_idx = 1:length(plot_coherences)
        c = plot_coherences(c_idx);
        ax1 = subplot(length(plot_coherences), 4, 1+4*(c_idx-1));
        hold on
        for i = 1:p+2
            shadedErrorBar(t, avg_correct_frs(c_idx, :, :, i), std_correct_frs(c_idx, :, :, i), ...
                'lineprops', {'color', default_colors(i, :)})
        end
        hold off
        if c_idx ~= length(plot_coherences)
            xticks([])
        else
            xlabel({"";"Correct"})
        end
        if c_idx ~= 2
            ylabel({sprintf("Coherence: %0.1f%%", c*100);""})
        else
            ylabel({sprintf("Coherence: %0.1f%%", c*100);"Firing Rate (Hz)"})
        end
        
        ax2 = subplot(length(plot_coherences), 4, 2+4*(c_idx-1));
        hold on
        for i = 1:p+2
            shadedErrorBar(t, avg_incorrect_frs(c_idx, :, :, i), std_incorrect_frs(c_idx, :, :, i), ...
                'lineprops', {'color', default_colors(i, :)})
        end
        hold off
        yticks([])
        if c_idx ~= length(plot_coherences)
            xticks([])
        else
            xlabel({"Time (s)";"Incorrect"})
        end
        if c_idx == 1
            title(sprintf("I-dc = %0.0fpA", I_dc(1)*1e12))
        end
        
        ax3 = subplot(length(plot_coherences), 4, 3+4*(c_idx-1));
        hold on
        for i = 1:p+2
            shadedErrorBar(t, avg_nodec_frs(c_idx, :, :, i), std_nodec_frs(c_idx, :, :, i), ...
                'lineprops', {'color', default_colors(i, :)})
        end
        hold off
        yticks([])
        if c_idx ~= length(plot_coherences)
            xticks([])
        else
            xlabel({""; "No-Decision"})
        end
        
        ax4 = subplot(length(plot_coherences), 4, 4+4*(c_idx-1));
        hold on
        for i = 1:p+2
            shadedErrorBar(t, avg_earlydec_frs(c_idx, :, :, i), std_earlydec_frs(c_idx, :, :, i), ...
                'lineprops', {'color', default_colors(i, :)})
        end
        hold off
        yticks([])
        if c_idx ~= length(plot_coherences)
            xticks([])
        else
            xlabel({""; "Early-Decision"})
        end
        linkaxes([ax1, ax2, ax3, ax4], "y")
    end
    legend([compose("Selective Population %0.0f", 1:p), "Non-Selective", "Inhibitory"])
end
%}

%{
%Decision Times
figure;
hAx = axes;
hAx.XScale = 'log';
hold on
for I_dc = I_dcs
    load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12))
    errorbar(coherences, avg_dts, std_dts)
end
hold off
xlabel("Coherence")
ylabel("Decision Time (s)")
legend(compose("I-dc=%0.0fpA", I_dcs(1, :)*1e12))
%}

%{
%Decision Times by outcome
figure;
ax1 = subplot(3, 1, 1);
hold on
ax2 = subplot(3, 1, 2);
hold on
ax3 = subplot(3, 1, 3);
hold on
ax1.XScale = 'log';
ax2.XScale = 'log';
ax3.XScale = 'log';
for I_dc = I_dcs
    load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12))
    errorbar(coherences, avg_dts, std_dts, "Parent", ax1)
    errorbar(coherences, avg_correct_dts, std_correct_dts, "Parent", ax2)
    errorbar(coherences, avg_incorrect_dts, std_incorrect_dts, "Parent", ax3)
end
hold off
xlabel("Coherence")
ylabel("Decision Time (s)")
legend(ax1, compose("I-dc=%0.0fpA", I_dcs(1, :)*1e12))
title(ax1, "Total Trials")
title(ax2, "Correct Trials")
title(ax3, "Incorrect Trials")
%}

%{
%Accuracies
c = 0:0.01:1;
figure;
hAx = axes;
hAx.XScale = 'log';
hold on
dc_colors = ["g", "k", "r"];
idx = 1;
for I_dc = I_dcs
    load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12))
    %coeffs
    %scatter(coherences, avg_acc, dc_colors(idx))
    %plot(c, weibull(coeffs, c), dc_colors(idx))
    plot(coherences, avg_acc, strcat(dc_colors(idx), 'o-'))
    idx = idx + 1;
end
hold off
xlabel("Coherence")
ylabel("Accuracy")
%f = flipud(get(gca, 'Children'));
%legend([f(2), f(4), f(6)], "I-dc=-4pA", "I-dc=0pA", "I-dc=4pA")
legend(compose("I-dc=%0.0fpA", I_dcs(1, :)*1e12))
%}

%{
%Batchwise Accuracies
figure;
hAx = axes;
hAx.XScale = 'log';
hold on
idx = 1;
for I_dc = I_dcs
    load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12))
    errorbar(coherences, batch_avg_acc, batch_std_acc, dc_colors(idx))
    idx = idx + 1;
end
hold off
xlabel("Coherence")
ylabel("Accuracy")
legend(compose("I-dc=%0.0fpA", I_dcs(1, :)*1e12))
%}

%{
%Percent No-decision
figure;
hAx = axes;
hAx.XScale = 'log';
hold on
idx = 1;
for I_dc = I_dcs
    load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12))
    plot(coherences, percent_nodec, strcat(dc_colors(idx), 'o-'))
    idx = idx + 1;
end
hold off
xlabel("Coherence")
ylabel("Percent No-Decisions")
legend(compose("I-dc=%0.0fpA", I_dcs(1, :)*1e12))
%}

%{
%Batchwise No-decision
figure;
hAx = axes;
hAx.XScale = 'log';
hold on
idx = 1;
for I_dc = I_dcs
    load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12))
    errorbar(coherences, batch_avg_nodec, batch_std_nodec, dc_colors(idx))
    idx = idx + 1;
end
hold off
xlabel("Coherence")
ylabel("Accuracy")
legend(compose("I-dc=%0.0fpA", I_dcs(1, :)*1e12))
%}

%{
%Percent Early-decision
figure;
hAx = axes;
hAx.XScale = 'log';
hold on
for I_dc = I_dcs
    load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12))
    %scatter(coherences, avg_acc)
    %plot(c, weibull(coeffs, c))
    percent_earlydec
    plot(coherences, percent_earlydec, 'o-')
end
hold off
xlabel("Coherence")
ylabel("Percent Early-Decisions")
legend(compose("I-dc=%0.0fpA", I_dcs(1, :)*1e12))
%}

%{
%Stacked Barplot of outcomes vs. coherence & Idc
%Percent No-decision
bar_data = zeros(length(coherences), 3, 3);
for I_dc_idx = 1:size(I_dcs, 2)
    I_dc = I_dcs(:, I_dc_idx);
    load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12))
    percent_truedec = (1-percent_nodec-percent_earlydec);
    bar_data(:, I_dc_idx, 1) = percent_truedec.*avg_acc;
    bar_data(:, I_dc_idx, 2) =  percent_truedec.*(1 - avg_acc);
    bar_data(:, I_dc_idx, 3) = percent_nodec;
end
h = plotBarStackGroups(bar_data, num2cell(compose("%0.1f%%", coherences*100)));
green = [0, 1, 0];
black = [0, 0, 0];
red = [1, 0, 0];
set(h(:, 1), "FaceColor", green)
set(h(:, 2), "FaceColor", red)
set(h(:, 3), "FaceColor", black)
xlabel("Coherence (-4pA | 0pA | 4pA)")
ylabel("Fraction of outcomes")
legend(["Correct", "Incorrect", "No-Decision"])
%}

%{
%Amplitude of peak FR
for i = 1:p+2
    figure;
    ax1 = subplot(3, 1, 1);
    hold on
    ax2 = subplot(3, 1, 2);
    hold on
    ax3 = subplot(3, 1, 3);
    hold on
    ax1.XScale = 'log';
    ax2.XScale = 'log';
    ax3.XScale = 'log';
    for I_dc = I_dcs
        load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12))
        errorbar(coherences, avg_correct_peaks(:, i), std_correct_peaks(:, i), "Parent", ax1)
        errorbar(coherences, avg_incorrect_peaks(:, i), std_incorrect_peaks(: , i), "Parent", ax2)
        errorbar(coherences, avg_nodec_peaks(:, i), std_nodec_peaks(:, i), "Parent", ax3)
    end
    hold off
    xlabel(sprintf("Coherence, Population %0.0f", i))
    ylabel("Peak Firing Rate (Hz)")
    legend(ax1, compose("I-dc=%0.0fpA", I_dcs(1, :)*1e12))
    title(ax1, "Correct Trials")
    title(ax2, "Incorrect Trials")
    title(ax3, "No-decision Trials")
end
%}


%Pre-Stimulus Bias
fig = figure;
hold on
invisible_fig = figure('Visible', 'Off');
hold on
for I_dc_idx = 1:size(I_dcs, 2)
    I_dc = I_dcs(:, I_dc_idx);
    load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12), ...
        "prestim_bias", "decisions")
    prestim_bias = prestim_bias(1, decisions(:, 1)~=0);
    decisions = decisions(decisions(:, 1)~=0, 1); %omit no-dec trials
    set(0, 'CurrentFigure', invisible_fig)
    h_tot = histogram(prestim_bias, 9);
    h1 = histogram(prestim_bias(decisions==1), h_tot.BinEdges);
    prestim_plot = (h1.BinEdges(1:end-1) + h1.BinEdges(2:end))/2;
    figure(fig)
    plot(prestim_plot, h1.Values ./ h_tot.Values, 'o-')
end
hold off
xlabel("Pre-Stimulus Bias (Hz)")
ylabel("Trials when P1 is selected (%)")
legend(compose("I-dc=%0.0fpA", I_dcs(1, :)*1e12))
title("Coherence 0%, t-I-dc = (0, 1), 12 Trials")
set(0, 'CurrentFigure', invisible_fig)
hold off
%}

%{
%Threshold
figure;
hold on
for I_dc_idx = 1:size(I_dcs, 2)
    I_dc = I_dcs(:, I_dc_idx);
    load(sprintf("data/I_dc=%0.0fpA/decisions.mat", I_dc(1)*1e12), ...
        "dec_thresh", "dec_thresh_idx", "decisions")
    bar(I_dc_idx, mean(dec_thresh), dc_colors(I_dc_idx))
    errorbar(I_dc_idx, mean(dec_thresh), diff(dec_thresh)/2, "k")
    dec_thresh_idx
    dec_thresh
end
%}