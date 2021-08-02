clear;
load("tdcs_constants.mat")
ei = 0;
dc_type = 4;
fr_win = 5;
mavgs = [5, 50];
srf = 3;
conds = 1:6;

fr1_avgs = zeros(length(t), length(conds));
fr2_avgs = zeros(length(t), length(conds));
fr3_avgs = zeros(length(t), length(conds));

for i = conds
    if i == 1
        trialname = sprintf("data/brain=1 ei=%0.0f k=3 dc_type=%0.0f/trial1fr_srf%0.0f_frwin%0.0f.mat", [ei, dc_type, srf, fr_win]);
    elseif i == 2
        trialname = sprintf("data/brain=1 ei=%0.0f k=3 dc_type=%0.0f/trial1fr_srf%0.0f.mat", [ei, dc_type, srf]);
    elseif i == 3
        trialname = sprintf("data/brain=1 ei=%0.0f k=3 dc_type=%0.0f/trial1fr_srf%0.0f_gausswin%0.0f.mat", [ei, dc_type, srf, fr_win]);
    elseif i == 4
        trialname = sprintf("data/brain=1 ei=%0.0f k=3 dc_type=%0.0f/trial1fr_srf%0.0f_mavg%0.0f.mat", [ei, dc_type, srf, mavgs(1)]);
    elseif i == 5
        trialname = sprintf("data/brain=1 ei=%0.0f k=3 dc_type=%0.0f/trial1fr_srf%0.0f_frwin%0.0f_mavg%0.0f.mat", [ei, dc_type, srf, fr_win, mavgs(2)]);
    else
        trialname = sprintf("data/brain=1 ei=%0.0f k=3 dc_type=%0.0f/trial1fr_srf%0.0f_mavg%0.0f.mat", [ei, dc_type, srf, mavgs(2)]);
        disp(trialname)
    end
    load(trialname)
    fr1_avgs(:, i) = fr1_avg;
    fr2_avgs(:, i) = fr2_avg;
    fr3_avgs(:, i) = fr3_avg;
end
fr_avgs = zeros(length(t), length(conds), 3);
fr_avgs(:, :, 1) = fr1_avgs;
fr_avgs(:, :, 2) = fr2_avgs;
fr_avgs(:, :, 3) = fr3_avgs;

for g = 1:3
    figure;
    plot(t, fr_avgs(:, :, g))
    legend(["5ms Window", "Linear Interpolation", "5ms Window + 5ms Gaussian Filter", "5ms Window + 5ms Moving Avg Filter", "5ms Window + 50ms Moving Avg", "50ms Window + 50ms Moving Avg"])
    title(sprintf("ei=%0.1f, dc-type=%0.0f, k=3, group=%0.0f, srf=%0.0f", [ei, dc_type, g, srf]))
    xlabel("Time (s)")
    ylabel("Average Firing Rate (Hz)")
end