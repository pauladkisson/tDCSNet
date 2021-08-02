clear;
load("tdcs_constants.mat")
ei = 0;
dc_type = 4;

srfs = [1, 2, 3, 4, 5];
fr1_avgs = zeros(length(t), length(srfs));
fr2_avgs = zeros(length(t), length(srfs));
fr3_avgs = zeros(length(t), length(srfs));

for i = 1:length(srfs)
    srf = srfs(i);
    trialname = sprintf("data/brain=1 ei=%0.0f k=3 dc_type=%0.0f/trial1fr_srf%0.0f.mat", [ei, dc_type, srf]);
    load(trialname)
    fr1_avgs(:, i) = fr1_avg;
    fr2_avgs(:, i) = fr2_avg;
    fr3_avgs(:, i) = fr3_avg;
end
fr_avgs = zeros(length(t), length(srfs), 3);
fr_avgs(:, :, 1) = fr1_avgs;
fr_avgs(:, :, 2) = fr2_avgs;
fr_avgs(:, :, 3) = fr3_avgs;

for g = 1:3
    figure;
    plot(t, (fr_avgs(:, :, g)))
    legend(["srf1", "srf2", "srf3-alt", "srf4", "srf5"])
    title(sprintf("ei=%0.1f, dc-type=%0.0f, k=3, group=%0.0f", [ei, dc_type, g]))
    xlabel("Time (s)")
    ylabel("Average Firing Rate (Hz)")
end

mape1 = get_mape(fr1_avgs)*100;
mape2 = get_mape(fr2_avgs)*100;
mape3 = get_mape(fr3_avgs)*100;
runtimes = [385.9, 82.1, 37.4, 24.7, 17.8];
figure;
hold on
plot(1./srfs, mape1)
plot(1./srfs, mape2)
plot(1./srfs, mape3)
ylabel("Mean Percent Error (%)")
yyaxis right
plot(1./srfs, runtimes)
ylabel("1-trial Runtime (s)")
hold off
legend(["C1", "C2", "Int", "Runtime"])
xlabel("Normalized Network Size")
title("brain=1, ei=0, k=3, dc-type=4")


function mape = get_mape(fr_avgs)
    ctrl_fr = fr_avgs(:, 1);
    non_zero_fr = ctrl_fr~=0;
    ctrl_fr = ctrl_fr(non_zero_fr);
    mape = mean(abs(fr_avgs(non_zero_fr, :) - ctrl_fr) ./ ctrl_fr, 1);
end