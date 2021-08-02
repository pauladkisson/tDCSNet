%Plot dt vs. ei results
%load("dt_varying_ei.mat", "dt_avgs", "dt_stds", "decisions")
%eis = 0:0.2:0.8;
figure;
hold on 
errorbar(eis, dt_avgs(:, 1, 1), dt_stds(:, 1, 1))
xlabel("E-I Balance")
ylabel("Decision Time (s)")
legend(["Hyperpolarizing DC Stim"])
title("Exponential-based Decision Time (5 trials)")

figure;
size(decisions)
acc = reshape(decisions(:, 1, : ), 5, 3);
acc = mean(acc, 2);
plot(eis, acc)
xlabel("ei balance")
ylabel("Accuracy")
title("")