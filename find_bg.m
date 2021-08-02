%%% Test 0 input
load("data/brain=1 ei=0 k=0 dc_type=4/trial1fr.mat")
load("tdcs_constants.mat")

figure; 
hold on
plot(t, fr1_avg)
plot(t, fr2_avg)
plot(t, fr3_avg)
hold off
legend(["C1", "C2", "Int"])
title("Background input=0, task input=0, srf=3")
ylabel("Firing Rate (Hz)")
xlabel("Time (s)")

load("data/brain=1 ei=0 k=0 dc_type=4/trial1.mat")
figure;
plot(t, Vm*1000)
title("Background input=0, task input=0, srf=3")
ylabel("Membrane Potential (mV)")
xlabel("Time (s)")