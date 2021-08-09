load("tdcs_constants.mat")
for i = 1:num_trials
    load(sprintf("data/brain=1 ei=0 k=4 dc_type=4/trial%0.0ffr_srf5.mat", i))
    figure;
    hold on 
    plot(t, fr1_avg)
    plot(t, fr2_avg)
    plot(t, fr3_avg)
    hold off
    legend(["C1", "C2", "Int"])
    title("875Hz bgp, 875Hz bgi, k=4")
end

load("data/brain=1 ei=0 k=4 dc_type=4/trial1.mat")
for i = 1:4
    figure;
    plot(t, Vm(:, i)*1000)
    xlabel("Time (s)")
    ylabel("Membrane Potential (mV)")
end

figure;
plot(t, Vm(:, end)*1000)
xlabel("Time (s)")
ylabel("Membrane Potential (mV)")
title("Inhibitory Neuron")

figure;
plot(t, Vm(:, end-1)*1000)
xlabel("Time (s)")
ylabel("Membrane Potential (mV)")
title("Inhibitory Neuron")