load("tdcs_constants.mat")
I = zeros(length(t), 3);
I_corrected = zeros(length(t), 3);
for i = 1:length(t)
    I(i, :) = synapse_current(-0.060, t(i));
    I_corrected(i, :) = synapse_current_corrected(-0.060, t(i));
end

figure;
plot(t, I)
title("Original Synaptic Currents (voltage clamp: -60mV)")
legend(["AMPA", "NMDA", "GABA"])

figure;
plot(t, I_corrected)
title("Corrected Synaptic Currents (voltage clamp: -60mV)")
legend(["AMPA", "NMDA", "GABA"])


function I = synapse_current(V, t)
    tau_AMPA = 2e-3;
    tau_NMDA_1 = 2e-3;
    tau_NMDA_2 = 100e-3;
    tau_GABA = 5e-3;
    Mg = 1e-3;
    E_AMPA = -70e-3;
    E_NMDA = 0;
    E_GABA = 0;
    I = zeros(length(V), 3);
    I(:, 1) = exp(-t/tau_AMPA).*(V-E_AMPA);
    g_NMDA = tau_NMDA_2/(tau_NMDA_2-tau_NMDA_1) ...
        *(exp(-t/tau_NMDA_1)-exp(-t/tau_NMDA_2));
    NMDA = g_NMDA.*(V-E_NMDA)./(1+Mg*exp(-0.062*V)/3.57);
    I(:, 2) = NMDA;
    I(:, 3) = exp(-t/tau_GABA).*(V-E_GABA);
    I = I';
end

function I = synapse_current_corrected(V, t)
    tau_AMPA = 2e-3;
    tau_NMDA_1 = 2e-3;
    tau_NMDA_2 = 100e-3;
    tau_GABA = 5e-3;
    Mg = 1e-3;
    E_AMPA = 0;
    E_NMDA = 0;
    E_GABA = -70e-3;
    I = zeros(length(V), 3);
    I(:, 1) = exp(-t/tau_AMPA).*(E_AMPA-V);
    g_NMDA = tau_NMDA_2/(tau_NMDA_2-tau_NMDA_1) ...
        *(exp(-t/tau_NMDA_2)-exp(-t/tau_NMDA_1));
    NMDA = g_NMDA.*(E_NMDA-V)./(1+Mg*exp(-0.062*V)/3.57);
    I(:, 2) = NMDA;
    I(:, 3) = exp(-t/tau_GABA).*(E_GABA-V);
    I = I';
    end