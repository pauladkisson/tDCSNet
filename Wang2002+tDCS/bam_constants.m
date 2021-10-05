%%% Paul Adkisson 
%%% 9.2.21
%%% Purpose Define Constants for Wang 2002 Biophysical Attractor
%%% Model (BAM)
clear;
tic;

%% Network Parameters
percent_size = 1;
N_E = floor(1600 * percent_size);
N_I = floor(400 * percent_size);
N = N_E + N_I;
num_sub = 1;
f = 0.15; % Fraction of cortical neurons asctivated by one type of stimulus
p = 2; % Number of different types of stimuli
num_selective = floor(p*f*N_E);
num_group = floor(f*N_E);
w_plus = 1.7; % Strength of "strong" synapses in the BAM network
w_minus = 1 - f*(w_plus - 1)/(1-f); %Strength of "weak" synapses in BAM
w = 1; %Strength of normal synapses in BAM
GenerateBAM(N_E, N_I, f, p, w_plus, w_minus, w)
GenerateConductances(N_E, N_I)
pop_type = ones(N, 1);
pop_type(N_E+1:end) = 2; % population_type = 1 for pyr, 2 for int

%% Simulation Parameters
dt = 0.02e-3; %ms
t_span = 3;
t = 0:dt:t_span;
start_trial = 73;
end_trial = 84;

%% Input Parameters
fr_bg = 2400; %Hz
% Synaptic Conductance = [pyramidal, interneuron]
G_ampa_ext = [2.1, 1.62]*1e-9; %nS
%coherences = [3.2, 6.4, 12.8, 25.6, 51.2, 100] / 100;
coherences = [0] / 100;
max_fr_task = 80;
t_task = 1;
t_dc_on = 0;
t_dc_off = 1;
GenerateSpikes(fr_bg, max_fr_task, coherences, f, N_E, N_I, t_task, ...
    t, start_trial, end_trial);

%% LIF Parameters
% Parameter = [pyramidal, interneuron]
C = [0.5, 0.2]*1e-9; %nF
gL = [25, 20]*1e-9; %nS
EL = -70e-3; %mV
Vs = -50e-3; %mV
Vr = -55e-3; %mV
tau_r = [2, 1]*1e-3; %ms
syn_delay = 0.5e-3; %ms
refract_ind = floor(tau_r/dt);
delay_ind = floor(syn_delay/dt);
tau_AMPA = 2e-3; %ms
tau_NMDA_1 = 2e-3; %ms
tau_NMDA_2 = 100e-3; %ms
tau_GABA = 5e-3; %ms
alpha = 500; %Hz

%% tDCS Parameters
I_dcs_pyr = [-4, 0, 4]*1e-12; %pA
I_dcs_int = I_dcs_pyr / (-2);
I_dcs = [I_dcs_pyr; I_dcs_int];


%% Firing Rate Parameters
win_size = 5e-3;
win_index = floor(win_size / dt);
avg_win_size = 50e-3;
avg_win_index = floor(avg_win_size / dt);
decision_thresh = 15; %Hz
earlystop_thresh = 30; %Hz
t_earlystop = t_span;
 
%% Save
save("bam_constants.mat")
toc


