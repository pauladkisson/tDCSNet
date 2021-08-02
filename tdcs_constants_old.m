%%% Define tDCS Constants
%%% Paul Adkisson
%%% 7.20.21
clear;
%% functions that only need to be run once
size_reduction_factor = 3;
num_cor = floor(400 / size_reduction_factor);
num_inter = floor(200 / size_reduction_factor);
num = 2*num_cor+num_inter;
num_sub = 12;
% the connectivity rates. [within cortical 1, within cortical 2, 
% cortical to interneuron, interneuron to cortical, within interneuron]
srf = size_reduction_factor;
srf_connectivity_corrections = [sqrt(srf), sqrt(srf), srf, srf, sqrt(srf)];
connectivity = [0.08, 0.08, 0.1, 0.1, 0.2].*srf_connectivity_corrections;
GenerateBrains(num_cor, num_inter, num_sub, connectivity);
%% basic parameters that may change across experiments
% parameters about the simulation time
dt = 1e-5;
tspan = 1; % the whole timespan
t_task = 0; % the time when task-related input comes in, also note that 
            % background input and DC input is present all the time
%number of trials 
num_trials = 1;
% parameters about the task-related input
rand_on = 0; % if random factor is turned on in the frequency of 
             % task-related input
step_freq = 20; % the difference of two task-related inputs will be 
                % 2*k*step_freq
                
std_freq = 10; % only effective when the rand option is tured on
% list of E/I balance which is the percentage of inhibitory neurons in
% cortical population
ei_balance = [0, 0.2, 0.4, 0.6, 0.8];
% generate max conductance map
GenerateConductanceMaps(num_cor, num_inter, ei_balance);

%% parameters that should not need change
% a pre-calculated excitatory stimulation (I_sti)
load ext.mat
% timespan
t = 0 : dt : tspan;
figure;
plot(t, I_sti)
% parameters of the neuron model
C = [0.5, 0.2]*1e-9; % the structure is [pyramidal neurons, interneurons]
gL = [25, 20]*1e-9; 
EL = -70e-3;
dT = 3e-3;
VT = -55e-3;
VS = -20e-3;
Vr = -53e-3;
tau_r = [2, 1]*1e-3; % refractory period
delay = 0.5e-3; % synapse delay
% count of refractory period and delay
RP = fix(tau_r/dt);
delay_ind = fix(delay/dt);
% marker of the populations
population_type = [ones(1, num_cor), 2*ones(1, num_cor), 3*ones(1, num_inter)];

save("tdcs_constants.mat")