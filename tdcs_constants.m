%%% Define tDCS Constants
%%% Paul Adkisson
%%% 7.20.21
clear;
tic;
%% functions that only need to be run once
size_reduction_factor = 400;
num_cor = floor(800 / size_reduction_factor);
num_inter = floor(400 / size_reduction_factor);
num = 2*num_cor+num_inter;
num_sub = 12;
% the connectivity rates:
% [within cortical 1, within cortical 2, within interneuron, cortical to interneuron, interneuron to cortical]
srf = size_reduction_factor;
srf_connectivity_corrections = [sqrt(srf), sqrt(srf), sqrt(srf), srf, srf];
%srf_connectivity_corrections = srf;
%connectivity = [0.08, 0.08, 0.1, 0.1, 0.2].*srf_connectivity_corrections;
connectivity = [0.08, 0.08, 0.1, 0.1, 0.2].*srf_connectivity_corrections;
GenerateBrains(num_cor, num_inter, num_sub, connectivity);
%% basic parameters that may change across experiments
% parameters about the simulation time
dt = 1e-5;
tspan = 4; % the whole timespan
t_task = 1; % the time when task-related input comes in, also note that 
            % background input and DC input is present all the time
t_taskoff = 3;
%number of trials
num_trials = 1;
% parameters about the task-related input
rand_on = 0; % if random factor is turned on in the frequency of 
             % task-related input
step_freq = 10; % the difference of two task-related inputs will be 
                % 2*k*step_freq
                
std_freq = 10; % only effective when the rand option is tured on
% list of E/I balance which is the percentage of inhibitory neurons in
% cortical population
ei_balance = [0, 0.2, 0.4, 0.6, 0.8];
% generate max conductance map
G_AMPA_ext = [2.1, 1.62]*1e-9;
G_AMPA_rec = [0.05, 0.04]*1e-9;
G_NMDA = [0.165, 0.13]*1e-9;
G_GABA = [1.3, 1.0]*1e-9;
GenerateConductanceMaps(num_cor, num_inter, ei_balance);

%% setups for different experiment
% fixed k
fixed_k = 4;
if ismember(fixed_k, 0:4)
    loop_k = fixed_k;
else
    loop_k = 0:4;
end

% fixed E/I balance
% index of the E/I balance wanted in the list above
% anything else put the E/I balance loop through the list
fixed_ei = 1;
if ismember(fixed_ei, 1:length(ei_balance))
    loop_ei = fixed_ei;
else
    loop_ei = 1 : length(ei_balance);
end

% this part determines whether or not there is DC in the experiment, and
% to which part DC is applied, options are +/-1) whole cortical population; 
% +/-2) excitatory neurons in cortical population; +/-3) inhibitory neurons in 
% cortical population; 4) no DC
dc_type = 4;
dc_amp = 2e-12;
% marker of the populations
population_type = [ones(1, num_cor), 2*ones(1, num_cor), 3*ones(1, num_inter)];
switch dc_type
    case 1
        dc = dc_amp*repmat((population_type==1), length(loop_ei), 1);
    case -1
        dc = -dc_amp*repmat((population_type==1), length(loop_ei), 1);
    case 2
        dc = dc_amp*repmat((population_type==1), length(loop_ei), 1).*(neuron_type)';
    case -2
        dc = -dc_amp*repmat((population_type==1), length(loop_ei), 1).*(neuron_type)';
    case 3
        dc = dc_amp*repmat((population_type==1), length(loop_ei), 1).*(~neuron_type)';
    case -3
        dc = -dc_amp*repmat((population_type==1), length(loop_ei), 1).*(~neuron_type)';
    case 4
        dc = zeros(length(loop_ei), num);
end

% the subjects
% you can choose a certain subject, or the code will loop through all of
% them
fixed_brains = 1; 
if ismember(fixed_brains, 1:num_sub)
    loop_brains = fixed_brains;
else
    loop_brains = 1 : num_sub;
end
%% parameters that should not need change
% timespan
t = 0 : dt : tspan;
% Generate Poisson Spikes
for k = loop_k
    frs = mapping(step_freq, std_freq, k, rand_on);
    for fr = frs
        GenerateSpikes(fr, t, num_cor, num_trials, G_AMPA_ext(1), G_NMDA(1));
    end
end
bgp_fr = 725;
bgi_fr = 0;
GenerateSpikes(bgp_fr, t, num_cor*2, num_trials, G_AMPA_ext(1), G_NMDA(1));
GenerateSpikes(bgi_fr, t, num_inter, num_trials, G_AMPA_ext(2), G_NMDA(2));

% parameters of the neuron model
C = [0.5, 0.2]*1e-9; % the structure is [pyramidal neurons, interneurons]
gL = [25, 20]*1e-9;
EL = -70e-3;
dT = 3e-3;
VT = -55e-3;
VS = -20e-3;
Vr = -53e-3;
tau_r = [2, 1]*1e-3; % refractory period
tau_AMPA = 2e-3;
tau_NMDA_1 = 2e-3;
tau_NMDA_2 = 100e-3;
tau_GABA = 5e-3;
delay = 0.5e-3; % synapse delay
% count of refractory period and delay
RP = fix(tau_r/dt);
delay_ind = fix(delay/dt);

save("tdcs_constants.mat")
toc
%%
function [freq] = mapping(step, std, k, rand_on)
    freq = zeros(1, 2);
    freq(1) = step*(4+k)+rand_on*std*(rand-0.5);
    freq(2) = step*(4-k)+rand_on*std*(rand-0.5);
    freq(freq<0) = 0;
end