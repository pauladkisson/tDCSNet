clear;
dt = 1e-5;
t = 0:dt:1;
num_cor = 400;
num_inter = 200;
num = 2*num_cor+num_inter;
type = [ones(1, num_cor), 2*ones(1, num_cor), 3*ones(1, num_inter)];
fold = dir('data/');

decision_time = zeros(1, length(fold)-2);
for setting = 3 : length(fold)
    trials = dir(['data/', fold(setting).name]);
    for trial =  3 : length(trials);
        dt_avg = 0;
        fr_avg = zeros(length(t), 3);
        load(['data/', fold(setting).name, '/', trials(trial).name]);
        for i = 1 : num
            [pks, locs] = findpeaks(Vm(:, i), 'MinPeakHeight', -25e-3);
            if isempty(locs)
                locs = [locs; length(t)];
            end
            fr_ins = 1./(diff([0; locs])*dt);
            fr = interp1([0; locs*dt], [0; fr_ins], t, 'linear', 0);
            fr_avg(:, type(i)) = fr_avg(:, type(i))+fr';
        end
%         dt_avg = dt_avg+t(find(fr_avg(:, 1)>0.9*max(fr_avg(:, 1)), 1));
        dt_avg = dt_avg+t(find(fr_avg(:, 1)/num_cor>20, 1)); %Finds the time at which firing rate exceeds 20Hz
    end
    if isempty(dt_avg)
        decision_time(setting-2) = Inf;
    else
        decision_time(setting-2) = dt_avg/(length(trials)-2);
    end
    break
end