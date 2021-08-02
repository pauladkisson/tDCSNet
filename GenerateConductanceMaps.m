function [] = GenerateConductanceMaps(num_cor, num_inter, ei_balance)
    num = 2*num_cor+num_inter;
    G_AMPA_ext = [2.1, 1.62]*1e-9;
    G_AMPA_rec = [0.05, 0.04]*1e-9;
    G_NMDA = [0.165, 0.13]*1e-9;
    G_GABA = [1.3, 1.0]*1e-9;
    for ei = 1 : length(ei_balance)
        AMPA = zeros(num, num);
        NMDA = zeros(num, num);
        GABA = zeros(num, num);
        eir = ei_balance(ei);
        num_i = num_cor*eir;
        AMPA(:, :) = G_AMPA_ext(1);
        AMPA(:, num-num_inter+1:num) = G_AMPA_ext(2);
        for i = 1 : num
            if i<=2*num_cor
                AMPA(i, i) = G_AMPA_rec(1);
            else
                AMPA(i, i) = G_AMPA_rec(2);
            end
        end
        AMPA(num-num_inter+1:num, :) = 0;
        AMPA(num_cor-num_i+1:num_cor, 1:num_cor) = 0;
        AMPA(2*num_cor-num_i+1:2*num_cor, num_cor+1:2*num_cor) = 0;
        NMDA(:, :) = G_NMDA(1);
        NMDA(:, num-num_inter+1:num) = G_NMDA(2);
        NMDA(num-num_inter+1:num, :) = 0;
        NMDA(num_cor-num_i+1:num_cor, 1:num_cor) = 0;
        NMDA(2*num_cor-num_i+1:2*num_cor, num_cor+1:2*num_cor) = 0;
        GABA(:, :) = G_GABA(1);
        GABA(:, num-num_inter+1:num) = G_GABA(2);
        GABA(1:2*num_cor, :) = 0;
        GABA(num_cor-num_i+1:num_cor, 1:num_cor-num_i) = G_GABA(1);
        GABA(num_cor-num_i+1:num_cor, num_cor-num_i+1:num_cor) = G_GABA(2);
        GABA(2*num_cor-num_i+1:2*num_cor, num_cor+1:2*num_cor-num_i) = G_GABA(1);
        GABA(2*num_cor-num_i+1:2*num_cor, 2*num_cor-num_i+1:2*num_cor) = G_GABA(2);
        neuron_type = diag(AMPA)>0;
        mkdir("conductance_maps")
        save(sprintf('conductance_maps/ei%0.1f.mat', eir), 'AMPA', 'NMDA', 'GABA', 'neuron_type');
    end
end

