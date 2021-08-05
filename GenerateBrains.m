function [] = GenerateBrains(num_cor, num_inter, num_sub, connectivity)
    num = 2*num_cor+num_inter;
    connectivity = num2cell(connectivity);
    [c_cor_l, c_cor_r, c_inter, c_cor_inter, c_inter_cor] = connectivity{:};
    C_thres = zeros(num, num);
    C_thres(1:num_cor, 1:num_cor) = c_cor_l;
    C_thres(num_cor+1:2*num_cor, num_cor+1:2*num_cor) = c_cor_r;
    C_thres(num-num_inter+1:num, num-num_inter+1:num) = c_inter;
    C_thres(1:num_cor, num_cor+1:2*num_cor) = -1;
    C_thres(num_cor+1:2*num_cor, 1:num_cor) = -1;
    C_thres(num-num_inter+1:num, 1:2*num_cor) = c_inter_cor;
    C_thres(1:2*num_cor, num-num_inter+1:num) = c_cor_inter;
    diag_idx = logical(eye(size(C_thres)));
    C_thres(diag_idx) = -1; %disallow self-connections
    for i = 1 : num_sub
        C_temp = rand(num, num);
        adja = zeros(num, num);
        adja(C_temp<=C_thres) = 1;
        mkdir("brains")
        save(sprintf('brains/%0.0f.mat', i), "adja");
    end
end

