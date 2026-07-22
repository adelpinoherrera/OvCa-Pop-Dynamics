function [res, BIC_conditions] = coculture_residuals_all_sharedK_m14(p, rS, alphaS, IC50S, rR, alphaR, IC50R, ...
    t_data, S_raw_all, R_raw_all, ce_all)
    
    K = p(1);
    res = [];
    BIC_conditions = zeros(6, 1); %for "6 conditions", 3 of 0uM and 3 of 1uM
    condition_names = cell(6,1);
    ratio_labels = {'75/25 S/R', '50/50 S/R', '25/75 S/R'};
    
    for i = 1:6
        S_raw = S_raw_all{i}; R_raw = R_raw_all{i};
        [~, nReps] = size(S_raw);
        
        res_S = []; res_R = [];
        for j = 1:nReps
            S0 = S_raw(1,j); R0 = R_raw(1,j);
            [S_mod, R_mod] = simulate_coculture_sharedK_m14(t_data, S0, R0, rS, alphaS, IC50S, ...
                rR, K, alphaR, IC50R, ce_all(i));
            res_S = [res_S; S_raw(:,j) - S_mod];
            res_R = [res_R; R_raw(:,j) - R_mod];
        end
        
        res_condition = [res_S; res_R];
        res = [res; res_condition];
        
        % BIC per condition (N=15×3×2=90 points, k=2 parameters)
        N_cond = length(res_condition);
        k = 1;  
        SSE_cond = sum(res_condition.^2);
        BIC_conditions(i) = N_cond * log(SSE_cond / N_cond) + k * log(N_cond);
        
        ce_vals = [0, 1.00];
        ratio_idx = mod(i-1,3) + 1;
        condition_names{i} = sprintf('Ce=%.2f_%s', ce_vals(ceil(i/3)), ratio_labels{ratio_idx});
    end
end