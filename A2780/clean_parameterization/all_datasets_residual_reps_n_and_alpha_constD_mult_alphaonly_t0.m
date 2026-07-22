function [res, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_mult_alphaonly_t0(p, rS, KS, IC50, ...
                                         t1,S1_raw,ce1, ...
                                         t2,S2_raw,ce2, ...
                                         t3,S3_raw,ce3)
    alpha = p(1);
    res = [];
    
    % IC50 dataset
    [n1, nReps1] = size(S1_raw);
    res_IC50 = [];
    for j = 1:nReps1
        S1 = S1_raw(:,j); S01 = S1(1);
        S1_model = simulate_treat_for_dataset_withConstD_mult_alphaonly_t0(alpha, t1, S01, ce1, rS, KS, IC50);
        res_IC50 = [res_IC50; S1 - S1_model];
    end
    res = [res; res_IC50];
    
    % BIC IC50
    N_IC50 = numel(res_IC50);
    k = 1;  % parameters (n, and alpha)
    SSE_IC50 = sum(res_IC50.^2);
    BIC_IC50 = N_IC50 * log(SSE_IC50 / N_IC50) + k * log(N_IC50);

    % IC25 dataset  
    [n2, nReps2] = size(S2_raw);
    res_IC25 = [];
    for j = 1:nReps2
        S2 = S2_raw(:,j); S02 = S2(1);
        S2_model = simulate_treat_for_dataset_withConstD_mult_alphaonly_t0(alpha, t2, S02, ce2, rS, KS, IC50);
        res_IC25 = [res_IC25; S2 - S2_model];
    end
    res = [res; res_IC25];
    
    % BIC IC25
    N_IC25 = numel(res_IC25);
    SSE_IC25 = sum(res_IC25.^2);
    BIC_IC25 = N_IC25 * log(SSE_IC25 / N_IC25) + k * log(N_IC25);
    
    % IC75 dataset
    [n3, nReps3] = size(S3_raw);
    res_IC75 = [];
    for j = 1:nReps3
        S3 = S3_raw(:,j); S03 = S3(1);
        S3_model = simulate_treat_for_dataset_withConstD_mult_alphaonly_t0(alpha, t3, S03, ce3, rS, KS, IC50);
        res_IC75 = [res_IC75; S3 - S3_model];
    end
    res = [res; res_IC75];
    
    % BIC IC75
    N_IC75 = numel(res_IC75);
    SSE_IC75 = sum(res_IC75.^2);
    BIC_IC75 = N_IC75 * log(SSE_IC75 / N_IC75) + k * log(N_IC75);
end