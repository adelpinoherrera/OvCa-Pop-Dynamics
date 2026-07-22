function [res, BIC_info] = joint_residuals_coculture_mono_sens2_m13( ...
    p, nS, IC50S, alphaR, nR, IC50R, ...
    t_co, S_co_all, R_co_all, ce_co, ...
    t_mono_S, S_mono_all, ce_mono_S, ...
    t_mono_R, R_mono_all, ce_mono_R)

    rS = p(1);
    K = p(2);
    alphaS = p(3);
    rR = p(4);

    res = [];

    res_co = [];
    for i = 1:numel(S_co_all)
        S_raw = S_co_all{i};
        R_raw = R_co_all{i};
        [~, nReps] = size(S_raw);

        for j = 1:nReps
            S0 = S_raw(1,j);
            R0 = R_raw(1,j);

            [S_mod, R_mod] = simulate_coculture_all_sens2_m13( ...
                t_co, S0, R0, rS, K, alphaS, nS, IC50S, ...
                rR, alphaR, nR, IC50R, ce_co(i));

            res_co = [res_co; S_raw(:,j) - S_mod; R_raw(:,j) - R_mod];
        end
    end

    res_mono_S = [];
    for i = 1:numel(S_mono_all)
        S_raw = S_mono_all{i};
        [~, nReps] = size(S_raw);

        for j = 1:nReps
            S0 = S_raw(1,j);
            S_mod = simulate_mono_sens2_m13( ...
                t_mono_S, S0, rS, K, alphaS, nS, IC50S, ce_mono_S(i));

            res_mono_S = [res_mono_S; S_raw(:,j) - S_mod];
        end
    end

    res_mono_R = [];
    for i = 1:numel(R_mono_all)
        R_raw = R_mono_all{i};
        [~, nReps] = size(R_raw);

        for j = 1:nReps
            R0 = R_raw(1,j);
            R_mod = simulate_mono_res_sens2_m13( ...
                t_mono_R, R0, rR, K, alphaR, nR, IC50R, ce_mono_R(i));

            res_mono_R = [res_mono_R; R_raw(:,j) - R_mod];
        end
    end

    res = [res_co; res_mono_S; res_mono_R];

    k = 4;
    BIC_info.co = calc_bic_from_residuals(res_co, k);
    BIC_info.monoS = calc_bic_from_residuals(res_mono_S, k);
    BIC_info.monoR = calc_bic_from_residuals(res_mono_R, k);
    BIC_info.total = calc_bic_from_residuals(res, k);
end
