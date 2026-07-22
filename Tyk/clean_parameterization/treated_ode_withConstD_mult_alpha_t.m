function dSdt = treated_ode_withConstD_mult_alpha_t(t, S, n, alpha, beta, ce, rS, KS, IC50)
    % t in DAYS, ce in uM
    D_t = constant_D_schedule(t, ce);  % always = Ce
    
    % Hill killing function
    kill = (D_t.^n) ./ (IC50.^n + D_t.^n);
    alpha_t = alpha * (1 + beta * t); %beta is delayed toxicity
    survival = 1 - 2 * alpha_t * kill;    % (1 - 2α(D,t))
    
    dSdt = rS * S * (1 - S / KS) * survival;
end