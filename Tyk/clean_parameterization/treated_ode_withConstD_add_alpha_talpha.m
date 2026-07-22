function dSdt = treated_ode_withConstD_add_alpha_talpha(t, S, n, alpha, ce, rS, KS, IC50)
    % t in DAYS, ce in uM
    D_t = constant_D_schedule(t, ce);  % always = Ce
    
    % Hill killing function
    kill = (D_t.^n) ./ (IC50.^n + D_t.^n);
    alpha_t = alpha * (1 + t); %beta is delayed toxicity
    survival = alpha_t * kill;    % (α(D,t))
    
    dSdt = rS * S * (1 - S / KS) - survival * S;
end