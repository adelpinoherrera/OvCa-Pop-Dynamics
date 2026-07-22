function dSdt = treated_ode_withConstantD_death(t, S, n, alpha, d, ce, rS, KS, IC50)
    % t in DAYS, ce in uM
    D_t = constant_D_schedule(t, ce);  % always = Ce
    
    % Hill killing function
    kill = (D_t.^n) ./ (IC50.^n + D_t.^n);
    survival = 1 - 2 * alpha * kill;    % (1 - 2α(D))
    
    dSdt = rS * S * (1 - S / KS) * survival - d * S;
end