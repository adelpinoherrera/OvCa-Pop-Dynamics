function dSdt = treated_ode_withConstD_add_alpha_t0(t, S, n, alpha, ce, rS, KS, IC50)
    % t in DAYS, ce in uM
    D_t = constant_D_schedule(t, ce);  % always = Ce
    
    % Hill killing function
    kill = (D_t.^n) ./ (IC50.^n + D_t.^n);
    alpha_t = alpha * t; %no beta
    survival = alpha_t * kill;    
    
    dSdt = rS * S * (1 - S / KS) - survival * S;
end