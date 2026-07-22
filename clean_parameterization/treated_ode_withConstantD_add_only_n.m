function dSdt = treated_ode_withConstantD_add_only_n(t, S, n, ce, rS, KS, IC50)
    % t in DAYS, ce in uM
    D_t = constant_D_schedule(t, ce);  % always = Ce
    
    % Hill killing function
    kill = (D_t.^n) ./ (IC50.^n + D_t.^n);
    survival = 2 * kill;    % (2(D))
    
    dSdt = rS * S * (1 - S / KS) - survival * S;
end