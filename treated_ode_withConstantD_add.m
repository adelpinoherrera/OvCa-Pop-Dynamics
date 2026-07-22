function dSdt = treated_ode_withConstantD_add(t, S, n, alpha, ce, rS, KS, IC50)
    % t in DAYS, ce in uM
    D_t = constant_D_schedule(t, ce);  % always = Ce
    
    % Hill killing function
    kill = (D_t.^n) ./ (IC50.^n + D_t.^n);
    
    growth = rS * S * (1 - S / KS);
    death_rate = alpha * kill * S; 
    dSdt = growth - death_rate;
end