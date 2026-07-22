function dSdt = treated_ode_withConstD_mult_alphaonly_t0(t, S, alpha, ce, rS, KS, IC50)
    % t in DAYS, ce in uM
    D_t = constant_D_schedule(t, ce);  % always = Ce
    
    % Hill killing function
    kill = (D_t) ./ (IC50 + D_t);
    alpha_t = alpha * t; %no beta
    survival = (1 - alpha_t * kill);    
    
    dSdt = rS * S * (1 - S / KS) * survival;
end