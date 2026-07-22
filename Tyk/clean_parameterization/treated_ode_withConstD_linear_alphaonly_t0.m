function dSdt = treated_ode_withConstD_linear_alphaonly_t0(t, S, alpha, ce, rS, KS)
    % t in DAYS, ce in uM
    D_t = constant_D_schedule(t, ce);  % always = Ce
    
    % linear killing function
    alpha_t = alpha * t; %no beta
    survival = alpha_t * D_t;    
    
    dSdt = rS * S * (1 - S / KS) - survival * S;
end