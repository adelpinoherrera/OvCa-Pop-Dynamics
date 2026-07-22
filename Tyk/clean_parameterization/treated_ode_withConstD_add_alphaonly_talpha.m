function dSdt = treated_ode_withConstD_add_alphaonly_talpha(t, S, alpha, ce, rS, KS, IC50)
    % t in DAYS, ce in uM
    D_t = constant_D_schedule(t, ce);  % always = Ce
    
    % Hill killing function
    kill = (D_t) ./ (IC50 + D_t);
    alpha_t = alpha * (1 + t); %no beta
    survival = alpha_t * kill;    
    
    dSdt = rS * S * (1 - S / KS) - survival * S;
end