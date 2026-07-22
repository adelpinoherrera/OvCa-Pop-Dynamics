function dydt = coculture_ode_all_sens2_m14(t, y, rS, K, alphaS, IC50S, rR, alphaR, IC50R, ce)
    S = y(1); R = y(2);
    D_t = constant_D_schedule(t, ce);  % always = Ce
    %hill kill terms, different for S (alpha is time dependent) and R
    kill_S = (D_t) ./ (IC50S + D_t);
    alpha_tS = alphaS * t;
    kill_R = (D_t) ./ (IC50R + D_t);
    alpha_tR = alphaR * t; 

    dSdt = rS * S * (1 - (S + R)/K) - alpha_tS * kill_S * S;
    dRdt = rR * R * (1 - (R + S)/K) - alpha_tR * kill_R * R;
    dydt = [dSdt; dRdt];
end