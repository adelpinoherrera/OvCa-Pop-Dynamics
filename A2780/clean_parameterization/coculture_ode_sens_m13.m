function dydt = coculture_ode_sens_m13(t, y, rS, KS, delta, alphaS, nS, IC50S, rR, KR, gamma, alphaR, nR, IC50R, ce)
    S = y(1); R = y(2);
    D_t = constant_D_schedule(t, ce);  % always = Ce
    %hill kill terms, different for S (alpha is time dependent) and R
    kill_S = (D_t.^nS) ./ (IC50S.^nS + D_t.^nS);
    alpha_tS = alphaS * t;
    kill_R = (D_t.^nR) ./ (IC50R.^nR + D_t.^nR);
    alpha_tR = alphaR * t; 

    dSdt = rS * S * (1 - (S + delta*R)/KS) - alpha_tS * kill_S * S;
    dRdt = rR * R * (1 - (R + gamma*S)/KR) - alpha_tR * kill_R * R;
    dydt = [dSdt; dRdt];
end