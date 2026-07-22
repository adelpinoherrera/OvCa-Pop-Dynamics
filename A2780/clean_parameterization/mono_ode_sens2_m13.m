function dSdt = mono_ode_sens2_m13(t, S, rS, K, alphaS, nS, IC50S, ce)
    D_t = constant_D_schedule(t, ce);
    kill_S = (D_t.^nS) ./ (IC50S.^nS + D_t.^nS);
    alpha_tS = alphaS * t;
    dSdt = rS * S * (1 - S / K) - alpha_tS * kill_S * S;
end