function dSdt = mono_ode_sens2_m14(t, S, rS, K, alphaS, IC50S, ce)
    D_t = constant_D_schedule(t, ce);
    kill_S = (D_t) ./ (IC50S + D_t);
    alpha_tS = alphaS * t;
    dSdt = rS * S * (1 - S / K) - alpha_tS * kill_S * S;
end