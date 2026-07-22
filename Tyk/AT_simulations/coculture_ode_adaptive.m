function dydt = coculture_ode_adaptive(t, y, p, D_t)
    S = y(1);
    R = y(2);

    kill_S = (D_t) ./ (p.IC50S + D_t);
    kill_R = (D_t) ./ (p.IC50R + D_t);

    alpha_tS = p.alphaS * t;
    alpha_tR = p.alphaR * t;

    dSdt = p.rS * S * (1 - (S + R)/p.K) - alpha_tS * kill_S * S;
    dRdt = p.rR * R * (1 - (R + S)/p.K) - alpha_tR * kill_R * R;

    dydt = [dSdt; dRdt];
end