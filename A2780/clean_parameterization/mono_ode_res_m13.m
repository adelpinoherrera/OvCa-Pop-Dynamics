function dRdt = mono_ode_res_m13(t, R, rR, K, alphaR, nR, IC50R, ce)
    D_t = constant_D_schedule(t, ce);
    kill_R = (D_t.^nR) ./ (IC50R.^nR + D_t.^nR);
    alpha_tR = alphaR * t;
    dRdt = rR * R * (1 - R / K) - alpha_tR * kill_R * R;
end