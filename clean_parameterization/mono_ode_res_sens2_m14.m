function dRdt = mono_ode_res_sens2_m14(t, R, rR, K, alphaR, IC50R, ce)
    D_t = constant_D_schedule(t, ce);
    kill_R = (D_t) ./ (IC50R + D_t);
    alpha_tR = alphaR * t;
    dRdt = rR * R * (1 - R / K) - alpha_tR * kill_R * R;
end