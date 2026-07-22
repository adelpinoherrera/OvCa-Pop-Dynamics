function [S_model, R_model] = simulate_coculture_sharedK_m14(t_data, S0, R0, rS, alphaS, IC50S, rR, K, alphaR, IC50R, ce)
    tspan = [min(t_data), max(t_data)];
    ode_fun = @(t,y) coculture_ode_sharedK_m14(t, y, rS, alphaS, IC50S, rR, K, alphaR, IC50R, ce);
    opts = odeset('RelTol',1e-8, 'AbsTol',1e-10);
    [t_sol, y_sol] = ode15s(ode_fun, tspan, [S0, R0], opts);
    S_model = interp1(t_sol, y_sol(:,1), t_data, 'pchip');
    R_model = interp1(t_sol, y_sol(:,2), t_data, 'pchip');
end
