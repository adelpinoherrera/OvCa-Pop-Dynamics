function IC = setInitialConditions(p) %set initial conditions of state variables
    S0 = p.S0;
    R0 = p.R0;
    T0 = S0 + R0;

    IC.S0 = S0;
    IC.R0 = R0;
    IC.T0 = T0; 
end
