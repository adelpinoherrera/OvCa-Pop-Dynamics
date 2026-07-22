function dcdt = odefun(t, c, p) %QC'd
    % Inputs:
    %   ~ - time (called t in "Dynamics" in step 3), changed to t because
    %   some of the variables depend on t
    %   c - state variables in this case sensitive and resistant
    %   populations 
    %   p - parameter structure of 35 parameters
    % Outputs:
    % dcdt - Derivatives of the state variables
 
    % ###Step 1
    % # c : levels of 2 biological factors at t. 
    S  = c(1); % Sensitive cells
    R  = c(2); % Resistant cells
    T  = c(3); % total cells

    % ###Step 2
        %% -- The rate change of the 2 populations (Dynamics of the system) -- %%
        % Pathways for the dynamics %
    ce = 1.0;
    D_t = ce;  % always = Ce
    %hill kill terms, different for S (alpha is time dependent) and R
    kill_S = (D_t.^p.nS) ./ (p.IC50S.^p.nS + D_t.^p.nS);
    alpha_tS = p.alphaS * t; 
    kill_R = (D_t.^p.nR) ./ (p.IC50R.^p.nR + D_t.^p.nR);
    alpha_tR = p.alphaR * t;

    % ###Step 3
    % Dynamics 
    dy(1)  = p.rS * S * (1 - (S + p.delta*R)/p.KS) - alpha_tS * kill_S * S; %dSdt
    dy(2)  = p.rR * R * (1 - (R + p.gamma*S)/p.KR) - alpha_tR * kill_R * R; %dRdt
    dy(3) = dy(1) + dy(2); %total population 

    dcdt = [dy(1), dy(2), dy(3)]';
end