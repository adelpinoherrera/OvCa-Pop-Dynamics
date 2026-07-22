function model = TykcocultureModel_top4_noT()
    syms rS rR K alphaS %no ce, IC50, n or alphaR, no t
    syms S R
    ce = 0.66; %same value as the sensitivity analysis
    IC50S = 1.974574;
    IC50R = 4.210567;
    alphaR = 0.0006;


    % Parameters (no initial conditions here)
    model.sym.p = [rS; rR; K; alphaS]; %no ce nor IC50s, n or alphaR, no t
    model.sym.x = [S; R];

    model.sym.g = [0; 0];

    D_t = ce;
    kill_S = (D_t) ./ (IC50S + D_t);
    kill_R = (D_t) ./ (IC50R + D_t);

    model.sym.xdot = [
        rS * S * (1 - (S + R)/K) - (alphaS) * kill_S * S %no t
        rR * R * (1 - (R + S)/K) - (alphaR) * kill_R * R %no t
    ];

    % Known initial conditions (set your actual values)
    S0 = 15000;
    R0 = 15000;
    model.sym.x0 = [S0; R0];

    % Observables
    model.sym.y = [S; R];
end

%ALL identifiable parameters genssiMain('TykcocultureModel_top4_noT',6);
%  ---> THE RESULTS ARE STORED IN: 
% /blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/identifiability/GenSSI/Examples/TykcocultureModel_top4_noT/run1
% => THE MODEL IS STRUCTURALLY GLOBALLY IDENTIFIABLE 
% Structurally globally identifiable parameters: 
%     rS
%     rR
%      K
% alphaS

%no identifiable parameters genssiMain('TykcocultureModel_top4_noT',2);
%  ---> THE RESULTS ARE STORED IN: 
% /blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/identifiability/GenSSI/Examples/TykcocultureModel_top4_noT/run2
% => THE MODEL IS STRUCTURALLY GLOBALLY IDENTIFIABLE 
% Structurally globally identifiable parameters: 
%     rS
%     rR
%      K
% alphaS