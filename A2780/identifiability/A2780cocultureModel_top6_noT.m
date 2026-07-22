function model = A2780cocultureModel_top6_noT()
    syms rS rR K alphaS nS nR %no ce, IC50, n or alphaR, no t
    syms S R
    ce = 1.0; %same value as the sensitivity analysis
    IC50S = 0.9960177;
    IC50R = 8.338535;
    alphaR = 1.7982;
    %nS = 3.7779;
    %nR = 2.5166; 


    % Parameters (no initial conditions here)
    model.sym.p = [rS; rR; K; alphaS; nS; nR]; %no ce nor IC50s, n or alphaR, no t
    model.sym.x = [S; R];

    model.sym.g = [0; 0];

    D_t = ce;
    kill_S = (D_t.^nS) ./ (IC50S.^nS + D_t.^nS);
    kill_R = (D_t.^nR) ./ (IC50R.^nR + D_t.^nR);

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

%no identifiable parameters genssiMain('A2780cocultureModel_all',2);
%  ---> THE RESULTS ARE STORED IN: 
% /blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/identifiability/GenSSI/Examples/A2780cocultureModel_top6_noT/run1

%no identifiable parameters genssiMain('A2780cocultureModel_all',6);
%  ---> THE RESULTS ARE STORED IN: 
% /blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/identifiability/GenSSI/Examples/A2780cocultureModel_top4_noT/run2

%no identifiable parameters genssiMain('A2780cocultureModel_all',8);
%  ---> THE RESULTS ARE STORED IN: 
% /blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/identifiability/GenSSI/Examples/A2780cocultureModel_top4_noT/run3