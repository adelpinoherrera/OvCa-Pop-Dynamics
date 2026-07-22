function p = setParameters()
%define parameters already estimated, in total 12 parameters
    
    % parameters: 
    %Sensitive data
    p.rS = 0.8547; %1
    %p.KS = 1114404.9924; %2
    p.IC50S = 0.9960177; %3
    p.nS = 3.7779; %4
    p.alphaS = 0.129; %5

    %Resistant data
    p.rR = 0.5553; %6
    %p.KR = 1356851.059; %7
    p.IC50R = 8.338535; %8
    p.nR = 2.5166; %9
    p.alphaR = 1.7982; %10

    %co-culture interaction parameters 
    %p.delta = 1.2244; %11
    %p.gamma = 0.9322; %12
    p.K= 990795.311; %

    %set the proportion of sensitive and resistant cells as a parameter too
    p.S0 = 50; %set to an initial proportion of 50/50 to not favor any population more than others
    p.R0 = 50;

end
