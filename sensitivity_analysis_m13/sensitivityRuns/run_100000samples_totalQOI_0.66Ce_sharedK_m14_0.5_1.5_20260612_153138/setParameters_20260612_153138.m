function p = setParameters()
%define parameters already estimated, in total 12 parameters

    % parameters: 
    %Sensitive data
    p.rS = 0.6433; %1
    %p.KS = 545681.7362; %2
    p.IC50S = 1.974574; %3
    p.alphaS = 0.1666; %4

    %Resistant data
    p.rR = 0.5339; %6
    %p.KR = 421668.8619; %7
    p.IC50R = 4.210567; %8
    p.alphaR = 0.0006; %9

    %co-culture interaction parameters 
    %p.delta = 1.2244; %11
    %p.gamma = 0.9322; %12
    p.K= 557834.9748; %10

    %set the proportion of sensitive and resistant cells as a parameter too
    p.S0 = 50; %set to an initial proportion of 50/50 to not favor any population more than others
    p.R0 = 50;

end
