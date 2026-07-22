%This file calculates the S1, S2, and ST values using the 
% samples generated in getSamplesSobol.
%It can also calculate the confidence intervals for these values, but that is commented out. 
%Author: Jaimit Parikh

function results = calculateSobolIndices(Y, k, N, numResamples, confLevel)
arguments
    Y % estimated output value at the samples
    k % dimension of input parameters
    N % number of base samples
    numResamples = 1; % # of samples used to calculate confidence intervals for S1, S2, ST
    confLevel = 0.95;
end

% Normalize the model output
% Estimates of the Sobol' indices can be 
% biased for non-centered outputs
Y = (Y - mean(Y))./std(Y);

r = randi(N, N, numResamples);
Z = norminv(0.5 + confLevel / 2);

[A, B, AB, BA] = splitOutputsIntoMatrices(Y,N,k);

ABall = mat2cell(AB, repmat(N, 1, k), size(Y, 2));
f = @(x) firstOrder(A, B, x);
results.S1 = cellfun(f, ABall, 'UniformOutput',false);

f = @(x)Z*std(firstOrder(A(r), B(r), x(r)));
results.S1conf = cellfun(f, ABall, 'UniformOutput', false);

% this block calculates the confidence interval for S1
% results.S1conf = cell(k, 1);
% for ii = 1:k
%     ABi = ABall{ii};
%     results.S1conf{ii} = std(Z * firstOrder(A(r), B(r), ABi(r)));
% end

f = @(x) totalOrder(A,B, x);
results.ST = cellfun(f, ABall, 'UniformOutput',false);


f = @(x)Z*std(totalOrder(A(r), B(r), x(r)));
results.STconf = cellfun(f, ABall, 'UniformOutput', false);

if ~isempty(BA)
    Sij = zeros(k, k); 
    Sijconf = zeros(k, k);
    ABall = mat2cell(AB, repmat(N, 1, k), ...
        size(Y, 2));
    BAall = mat2cell(BA, repmat(N, 1, k),...
        size(Y, 2));
    for i = 1:k
        for j = i+1:k
            Sij(i, j) = secondOrder(A, B, ...
                ABall{i}, ABall{j}, BAall{i});

            Sijconf(i, j) = Z*std(secondOrder(A(r), B(r), ...
                ABall{i}(r), ABall{j}(r), BAall{i}(r)));
        end
    end
    % this block calculates the confidence interval for S2
    % ij = nchoosek(1:k, 2); % all ij combinations
    % f = @(x, y, z)secondOrder(A, B, x, y, z);
    % Sij = cellfun(f,  ABall(ij(:, 1)),...
    %     ABall(ij(:, 2)), BAall(ij(:, 1)), 'UniformOutput',false);
    % results.Sij = Sij;
    % Sij = secondOrder(A, B, ABi, ABj, BAi);
    results.S2 = Sij;
    results.S2conf = Sijconf;
end


end

function [A, B, AB, BA] = splitOutputsIntoMatrices(modelOutputs, N,k)
A = modelOutputs(1:N, :);
B = modelOutputs(N+1:2*N, :);
AB = modelOutputs(2*N+1:N*(k+2), :);
if length(modelOutputs) > N*(k+2)
    BA = modelOutputs(N*(k+2)+1:2*N*(k+1), :);
else
    BA = [];
end
end

function S1i = firstOrder(A, B, ABi)
V = var([A; B]);
S1i = mean(B .* (ABi - A)) ./ V;
end

function STi = totalOrder(A, B, ABi)
V = var([A; B]);
STi = mean(A .* (A - ABi)) ./ V;
end

%
function Sij = secondOrder(A, B, ABi, ABj, BAi)
V = var([A; B]);
Vij = mean(BAi .* ABj - A .* B) ./ V;
Si = firstOrder(A, B, ABi);
Sj = firstOrder(A, B, ABj);
Sij =  Vij - Si - Sj;
end
