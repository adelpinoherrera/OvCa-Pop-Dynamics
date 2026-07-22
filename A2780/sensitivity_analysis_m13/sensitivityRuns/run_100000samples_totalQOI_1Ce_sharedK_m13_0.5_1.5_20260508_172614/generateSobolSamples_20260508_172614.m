%This file creates the sample sets needed to compute Sobol indices
function [samples, sampleDist] = generateSobolSamples(paramObj, ...
    calcSecondOrder)
arguments
    paramObj;%note: this is called parsObj in SobolMain
    calcSecondOrder logical; % 
end
% getSamplesSobol provides the sobol samples using structured sampling
% Inputs:
%   paramObj; structure with parameter information
% Outputs:
%   samples: the number of returned samples = N(k+2) if only calculating
%   first order index or 2N(k+1) if also calculating second order index
% Author: Jaimit Parikh
% Last modified 07-10-2023 by Jaimit Parikh
% see also sobolset sgenerator makedist

% comments by Skylar Grey and Kyle Adams

N = paramObj.N; % Number of base samples - this will get multiplied by (K+2) for first-order S1
k = length(paramObj.name);% number of parameters
sgenerator = sobolset(k * 2); % double the number of samples to create 2 matrices (A & B) below
%the line above creates a higher-dimensional object of dimension k*2
%"p = sobolset(d) constructs a d-dimensional  point set p, which is a
%sobolset object with default property settings"
samples = sgenerator(1:N, :); %gets N rows of k*2 columns


samples = mat2cell(samples, N, [k, k]);%separates the matrix into a row of N kxk cells
%separates matrix into 2 unique cells, which look like a smaller matrix
%that is Nxk
samples = cell2mat(samples');%stacks cells into a N*2 x k matrix, but in reverse
%looks like cell(2) on top of cell(1)
%calls the function below beginning on line 44 to take the samples and generate the matrices and their products
[A, B, AB, BA] = createMatricesForSobolIndices(samples);

if ~calcSecondOrder %(if not calcSecondOrder)
    samples = [A; B; AB];%puts matrices into desired object for second order
else
    samples = [A;B;AB;BA];%puts matrices into desired object for first order
end

sampleDist = cell(1, k);%makes an object of type cell with dimension 1xk
[sampleDist{:}] = deal(makedist('Uniform'));%makes a standard uniform distribution [0,1] and fills each cell with that
%Calls a function to take the parameter object along with the sobol matrices and the
%standard uniform distributions and return sobol matrices that are in the parameter space
samples = transformSobolSamplesDistributions(paramObj, samples, sampleDist);
end
%This function takes the NxKxK object and transforms it into sobol matrices
function [A, B, AB, BA] = createMatricesForSobolIndices(samples)
arguments
    samples = reshape(1:8, 2, [])'; % default matrix for test
end


N = length(samples) / 2; %sets N to be what it originally was
nInputs = size(samples, 2); %gets number of parameters

% randomly permute the rows of the samples
% samples = samples(randperm(2*N), :);

% create matrix A, B and AB by splitting the samples
A = samples(1:N, :);
B = samples(N + 1:end, :);

Arepeat = repmat(A, nInputs, 1); %matrix A is duplicated nInputs times

Bcell = num2cell(B, 1); %splits cols of B into cells
Bdiag = blkdiag(Bcell{:}); %makes a block diagonal matrix, each block is a column from B
AB = Bdiag + Arepeat .* not(Bdiag);
%not(Bdiag) seems to return 1 in element if Bdiag is 0, return 0 if element
%in Bdiag is NOT 0

Brepeat = repmat(B, nInputs, 1); %matrix  B is duplicated nInputs times
Acell = num2cell(A, 1); %splits cols of A into cells
Adiag = blkdiag(Acell{:}); %makes a block diagonal matrix, each block is a column from A
BA = Adiag + Brepeat .*not(Adiag); 


end
