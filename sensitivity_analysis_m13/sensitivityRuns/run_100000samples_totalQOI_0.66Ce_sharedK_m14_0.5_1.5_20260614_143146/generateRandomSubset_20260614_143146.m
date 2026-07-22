%This function selects a random subset of specified size from any 
% bigger set. It is used to compare QOI distributions as histograms, 
% as the Sobol main file ends with a sample of (2N+1)*base samples 
% QOI values, where N is # of parameters. So it was necessary to 
% restrict the QOI distribution to the same number of samples as used
% for varying different parameters (most- and least-influential).
%Code written by Kyle Adams

function random_subset = generateRandomSubset(data, subset_size, seed)
    rng(seed); %set random seed
    random_indices = randperm(length(data), subset_size); %get random indices
    random_subset = data(random_indices); %return the random subset
end
