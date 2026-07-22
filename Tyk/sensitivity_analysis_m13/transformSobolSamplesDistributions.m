function samples = transformSobolSamplesDistributions(paramObj, samples, sampleDist)
%Author: Jaimit Parikh
%Comments by Skylar Grey
%This file specifies the exact distribution for each parameter

%get the parameter distribution types from parameter object
paramDistTypes = paramObj.dist;
%get the lower bounds for truncation function
lb = paramObj.lb;
%get the upper bounds for truncation function
ub = paramObj.ub;
%get the parameter names
distParamNamesValues = paramObj.parameters;

%cellfun applies the function to each individual cell of an array
%makedist makes the distributions of the specified types with the 
%specified properties, while truncate cuts those distributions off 
%at the specified upper and lower bounds (not needed if only using 
%uniform distributions, however this code works for other distributions)
paramDists = cellfun(@(x, y, l, u)truncate(makedist(x, y{:}), l, u), ...
    paramDistTypes, ...
    distParamNamesValues, lb, ub, ...
    'UniformOutput',false);
%splitapply splits each object into groups (in this case along columns)
%and applies the function to each group one at a time
%In this case, it applies the user-written function transformedSamples
%which transforms the Sobol matrices into the parameter space
samples = splitapply(...
    @(a,b,c)inverseTransformSampling(a, b{:}, c{:}), ...
    samples, sampleDist, paramDists,...
    1:size(samples, 2));
end

%inverseTransformSampling passes elements of samples, and 
% unpacks sampleDist and paramDists. For each column of samples, 
% inverseTransform sampling transforms the distribution
% of the column from a uniform distribution to a desired 
% distribution, specified by paramDists

%this function is called by splitapply above, and is described below
function transformedSamples = inverseTransformSampling(samples, ...
    sampleDist, desiredDist)
arguments
    samples;
    sampleDist;
    desiredDist;
end
%% transformedSamples = inverseTransformSampling(samples, sampleDist,
% desiredDist) returns the transformedSamples with the desired
% distribution given orginal samples sampled from sampleDist
%samples contains the Sobol matrices and their products
%this function takes a uniform distribution with min 0 and max 1
%and applies its cdf to the Sobol matrices, then it creates an inverse
%cumulative distribution function with the parameter distributions. 
%This transforms the samples into the parameter space.
transformedSamples = desiredDist.icdf(sampleDist.cdf(samples));

end
