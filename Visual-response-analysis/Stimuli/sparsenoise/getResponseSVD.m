% Function to return SVD
% Inputs:
%   responseMatrix: [N x t] where N is number of stimuli (or trials)
%
% Outputs:
%    stimulusWeights - N x 1 of responses to each stimulus
%    overallTimeCourse -1 x t time course of response
%
% History
% SGS 13th April 2020. Wrote it.
function [stimulusWeights,overallTimeCourse] = getResponseSVD(responseMatrix)
% Get mean LFP
meanLFP = mean(responseMatrix,1);
[maxF,maxI] = max(abs(meanLFP));
signData = sign(meanLFP(maxI));

% Get SVD 
[U, S, V] = svd(responseMatrix);
overallTimeCourse = V(:,1)';
stimulusWeights = U(:,1);

% Correct for any possible flipping
[maxF,maxI] = max(abs(overallTimeCourse(:)));
signU = sign(overallTimeCourse(maxI));


if signU ~= signData
    overallTimeCourse = -1*overallTimeCourse;
    stimulusWeights = -1*stimulusWeights;
end

