% Function to contruct RF map from minimum inputs
%   Intended to support the generation of a forward receptive field map (ie. as a function of time after onset of stimulus) 
% Inputs:
%   stimulusMatrix: m x n x t matrix where m is number of rows in spatial
%                   array, n is number of columns, and t in number of timebins
%
%   responseVector: [1 x t] or [t x 1] vector of response (whatever metric
%                    that is)
%
%   delaysInBins: [default 0] offset applied to stimulus before calculating map, in bins.
%                   -3 indicates map if response is correlated with a stimulus 3 time bins
%                   before response
%   method:[default mean] mean or fitlm
%
% Outputs:
%   thisAvgMap - n x m x N receptive field map (average, fitlm)
%   thisVarMap - n x m x N receptive field map (variance, p value)
%
% History
% SGS 13th April 2020. Wrote it.
%
function [thisAvgMap,thisVarMap] = makeRFmap(stimulusMatrix,responseVector,delaysInBins,method)
if ~exist('method','var')
    method = 'mean';   % If no method specified, use the mean
end
if ~exist('delaysInBins','var')
    delaysInBins = 0;   % If no delays asked for, calculate only for aligned
end
if ~exist('responseVector','var')
    error('Both stimulus matrix and response vector needed')
end

% Make sure response vector is a column vector
responseVector  = responseVector(:);

% For each requested delay, estimate the receptive field from the combination of
% response vector and stimulus matrix
nDelaysRequested = length(delaysInBins);

% Reshape stimulus matrix into row vector (N times  x (m x n) positions matrix)
rowVec_sm = reshape(permute(stimulusMatrix,[3 1 2]),size(stimulusMatrix,3),size(stimulusMatrix,1)*size(stimulusMatrix,2));

% For N delays....
allResponses = NaN(length(delaysInBins),size(rowVec_sm,2)); % N delays x (m x n) positions matrix
allVariances = NaN(length(delaysInBins),size(rowVec_sm,2)); % N delays x (m x n) positions matrix
for thisDelay = 1:length(delaysInBins)
    % Shift the stimulus matrixby appropriate n
    delayBin = delaysInBins(thisDelay);
    % Get weighted mean of stimulus matrix (ie amplitude multiplied by stimulus)
    if delaysInBins < 0 % negative means using the stimuli at X bins before the response
        thisResponseVector = responseVector(1:end-(-delayBin));
        thisStimulusMatrix = rowVec_sm(-delayBin:end,:);
    elseif delaysInBins > 0 % positive means using the stimuli at X bins after the response (non causal!)
        thisResponseVector = responseVector(delayBin:end);
        thisStimulusMatrix = rowVec_sm(1:end-delayBin,:);
    else % if no offset
        thisResponseVector = responseVector;
        thisStimulusMatrix = rowVec_sm;
    end
    
    % Depending on request estimate in different ways
    switch(method)
        case 'fitlm'
            td = fitlm(thisStimulusMatrix,thisResponseVector);
            theseResponseWeights = td.Coefficients.Estimate(2:end);
            theseVariances = td.Coefficients.pValue(2:end);
            allResponses(thisDelay,:) = theseResponseWeights;
            allVariances(thisDelay,:) = theseVariances;
        case 'mean'
            theseResponseWeights = thisResponseVector'*thisStimulusMatrix/size(thisStimulusMatrix,2);
            allResponses(thisDelay,:) = theseResponseWeights;
    end
    
    
end

% Move back into 3D state
thisAvgMap = reshape(permute(allResponses,[2 1]),[size(stimulusMatrix,1),size(stimulusMatrix,2),length(delaysInBins)]);
thisVarMap = reshape(permute(allVariances,[2 1]),[size(stimulusMatrix,1),size(stimulusMatrix,2),length(delaysInBins)]);