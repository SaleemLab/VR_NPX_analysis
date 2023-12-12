% Programme to estiamte receptive fields from sparse noise in NP1 / NP2
% data
% Dependencies: NPXAnalysis2022;
%               Visual-response-analysis/Stimuli/sparsenoise
% FR used the original Script to make it easier to call this script on
% another function
% 09/06/22

function [initMap,tMap, scal_f, sigma, framesToShow] = getSparseNoise_NPX(resps,otherData,stimData,nChannelsToBin,channelRange,sn_options)

%%
% Calculate sparsenoise RFs
stim_matrix = cat(3,stimData.stim_matrix{:}); % N rows x M cols x nFrames
sn_options.grid_size = [size(stim_matrix,1) size(stim_matrix,2)];
sn_options.mapSampleRate = 60; % Hz
sn_options.mapsToShow = {'linear','black','white','contrast'};
sn_options.mapMethod = 'fitlm'; % fitlm mean
sn_options.framesToShow = [1 4 6 8 10 12];  % at 60 Hz would be 8 ms, 40, 72 etc
sn_options.plotflag = 1;

% Get the data for the specific channels
newresps = resps(channelRange,:,:);

% Extract wheel and eye data
wheel_data = squeeze(nanmean(otherData(1,:,:),2))';
eye_data = squeeze(nanmean(otherData(2,:,:),2))';
% Choose a channel...or two
channelBins = 1:nChannelsToBin:length(channelRange);
itCount = 0;
for theseChannels = 1:length(channelBins)-1
    itCount = itCount+1;
    % Get average across channels (zscore?)
    switch(sn_options.importMode)
        case 'KS'
            tp = sn_options.peakChannel >= channelBins(theseChannels) & sn_options.peakChannel < channelBins(theseChannels+1);
        otherwise
            tp = channelBins(theseChannels):channelBins(theseChannels+1);
    end
    lfp_data = mean(newresps(tp,:,:),1); % returns 1 x N time bins x nFrames
    lfp_data = squeeze(lfp_data)';   % sparsenoise needs frames x N time bins

    [initMap(:,:,:,:,itCount), tMap, scal_f, sigma, framesToShow] = sparseNoiseAnalysis(stim_matrix,lfp_data,wheel_data,eye_data,sn_options); % check to make sure this is the correct one - should be in Visua;-response-analysis/Stimuli/sparsenoise

end
end