function [LFP,config,info] = getNPLFP(binName,binPath,windowS)
%GETNPLFP gets LFP from bin file and sorts in depth order

meta = NPadmin.ReadMeta(binName, binPath); %Read meta
% Read in file
fid = fopen(fullfile(binPath, binName), 'rb');
% current just have a short window but if want to load in the
% whole thing then use: meta.fileTimeSecs

if isempty(windowS) || windowS > meta.fileTimeSecs
    windowS = meta.fileTimeSecs;
end
inData = fread(fid, [meta.nSavedChans,round(windowS*meta.imSampRate)], 'int16=>double'); % read LFP data, should be in imro channel order?
fclose(fid);
% Check probe type and resample if needed
probeType = str2double(meta.imDatPrb_type);
if probeType == 24 && contains(binName,'.ap.bin')
    % Need to resample down to get to 1kHz LFP
    disp('Need to resample, might take a while...')
    rsmpfactor = round(meta.imSampRate/1000);
    fs = meta.imSampRate/rsmpfactor;
    for currChan = 1:size(inData,1)
        LFP(currChan,:) = resample(inData(currChan,:), 1, rsmpfactor);
    end
    disp('Finished resampling')
else
    LFP = inData;
    fs = meta.imSampRate;
end
% syncArray = LFP(end,:);
% syncPulse = logical(syncArray);
% NSyncSamples = strfind(syncPulse,[0 1])+1;
% SyncTimes = (NSyncSamples/fs);

LFP(end,:) = []; % last channel is sync square wave, remove.
[config,type] = NPadmin.getNPChannelConfig(meta);

% No gain correction to be performed on type 24
if probeType ~= 24
    LFP = NPadmin.GainCorrectIM(LFP,config.SpikeGLXchan0,meta);
end

LFP(config.Internal_reference,:) = []; % remove internal reference
config(config.Internal_reference,:) = []; %  now indexes data


% Sort in order
[config, electrode_order] = sortrows(config,{'Shank','Electrode'},'ascend');
LFP = transpose(LFP(electrode_order,:));  % sort in depth order
info.probeType = probeType;
info.fs = fs;
end

