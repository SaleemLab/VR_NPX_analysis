% function [thisOEData, ttl, wheelData, fileinfo, timestamps] = readSpikesAndPDandWheel(thisOEfileName,options)
% Function to load the open ephys data, photodiode pulses and wheel data.
% AP 10/07/19

function [thisOEData, ttl, wheelData, fileinfo, timestamps] = readSpikesAndPDandWheel(thisOEfileName,options)

if ~isfield(options,'chanNumbersToInclude') || isempty(options.chanNumbersToInclude)
    channels = -1;    % Channels number associated with the stimulus, to include in the export (-1 = all) or [1 2 3 4]
else
    channels = options.chanNumbersToInclude;
end
if length(channels) == 1 && channels == -1
    channels = 1:32;
end
if ~isfield(options,'fileID')
    options.fileID = 1; % Default log to command prompt
end

%%%%%%%%%%%%%%%%---------------------------------%%%%%%%%%%%%%%%
% Fix sample rates
mainSampleRate = 30000;  % Main sample rate for OE files is 30 kHz

%%%%%%%%%%%%%%%%---------------------------------%%%%%%%%%%%%%%%
% PD options (built in)
if ~isfield(options,'baselineTimeForPDInS') || isempty(options.baselineTimeForPDInS)
    options.baselineTimeForPDInS = 5; % By default take the first 5s of recording to estimate PD baseline
end

%%%%
% Get spike data


% First check to see if already exists in processed data

clear SpikeData
if kilosort == 1
    SpikeData = get_ks_Spikes_new(options.spikesToppath);
else
    SpikeData = getSpikes(options.spikesToppath);
end
OE.MetaData = SpikeData.MetaData;

r=0; i=1;
while r==0
    if strcmp(OE.MetaData.FoldersList{i}(end-length(OE.ACInfo.ExpDate)+1:end),OE.ACInfo.ExpDate)
        r=i;
    end
    i=i+1;
end
OE.ROI = r;

OE.SpikeInfo = SpikeData.SpikeInfo;

if r==1
    stim_lims=[0 OE.MetaData.lims(r)-1]/OE.ACInfo.SamplingRateOE;
else
    stim_lims=[sum(OE.MetaData.lims(1:r-1)) sum(OE.MetaData.lims(1:r))]/OE.ACInfo.SamplingRateOE;
end

k=1;
for icell=2:size(OE.SpikeInfo,2)
    st=OE.SpikeInfo{1,icell};
    OE.StimSpiketimes{k}=st(st>stim_lims(1) & st<=stim_lims(2))-stim_lims(1);
    k=k+1;
end

%%
% 5) Count spikes per frame
numframes=length(OE.SN_onsets);
delays=-10:10:120; % -10 to 90ms delays

for icell = 1:length(OE.StimSpiketimes)
    st = OE.StimSpiketimes{icell};
    
    for delayInd=1:length(delays)
        for iframe = 1:numframes
            startTime=OE.SN_onsets(iframe);
            if iframe==numframes
                endTime=OE.SN_onsets(iframe)+OE.SN_singleI_dur;
            else
                endTime=OE.SN_onsets(iframe+1);
            end
            OE.frameSpikeCount{icell}(iframe, delayInd) = sum(st>(startTime + delays(delayInd)) & st<=(endTime + delays(delayInd))); %spike count per frame, per delay
        end
    end
end















switch (checkNew)
    case 0 % Old style
        oe_ch = cellfun(@any, regexp({oe_dir.name},'CH\d.*?.continuous'));
        oe_ch_filenames = natsortfiles(cellfun(@(x) [oe_path filesep x], ...
            {oe_dir(oe_ch).name},'uni',false));
        
        % dmem = memory;
        % memToLeaveFree = 4 * 2^30; % num of GB to keep free
        % memToAllocate = dmem.MemAvailableAllArrays - memToLeaveFree;
        % memToAllocate = max(0, memToAllocate);
        % nint16s = memToAllocate/2;
        
        %% Load the data
        chs = oe_ch_filenames(channels);
        nCH = length(chs);
        % chunkSizeSamps = nint16s/nCH;
        % nChunks = ceil(totalSamps/chunkSizeSamps);
        
        % Load in first channel to set memory restrictions
        fprintf(options.fileID, 'Reading channel 1 for initialization\n');
        [data, timestamps, data_info] = load_open_ephys_data_faster(chs{1});
        %[data, timestamps, data_info] = load_open_ephys_data(chs{1});
        totalSamps = length(data);
        
        % OEDataAllocate = [];
        OEDataAllocate = zeros(nCH, totalSamps, 'int16');
        
        for n = 1:nCH
            fprintf(options.fileID, 'reading CH %d/%d\n', n, nCH);
            if n>1 % can skip this for the first channel, it was loaded above
                [data, timestamps, info] = load_open_ephys_data_faster(chs{n});
            end
            OEDataAllocate(n,:) = int16(data);
        end
        
    case 1 % New style?
        OEpath = fullfile(oe_path,'experiment1','recording1','structure.oebin');
        OEdat = load_open_ephys_binary(OEpath, 'continuous',1);
        % Reassign to variables for consistency
        timestamps = OEdat.Timestamps;
        % Work out which channels are present
        tch = {}; [tch{1:length(OEdat.Header.channels)}] = deal(OEdat.Header.channels(:).description);
        % Find std channels
        stdCh = find(contains(tch,'headstage','ignorecase',true));
        OEDataAllocate = OEdat.Data(stdCh(channels),:);
        % Find ADC channels
        adcCh = find(contains(tch,'adc','ignorecase',true));
        dataForTTL = OEdat.Data(adcCh(options.sync_channel),:);
        dataForTTL = dataForTTL-min(dataForTTL);% rescale data to between 0 and 1
        dataForTTL = dataForTTL./max(dataForTTL);
        dataForTTL = dataForTTL'; % seems to think it was a column vector before

        if options.includeWheel
            dataForWheel = int16(OEdat.Data(adcCh(options.wheel_channels),:));
        end
        data_info.header = OEdat.Header;
        data_info.header.sampleRate = data_info.header.sample_rate;
        % Clear full variable
        clear OEdat;
end
% Invert if necessary
if options.multiplierDat ~= 1
    meanOEvals = int16(mean(OEDataAllocate,2));
    
    tmean = repmat(meanOEvals,1,size(OEDataAllocate,2));
    OEDataAllocate = OEDataAllocate-tmean;
    OEDataAllocate = OEDataAllocate.*options.multiplierDat;
    OEDataAllocate = OEDataAllocate+tmean;
end

% Process data if requested (subtract average of other electrodes)
if options.referenceDat
    thisDat = zeros(length(channels),size(OEDataAllocate,2),'int16');
    for thisChannel = 1:length(channels)
        
        % And not include 'dead' or otherwise corrupted channels
        [~, ignElectrodes] = setdiff(channels, options.ignoreElectrodes);
        
        % Comparison electrodes should only be those outside this
        % tetrode...
        refElectrodes = setdiff(1:length(channels),thisChannel);
        refElectrodes = intersect(refElectrodes, ignElectrodes);
        % Debug report
        fprintf(options.fileID,'\nCh %01d: referenced to Chs',thisChannel);
        fprintf(options.fileID,' %01d',refElectrodes);
        
        % Create average of z-scores across included electrodes
        zscored= zscore(double(OEDataAllocate(refElectrodes,:))')';
        meanZscore = mean(zscored,1);
        mChDat = mean(OEDataAllocate(thisChannel,:));
        thisChDat = double(OEDataAllocate(thisChannel,:))-mChDat;
        % Regress against target channel
        weight = meanZscore'\thisChDat';
        % Subtract weighted version of mean
        thisDat(thisChannel,:) = int16((thisChDat - weight*meanZscore)+mChDat);
    end
    OEDataAllocate = thisDat;
end

thisOEData = OEDataAllocate;
if options.logErrorsFlag
    fprintf(options.errorlogfile, 'done\n');
end

%% Load photodiote or TTL data
switch options.sync_input
    % Not sure if this code works but I keep it here in case it is ever needed
    % in the future
    case 'ttl'
        
        % Get events filename
        ttl_file = cellfun(@any, regexp({oe_dir.name},'all_channels.*?.events'));
        ttl_filename = [oe_path filesep oe_dir(ttl_file).name];
        [data, timestamps, info] = load_open_ephys_data_faster(ttl_filename);
        
        % Structure of file: data = channel
        ttl_timestamps = timestamps(data == (sync_channel - 1));
        
        % First timestamp is arbitrary with regards to when open ephys was
        % started, so subtract the first timestamp (should be the start of
        % recording timestamp, and be equal to the continuous file
        % timestamps(1)/sample rate)
        ttl_timestamps = ttl_timestamps - timestamps(1);
        
    case 'adc'
        switch (checkNew)
            case 0 % old style
                % Get ADC sync channel filename
                oe_adc = cellfun(@any, strfind({oe_dir.name},['ADC' num2str(options.sync_channel)]));
                oe_adc_filename = [oe_path filesep oe_dir(oe_adc).name];
                nADC = 1;
                ttlData = zeros(nADC, totalSamps);
                fprintf(options.fileID, 'reading PD ADC %d/%d\n', options.sync_channel, nADC);
                [data, timestamps, sync_info] = load_open_ephys_data_faster(oe_adc_filename);
                ttlData = data;
                
            case 1 % new style
                ttlData = dataForTTL; % Kept over from above
                sync_info = data_info;
        end
        
    case 'ch'
        switch (checkNew)
            case 0 % old style
                % Get ADC sync channel filename
                oe_adc_filename = oe_ch_filenames(options.sync_channel);
                nADC = 1;
                fprintf(options.fileID, 'reading PD ADC %d/%d\n', options.sync_channel, nADC);
                [ttlData, timestamps, sync_info] = load_open_ephys_data_faster(oe_adc_filename{1});
                
            case 1 % new style
                ttlData = dataForTTL; % Kept over from above
                sync_info = data_info;
        end
end

%% Process sync/PD data
switch options.sync_input
    case 'ttl'
        ; % Unnecessary
    otherwise
        switch(options.photo_br)
            case 1
                % Find the baseline within the first bit of recording
                base = mean(ttlData(1:mainSampleRate*options.baselineTimeForPDInS));
                base_std = std(ttlData(1:mainSampleRate*options.baselineTimeForPDInS));
                if 3*base_std>0.02
                    thresh_pos = base+3*base_std;
                    thresh_neg = base-3*base_std;
                else
                    thresh_pos = base+0.02;
                    thresh_neg = base-0.02;
                end
            otherwise
                thresh_pos = 1.5;
                thresh_neg = 0.5;
        end
        if isfield(options,'OE_pD_thresholds') && ~isempty(options.OE_pD_thresholds)
            thresh_neg = options.OE_pD_thresholds(1);
            thresh_pos = options.OE_pD_thresholds(2);
        end
        ttl = convertPDtoTimeStamps(ttlData,thresh_pos,thresh_neg,mainSampleRate,options);
        ttl.ttlData = ttlData;
end

%% Load wheel data
if options.includeWheel
    switch (checkNew)
        case 0 % old style
            switch options.wheel_input
                case 'adc'
                    nADC = length(options.wheel_channels);
                    wheelData = zeros(nADC, totalSamps, 'int16');
                    for i=1:nADC
                        % Get ADC wheel channel filename
                        oe_adc = cellfun(@any, strfind({oe_dir.name},['ADC' num2str(options.wheel_channels(i))]));
                        oe_adc_filename = [oe_path filesep oe_dir(oe_adc).name];
                        fprintf(options.fileID, 'reading Wheel ADC %d/%d\n', i, nADC);
                        [data, timestamps, wheel_info] = load_open_ephys_data_faster(oe_adc_filename);
                        %[data, timestamps, wheel_info] = load_open_ephys_data(oe_adc_filename);
                        wheelData(i,:) = int16(data);
                    end
                case 'ch'
                    nADC = length(options.wheel_channels);
                    wheelData = zeros(nADC, totalSamps, 'int16');
                    for i=1:nADC
                        oe_adc_filename = oe_ch_filenames(options.wheel_channels(i));
                        fprintf(options.fileID, 'reading Wheel Channel %d\n', options.wheel_channels(i));
                        [data, timestamps, wheel_info] = load_open_ephys_data_faster(oe_adc_filename{1});
                        wheelData(i,:) = int16(data);
                    end
            end
        case 1
            wheelData = dataForWheel;
            wheel_info = data_info;
    end
    
else
    wheelData = [];
end

% Save some info for later
fileinfo.data_info = data_info;
fileinfo.sync_info = sync_info;
if options.includeWheel
    fileinfo.wheel_info = wheel_info;
end

