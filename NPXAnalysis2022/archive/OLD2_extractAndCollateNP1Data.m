% function [resps,otherData,stimData,eyeData,wheelData,photodiodeData] = extractAndCollateNPData(options)
%
% Programe to extract data from Neuropixel 1 probes and Bonsai/BonVision 
%   for subsequent analysis
%
% Input: options structure 
%   
%     Mandatory fields
%      importMode               % LF or MUA or KS
%      BinWidth                 % resolution (in s) of output resps (e.g. 1/60)
%      AnalysisTimeWindow       % two-element vector around stim-on (e.g. [-0.25 1.25])
%      EPHYS_DATAPATH           % full path to Imec Ephys recording
%      PERIPHERALS_DATAPATH     % full path to BehaviourData bonsai logfile
%      EYEDATA_DATAPATH         % full path to EyeTrackingData bonsai logfile
%      TRIALDATA_DATAPATH       % full path to TrialData bonsai logfile
%
%     Optional fields
%      verboseFlag              [Default 1]
%      MAP_FILE                 [Default neuropixPhase3A_kilosortChanMap.mat]
%      otherMetrics             [Default {'wheelspeed','pupilarea'}] that
%                                   are outputted in otherData matrix at
%                                   the same resolution and timebase as the ePhys data
%      KS_DATAPATH              % full path to Kilosort3 outputs for this session
%
% History: 
%   SGS: 24th Feb 2022: Wrote ExtractLFPAndMUAFromNP1 
%        27th Feb 2022: Updated to extractAndCollateNP1Data to be stand alone 
%                       and removed filepath declarations to caller 
%        4th April 2022:Allowed interface with KS files
%       10th April 2022:Allowed compatability with NP2 as well as NP1 (by including the
%               channel map
% TO DO: 
%        Allow user to select which time to sync to (e.g. photodiode up or
%           down, bonsai stimulus log etc)
%
% Dependencies
%   importAndAlignBonsaiLogs and associated dependencies
%   downsample
%   filtfilt
%   designfilt
%   https://github.com/djoshea/neuropixel-utils
%   https://github.com/cortex-lab/neuropixels
%
function [resps,otherData,stimData,eyeData,wheelData,photodiodeData,timeVector,options] = extractAndCollateNPData(options)

% Move mandatory info from options
BinWidth = options.BinWidth;
AnalysisTimeWindow = options.AnalysisTimeWindow;
importMode = options.importMode;
EPHYS_DATAPATH = options.EPHYS_DATAPATH;
PERIPHERALS_DATAPATH = options.PERIPHERALS_DATAPATH;
EYEDATA_DATAPATH = options.EYEDATA_DATAPATH;
TRIALDATA_DATAPATH = options.TRIALDATA_DATAPATH;
switch(options.importMode)
    case 'KS'
        KS_DATAPATH = options.KS_DATAPATH;
end

%%%%%
% Set defaults
if ~isfield(options,'verboseFlag') || isempty(options.verboseFlag)
    verboseFlag = 1;
else
    verboseFlag = options.verboseFlag;
end
if ~isfield(options,'otherMetrics') || isempty(options.otherMetrics)
    otherMetrics = {'wheelspeed','pupilarea'};
else
    otherMetrics = options.otherMetrics;
end
% Set up MAP_FILE 
if ~isfield(options,'MAP_FILE') || isempty(options.MAP_FILE)
    MAP_FILE = 'neuropixPhase3A_kilosortChanMap.mat'; %already added to path
else
    MAP_FILE = options.MAP_FILE;
end

%%%%%
% Check to see if MAP file on system
if isempty(which(MAP_FILE))
    error('MAP_FILE %s not on system; you may need to add ''neuropixels'' and/or ''neuropixels-utils'' to your file path',MAP_FILE)
else
    setenv('NEUROPIXEL_MAP_FILE',MAP_FILE)
end

%%%%%
% Step 1: Extract relevant data from wheel, pupil and stimulus files and
% sync different timebases
if verboseFlag
    fprintf('\n\tImporting stimulus, peripherals and eye data from \n\t\t%s and \n\t\t%s and \n\t\t%s ...',TRIALDATA_DATAPATH,PERIPHERALS_DATAPATH,EYEDATA_DATAPATH)
end
tstart = tic;
[stimData,eyeData,wheelData,photodiodeData,stimTimes] = importAndAlignBonsaiLogs(EPHYS_DATAPATH,TRIALDATA_DATAPATH,PERIPHERALS_DATAPATH,EYEDATA_DATAPATH);
basicDataLoadTime = toc(tstart);
if verboseFlag
    fprintf('\n\t\tTook %3.1fs to load Bonsai wheel and eye data',basicDataLoadTime)
end

%%%%%
% Step 2: Cycle through the trials and extract LFP or MUA around each time
if verboseFlag
    fprintf('\n\tImporting ephys data from %s ...',EPHYS_DATAPATH)
end
imec = Neuropixel.ImecDataset(EPHYS_DATAPATH);

tstart = tic;

%%%%%
% If this stimulus has a trial structure (ie. multiple rows in stimData),
% then process as trials, producing 3D matric of channels/units vs time vs trial), 
% otherwise return data as single matrix of channels/units vs time into recording
switch(importMode)
    case 'LF'
        % Design low pass filter (with corner frequency determined by the
        % desired output binwidth
        d1 = designfilt('lowpassiir','FilterOrder',12, ...
            'HalfPowerFrequency',(1/BinWidth)/2,'DesignMethod','butter','SampleRate', imec.fsLF);
        
        % Preallocate matrix for speed - so use the first stimulus to work
        % out how much space we will need for the resampled data
        first_stim = stimTimes(1);
        downSampleRate = fix(imec.fsLF*BinWidth); % the rate at which we downsample depends on the acquisition rate and target binwidth
        tresps = imec.readLF_timeWindow([first_stim+AnalysisTimeWindow(1) first_stim+AnalysisTimeWindow(2)]);
        % Downsample
        tresps = downsample(tresps',downSampleRate)';
        % Create vector for PSTH times at the downsampled rate (NB I am not
        % sure this is completely correct - it depends on how the imec
        % class does this [in the imec class, it does: 
        %   idxWindow = imec.closestSampleLFForTime(timeWindowSec);
        %   sampleIdx = idxWindow(1):idxWindow(2);
        timeVector = linspace(AnalysisTimeWindow(1),AnalysisTimeWindow(2),size(tresps,2));
        % Initialise matrix for ephys data
        resps = zeros([size(tresps,1) size(tresps,2) length(stimTimes)]);
        % For each trial
        H = waitbar(0,'Importing trials');
        for thisTrial = 1:length(stimTimes)
            waitbar(thisTrial/length(stimTimes),H)
            % Establish what the stimulus time would be in ephys
            stim_on = stimTimes(thisTrial);
            % Get the LFP from IMEC file
            tresps = imec.readLF_timeWindow([stim_on+AnalysisTimeWindow(1) stim_on+AnalysisTimeWindow(2)]);
            % Zero-mean the data ...
            tresps = tresps-repmat(mean(tresps,2),1,size(tresps,2));
            % ... and lowpass filter ...
            tresps = filtfilt(d1,double(tresps)')';
            % ... then downsample that data
            ttresps = downsample(tresps',downSampleRate)';
            % ... and save
            resps(:,:,thisTrial) = single(ttresps);
        end
        close(H)
        
        
    case 'MUA'
        % Design low pass filter (with corner frequency determined by the
        % desired output binwidth
        d1 = designfilt('lowpassiir','FilterOrder',12, ...
            'HalfPowerFrequency',(1/BinWidth)/2,'DesignMethod','butter','SampleRate', imec.fsAP);
        
        % Preallocate matrix for speed - so use the first stimulus to work
        % out how much space we will need for the resampled data
        first_stim = stimTimes(1);
        downSampleRate = fix(imec.fsAP*BinWidth); % the rate at which we downsample depends on the acquisition rate and target binwidth
        tresps = imec.readAP_timeWindow([first_stim+AnalysisTimeWindow(1) first_stim+AnalysisTimeWindow(2)]);
        % Downsample
        tresps = downsample(tresps',downSampleRate)';
        % Create vector for PSTH times at the downsampled rate (NB I am not
        % sure this is completely correct - it depends on how the imec
        % class does this [in the imec class, it does:
        %   idxWindow = imec.closestSampleLFForTime(timeWindowSec);
        %                   sampleIdx = round(timeSeconds * imec.fsLF)
        %   sampleIdx = idxWindow(1):idxWindow(2);
        timeVector = linspace(AnalysisTimeWindow(1),AnalysisTimeWindow(2),size(tresps,2));
        
        % Initialise
        resps = zeros([size(tresps,1) size(tresps,2) length(stimTimes)]);
        % For each trial
        H = waitbar(0,'Importing trials');
        for thisTrial = 1:length(stimTimes)
            waitbar(thisTrial/length(stimTimes),H)
            % Establish what the stimulus time would be in ephys
            stim_on = stimTimes(thisTrial);
            % Get the LFP from IMEC file
            tresps = imec.readAP_timeWindow([stim_on+AnalysisTimeWindow(1) stim_on+AnalysisTimeWindow(2)]);
            % Zero-mean the data ...
            tresps = tresps-repmat(mean(tresps,2),1,size(tresps,2));
            % ... fullwave rectify the signal ...
            tresps = abs(tresps);
            % ... and lowpass filter ...
            tresps = filtfilt(d1,double(tresps)')';
            % ... then downsample that data
            ttresps = downsample(tresps',downSampleRate)';
            % ... and save
            resps(:,:,thisTrial) = single(ttresps);
        end
        close(H)
        
    case 'KS'
        if ~isfield(options,'gFileNum') || isempty(options.gFileNum)
            error('Need to specify ''g'' file number to proceed with KS import')
        end
        if ~isfield(options,'KS_DATAPATH') || isempty(options.KS_DATAPATH)
            error('Need to specify KS_DATAPATH to proceed with KS import')
        end        
        if ~isfield(options,'ks_unitType') || isempty(options.ks_unitType)
            options.ks_unitType = '';
            fprintf('\n'); warning('Returning all units by default, both good and mua'); fprintf('\n');
        end     
        
        % Extract spike times for each unit
        SampleRate = imec.fsAP;
        [these_spike_times,nominal_KSLabel,cluster_id,peakChannel,maxSpkTime] = import_ks_spiketimes(options.KS_DATAPATH,options.gFileNum,[],SampleRate);

        % Keep only desired units
        switch(lower(options.ks_unitType))
            case 'good'
                % if user specifed only 'good'
                unitsToKeep = find(nominal_KSLabel == nominal('good'));
            case 'mua'
                % ... or 'mua'
                unitsToKeep = find(nominal_KSLabel == nominal('mua'));
            otherwise
                % ... or nothing
                unitsToKeep = 1:length(nominal_KSLabel);
        end
        these_spike_times = these_spike_times(unitsToKeep);
        options.nominal_KSLabel = nominal_KSLabel(unitsToKeep);
        options.peakChannel = peakChannel(unitsToKeep);
        options.cluster_id = cluster_id(unitsToKeep);
        
        % Initialise 
        if size(stimData,1) == 1
            stimTimes = 0;
            binEdges = stimTimes(1):BinWidth:maxSpkTime;
            timeVector = binEdges(1:end-1);
            resps = zeros([length(unitsToKeep) length(binEdges)-1]);
            % For each unit
            for thisUnit = 1:length(these_spike_times)
                tspk = these_spike_times{thisUnit};
                tspk = tspk(tspk >= 0 & tspk <= maxSpkTime);
                resps(thisUnit,:) = histcounts(tspk,binEdges)/BinWidth;
            end
        else
            binEdges = AnalysisTimeWindow(1):BinWidth:AnalysisTimeWindow(2);
            timeVector = binEdges(1:end-1);
            resps = zeros([length(unitsToKeep) length(binEdges)-1 length(stimTimes)]);
            % For each trial
            H = waitbar(0,'Importing trials');
            for thisTrial = 1:length(stimTimes)
                waitbar(thisTrial/length(stimTimes),H)
                % Establish what the stimulus time would be in ephys
                stim_on = stimTimes(thisTrial);
                % Get the spike times for each unit
                tresps = NaN*zeros([length(these_spike_times) length(binEdges)-1]);
                for thisUnit = 1:length(these_spike_times)
                    tspk = these_spike_times{thisUnit};
                    tspk = tspk(tspk >= stim_on+AnalysisTimeWindow(1) & tspk <= stim_on+AnalysisTimeWindow(2))-stim_on;
                    tresps(thisUnit,:) = histcounts(tspk,binEdges)/BinWidth;
                end
                % ... and save
                resps(:,:,thisTrial) = tresps;
            end
            close(H)
        end
end
ePhysDataLoadTime = toc(tstart);
if verboseFlag
    fprintf('\n\t\tTook %3.1fs to load ePhys data',ePhysDataLoadTime)
end

% Step 3: extract other data
% Obtain trial based dat for othermetrics
% Initialise matrix for other data (wheel, pupil etc)
tstart = tic;
otherData = [];
if ~isempty(otherMetrics)
    for thisMetric = 1:length(otherMetrics)
        % You can add in more cases here
        switch(lower(otherMetrics{thisMetric}))
            case 'wheelspeed'
                inData.sglxTime = wheelData.sglxTime;
                inData.Data = wheelData.smthSpeed;
            case 'pupilarea'
                inData.sglxTime = eyeData.sglxTime;
                inData.Data = eyeData.Area;
        end
        
        % Accumulate these data
        t_otherdata = getSampledOtherMetrics(inData,stimTimes,timeVector);
        if thisMetric == 1
            otherData = t_otherdata;
        else
            otherData(thisMetric,:,:) = t_otherdata;
        end
    end
end
otherDataLoadTime = toc(tstart);
if verboseFlag
    fprintf('\n\t\tTook %3.1fs to load other Bonsai data',otherDataLoadTime)
end
end

function otherdata = getSampledOtherMetrics(inData,stimTimes,analysisTimeWindowVector)
% For eye, wheel data etc
% Input inData must have fields sglxTime and Data
otherdata = zeros([1 length(analysisTimeWindowVector) length(stimTimes)]);
% For each trial
for thisTrial = 1:length(stimTimes)
    [~,ind1] = min(abs(inData.sglxTime - (stimTimes(thisTrial)+analysisTimeWindowVector(1)))); % find index where timestamp is closest to that desired
    [~,ind2] = min(abs(inData.sglxTime - (stimTimes(thisTrial)+analysisTimeWindowVector(end)))); % find index where timestamp is closest to that desired
    centredTime = inData.sglxTime(ind1:ind2)-stimTimes(thisTrial);
    td = inData.Data(ind1:ind2);
    tresps = interp1(centredTime(:)', td(:)', analysisTimeWindowVector(:)');
    otherdata(1,:,thisTrial) = tresps;
end
end