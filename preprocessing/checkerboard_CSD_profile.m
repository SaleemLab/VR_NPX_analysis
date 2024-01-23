function [lfpAvg,csd,power,best_channels] = checkerboard_CSD_profile(options)

load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'))
load(fullfile(erase(options.ANALYSIS_DATAPATH,['\','Checkerboard']),"best_channels.mat"))
load(fullfile(erase(options.ANALYSIS_DATAPATH,['\','Checkerboard']),"extracted_PSD.mat"),'power') % Load only power variable

BinWidth = 1/1250;
options.importMode = 'LF'; % LF or MUA or KS
% options.importMode = 'LF'; % LF or MUA or KS
AnalysisTimeWindow = [-0.1 0.5];% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])
%             options.ks_unitType = 'good'; % 'mua', 'good' or ''
nprobe = options.probe_id+1; % probe_no is [1,2] based on options.probe_id (0 and 1)

column = 1;
[FILE_TO_USE,~,chan_config,sorted_config] = extract_NPX_channel_config(options,column);


% Extract LFP data relative to the stimulu onset time
% Preallocate matrix for speed - so use the first stimulus to work
% out how much space we will need for the resampled data
stimTimes = Task_info.stim_onset;
first_stim = stimTimes(1);

%%%% OLD - use Neuropixels library
%         downSampleRate = fix(imec.fsLF*BinWidth); % the rate at which we downsample depends on the acquisition rate and target binwidth
%         tresps = imec.readLF_timeWindow([first_stim+AnalysisTimeWindow(1) first_stim+AnalysisTimeWindow(2)]);
%         % Design low pass filter (with corner frequency determined by the
%         % desired output binwidth
%         d1 = designfilt('lowpassiir','FilterOrder',12, ...
%             'HalfPowerFrequency',(1/BinWidth)/2,'DesignMethod','butter','SampleRate', imec.fsLF);
%%%% NEW - use SpikeGLX_Datafile_Tools

imecMeta = ReadMeta(fullfile(options.EPHYS_DATAPATH,FILE_TO_USE));
nSamp = fix(SampRate(imecMeta)*range(AnalysisTimeWindow));
chanTypes = str2double(strsplit(imecMeta.acqApLfSy, ','));
nEPhysChan = chanTypes(1);
downSampleRate = fix(SampRate(imecMeta)*BinWidth); % the rate at which we downsample depends on the acquisition rate and target binwidth
% Design low pass filter (with corner frequency determined by the
% desired output binwidth)
d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',(1/BinWidth)/2,'DesignMethod','butter','SampleRate',SampRate(imecMeta));

% Read the data
tresps = ReadBin(fix(SampRate(imecMeta)*first_stim), nSamp, imecMeta, FILE_TO_USE, options.EPHYS_DATAPATH);
tresps(nEPhysChan+1:end,:) = []; % Get rid of sync channel
if strcmp(imecMeta.typeThis, 'imec')
    tresps = GainCorrectIM(tresps, 1:nEPhysChan, imecMeta);
else
    tresps = GainCorrectNI(tresps, 1:nEPhysChan, imecMeta);
end

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
    %%%% OLD - use Neuropixels library
    %             tresps = imec.readLF_timeWindow([stim_on+AnalysisTimeWindow(1) stim_on+AnalysisTimeWindow(2)]);
    %%%% NEW - use SpikeGLX_Datafile_Tools
    tresps = ReadBin(fix(SampRate(imecMeta)*(stim_on+AnalysisTimeWindow(1))), nSamp, imecMeta, FILE_TO_USE, options.EPHYS_DATAPATH);
    tresps(nEPhysChan+1:end,:) = []; % Get rid of sync channel
    if strcmp(imecMeta.typeThis, 'imec')
        tresps = GainCorrectIM(tresps, 1:nEPhysChan, imecMeta);
    else
        tresps = GainCorrectNI(tresps, 1:nEPhysChan, imecMeta);
    end

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


timestamps = linspace(AnalysisTimeWindow(1),AnalysisTimeWindow(2),size(resps,2));

lfpAvg = [];
csd = [];
[csd.all, lfpAvg.all ] = perievent_CSD_LFP_amplitude_phase(permute(resps(sorted_config.Channel,:,:), [2, 1, 3]), timestamps,[],'twin',[0.1 0.5]);

lfpAvg.filter_type = {'all'};
lfpAvg.event_group = {'Checkerboard'};
lfpAvg.probe_no = nprobe;
% cd(fullfile(options.ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
%             plot_perievent_CSD_LFP_amplitude_phase(lfpAvg,csd,power{nprobe},chan_config,sorted_config,best_channels{nprobe})

plot_perievent_CSD_LFP(lfpAvg,csd,power{nprobe},chan_config,sorted_config,best_channels{nprobe},options)
