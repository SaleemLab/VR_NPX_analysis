function [lfpAvg,csd,PSD,best_channels] = checkerboard_CSD_profile(options)

BinWidth = 1/1250;
% options.importMode = 'LF'; % LF or MUA or KS
AnalysisTimeWindow = [-0.1 0.5];% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])
%             options.ks_unitType = 'good'; % 'mua', 'good' or ''
nprobe = options.probe_id+1; % probe_no is [1,2] based on options.probe_id (0 and 1)

power = [];
load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'))
DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'..','best_channels.mat'));
if ~isempty(DIR)
    load(fullfile(options.ANALYSIS_DATAPATH,'..','best_channels.mat'))
else
    best_channels{nprobe} = [];
end

clusters = [];
DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));
if ~isempty(DIR)
    if sum(contains({DIR.name},'ks3')) > 0
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'))
        clusters = clusters_ks3;
    else
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters.mat'))
    end
else
    clusters{nprobe} = [];
end

load(fullfile(options.ANALYSIS_DATAPATH,'..','extracted_PSD.mat'))



column = 1;
options.importMode = 'KS';
[file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column);% Since it is LF

if ~contains(imecMeta.acqApLfSy,'384,0') % NPX2 only has AP but NPX1 has AP and LF
    options.importMode = 'LF';
    [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column);% Since it is LF
    probe_type = 1;
else
    probe_type = 2;
end


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

imecMeta = ReadMeta(fullfile(options.EPHYS_DATAPATH,file_to_use));
nSamp = fix(SampRate(imecMeta)*range(AnalysisTimeWindow));
chanTypes = str2double(strsplit(imecMeta.acqApLfSy, ','));
nEPhysChan = chanTypes(1);
downSampleRate = fix(SampRate(imecMeta)*BinWidth); % the rate at which we downsample depends on the acquisition rate and target binwidth
% Design low pass filter (with corner frequency determined by the
% desired output binwidth)
d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',(1/BinWidth)/2,'DesignMethod','butter','SampleRate',SampRate(imecMeta));

% Read the data
tresps = ReadBin(fix(SampRate(imecMeta)*first_stim), nSamp, imecMeta, file_to_use, options.EPHYS_DATAPATH);
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
    tresps = ReadBin(fix(SampRate(imecMeta)*(stim_on+AnalysisTimeWindow(1))), nSamp, imecMeta, file_to_use, options.EPHYS_DATAPATH);
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

for nchannel = 1:size(chan_config,1)
%     power(nchannel,:) = PSD(nchannel).mean_power;
    xcoord(nchannel) = PSD{nprobe}(nchannel).xcoord;
    ycoord(nchannel) = PSD{nprobe}(nchannel).ycoord;
end
xcoord_avaliable = unique(xcoord);

% sort channel according to y coordinate
[ycoord idx] = sort(ycoord,'ascend');
xcoord = xcoord(idx);
power{nprobe} = power{nprobe}(idx,:);
chan_config = chan_config(idx,:);
resps = resps(idx,:,:);

if isempty(best_channels{nprobe})
    [~,ycoord_max] = max(ycoord);
    best_channels{nprobe}.surface_channel = ycoord_max;
    best_channels{nprobe}.surface_depth = max(ycoord);
end

lfpAvg = [];
csd = [];

for col = 1:length(xcoord_avaliable)
    lfpAvg(col).filter_type = {'all'};
    lfpAvg(col).event_group = {'Checkerboard'};
    lfpAvg(col).probe_no = nprobe;
%     lfpAvg(col).column = col;
    lfpAvg(col).xcoord = xcoord_avaliable(col);

    [csd(col).all, lfpAvg(col).all ] = perievent_CSD_LFP_amplitude_phase(permute(resps(xcoord == xcoord_avaliable(col),:,:), [2, 1, 3]), timestamps,[],'twin',[0.1 0.5]);


    % cd(fullfile(options.ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
    %             plot_perievent_CSD_LFP_amplitude_phase(lfpAvg,csd,power{nprobe},chan_config,sorted_config,best_channels{nprobe})

    plot_perievent_CSD_LFP(lfpAvg(col),csd(col),power{nprobe}(xcoord == xcoord_avaliable(col),:),chan_config,chan_config(xcoord == xcoord_avaliable(col),:),best_channels{nprobe},options)
    
    if ~isempty(clusters)
        [~,cluster_id] = intersect(clusters(nprobe).peak_channel,find(xcoord == xcoord_avaliable(col)));
        all_fields = fieldnames(clusters);
        clusters_info = struct();

        for nfield = 1:length(all_fields)
            if length(clusters(nprobe).(all_fields{nfield})) > length(cluster_id)
                clusters_info.(all_fields{nfield}) = clusters(nprobe).(all_fields{nfield})(cluster_id);
            end
        end

        plot_cluster_density_profile(power{nprobe}(xcoord == xcoord_avaliable(col),:),chan_config,chan_config(xcoord == xcoord_avaliable(col),:),best_channels{nprobe},clusters_info,options);
    end
end