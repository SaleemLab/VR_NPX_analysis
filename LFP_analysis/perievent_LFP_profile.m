function [lfpAvg,csd,PSD] = perievent_LFP_profile(event_name,stimTimes,AnalysisTimeWindow,PSD,clusters,options,varargin)

p = inputParser;
addParameter(p,'place_fields',[],@isstruct);
% addParameter(p,'doDetrend',false,@islogical);
addParameter(p,'filter_type',[],@iscell);
addParameter(p,'filter_freq',[0.5 3;4 12;9 17;30 60;60 100;125 300; 300 600],@isnumeric) % Frequency range: [freq1(1) freq1;freq2 freq2;......]
addParameter(p,'CSD_V1_CA1_normalisation',0,@isnumeric) % Normalised CSD within region or not
addParameter(p,'x_col',1,@isnumeric) % which x_col to use


parse(p,varargin{:})
filter_type = p.Results.filter_type;
filter_freq = p.Results.filter_freq;
place_fields = p.Results.place_fields;
CSD_V1_CA1_normalisation = p.Results.CSD_V1_CA1_normalisation;
x_col  = p.Results.x_col;

extractedTimeWindow = AnalysisTimeWindow;
extractedTimeWindow(1) = extractedTimeWindow(1)-5;
extractedTimeWindow(2)= extractedTimeWindow(2)+5;

if ~isempty(place_fields)

end

% 
% if isfield(clusters,'merged_spike_id')
%     clusters.cluster_id = unique(clusters.merged_cluster_id);
%     for ncell = 1:length(unique(clusters.merged_cluster_id))
%         tempt_peak_channel = clusters.peak_channel(clusters.merged_cluster_id == clusters.cluster_id(ncell));
%         tempt_peak_depth = clusters.peak_depth(clusters.merged_cluster_id == clusters.cluster_id(ncell));
%         tempt_peak_waveform = clusters.peak_channel_waveforms(clusters.merged_cluster_id == clusters.cluster_id(ncell),:);
%         tempt_cell_type = clusters.cell_type(clusters.merged_cluster_id == clusters.cluster_id(ncell));
% 
%         if length(tempt_peak_depth)== 2
%             tempt_peak_channel = tempt_peak_channel(1);
%             tempt_peak_depth = tempt_peak_depth(1);
%             tempt_peak_waveform = tempt_peak_waveform(1,:);
%             tempt_cell_type = tempt_cell_type(1);
%         else % find median peak depth assign that value to the unit
%             [~,index]= min(tempt_peak_depth - median(tempt_peak_depth));
%             tempt_peak_channel = tempt_peak_channel(index);
%             tempt_peak_depth = tempt_peak_depth(index);
%             tempt_peak_waveform = tempt_peak_waveform(index,:);
%             tempt_cell_type = tempt_cell_type(index);
%         end
% 
%         merged_peak_channel(ncell) = tempt_peak_channel;
%         merged_peak_depth(ncell) = tempt_peak_depth;
%         merged_peak_waveform(ncell,:) = tempt_peak_waveform;
%         merged_cell_type(ncell,:) = tempt_cell_type;
%     end
% 
%     clusters.peak_channel = merged_peak_channel;
%     clusters.peak_depth = merged_peak_depth;
%     clusters.peak_channel_waveforms = merged_peak_waveform;
%     clusters.cell_type = merged_cell_type;
% 
%     clusters.spike_id = clusters.merged_spike_id;
%     clusters.cluster_id = unique(clusters.merged_cluster_id);
% end

BinWidth = 1/1250;
% options.importMode = 'LF'; % LF or MUA or KS
% AnalysisTimeWindow two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])

%             options.ks_unitType = 'good'; % 'mua', 'good' or ''
nprobe = options.probe_id+1; % probe_no is [1,2] based on options.probe_id (0 and 1)

power = [];

DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'..','best_channels.mat'));
if ~isempty(DIR)
    load(fullfile(options.ANALYSIS_DATAPATH,'..','best_channels.mat'))

    if nprobe > 1 & length(best_channels) < nprobe
        best_channels{nprobe} = [];
    end
else
    best_channels{nprobe} = [];
end

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


% Extract LFP data relative to the EVENT onset time
% Preallocate matrix for speed - so use the first stimulus to work
% out how much space we will need for the resampled data
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
nSamp = fix(SampRate(imecMeta)*range(extractedTimeWindow));
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
timeVector = linspace(extractedTimeWindow(1),extractedTimeWindow(2),size(tresps,2));

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
    tresps = ReadBin(fix(SampRate(imecMeta)*(stim_on+extractedTimeWindow(1))), nSamp, imecMeta, file_to_use, options.EPHYS_DATAPATH);
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


timestamps = linspace(extractedTimeWindow(1),extractedTimeWindow(2),size(resps,2));

for nchannel = 1:size(chan_config,1)
    power(nchannel,:) = PSD{nprobe}(nchannel).mean_power;
    xcoord(nchannel) = PSD{nprobe}(nchannel).xcoord;
    ycoord(nchannel) = PSD{nprobe}(nchannel).ycoord;
    shanks(nchannel) = PSD{nprobe}(nchannel).shank;
end
xcoord_avaliable = unique(xcoord);

% sort channel according to y coordinate
[ycoord idx] = sort(ycoord,'ascend');
xcoord = xcoord(idx);
power = power(idx,:);
chan_config = chan_config(idx,:);
resps = resps(idx,:,:);
% resps1 = resps;

resps = permute(resps(xcoord == xcoord_avaliable(x_col),:,:), [2, 1, 3]);

lfpAvg = [];
csd = [];
% lfpAvg.filter_type = {'all'};
lfpAvg.filter_type = filter_type;
lfpAvg.event_group = {event_name};
lfpAvg.probe_id = options.probe_id;
%     lfpAvg(col).column = col;
lfpAvg.xcoord = xcoord_avaliable(1);
AnalysisTimeWindow(1) = -AnalysisTimeWindow(1);

for type = 1:length(filter_type)
    [csd.(filter_type{type}), lfpAvg.(filter_type{type}) ] = perievent_CSD_LFP_amplitude_phase(resps, timestamps,[],'twin',AnalysisTimeWindow,'filter',filter_freq(type,:));


    % cd(fullfile(options.ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
    %             plot_perievent_CSD_LFP_amplitude_phase(lfpAvg,csd,power{nprobe},chan_config,sorted_config,best_channels{nprobe})

end

lfpAvg.filter_type(length(filter_type)+1) = {'all'};
[csd.all, lfpAvg.all ] = perievent_CSD_LFP_amplitude_phase(resps, timestamps,[],'twin',AnalysisTimeWindow);

if CSD_V1_CA1_normalisation==1 %plot both versions
    plot_perievent_CSD_LFP_amplitude_phase(lfpAvg,csd,power(xcoord == xcoord_avaliable(x_col),:),chan_config,chan_config(xcoord == xcoord_avaliable(x_col),:),best_channels{nprobe},options)
    options.CSD_V1_CA1_normalisation = 1;
    plot_perievent_CSD_LFP_amplitude_phase(lfpAvg,csd,power(xcoord == xcoord_avaliable(x_col),:),chan_config,chan_config(xcoord == xcoord_avaliable(x_col),:),best_channels{nprobe},options)
else
    plot_perievent_CSD_LFP_amplitude_phase(lfpAvg,csd,power(xcoord == xcoord_avaliable(x_col),:),chan_config,chan_config(xcoord == xcoord_avaliable(x_col),:),best_channels{nprobe},options)
end

% if ~isempty(clusters)
%     [~,cluster_id] = intersect(clusters(nprobe).peak_channel,find(xcoord == xcoord_avaliable(col)));
%     all_fields = fieldnames(clusters);
%     clusters_info = struct();
% 
%     for nfield = 1:length(all_fields)
%         if length(clusters(nprobe).(all_fields{nfield})) > length(cluster_id)
%             clusters_info.(all_fields{nfield}) = clusters(nprobe).(all_fields{nfield})(cluster_id);
%         end
%     end
% 
%     plot_cluster_density_profile(power{nprobe}(xcoord == xcoord_avaliable(col),:),chan_config,chan_config(xcoord == xcoord_avaliable(col),:),best_channels{nprobe},clusters_info,options);
% end
