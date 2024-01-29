% function [these_spike_times,nominal_KSLabel,cluster_id,peakChannel,maxSpkTime] = import_spikeinterface_spiketimes(options,SampleRate)
% Programme to extract spikeinterface sorter spiketimes for relevant file
%  Inputs
%       Path to kilosort folder for this session
%       'g' file number associated with this recording
%       [optional] file name where sample indices in concatenated file are
%           stored (default: 'catGT_number.mat')
%       [optional] AP sample rate for this recording
%           (default: 30000)
% SGS 5/4/2022 adapted from AT post_kilosort_processing
%
function [these_spike_times,cluster_id,peakChannel] = import_spikeinterface_spiketimes(options,SampleRate)

%%%%%%%
options.sorter

% Set default sampling rate for AP files
if ~exist('SampleRate','var') || isempty(SampleRate)
    SampleRate = 30000;
end
% If not defined use default name for catGT file
if ~exist('segment_frames','var') || isempty(KS_CATGT_FNAME)
    KS_CATGT_FNAME = 'catGT_number.mat';
end

SORTER_DATAPATH = fullfile(options.SORTER_DATAPATH,'sorters',options.sorter_folder,'sorter_output');
gfileNum = str2num(cell2mat(extractBetween(options.EPHYS_DATAPATH,'_g','\')));% g file number
folder_names = cell2mat(extractBetween(options.EPHYS_DATAPATH,'ephys\','\'));% to search for the session (e.g. '20240231_1')
segment_frames = readtable(options.segment_frames);% corresponding start and end sample point after concatenation via spike interface

%%%%%%%
% Get spike times
ks_spike_times = readNPY(fullfile(SORTER_DATAPATH,'spike_times.npy'));
% Get templates corresponding to each spike time
% spike_templates = readNPY(fullfile(SORTER_DATAPATH,'spike_templates.npy'));
spike_clusters = readNPY(fullfile(SORTER_DATAPATH,'spike_clusters.npy'));
similar_templates = readNPY(fullfile(SORTER_DATAPATH,'similar_templates.npy'));

spike_clusters = readNPY(fullfile(SORTER_DATAPATH,'templates.npy'));
tem1 = readNPY(fullfile(options.KS_DATAPATH,'templates.npy'));

% Get types of unit cluster Id
% cd('X:\ibn-vision\DATA\SUBJECTS\M22069\ephys\20221201\kilosort')
% cluster_group = tdfread(fullfile(KS_DATAPATH,'cluster_group.tsv'));
cluster_group = tdfread(fullfile(SORTER_DATAPATH,'cluster_KSLabel.tsv'));
%
% amplitudes = readNPY(fullfile(KS_DATAPATH,'amplitudes.npy'));


% cluster_ContamPct = tdfread(fullfile(KS_DATAPATH,'cluster_contam_rate.tsv'));
% cluster_ContamPct.cluster_id(cluster_ContamPct.contam_rate < 1);

%

% KS_DATAPATH = 'X:\ibn-vision\DATA\SUBJECTS\M22021\ephys\20220426\kilosort'

%%%%%%%
% Work out where your spike times will be
% Load the segment frame file for this sorting: check to see if it is a csv file
% (ie. already parsed) or not
if contains(segment_frames.Properties.VariableNames,'session') & contains(segment_frames.Properties.VariableNames,'gfileNum')
    this_recording = segment_frames.session == gfileNum & strcmp(segment_frames.session,folder_names);
    segmentStartFrame = segment_frames.segmentStartFrame(this_recording);
    segmentEndFrame = segment_frames.segmentEndFrame(this_recording);
end

% In both cases we are assuming that the first file in the list was ignored
% in the output log; and further that this was one g file less than the
% first in the list
if ~isempty(catGT_table) % If catGT file exists (with more than one recording)


    g_number = [g_number(1)-1 g_number(:)'];
    out_start_smp = [NaN out_start_smp(:)'];
    out_zeros_smp = [NaN out_zeros_smp(:)'];

    % Find where this file is in the concatenated list
    thisGFile = find(g_number == gfileNum);
    if thisGFile == 1 % First file in the concatenated data
        sampleStart = 1;
        sampleEnd = out_start_smp(2)-out_zeros_smp(2);
    elseif thisGFile == length(g_number) % Last file
        sampleStart = out_start_smp(end)+out_zeros_smp(end);
        sampleEnd = max(ks_spike_times);
    else
        sampleStart = out_start_smp(thisGFile)+out_zeros_smp(thisGFile);
        sampleEnd = out_start_smp(thisGFile+1)-1;
    end

else % Just one recording (grab whole things)
    sampleStart = 1;
    sampleEnd = max(ks_spike_times);
end
% Get the spike times for this file, referenced to start of the recording
idx = ks_spike_times >= sampleStart & ks_spike_times <= sampleEnd;
ks_spike_times = double(ks_spike_times(idx)) - sampleStart;


if options.probe_id ~= 0 % if probe 2 convert to probe 1 based time
    [~,fname] = fileparts(options.EPHYS_DATAPATH);
    if exist(fullfile(options.EPHYS_DATAPATH,[fname,'_aligned_AP_sample_number.mat'])) ~= 0
        load(fullfile(options.EPHYS_DATAPATH,[fname,'_aligned_AP_sample_number.mat'])); % For some sessions, this is not generated due to probe brokage
        ks_spike_times = aligned_AP_sample_number(ks_spike_times)';
    end
end

ks_spike_templates = spike_templates(idx);
ks_spike_times = ks_spike_times/SampleRate; % Convert to time

% If there are spikes that are negative timestamped due to probe alignment,
% remove them
ks_spike_templates(ks_spike_times<0) = [];
ks_spike_times(ks_spike_times<0) = [];

% For each unit present in (entire) recording [for consistency - ie, some
% may not be responsive during this recording]
nUnit = length(cluster_group.cluster_id);
these_spike_times = cell(1,nUnit);
for thisUnit = 1:nUnit
    these_spike_times{thisUnit} = ks_spike_times(ks_spike_templates == cluster_group.cluster_id(thisUnit));
end

%%%%%%%
% For each unit find central channel
% Load the templates to define where on the array spikes are
templates = readNPY(fullfile(options.KS_DATAPATH,'templates.npy')); % returns nUnits x nTimepoints x nChannels array
[~,peakChannel] = max(squeeze(mean(abs(templates),2)),[],2); % channel index providing the maximum of the mean of the absolute

%%%%% OTHER Features available in KS stored files
% % channel_map=readNPY([GLXDir,'\kilosort3\channel_map.npy']);
% % channel_positions=readNPY([GLXDir,'\kilosort3\channel_positions.npy']);
% % similar_templates = readNPY([GLXDir,'\kilosort3\similar_templates.npy']);
% % templates = readNPY([GLXDir,'\kilosort3\templates.npy']);
% % templates_ind = readNPY([GLXDir,'\kilosort3\templates_ind.npy']);
% % whitening_mat = readNPY([GLXDir,'\kilosort3\whitening_mat.npy']);
% % whitening_mat_inv = readNPY([GLXDir,'\kilosort3\whitening_mat_inv.npy']);
% % amplitudes = readNPY([GLXDir,'\kilosort3\amplitudes.npy']);
% % spike_clusters = readNPY([GLXDir,'\kilosort3\spike_clusters.npy']);
% % spike_templates = readNPY([GLXDir,'\kilosort3\spike_templates.npy']);
% % spike_times = readNPY([GLXDir,'\kilosort3\spike_times.npy']);
% % cluster_Amplitude = tdfread([GLXDir,'\kilosort3\cluster_Amplitude.tsv']);
% % cluster_ContamPct = tdfread([GLXDir,'\kilosort3\cluster_ContamPct.tsv']);
% % cluster_group = tdfread([GLXDir,'\kilosort3\cluster_group.tsv']);
% % cluster_KSLabel = tdfread([GLXDir,'\kilosort3\cluster_KSLabel.tsv']);



