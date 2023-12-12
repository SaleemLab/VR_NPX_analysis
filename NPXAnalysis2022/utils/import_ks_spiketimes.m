% Programme to extract KS spiketimes for relevant file
%  Inputs
%       Path to kilosort folder for this session
%       'g' file number associated with this recording
%       [optional] file name where sample indices in concatenated file are
%           stored (default: 'catGT_number.mat')
%       [optional] AP sample rate for this recording
%           (default: 30000)
% SGS 5/4/2022 adapted from AT post_kilosort_processing
%
function [these_spike_times,nominal_KSLabel,cluster_id,peakChannel,maxSpkTime] = import_ks_spiketimes(options,gfileNum,KS_CATGT_FNAME,SampleRate)

%%%%%%%
% Set default sampling rate for AP files
if ~exist('SampleRate','var') || isempty(SampleRate)
    SampleRate = 30000;
end 
% If not defined use default name for catGT file
if ~exist('KS_CATGT_FNAME','var') || isempty(KS_CATGT_FNAME)
    KS_CATGT_FNAME = 'catGT_number.mat';
end
KS_DATAPATH = options.KS_DATAPATH;
gfileNum = str2num(cell2mat(extractBetween(options.EPHYS_DATAPATH,'_g','\')));

%%%%%%%
% Get spike times
ks_spike_times = readNPY(fullfile(KS_DATAPATH,'spike_times.npy'));
% Get templates corresponding to each spike time
spike_templates = readNPY(fullfile(KS_DATAPATH,'spike_templates.npy'));
spike_clusters = readNPY(fullfile(KS_DATAPATH,'spike_clusters.npy'));
% Get types of unit cluster Id
% cd('X:\ibn-vision\DATA\SUBJECTS\M22069\ephys\20221201\kilosort')
% cluster_group = tdfread(fullfile(KS_DATAPATH,'cluster_group.tsv'));
cluster_group = tdfread(fullfile(KS_DATAPATH,'cluster_KSLabel.tsv'));
% 
% amplitudes = readNPY(fullfile(KS_DATAPATH,'amplitudes.npy'));


% cluster_ContamPct = tdfread(fullfile(KS_DATAPATH,'cluster_contam_rate.tsv'));
% cluster_ContamPct.cluster_id(cluster_ContamPct.contam_rate < 1);

% 

% KS_DATAPATH = 'X:\ibn-vision\DATA\SUBJECTS\M22021\ephys\20220426\kilosort'

if ~isfield(cluster_group,'KSLabel') & isfield(cluster_group,'group')
    cluster_group.KSLabel = cluster_group.group;
end

%%%%%%%
% Work out where your spike times will be
% Load the CatGT file for this sorting: check to see if it is a mat file
% (ie. already parsed) or not
if sum(contains(KS_CATGT_FNAME,'.mat'))
   load(fullfile(KS_DATAPATH,KS_CATGT_FNAME),'g_number','out_start_smp','out_zeros_smp');
   if length(g_number) == (length(out_start_smp)+1) % in this case it is probable that the first was ignored
       g_number(1) = [];
   end
else
    catGT_table = readCatGTlog(fullfile(KS_DATAPATH,KS_CATGT_FNAME));
    if length(catGT_table.g_idx) ~= length(unique(catGT_table.g_idx)) % If not one probe
        if contains(KS_DATAPATH,'probe_2')
            catGT_table = catGT_table(length(unique(catGT_table.g_idx))+1:length(catGT_table.g_idx),:);
        else
            catGT_table = catGT_table(1:length(unique(catGT_table.g_idx)),:);
        end
    end
        g_number = catGT_table.g_idx;
        out_start_smp = catGT_table.out_start_smp;
        out_zeros_smp = catGT_table.out_zeros_smp;
end

% In both cases we are assuming that the first file in the list was ignored
% in the output log; and further that this was one g file less than the
% first in the list
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

% Get the spike times for this file, referenced to start of the recording
idx = ks_spike_times >= sampleStart & ks_spike_times <= sampleEnd;
ks_spike_times = double(ks_spike_times(idx)) - sampleStart;


if contains(KS_DATAPATH,'probe_2') % if probe 2 convert to probe 1 based time
    [~,fname] = fileparts(options.EPHYS_DATAPATH);
    if exist(fullfile(options.EPHYS_DATAPATH,[fname,'_aligned_AP_sample_number.mat'])) ~= 0
        load(fullfile(options.EPHYS_DATAPATH,[fname,'_aligned_AP_sample_number.mat'])); % For some sessions, this is not generated due to probe brokage
        ks_spike_times = aligned_AP_sample_number(ks_spike_times)';
    end
end

ks_spike_templates = spike_templates(idx);
maxSpkTime = max(ks_spike_times);
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
templates = readNPY(fullfile(KS_DATAPATH,'templates.npy')); % returns nUnits x nTimepoints x nChannels array
[~,peakChannel] = max(squeeze(mean(abs(templates),2)),[],2); % channel index providing the maximum of the mean of the absolute
% And retrieve class of unit 
FieldNames = fields(cluster_group);                    % FR: Added this if statement to deal with this version of KS
if ismember('KSLabel',FieldNames)
    nominal_KSLabel = nominal(cluster_group.KSLabel); % should have nominals of 'mua', 'good'
else
    nominal_KSLabel = nominal(cluster_group.group);
end
cluster_id = cluster_group.cluster_id;

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



