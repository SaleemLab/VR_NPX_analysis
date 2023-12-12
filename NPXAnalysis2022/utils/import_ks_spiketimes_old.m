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
function [these_spike_times,nominal_KSLabel,cluster_id,peakChannel,maxSpkTime] = import_ks_spiketimes_old(KS_DATAPATH,gfileNum,KS_CATGT_FNAME,SampleRate)

%%%%%%%
% Set default sampling rate for AP files
if ~exist('SampleRate','var') || isempty(SampleRate)
    SampleRate = 30000;
end 
% If not defined use default name for catGT file
if ~exist('KS_CATGT_FNAME','var') || isempty(KS_CATGT_FNAME)
    KS_CATGT_FNAME = 'catGT_number.mat';
end


%%%%%%%
% Get spike times
ks_spike_times = readNPY(fullfile(KS_DATAPATH,'spike_times.npy'));
% Get templates corresponding to each spike time
spike_templates = readNPY(fullfile(KS_DATAPATH,'spike_templates.npy'));
% Get types of unit cluster Id
cluster_group = tdfread(fullfile(KS_DATAPATH,'cluster_group.tsv'));

%%%%%%%
% Work out where your spike times will be
% Load the CatGT file for this sorting
%load(fullfile(KS_DATAPATH,KS_CATGT_FNAME),'g_number','out_start_smp','out_zeros_smp','inp_gamp_smp');

catGT_table = readCatGTlog(fullfile(KS_DATAPATH,KS_CATGT_FNAME));
g_number = catGT_table.g_idx;
out_start_smp = catGT_table.out_start_smp;
out_zeros_smp = catGT_table.out_zeros_smp;
% Find where this file is in the concatenated list
thisGFile = find(g_number == gfileNum);
if thisGFile == 1 % First file in the concatenated data
    sampleStart = 1;
    sampleEnd = out_start_smp(1)-1;
elseif thisGFile == length(g_number) % Last file
    sampleStart = out_start_smp(end)+out_zeros_smp(end);
    sampleEnd = max(ks_spike_times);
else
    sampleStart = out_start_smp(thisGFile)+out_zeros_smp(thisGFile);
    sampleEnd = out_start_smp(thisGFile+1)-1;
end

% Get the spike times for this file, referenced to start of the recording
idx = ks_spike_times >= sampleStart & ks_spike_times <= sampleEnd;
ks_spike_times = double(ks_spike_times(idx) - sampleStart)/SampleRate;
ks_spike_templates = spike_templates(idx);
maxSpkTime = max(ks_spike_times);

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



