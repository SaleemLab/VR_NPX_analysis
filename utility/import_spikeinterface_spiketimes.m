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
function [these_spike_times,cluster_id,peak_channel,peak_depth,peak_channel_waveforms] = import_spikeinterface_spiketimes(options,SampleRate)

%%%%%%%

% Set default sampling rate for AP files
if ~exist('SampleRate','var') || isempty(SampleRate)
    SampleRate = 30000;
end

SORTER_DATAPATH = fullfile(options.SORTER_DATAPATH,'sorters',options.sorter_folder,'sorter_output');
gfileNum = str2num(cell2mat(extractBetween(options.EPHYS_DATAPATH,'_g','\')));% g file number
folder_names = cell2mat(extractBetween(options.EPHYS_DATAPATH,['_g',num2str(gfileNum),'\'],'_imec'));% to search for the session (e.g. 'M24010_20240231_1_g0')
segment_frames = readtable(options.segment_frames);% corresponding start and end sample point after concatenation via spike interface

%%%%%%%
% Get KS spike times
if contains(options.sorter_folder,'kilosort')
    spike_times = readNPY(fullfile(SORTER_DATAPATH,'spike_times.npy'));
    spike_clusters = readNPY(fullfile(SORTER_DATAPATH,'spike_clusters.npy'));

    %     cluster_group= tdfread(fullfile(SORTER_DATAPATH,'cluster_group.tsv')); % KS output cluster ID (before postprocessing)
    cluster_group = readtable(fullfile(options.SORTER_DATAPATH,'waveform',options.sorter_folder,'quality_metrics','metrics.csv')); % KS output cluster ID (before postprocessing)
    cluster_id = table2array(cluster_group(:,1));% 0 based original cluster id

    %%%%%%%
    % Work out where your spike times will be
    % Load the segment frame file for this sorting: check to see if it is a csv file
    % (ie. already parsed) or not
    if sum(contains(segment_frames.Properties.VariableNames,'segment_info')) > 0
        this_segment = strcmp(segment_frames.segment_info,folder_names);
        sampleStart = segment_frames.segmentStartFrame(this_segment);
        sampleEnd = segment_frames.segmentEndFrame(this_segment);
    end

    % Get the spike times for this file, referenced to start of the recording
    idx = spike_times > sampleStart & spike_times <= sampleEnd;
    spike_times = double(spike_times(idx)) - sampleStart;

    if options.probe_id ~= 0 % if probe 2 convert to probe 1 based time
        [~,fname] = fileparts(options.EPHYS_DATAPATH);
        if exist(fullfile(options.EPHYS_DATAPATH,[fname,'_aligned_AP_sample_number.mat'])) ~= 0
            load(fullfile(options.EPHYS_DATAPATH,[fname,'_aligned_AP_sample_number.mat'])); % For some sessions, this is not generated due to probe brokage
            spike_times = aligned_AP_sample_number(spike_times)';
        end
    end

    spike_clusters = spike_clusters(idx);
    spike_times = spike_times/SampleRate; % Convert to time

    % If there are spikes that are negative timestamped due to probe alignment,
    % remove them
    spike_clusters(spike_times<0) = [];
    spike_times(spike_times<0) = [];


    %%%%%%%
    % For each unit find central channel
    % Load the templates to define where on the array spikes are
    options.importMode = 'KS';
    [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,1);
    cluster_coordiantes = readNPY(fullfile(SORTER_DATAPATH,'channel_positions.npy')); % load all good channels (not noise channels) coordinate used for spike sorting
%     templates= readNPY(fullfile(SORTER_DATAPATH,'templates.npy'));
    templates = readNPY(fullfile(options.SORTER_DATAPATH,'waveform',options.sorter_folder,'templates_average.npy')); % Information about template waveform (for postprocessed clusters);

    for nchannel = 1:size(cluster_coordiantes,1)
        % for these good channels used for spike sorting, what is the actual channel number 
        all_good_channels(nchannel) = chan_config.Channel(chan_config.Ks_xcoord == cluster_coordiantes(nchannel,1)+11 & chan_config.Ks_ycoord == cluster_coordiantes(nchannel,2));
    end
    
    [~,peak_channel] = max(squeeze(mean(abs(templates),2)),[],2); % channel index providing the maximum of the mean of the absolute
   
    % For each unit present in (entire) recording [for consistency - ie, some
    % may not be responsive during this recording]
    nUnit = length(cluster_id);
    these_spike_times = cell(1,nUnit);
    for thisUnit = 1:nUnit
        these_spike_times{thisUnit} = spike_times(spike_clusters == cluster_id(thisUnit));
        peak_channel_waveforms(thisUnit,:) = templates(thisUnit,:,peak_channel(thisUnit));
    end
    
    peak_channel = all_good_channels(peak_channel); % peak channel for each cluster in terms of all channel index (not just good channel index)
    peak_depth = chan_config.Ks_ycoord(peak_channel); % peak depth for each cluster
    
else % reserved for future non-ks spike sorter

end

end
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


