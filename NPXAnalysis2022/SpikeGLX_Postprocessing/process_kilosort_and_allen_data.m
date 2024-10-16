function cluster_table = process_kilosort_and_allen_data(kilosort_directory)

% Written by EH, first push (xx/05/2022)

%% get required files
% To do:
% check if any files are missing
% read sample_rate from meta file.

meta_file = dir(fullfile(kilosort_directory,'*tcat.imec0.ap.meta'));

st_file = dir(fullfile(kilosort_directory,['spike_times.npy']));
clusterID_file = dir(fullfile(kilosort_directory,['spike_clusters.npy']));
chanPos_file = dir(fullfile(kilosort_directory,['channel_positions.npy']));
clusterTable_file = dir(fullfile(kilosort_directory,['clus_Table*.npy']));
clusterTable_file = clusterTable_file(end); % get most recent

clusterNoise_file = dir(fullfile(kilosort_directory,['cluster_group.tsv'])); % noise not noise
clusterKSlabel_file = dir(fullfile(kilosort_directory,['cluster_KSLabel.tsv']));

allenMetrics_file = dir(fullfile(kilosort_directory,['metrics*.csv']));
allenMetrics_file = allenMetrics_file(end); % get most recent.

allenWaveforms_file = dir(fullfile(kilosort_directory,['waveform_metrics*.csv']));
allenWaveforms_file = allenWaveforms_file(end);

meanWaveforms_file = dir(fullfile(kilosort_directory,['mean_waveforms*.npy']));
meanWaveforms_file = meanWaveforms_file(end); % get most recent

%% create the cluster table


%wf_temp = readtable(fullfile(allenWaveforms_file.folder, allenWaveforms_file.name)); % quality metrics table

% read in the metrics table
ct_temp = readtable(fullfile(allenMetrics_file.folder, allenMetrics_file.name)); % quality metrics table

% remove these variables that are not computed using ks3/pipeline
cluster_table = removevars(ct_temp, {'isolation_distance', 'l_ratio', 'd_prime', 'nn_hit_rate',...
    'nn_miss_rate', 'silhouette_score', 'max_drift', 'cumulative_drift', 'epoch_name_quality_metrics', ...
    'epoch_name_waveform_metrics'});

%% get spike times for each cluster

sample_rate = 30000; % change to read from .meta file
% load spike times in samples and convert to seconds
spike_times_samples = double(readNPY(fullfile(kilosort_directory, st_file.name)));
spike_times = spike_times_samples/sample_rate;
% get cluster_id for each spike time
cluster_id = double(readNPY(fullfile(kilosort_directory,  clusterID_file.name)));
% add table var for spike times
cluster_table = addvars(cluster_table, cell(height(cluster_table),1), 'NewVariableNames', 'spike_times');
% unique cluster ids in the cluster_table
uniqueClusters = unique(cluster_table.cluster_id);

for icluster = 1:numel(uniqueClusters)
    table_row = find(cluster_table.cluster_id==uniqueClusters(icluster)); % what is the table row for this cluster_id?
    tempSpikeTimes = spike_times(cluster_id==uniqueClusters(icluster)); % get spike_times that correspond to this cluster_id
    cluster_table.spike_times(table_row) = {tempSpikeTimes}; % add these spike_times to the correct table_row
end


%% get anatomical location of clusters

% x,y locations of each probe channel
[channel_locs] = readNPY(fullfile(kilosort_directory, chanPos_file.name));

cluster_locs = channel_locs(cluster_table.peak_channel+1,:); % get (x,y) channel location of each clusters peak_channel
shank_id = round(cluster_locs(:,1)/250); % get shank id from x-coord using fact that shanks are 250um apart.

% add these new vars to cluster_table
cluster_table.shank = shank_id;
cluster_table = addvars(cluster_table, cluster_locs, 'NewVariableNames', 'probe_location_xy');

%% get cluster labels (good, mua, noise)

% load cluster_KSlabel.tsv and cluster_group.tsv
% these files seem to have labels for all cluster ids (i.e. before
% kilosorts automerging. I'm not sure how this is possible, so proceed with
% caution using these labels for anything).
kslabels = tdfread(fullfile(kilosort_directory,clusterKSlabel_file.name));
noiselabels = tdfread(fullfile(kilosort_directory,clusterNoise_file.name));

% assign clusters as noise, or default ks info
for irow = 1:height(cluster_table)
    idx = find([kslabels.cluster_id] == cluster_table.cluster_id(irow));
    if strcmp(noiselabels.group(idx,:), 'noise')
        cluster_label{irow} = noiselabels.group(idx,:);
    else
        A = kslabels.KSLabel(idx,:);
        A = strrep(A,' ','');
        cluster_label{irow} = A;
    end
end

cluster_table = addvars(cluster_table, cluster_label', 'NewVariableNames', 'label');

%% get the mean waveforms for each cluster

mean_waveforms =  double(readNPY(fullfile(kilosort_directory, meanWaveforms_file.name)));
% cluster(0-index),channel,samples

for irow = 1:height(cluster_table)
    cid = cluster_table.cluster_id(irow);
    cpc = cluster_table.peak_channel(irow);
    waveforms(irow,:) = mean_waveforms(cid+1, cpc+1, :);
end

cluster_table = addvars(cluster_table, waveforms, 'NewVariableNames', 'mean_waveform');

