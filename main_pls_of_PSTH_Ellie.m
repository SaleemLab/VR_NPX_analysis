%%% For PLS analysis of spiking in response to visual stimuli ERB 2025
% PLS regression is run to obtain PLS latent components that are maximised for covariance between stimulus identity (A, B, C
% or D) and spiking data during stimulus presentation.

addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))
rmpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis\toolboxes\LMT\drtoolbox\techniques') % to remove the LMT version of pca

%% setting metrics to screen good clusters
clear all
% 1/5 Choose your probe depth of interest
depth_for_analysis = 'V1'; % choose 'L4' or 'V1' or 'CA1' or 'Sub_CA1'
SUBJECTS = {'M00069'};
params = create_cluster_selection_params('sorting_option','ellie');
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);

%%% 2/5 choose mode for PLS - cross-validation or using all trials
pls_mode = 'cross_validation'; % cross_validation or 'full'
train_fraction = 0.8;
n_pls_components = 3;
rng(1); % for reproducibility, fix the random number generator for the cross-validation split

%%% 3/5
Stimulus_type = 'TRAIN_2'; % OMIT, 'TRAIN_1' 'GAVNIK_ABCD_1' MUST specify which recording _#
temporal_structure = 'trial spikecounts'; % 'trial spikecounts' or 'timebinned spikes' to analyse spiking per neuron per trial across timebins, 'trial spikecounts' to just consider mean spiking per timebin during each trial
% simple spikecounts gives better PCA silhouette scores for my TRAIN protocol 20250203
% timebinned spikes gives better PCA silhouette scores for GAVNIK protocol 20250217
z_score_period = 'entire_session'; % z score either over 'entire_session' or 'first30secs' or 'none' (for every stimulus recording
% session from 20250205 onward, I presented grey screen to the mouse for at least 30s before starting the stimulus. 'none' may be useful 
% to try for the aggregate TRAIN case across days)
psthBinSize = 0.01; % for GAVNIK protocol 20250217, using 1ms worsens cluster separation a lot. 0.02 worsens separation somewhat.
stim_window = [0.03 0.16];  % in seconds: [start, end] For GAVNIK protocol 20250217, [0.03 0.16] is superior  to a longer or shorter window, in terms of mean silhouette score and centroid separation
                            % For my TRAIN protocol 20250203, [0.02 0.21] seems best. But across training days an uptick in spiking develops in the greyscreen period following stimulus offset

cd('V:\Ellie\DATA\SUBJECTS\M00069\analysis\20251006\TRAIN_2') % 3/4 files will be saved here in the cd

%% SET THIS 4/5***** For NPX2.0 you will use a different L4 channel for each shank. Use CSD to estimate the best channel to use in L4
probe_type = 1; % NPX1.0 is type 0, NPX2.0 is type 1.
csd_file_path = fullfile(pwd, '..', 'earliest_V1sink_CSD.mat');
loaded_data = load(csd_file_path);
if probe_type == 0
    layerfour_channels = loaded_data.earliest_V1sink_CSD.overall_best_halfmax_channel;
    layerfour_channels = layerfour_channels(:);
    shank_ids = ones(size(layerfour_channels)); % dummy shank ID = 1
elseif probe_type == 1
    % NPX2.0 → multiple shanks
    % earliest_V1sink_CSD is a struct array, extract fields safely
    layerfour_channels = [loaded_data.earliest_V1sink_CSD.best_channel_this_shank]; 
    shank_ids = [loaded_data.earliest_V1sink_CSD.shank_id]; 
    % ensure column vectors
    layerfour_channels = layerfour_channels(:); 
    shank_ids = shank_ids(:);
end

for nsession = 9 %5/5 row number of recording date in "experiment_info" 
    session_info = experiment_info(nsession).session(strcmp(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(strcmp(experiment_info(nsession).StimulusName,Stimulus_type));
    
    for n = 1:length(session_info) % How many recording sessions 
        options = session_info(n).probe(1);
        subject_number = session_info(n).probe(1).SUBJECT;
        DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters = clusters_ks4;
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH, '..', 'earliest_V1sink_CSD.mat'))
        load(fullfile(options.ANALYSIS_DATAPATH, '..', 'depths_from_PSD.mat'))
        files = dir(fullfile(options.EPHYS_DATAPATH, '*ChanMap*.mat')); %the channel map has the y coordinate of each channel
        file_to_load = fullfile(options.EPHYS_DATAPATH, files(1).name); %load() does not accept wildcards like *
        load(file_to_load);
        if strcmp(subject_number, 'M00013') || strcmp(subject_number, 'M00014')
            load(fullfile(fullfile(fileparts(options.ANALYSIS_DATAPATH), 'OP_Tuning'), 'OP_tuning.mat'));
        end

        % --- Identify shanks and their layer 4 channels ---
        [unique_shanks, ia, ~] = unique(shank_ids, 'stable');
        best_channels_per_shank = layerfour_channels(ia); % one per shank
        
        % --- Preallocate shank-specific depth ranges ---
        L4_depth_range = cell(numel(unique_shanks), 1);
        V1_depth_range = cell(numel(unique_shanks), 1);
        CA1_depth_range = cell(numel(unique_shanks), 1);
        Sub_CA1_depth_range = cell(numel(unique_shanks), 1);

        % --- Compute per-shank depth ranges using CSD and PSD data ---
        for iShank = 1:numel(unique_shanks)
            this_shank   = unique_shanks(iShank);
            this_channel = best_channels_per_shank(iShank);
            
            if strcmp(subject_number, 'M00013') || strcmp(subject_number, 'M00014')
                Brain_surface_depth = depths_from_PSD.surface_depth_PSD;
                L4_channel_depth    = earliest_V1sink_CSD.overall_best_halfmax_depth;
                L5_depth            = depths_from_PSD.L5_depth_PSD;
                CA1_depth           = depths_from_PSD.CA1_depth_PSD;
            else    
                this_shank_name = ['shank', num2str(this_shank)];  % creates e.g. 'shank2'
                Brain_surface_depth = depths_from_PSD.(this_shank_name).surface_depth_PSD;
                L4_channel_depth    = earliest_V1sink_CSD(find([earliest_V1sink_CSD.shank_id] == this_shank, 1)).best_depth_this_shank;
                L5_depth            = depths_from_PSD.(this_shank_name).L5_depth_PSD;
                CA1_depth           = depths_from_PSD.(this_shank_name).CA1_depth_PSD;
            end 

            L4_depth_range{iShank}   = [L4_channel_depth - 60, L4_channel_depth + 60]; % giving electrodes the full extent of 120um inclusive 
            % (hence measuring spiking over depth range greater than 120um. As 60 is divisible by 15 and 10, this will give the same effective range for NPX1.0 (which has staggered electrodes every 10um down the shank) and NPX2.0 (which has electrode rows every 15um down each shank)
            V1_depth_range{iShank}   = [L5_depth - 330, L5_depth + 700];
            CA1_depth_range{iShank}  = [CA1_depth - 150, CA1_depth + 150];
            Sub_CA1_depth_range{iShank} = [min(CA1_depth_range{iShank}) - 1000, min(CA1_depth_range{iShank})];  
        end

        ordered_oris = unique(Task_info.stim_orientation, 'stable'); 
               
        params = create_cluster_selection_params('sorting_option','ellie');
        %params.orientation_tuned = ...
        
        switch depth_for_analysis
            case 'L4' 
                depth_ranges = L4_depth_range;
            case 'V1'
                depth_ranges = V1_depth_range;
            case 'CA1' 
                depth_ranges = CA1_depth_range;
            case 'Sub_CA1'
                depth_ranges = Sub_CA1_depth_range;
        end

        

        if (contains(Stimulus_type, 'GAVNIK_ABCD')) || (contains(Stimulus_type, 'TRAIN'))
            for nprobe = 1:length(clusters)
                selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters, which pass the set parameters
                
                depth_selected_clusters = selected_clusters; % initialize with same structure
                for np = 1:length(clusters)
                    sc = selected_clusters(np);
                    cluster_channels = clusters(np).peak_channel(sc.cluster_id);
                    cluster_depths = ycoords(cluster_channels); % get depths of selected clusters from ycoords (peak_depths from SI are 15 microns different..)
                    %peak_depths = clusters(np).peak_depth(sc.cluster_id); % get depths of selected clusters

                    % Map each cluster to its shank
                    cluster_shanks = kcoords(cluster_channels);
    
                    % Initialize container for depth-filtered cluster IDs
                    depth_cluster_ids = [];
                    
                    % Loop over each shank
                    for iShank = 1:numel(unique_shanks)
                        this_shank = unique_shanks(iShank);
                        this_range = depth_ranges{iShank};
            
                        % Logical mask for clusters on this shank and within this shank’s depth range
                        shank_mask = cluster_shanks == this_shank;
                        depth_mask = cluster_depths >= min(this_range) & cluster_depths <= max(this_range);
            
                        % Keep only clusters satisfying both conditions
                        keep_mask = shank_mask & depth_mask;
            
                        % Append these cluster IDs
                        depth_cluster_ids = [depth_cluster_ids; sc.cluster_id(keep_mask(:))];
                    end

                    % Keep only spike times and IDs for clusters at selected depth
                    depth_selected_clusters(np).cluster_id = depth_cluster_ids;
                    depth_selected_clusters(np).spike_times = sc.spike_times(ismember(sc.spike_id, depth_cluster_ids));
                    depth_selected_clusters(np).spike_id = sc.spike_id(ismember(sc.spike_id, depth_cluster_ids));
                end
                cluster_id = depth_selected_clusters(nprobe).cluster_id; % cluster_ids of units which pass the set parameters [NB these are one count higher than per zero-based pythonic SI output cluster IDs...]
                %If you want to use manually chosen clusters within the depth range:
                %cluster_id = [435 438 439 440 445 458 460 464 469 474 471 483 482 484 488 490 493 494 496 497 498 500 502 504 508 510 511 513 515 517 518 522 520 525 526]; % manually selected (visually responsive) for 20250217
                
                               
                %%%% Filter out untuned neurons - OP tuning was only run daily for M00013 and M00014
                % if strcmp(subject_number, 'M00013') || strcmp(subject_number, 'M00014')
                %     [is_member, loc_in_OP] = ismember(cluster_id, OP_tuning.cluster_id); % Get index of cluster_id within OP_tuning
                %     totalResponse = sum(OP_tuning.summedNetSpikeCounts_per_ori(:, loc_in_OP(is_member)), 1);
                %     visuallyresponsive_mask = abs(totalResponse) >= 2; % very lenient; summed net spiking across all orientations of more than 1 spike or less than -1 spike
                % 
                %     gOSI_values = OP_tuning(nprobe).gOSI(loc_in_OP(is_member));
                %     op_tuned_mask = gOSI_values > gOSI_threshold;
                %     gDSI_values = OP_tuning(nprobe).gDSI(loc_in_OP(is_member));
                %     direction_tuned_mask = gDSI_values > gDSI_threshold;
                % 
                %     combined_mask = visuallyresponsive_mask & op_tuned_mask & direction_tuned_mask;
                % 
                %     cluster_id = cluster_id(combined_mask);
                %     % --- Filter all OP_tuning fields based on the updated cluster_id list (i.e. excluding non-visually-responsive and untuned clusters)
                %     % Keep only the entries corresponding to the updated cluster_id
                %     [~, loc_in_OP_final] = ismember(cluster_id, OP_tuning(nprobe).cluster_id);
                %     filtered_OP = struct();
                %     filtered_OP.cluster_id = OP_tuning(nprobe).cluster_id(loc_in_OP_final);
                %     filtered_OP.OriAngle = OP_tuning(nprobe).OriAngle(loc_in_OP_final);
                %     filtered_OP.DirAngle = OP_tuning(nprobe).DirAngle(loc_in_OP_final);
                %     filtered_OP.gOSI = OP_tuning(nprobe).gOSI(loc_in_OP_final);
                %     filtered_OP.gDSI = OP_tuning(nprobe).gDSI(loc_in_OP_final);
                %     filtered_OP.summedNetSpikeCounts_per_ori = OP_tuning(nprobe).summedNetSpikeCounts_per_ori(:, loc_in_OP_final);
                %     filtered_OP.all_orientations = OP_tuning(nprobe).all_orientations;
                % end

                baseline_window = [0 30]; % in seconds - first 30s of recording is grey screen - can use for z-scoring

                % Define time_edges depending on z_score_period
                if contains(z_score_period, 'entire_session')
                    time_edges = 0:psthBinSize:max(depth_selected_clusters(nprobe).spike_times); %time bin edges from 0 to the max time in steps of ptshBinSize
                elseif contains(z_score_period, 'first30secs')
                    time_edges = baseline_window(1):psthBinSize:baseline_window(2);
                end
                
               
                % Count total number of stimulus presentations
                total_presentations = length(Task_info.stim_onset);
                
                if contains(temporal_structure, 'timebinned spikes')
                    n_timebins = round(diff(stim_window)/psthBinSize);  % e.g., for a 30-150ms window with 10ms bins, this is 12 bins
                    z_trial_responses = nan(total_presentations, n_timebins, length(cluster_id));
                    trial_responses = NaN(total_presentations, n_timebins, length(cluster_id));  % stores raw spike counts
                else
                    z_trial_responses = zeros(total_presentations, length(cluster_id));
                    trial_responses = NaN(total_presentations, length(cluster_id));  % stores raw spike counts
                end
                
                trial_labels = zeros(total_presentations, 1);
                trial_counter = 0;
                baseline_variance = nan(1, length(cluster_id)); % neurons with zero baseline variance can generate NaNs
                
                for ori = 1:length(ordered_oris) % loop through each unique stimulus orientation shown to the mouse during this task
                    onset_times = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));
                    n_trials = length(onset_times);

                    for trial_idx = 1:n_trials
                        trial_counter = trial_counter + 1;
                        trial_labels(trial_counter) = ordered_oris(ori);
                        z_trial_vector = zeros(1, length(cluster_id)); % one value per neuron
                        trial_vector = zeros(1, length(cluster_id)); % one value per neuron
                        
                        for nCluster = 1:length(cluster_id) % loop through each good cluster
                            spike_times_this_cluster = depth_selected_clusters(nprobe).spike_times(depth_selected_clusters(nprobe).spike_id == cluster_id(nCluster)); % extract the spike times for this cluster    
                            % Compute appropriate histogram counts for z-scoring
                            if contains(z_score_period, 'entire_session')
                                zscore_counts = histcounts(spike_times_this_cluster, time_edges);
                            elseif contains(z_score_period, 'first30secs')
                                baseline_spikes = spike_times_this_cluster(spike_times_this_cluster >= baseline_window(1) & spike_times_this_cluster <= baseline_window(2));
                                zscore_counts = histcounts(baseline_spikes, time_edges);
                            end

                            % store baseline variance for this neuron (only once)
                            if trial_counter == 1
                                baseline_variance(nCluster) = var(zscore_counts);
                            end
                            
                            % Get spike per trial (1 row per trial)
                            onset = onset_times(trial_idx);
                            % in the [psth..] function below, psth and spikeCounts are derived from binnedArray (psth is binnedArray./psthBinsize i.e. converted to Hz, and spikeCounts is the sum
                            % of the columns of binnedArray (spikeCounts = sum(binnedArray,2)) i.e. each timebin in the window of interest) for each trial (row)
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  onset, stim_window, psthBinSize); % for this orientation, get the binnedArray (matrix of spike counts per timebin)
                            
                            % Now z-score. With z-scoring, all units contribute equally to the PCA, regardless of firing rate.
                            z_binnedArray = (binnedArray - mean(zscore_counts))./std(zscore_counts);
                            
                            if contains(temporal_structure, 'timebinned spikes')
                                z_trial_responses(trial_counter, :, nCluster) = z_binnedArray;  % store timebin-wise vector for this neuron
                                trial_responses(trial_counter, :, nCluster) = binnedArray;
                            else
                                z_trial_vector(nCluster) = mean(z_binnedArray, 2); % mean binned spike count during this trial window for this cluster
                                trial_vector(nCluster) = sum(spikeCounts);  % total spikes during trial window
                            end
                        end

                        
                        % Fill in the current row of trial_responses (i.e., current trial) with the activity from all selected clusters (neurons).
                        if contains(temporal_structure, 'trial spikecounts')
                            z_trial_responses(trial_counter, :) = z_trial_vector; % trial_counter is a running index that tracks which trial we're on (regardless of orientation), : refers to all columns
                            trial_responses(trial_counter, :) = trial_vector;
                        end
                    end
                end
                
                neurons_to_keep = baseline_variance > 0;
                fprintf('Removed %d units with zero baseline variance.\n', sum(~neurons_to_keep));
                if contains(temporal_structure, 'timebinned spikes')
                    z_trial_responses = z_trial_responses(:, :, neurons_to_keep);
                    trial_responses   = trial_responses(:, :, neurons_to_keep);
                else
                    z_trial_responses = z_trial_responses(:, neurons_to_keep);
                    trial_responses   = trial_responses(:, neurons_to_keep);
                end
                cluster_id = cluster_id(neurons_to_keep);

                
                % now build X for PLS
               
                if contains(temporal_structure, 'timebinned spikes')
                     [num_trials, num_bins, num_neurons] = size(z_trial_responses);
                     %%%% RESHAPE the data to trials (rows) x features (columns, being spiking per timebin for each neuron sequentially) - consider info (features) for each trial as a whole to analyse how trials differ from each other
                    X = reshape(z_trial_responses, num_trials, num_bins * num_neurons);                    
                else
                    X = z_trial_responses;  % or reshaped version if time-binned
                end
              
                %%%% Now run PLS
                
                
                % Map orientations to condition labels
                cond_labels = zeros(size(trial_labels));
                for i = 1:numel(ordered_oris)
                    cond_labels(trial_labels == ordered_oris(i)) = i;
                end
                
                % One-hot encode Y
                Y = zeros(numel(cond_labels), 4);
                for c = 1:4
                    Y(cond_labels == c, c) = 1;
                end
                
                % Run PLS (XS is trial scores in PLS space (like PCA scores, but supervised); 
                % XL is neural weights (which neurons define the axes). Y is stimulus identity
                %PCTVAR is how much variance in X and Y each component explains
                
                if strcmp(pls_mode, 'cross_validation')
                    num_trials = size(X,1);
                    idx = randperm(num_trials);
                    
                    n_train = floor(train_fraction * num_trials);
                    train_idx = idx(1:n_train);
                    test_idx  = idx(n_train+1:end);
                    
                    X_train = X(train_idx,:);
                    Y_train = Y(train_idx,:);
                    
                    X_test  = X(test_idx,:);
                    Y_test  = Y(test_idx,:);
                    true_labels = cond_labels(test_idx);
                    
                    [XL, YL, XS_train, YS_train, beta, PCTVAR, ~, stats] = ...
                        plsregress(X_train, Y_train, n_pls_components);

                    % project held-out trials   
                    X_test_aug = [ones(size(X_test,1),1) X_test];
                    Y_hat = X_test_aug * beta;   % predicted stimulus evidence
                    [~, decoded_labels] = max(Y_hat, [], 2);

                    cv_results = struct();
                    cv_results.true_labels = true_labels;
                    cv_results.decoded_labels = decoded_labels;
                    cv_results.Y_hat = Y_hat;
                    cv_results.PCTVAR = PCTVAR;
                    cv_results.beta = beta;

                elseif strcmp(pls_mode, 'full')         
                    [XL, YL, XS, YS, beta, PCTVAR, ~, stats] = plsregress(X, Y, n_pls_components); %  start with 3 PLS components
                    pls_model = struct();
                    pls_model.W = stats.W; % projection weights               
                    pls_model.XL = XL; % X loadings
                    pls_model.beta = beta; % regression coefficients
                    pls_model.PCTVAR = PCTVAR; % variance explained
                    pls_model.n_components = n_pls_components;
                    
                    pls_model.X = z_trial_responses;       % z-scored input for PLS
                    pls_model.Y = Y;                       % one-hot stimulus labels
                    pls_model.trial_responses = trial_responses; %raw spike counts
                    
                    pls_model.ordered_oris = ordered_oris;
                    pls_model.temporal_structure = temporal_structure;
                    pls_model.stim_window = stim_window;
                    pls_model.psthBinSize = psthBinSize;
                    pls_model.cluster_id = cluster_id;
                    pls_model.session_info = session_info;
                    pls_model.created_on = datestr(now);

                    pls_model.temporal_structure = temporal_structure;
                    pls_model.n_neurons = numel(cluster_id);
                    pls_model.z_score_period = z_score_period;

                    pls_filename = sprintf('pls_V1_%s_%s.mat', session_info.probe.SESSION, session_info.probe.StimulusName);
                    save(pls_filename, 'pls_model');
                end

                %%% ---- Centroid calculations based on 3 PLS components ----

                PLS_for_centroids = n_pls_components;
                
                % XS = trial scores from plsregress (trial × PLS components)
                %centroid_score_pls = XS(:, 1:PLS_for_centroids);

                if strcmp(pls_mode, 'cross_validation')
                    XS_train_plot = XS_train(:, 1:PLS_for_centroids);
                
                    % Project held-out trials into PLS space
                    Xmean_train = mean(X_train, 1);
                    X_test_centered = X_test - Xmean_train;
                    XS_test_plot = X_test_centered * stats.W(:, 1:PLS_for_centroids);
                
                    labels_train = cond_labels(train_idx);
                    labels_test  = cond_labels(test_idx);
                else
                    XS_train_plot = XS(:, 1:PLS_for_centroids);
                    labels_train = cond_labels;
                end
                
                centroid_score_pls = XS_train_plot;
                
                pls_centroids = zeros(length(ordered_oris), PLS_for_centroids);
                
                for i = 1:length(ordered_oris)
                    %idx = trial_labels == ordered_oris(i);
                    idx = labels_train == i;
                    pls_centroids(i, :) = mean(centroid_score_pls(idx, :), 1);
                end
                
                %% ---- Centroid-to-centroid distances ----
                centroid_distances = pdist(pls_centroids, 'euclidean');
                mean_centroid_distance = mean(centroid_distances);
                
                disp(['[PLS] Mean centroid-to-centroid distance: ', ...
                      num2str(mean_centroid_distance)]);
                
                %% ---- Within-cluster spread ----
                within_cluster_spread = zeros(length(ordered_oris), 1);
                
                for i = 1:length(ordered_oris)
                    %idx = trial_labels == ordered_oris(i);
                    idx = labels_train == i;
                    distances_to_centroid = sqrt(sum( ...
                        (centroid_score_pls(idx, :) - pls_centroids(i, :)).^2, 2));
                    within_cluster_spread(i) = mean(distances_to_centroid);
                end
                
                mean_within_cluster_spread = mean(within_cluster_spread);
                
                disp(['[PLS] Mean within-cluster spread: ', ...
                      num2str(mean_within_cluster_spread)]);
                
                %% ---- Global separation score ----
                separation_score = mean_centroid_distance / mean_within_cluster_spread;
                
                disp(['[PLS] Separation score (higher is better): ', ...
                      num2str(separation_score)]);
                
                %% ---- Per-class separation scores ----
                per_class_sep_scores = zeros(length(ordered_oris), 1);
                
                for i = 1:length(ordered_oris)
                    idx_i = labels_train == i;
                    centroid_i = pls_centroids(i, :);
                
                    % Spread of cluster i
                    spread_i = mean(sqrt(sum( ...
                        (centroid_score_pls(idx_i, :) - centroid_i).^2, 2)));
                
                    % Mean distance to all other centroids
                    other_centroids = pls_centroids(setdiff(1:length(ordered_oris), i), :);
                    dist_to_others = sqrt(sum((other_centroids - centroid_i).^2, 2));
                    mean_dist_to_others = mean(dist_to_others);
                
                    per_class_sep_scores(i) = mean_dist_to_others / spread_i;
                end


                % Define unique orientations and colors
                colors = lines(length(ordered_oris));  % One distinct color per orientation in order
                % Generate legend labels: A 30°, B 60°, etc.
                letters = 'A':'D';
                if contains(Stimulus_type, 'TRAIN')
                    ori_labels = arrayfun(@(i, ori) sprintf('%s %.0f°', letters(i), round(rad2deg(ori))), ...
                                      (1:length(ordered_oris))', ordered_oris, 'UniformOutput', false);
                else 
                    ori_labels = arrayfun(@(i, ori) sprintf('%s %.0f°', letters(i), ori), ...
                                      (1:length(ordered_oris))', ordered_oris, 'UniformOutput', false);
                end
                
                %%%% ---- Plot PLS (z-scored data) ----
                fig_pls = figure;
                cla;
                hold on;
                
                for i = 1:length(ordered_oris)
                    idx = labels_train == i;
                                
                    scatter3( ...
                        XS_train_plot(idx,1), ...
                        XS_train_plot(idx,2), ...
                        XS_train_plot(idx,3), ...
                        36, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.6);
                end
                hold on
                if strcmp(pls_mode, 'cross_validation') % outlined markers for test data
                    for i = 1:length(ordered_oris)
                        idx = labels_test == i;
                
                        scatter3( ...
                            XS_test_plot(idx,1), ...
                            XS_test_plot(idx,2), ...
                            XS_test_plot(idx,3), ...
                            36, colors(i,:), ...
                            'o', 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
                    end
                end    


                xlabel(sprintf('PLS 1 (%.1f%% Y-var)', PCTVAR(2,1)));
                ylabel(sprintf('PLS 2 (%.1f%% Y-var)', PCTVAR(2,2)));
                zlabel(sprintf('PLS 3 (%.1f%% Y-var)', PCTVAR(2,3)));
                
                grid on;
                view(3);
                axis equal;
                
                %% ---- plot Centroids ----
                hold on;                
                h_centroids = scatter3(pls_centroids(:,1), pls_centroids(:,2), pls_centroids(:,3), ...
                    120, colors, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                
                legend_entries = ori_labels;    
                h_legend = zeros(1,length(legend_entries));

                for i = 1:length(ordered_oris)
                    %idx = labels_train == i;
                    h_legend(i) = scatter3(NaN, NaN, NaN, 36, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.6);
                end

                hold on;

                if strcmp(pls_mode, 'cross_validation')
                    legend_entries{end+1} = 'Training centroids';

                    % Add "Test datapoints" entry manually
                    legend_entries{end+1} = 'Test datapoints';
                    
                    %% Training centroids handle (largest, 120)
                    h_legend(length(ordered_oris)+1) = scatter3(NaN, NaN, NaN, 120, colors(1, :), 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor','k','LineWidth',1.5);
                
                    % Test datapoints handle (medium, 36, hollow)
                    h_legend(length(ordered_oris)+2) = scatter3(NaN, NaN, NaN, 36, colors(1, :), 'o', 'LineWidth', 1.5, 'MarkerFaceColor', 'none');

                
                    legend(h_legend, legend_entries, 'Location', 'southeast');
                else
                    legend([ori_labels', {'Centroids'}], 'Location', 'southeast');
                end
                
                %% ---- Format per-class separation scores ----
                per_class_lines = cell(length(per_class_sep_scores), 1);
                for i = 1:length(per_class_sep_scores)
                    per_class_lines{i} = sprintf('%s separation score: %.4f', ...
                                                 ori_labels{i}, per_class_sep_scores(i));
                end
                
                per_class_text = strjoin(per_class_lines, '\n');
                
                %% ---- Combine metrics into annotation text ----
                if strcmp(pls_mode, 'cross_validation')
                    header_line = sprintf('PLS metrics derived using %.0f%% of trials:\n\n', train_fraction*100);
                else
                    header_line = '';
                end
                
                metrics_text = sprintf([ ...
                    '%s' ...
                    'Mean centroid-to-centroid distance: %.4f\n' ...
                    'Mean within-cluster spread: %.4f\n' ...
                    'Separation score: %.4f\n\n%s'], ...
                    header_line, mean_centroid_distance, mean_within_cluster_spread, ...
                    separation_score, per_class_text);
                
                %% ---- Annotation textbox (top-right, normalized coords) ----
                n_types = length(per_class_sep_scores);
                box_height = 0.15 + 0.03 * n_types;
                
                x_pos = 0.70;
                y_top = 0.90;
                y_pos = y_top - box_height;
                
                annotation('textbox', [x_pos, y_pos, 0.28, box_height], ...
                           'String', metrics_text, ...
                           'FitBoxToText', 'on', ...
                           'BackgroundColor', 'white', ...
                           'EdgeColor', 'black', ...
                           'FontSize', 10, ...
                           'VerticalAlignment', 'top', ...
                           'Interpreter', 'none');
                %%% ---- Metrics textbox for test trials (SW) ----
                if strcmp(pls_mode,'cross_validation')
                    % ---- Compute classification accuracy of test trials ----
                    pred_labels_test = zeros(size(labels_test));
                    for i = 1:length(labels_test)
                        % assign test trial to nearest training centroid
                        [~, pred_labels_test(i)] = min(vecnorm(pls_centroids - XS_test_plot(i,:), 2, 2));
                    end
                    accuracy_test = mean(pred_labels_test == labels_test);  % fraction correct
                
                    % ---- Text for SW box ----
                    metrics_text_test = sprintf('Test trial classification based on nearest\ntraining centroid: %.1f%% accuracy', ...
                                                accuracy_test*100);
                
                    % ---- Place box in SW ----
                    annotation('textbox', [0.72, 0.45, 0.20, 0.1], ... % smaller height
                               'String', metrics_text_test, ...
                               'FitBoxToText', 'on', ...
                               'BackgroundColor', 'white', ...
                               'EdgeColor', 'black', ...
                               'FontSize', 10, ...
                               'VerticalAlignment', 'top', ...
                               'Interpreter', 'none');   
                end    
                hold off;

                if strcmp(pls_mode, 'cross_validation')
                    pls_desc = sprintf('PLS (train %.0f%% / test %.0f%%)', ...
                                       train_fraction*100, (1-train_fraction)*100);
                else
                    pls_desc = 'PLS (all trials)';
                end

                
                sgtitle(sprintf(['%s - %s: %s of %s %s ' ...
                                 '(%.2f–%.2fs after stim onset) ' ...
                                 'Z-scored over %s, Day %d'], ...
                                 subject_number, Stimulus_type, ...
                                 pls_desc, ...
                                 depth_for_analysis, temporal_structure, ...
                                 stim_window(1), stim_window(2), ...
                                 z_score_period, ...
                                 experiment_info(nsession).date), ...
                        'Interpreter', 'none');
            
                
                if strcmp(pls_mode, 'cross_validation')
                    pls_tag = sprintf('PLS_Cross_Validation_train%.0f', train_fraction*100);
                else
                    pls_tag = 'PLS_full';
                end
                
                fig_filename1 = sprintf(['%s - %s - %s of %s %s ' ...
                                         '(%.2f–%.2fs after stim onset).fig'], ...
                                         subject_number, Stimulus_type, ...
                                         pls_tag, ...
                                         depth_for_analysis, temporal_structure, ...
                                         stim_window(1), stim_window(2));             
                
                savefig(fig_pls, fullfile(pwd, fig_filename1));
            end
        end
    end
end
