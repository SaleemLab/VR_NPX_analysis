%%% For PCA of spiking in response to visual stimuli

addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

%% setting metrics to screen good clusters
clear all
% Choose your probe depth of interest
L4_depth_range = 4440:4580; % 1/5. um. Set for each SESSION based on CSD +/- 70um
V1_depth_range = (min(L4_depth_range) - 400) : (max(L4_depth_range) + 500); 
CA1_depth_range = 3640:3940; % 2/5. um. Set for each SESSION based on PSD; ~300um around Ripple power "bump"
Sub_CA1_depth_range = 1550:(min(CA1_depth_range));
depth_for_analysis = 'V1'; % choose 'L4' or 'V1' or 'CA1' or 'Sub_CA1'

SUBJECTS = {'M00013'};

params = create_cluster_selection_params('sorting_option','ellie');
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);

%%% 3/5
Stimulus_type = 'TRAIN'; % OMIT
temporal_structure = 'trial spikecounts'; % 'trial spikecounts' or 'timebinned spikes' to analyse spiking per neuron per trial across timebins, 'trial spikecounts' to just consider mean spiking per timebin during each trial
% simple spikecounts gives better silhouette scores for my TRAIN protocol 20250203
% timebinned spikes gives better silhouette scores for GAVNIK protocol 20250217

z_score_period = 'entire_session'; % z score either over 'entire_session' or 'first30secs' or 'none' (for every stimulus recording
% session from 20250205 onward, I presented grey screen to the mouse for at least 30s before starting the stimulus. 'none' may be useful 
% to try for the aggregate TRAIN case across days)
psthBinSize = 0.01; % for GAVNIK protocol 20250217, using 1ms worsens cluster separation a lot. 0.02 worsens separation somewhat.
stim_window = [0.02, 0.21];  % in seconds: [start, end] For GAVNIK protocol 20250217, 0.03-0.16 is superior  to a longer or shorter window, in terms of mean silhouette score and centroid separation
                            % For my TRAIN protocol 20250203, 0.02-0.21 seems best. But across training days an uptick in spiking develops in the greyscreen period following stimulus offset
%nprobe = 1;
%base_folder='V:\Ellie\DATA\SUBJECTS';
cd('V:\Ellie\DATA\SUBJECTS\M00013\analysis\20250204\TRAIN') % 4/5 files will be saved here in the cd


for nsession = 5 %5/5 row number of recording date in "experiment_info" 
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    % load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));
    % SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT},option);
    % % find right date number based on all experiment dates of the subject
    % iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        subject_number = session_info(n).probe(1).SUBJECT;
        DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters = clusters_ks4;
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));

        ordered_oris = unique(Task_info.stim_orientation, 'stable'); 
        %Task_info.stim_onset
        
       
        params = create_cluster_selection_params('sorting_option','ellie');
        %params.orientation_tuned = ...
        
        switch depth_for_analysis
            case 'L4' 
                depth_range = L4_depth_range;
            case 'V1'
                depth_range = V1_depth_range;
            case 'CA1' 
                depth_range = CA1_depth_range;
            case 'Sub_CA1'
                depth_range = Sub_CA1_depth_range;
        end

        

        if (contains(Stimulus_type, 'GAVNIK_ABCD')) || (contains(Stimulus_type, 'TRAIN'))
            for nprobe = 1:length(clusters)
                selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters, which pass the set parameters
                
                depth_selected_clusters = selected_clusters; % initialize with same structure
                for np = 1:length(clusters)
                    sc = selected_clusters(np);
                    peak_depths = clusters(np).peak_depth(sc.cluster_id); % get depths of selected clusters

                    % Find cluster_ids within selected depth range
                    keep_mask = peak_depths >= min(depth_range) & peak_depths <= max(depth_range);
                    depth_cluster_ids = sc.cluster_id(keep_mask);

                    % Keep only spike times and IDs for clusters at selected depth
                    depth_selected_clusters(np).cluster_id = depth_cluster_ids;
                    depth_selected_clusters(np).spike_times = sc.spike_times(ismember(sc.spike_id, depth_cluster_ids));
                    depth_selected_clusters(np).spike_id = sc.spike_id(ismember(sc.spike_id, depth_cluster_ids));
                end
                cluster_id = depth_selected_clusters(nprobe).cluster_id; % cluster_ids of units which pass the set parameters [NB these are one count higher than per zero-based pythonic SI output cluster IDs...]
                %If you want to use manually chosen clusters within the depth range:
                %cluster_id = [435 438 439 440 445 458 460 464 469 474 471 483 482 484 488 490 493 494 496 497 498 500 502 504 508 510 511 513 515 517 518 522 520 525 526]; % manually selected (visually responsive) for 20250217
                

                % ------------------- REMOVE ZERO-SPIKING NEURONS -------------------
                spike_counts_per_cluster = zeros(size(cluster_id));
                for i = 1:length(cluster_id)
                    spike_counts_per_cluster(i) = sum(depth_selected_clusters(nprobe).spike_id == cluster_id(i));
                end
                keep_idx = spike_counts_per_cluster > 0;
                cluster_id = cluster_id(keep_idx); 
                % Update cluster_id list and corresponding spike_times and spike_ids
                kept_ids = cluster_id;
                is_kept = ismember(depth_selected_clusters(nprobe).spike_id, kept_ids);
                depth_selected_clusters(nprobe).cluster_id = kept_ids; 
                depth_selected_clusters(nprobe).spike_times = depth_selected_clusters(nprobe).spike_times(is_kept); 
                depth_selected_clusters(nprobe).spike_id = depth_selected_clusters(nprobe).spike_id(is_kept); 

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
                                % Get spike counts per trial (1 row per trial)
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
                
                %%%% Now run PCA
                if contains(temporal_structure, 'timebinned spikes')
                    [num_trials, num_bins, num_neurons] = size(z_trial_responses);
                    z_trial_reshaped = reshape(z_trial_responses, num_trials, num_bins * num_neurons);
                    [coeff_z, score_z, latent_z] = pca(z_trial_reshaped);
                else
                    [coeff_z, score_z, latent_z] = pca(z_trial_responses);
                end
                
                %%%% Make some plots to check the PCA analysis makes sense
                if contains(temporal_structure, 'trial spikecounts')
                    figure;
                    imagesc(coeff_z);
                    colorbar;
                    xlabel('Principal Component Number');    
                    xticks(1:length(cluster_id));
                    ylabel('Cluster')
                    title('Weights per coeff matrix')
                end    

                figure; 
                for i = 1:length(ordered_oris)
                    stim_temp = ordered_oris(i);
                    subplot(2,2,i);
                    hold on; 
                    for iPC = 1:13
                        %disp(mean(score_z(trial_labels==stim_temp,iPC)));
                        scatter(iPC,mean(score_z(trial_labels==stim_temp,iPC)));
                    end
                    xlabel('PC');
                    ylabel('Mean PCA Score') %%%% 
                    if contains(Stimulus_type, 'TRAIN')
                        title('Orientation', round(rad2deg(ordered_oris(i))));
                    else 
                        title('Orientation', ordered_oris(i));
                    end
                end
              
                
                % Calculate % variance explained
                var_exp_z = 100 * latent_z / sum(latent_z);
                
                % Cumulative explained variance
                cumulative_var_exp_z = cumsum(var_exp_z);
                % Find elbow point using the "maximum distance to line" method
                x = 1:length(cumulative_var_exp_z);
                y = cumulative_var_exp_z;
                
                % Vector between first and last point
                line_vec = [x(end) - x(1), y(end) - y(1)];
                line_vec = line_vec / norm(line_vec);  % Normalize
                
                % Compute perpendicular distance from each point to the line
                distances = zeros(size(x));
                for i = 1:length(x)
                    point_vec = [x(i) - x(1), y(i) - y(1)];
                    proj_len = dot(point_vec, line_vec);
                    proj_point = [x(1), y(1)] + proj_len * line_vec;
                    distances(i) = norm([x(i), y(i)] - proj_point);
                end
                
                [~, elbow_idx] = max(distances);  % Index of elbow point
                
                %%%% Elbow Point - Plot cumulative explained variance with elbow point i.e. where additional PCs yield diminishing returns
                figure;
                plot(x, y, 'o-', 'LineWidth', 2); hold on;
                plot(elbow_idx, y(elbow_idx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);  % Mark elbow
                xlabel('Number of Principal Components');
                ylabel('Cumulative Variance Explained (%)');
                title('Cumulative Explained Variance (Z-scored PCA)');
                yline(90, '--r', '90% Threshold');
                grid on;
                legend('Cumulative Variance', 'Elbow Point', '90% Threshold');
                
                % Optional: print elbow result
                fprintf('Elbow point detected at PC %d (%.2f%% variance explained).\n', ...
                    elbow_idx, y(elbow_idx));
                
                %%% Centroid calculations based on 3 PCs to characterise separation of trial-type clusters on the 3D scatter plot
                PCs_for_centroids = 3;
                centroid_score_z = score_z(:, 1:PCs_for_centroids);              
                centroids = zeros(length(ordered_oris), PCs_for_centroids); % n neuron dimensions for each orientation. If there are three neurons, there will be x, y, z dimensions for each orientation
                
                for i = 1:length(ordered_oris)
                    idx = trial_labels == ordered_oris(i);
                    centroids(i, :) = mean(centroid_score_z(idx, :), 1);
                end

                centroid_distances = pdist(centroids, 'euclidean');  % pdist gives all pairwise distances
                mean_centroid_distance = mean(centroid_distances);  % mean separation between clusters
                disp(['Mean centroid-to-centroid distance: ', num2str(mean_centroid_distance)]);

                within_cluster_spread = zeros(length(ordered_oris), 1);
                for i = 1:length(ordered_oris)
                    idx = trial_labels == ordered_oris(i);
                    distances_to_centroid = sqrt(sum((centroid_score_z(idx, :) - centroids(i, :)).^2, 2));
                    within_cluster_spread(i) = mean(distances_to_centroid);
                end
                mean_within_cluster_spread = mean(within_cluster_spread);
                disp(['Mean within-cluster spread: ', num2str(mean_within_cluster_spread)]);

                separation_score = mean_centroid_distance / mean_within_cluster_spread;
                disp(['Separation score (higher is better): ', num2str(separation_score)]);
                
                per_class_sep_scores = zeros(length(ordered_oris), 1);

                for i = 1:length(ordered_oris)
                    idx_i = trial_labels == ordered_oris(i);
                    centroid_i = centroids(i, :);
                
                    % Spread of cluster i
                    spread_i = mean(sqrt(sum((centroid_score_z(idx_i, :) - centroid_i).^2, 2)));
                
                    % Mean distance to all *other* centroids
                    other_centroids = centroids(setdiff(1:1:length(ordered_oris), i), :);
                    dist_to_others = sqrt(sum((other_centroids - centroid_i).^2, 2));
                    mean_dist_to_others = mean(dist_to_others);
                
                    % Separation score for this class
                    per_class_sep_scores(i) = mean_dist_to_others / spread_i;
                end

                % Define unique orientations and colors
                colors = lines(length(ordered_oris));  % One distinct color per orientation in order
                % Generate legend labels: A 30°, B 60°, etc.
                letters = 'A':'Z';
                if contains(Stimulus_type, 'TRAIN')
                    ori_labels = arrayfun(@(i, ori) sprintf('%s %.0f°', letters(i), round(rad2deg(ori))), ...
                                      (1:length(ordered_oris))', ordered_oris, 'UniformOutput', false);
                else 
                    ori_labels = arrayfun(@(i, ori) sprintf('%s %.0f°', letters(i), ori), ...
                                      (1:length(ordered_oris))', ordered_oris, 'UniformOutput', false);
                end
                
                
                            
                
                                              
                %%%% ---- Plot PCA (z-scored data) ----
                fig1 = figure;
                cla;  % Clear axes
                hold on;
                for i = 1:length(ordered_oris)
                    idx = trial_labels == ordered_oris(i);
                    % use a jitter as many dots overlap
                    jitter_x = 0.1;  
                    jitter_y = 0.1;   
                    jitter_z = 0.1;  

                    scatter3(score_z(idx,1) + randn(sum(idx),1)*jitter_x, ...
                        score_z(idx,2) + randn(sum(idx),1)*jitter_y, ...
                        score_z(idx,3) + randn(sum(idx),1)*jitter_z, ...
                        36, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.6);
                end
                
                xlabel(sprintf('PC 1 (%.1f%%)', var_exp_z(1)));
                ylabel(sprintf('PC 2 (%.1f%%)', var_exp_z(2)));
                zlabel(sprintf('PC 3 (%.1f%%)', var_exp_z(3)));
                lgd3 = legend(ori_labels);
                title(lgd3, 'Orientation');
                grid on;
                view(3);
                axis equal;
                
                hold on;
                scatter3(centroids(:,1), centroids(:,2), centroids(:,3), 100, colors, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                legend([ori_labels', {'Centroids'}], 'Location', 'southeast');
                % Format the per-class separation scores into lines of text
                per_class_lines = cell(length(per_class_sep_scores), 1);
                for i = 1:length(per_class_sep_scores)
                    per_class_lines{i} = sprintf('%s separation score: %.4f', ori_labels{i}, per_class_sep_scores(i));
                end
                
                % Combine into a single string with newlines
                per_class_text = strjoin(per_class_lines, '\n');
                
                % Combine with your original metrics
                metrics_text = sprintf(['Mean centroid-to-centroid distance: %.4f\n' ...
                                        'Mean within-cluster spread: %.4f\n' ...
                                        'Separation score: %.4f\n\n%s'], ...
                                        mean_centroid_distance, ...
                                        mean_within_cluster_spread, ...
                                        separation_score, ...
                                        per_class_text);
                
                % Add annotation textbox in top-right corner (normalized coordinates)
                % Calculate height based on number of lines
                n_types = length(per_class_sep_scores);
                box_height = 0.15 + 0.03 * n_types;
                
                % Place the box so its top aligns with y = 1 (top of figure)
                x_pos = 0.70;
                y_top = 0.9;
                y_pos = y_top - box_height;  % anchor from top
                
                annotation('textbox', [x_pos, y_pos, 0.28, box_height], ...
                           'String', metrics_text, ...
                           'FitBoxToText', 'on', ...
                           'BackgroundColor', 'white', ...
                           'EdgeColor', 'black', ...
                           'FontSize', 10, ...
                           'VerticalAlignment', 'top', ...
                           'Interpreter', 'none');  % Avoid italicizing if any LaTeX characters

                hold off;
                sgtitle(sprintf('%s - %s: PCA of %s %s (%.2f–%.2fs after stim onset) Z-scored over %s, Day %d', subject_number, Stimulus_type, depth_for_analysis, temporal_structure, stim_window(1), stim_window(2), z_score_period, experiment_info(nsession).date), 'Interpreter', 'none');
                fig_filename1 = sprintf('%s - %s - PCA of %s %s (%.2f–%.2fs after stim onset.fig)', subject_number, Stimulus_type, depth_for_analysis, temporal_structure, stim_window(1), stim_window(2));
                savefig(fig1, fullfile(pwd, fig_filename1));
                        
                

                %%%% Display the top 5 clusters contributing to PC1 
                if contains(temporal_structure, 'timebinned spikes')
                    num_timebins = round(diff(stim_window)/psthBinSize);
                    num_clusters = length(cluster_id);
                
                    % Reshape loadings to [timebins x clusters]
                    pc1_matrix = abs(reshape(coeff_z(:,1), num_timebins, num_clusters)); 
                    
                    % Aggregate across timebins (e.g., mean or max)
                    cluster_loadings = mean(pc1_matrix, 1); % one value per cluster
                
                else
                    cluster_loadings = abs(coeff_z(:,1)); % already one value per cluster
                end
                
                % Sort and display
                [sorted_loadings, sort_idx] = sort(cluster_loadings, 'descend');
                disp('Top contributing clusters to PC1:');
                for i = 1:min(5, length(cluster_id))
                    fprintf('Cluster %d (ID: %d): loading = %.3f\n', ...
                        i, cluster_id(sort_idx(i)), sorted_loadings(i));
                end


                %%%%% Create a figure for correlation matrices
                if contains(temporal_structure, 'trial spikecounts') % simpler case: trials x neurons
                    % Preallocate
                    stored_corrmats = cell(length(ordered_oris), 1);
                    stored_cluster_ids = cell(length(ordered_oris), 1);
                    stored_cluster_orders = cell(length(ordered_oris), 1);
                    stored_valid_cluster_id = cell(length(ordered_oris), 1);
                    all_corr_values = [];                    
                
                    for i = 1:length(ordered_oris)
                        ori = ordered_oris(i);
                        trial_mask = trial_labels == ori;
                        responses_this_ori = z_trial_responses(trial_mask, :);

                        %exclude neurons with extremely low variance.
                        neuron_var = var(responses_this_ori);  % variance across trials for each neuron
                        threshold_var = 1e-30;
                        valid_neurons = neuron_var > threshold_var;
                        fprintf('Excluding %d neurons with very low variance\n', sum(~valid_neurons));
                        valid_cluster_id = cluster_id(valid_neurons);
                        
                        responses_this_ori = responses_this_ori(:, valid_neurons);
                        valid_cluster_id_this_ori = cluster_id(valid_neurons);  % 
                        
                        corrmat = corr(responses_this_ori);
                        corrmat(1:size(corrmat,1)+1:end) = 1;
                        off_diagonal = corrmat(~eye(size(corrmat)));  % exclude autocorrelations
                        all_corr_values = [all_corr_values; off_diagonal(:)];
                        
                        % Clustering
                        Y = 1 - corrmat;
                        D = squareform(Y, 'tovector');
                        Z = linkage(D, 'average');
                        cluster_order = optimalleaforder(Z, D);
                    
                        % Store
                        stored_corrmats{i} = corrmat;
                        stored_valid_cluster_id{i} = valid_cluster_id;
                        stored_cluster_orders{i} = cluster_order;

                    end
                    % set global scale
                    vmin = min(all_corr_values);
                    vmax = max(all_corr_values);
                
                    fig2 = figure;
                    t = tiledlayout(ceil(sqrt(length(ordered_oris))), ceil(sqrt(length(ordered_oris))), ...
                            'TileSpacing', 'compact', 'Padding', 'compact');

                    for i = 1:length(ordered_oris)
                        nexttile;
                        corrmat = stored_corrmats{i};
                        valid_cluster_ids = stored_valid_cluster_id{i};
                        cluster_order = stored_cluster_orders{i};
                        
                        corrmat_reordered = corrmat(cluster_order, cluster_order);
                        corrmat_reordered(1:size(corrmat_reordered,1)+1:end) = NaN;  % set diagonal to NaN for plotting only
 
                        ori = ordered_oris(i);
                        
                                               
                        % Plot
                        imagesc(corrmat_reordered, 'AlphaData', ~isnan(corrmat_reordered), [vmin vmax]);
                        axis square;
                        colormap('parula');
                        cb = colorbar;
                        ylabel(cb, 'Pearson Correlation');
                        if contains(Stimulus_type, 'TRAIN')
                            title([letters(i), ': orientation ', num2str(round(rad2deg(ori))), '°']);
                        else 
                            title([letters(i), ': orientation ', num2str(ori), '°']);
                        end    
                        xticks(1:length(cluster_order));
                        yticks(1:length(cluster_order));
                        xticklabels(string(valid_cluster_ids(cluster_order)));
                        yticklabels(string(valid_cluster_ids(cluster_order)));
                        xtickangle(90);
                        set(gca, 'FontSize', 8);
                    end
                
                    sgtitle(sprintf('%s - %s: %s cluster correlation matrices of trial spikecounts (%.2f–%.2fs after stim onset) Z-scored over %s, Day %d', ...
                        subject_number, Stimulus_type, depth_for_analysis, stim_window(1), stim_window(2), z_score_period, experiment_info(nsession).date), ...
                        'Interpreter', 'none');
                    fig_filename2 = sprintf('%s - %s - %s cluster correlation matrices of %s (%.2f–%.2fs after stim onset).fig', subject_number, Stimulus_type, depth_for_analysis, temporal_structure, stim_window(1), stim_window(2));
                    savefig(fig2, fullfile(pwd, fig_filename2));
                
                elseif contains(temporal_structure, 'timebinned spikes') % trials x timebins x neurons
                    % Preallocate
                    stored_corrmats = cell(length(ordered_oris), 1);
                    stored_cluster_ids = cell(length(ordered_oris), 1);
                    stored_cluster_orders = cell(length(ordered_oris), 1);
                    stored_valid_cluster_id = cell(length(ordered_oris), 1);
                    all_corr_values = [];
                
                    for i = 1:length(ordered_oris)
                        ori = ordered_oris(i);
                        trial_mask = trial_labels == ori;
                        responses_this_ori = z_trial_responses(trial_mask, :, :);
                        
                        [num_trials, num_bins, num_neurons] = size(responses_this_ori);
                        reshaped = permute(responses_this_ori, [3, 1, 2]);  % neurons x trials x timebins
                        reshaped = reshape(reshaped, num_neurons, num_trials * num_bins)';  % (trials*timebins) x neurons
                
                        %exclude neurons with extremely low variance.
                        neuron_var = var(reshaped, 0, 1);  % 1 x neurons
                        threshold_var = 1e-30;
                        valid_neurons = neuron_var > threshold_var;
                        fprintf('Excluding %d neurons with very low variance\n', sum(~valid_neurons));
                        valid_cluster_id = cluster_id(valid_neurons);
                        % Apply filter to original 3D data
                        responses_this_ori = responses_this_ori(:, :, valid_neurons);
                        valid_cluster_id_this_ori = cluster_id(valid_neurons);    
                        stored_valid_cluster_id{i} = valid_cluster_id_this_ori;

                        % Recalculate reshaped after excluding invalid neurons
                        reshaped = permute(responses_this_ori, [3, 1, 2]);  
                        reshaped = reshape(reshaped, sum(valid_neurons), num_trials * num_bins)';

                        corrmat = corr(reshaped);
                        corrmat(1:size(corrmat,1)+1:end) = 1;
                        off_diagonal = corrmat(~eye(size(corrmat)));  % exclude autocorrelations
                        all_corr_values = [all_corr_values; off_diagonal(:)];
                        
                        % Hierarchical clustering
                        Y = 1 - corrmat;  % convert to distance matrix
                        D = squareform(Y, 'tovector');  % get condensed distance vector
                        Z = linkage(D, 'average');
                        cluster_order = optimalleaforder(Z, D);
                        % Store
                        stored_corrmats{i} = corrmat;
                        stored_cluster_orders{i} = cluster_order;
                    end
                
                    vmin = min(all_corr_values);
                    vmax = max(all_corr_values);
                
                    fig2 = figure;
                    t = tiledlayout(ceil(sqrt(length(ordered_oris))), ceil(sqrt(length(ordered_oris))), ...
                    'TileSpacing', 'compact', 'Padding', 'compact');

                    for i = 1:length(ordered_oris)
                        nexttile;
                        corrmat = stored_corrmats{i};
                        valid_cluster_ids = stored_valid_cluster_id{i}; 
                        cluster_order = stored_cluster_orders{i};
                        
                        corrmat_reordered = corrmat(cluster_order, cluster_order);
                        corrmat_reordered(1:size(corrmat_reordered,1)+1:end) = NaN;  % hide diagonal in plot
                        ori = ordered_oris(i);

                        
                        % Plot
                        imagesc(corrmat_reordered, 'AlphaData', ~isnan(corrmat_reordered), [vmin vmax]);
                        axis square;
                        colormap('parula');
                        cb = colorbar;
                        ylabel(cb, 'Pearson Correlation');
                        if contains(Stimulus_type, 'TRAIN')
                            title([letters(i), ': orientation ', num2str(round(rad2deg(ori))), '°']);
                        else 
                            title([letters(i), ': orientation ', num2str(ori), '°']);
                        end
                        xticks(1:length(cluster_order));
                        yticks(1:length(cluster_order));
                        xticklabels(string(valid_cluster_ids(cluster_order)));
                        yticklabels(string(valid_cluster_ids(cluster_order)));
                        xtickangle(90);
                        set(gca, 'FontSize', 8);
                    end
                
                    sgtitle(sprintf('%s - %s: %s cluster correlation matrices of spiking across time (%.2f–%.2fs after stim onset) Z-scored over %s, Day %d', ...
                        subject_number, Stimulus_type, depth_for_analysis, stim_window(1), stim_window(2), z_score_period, experiment_info(nsession).date), ...
                        'Interpreter', 'none');
                    fig_filename2 = sprintf('%s - %s - %s cluster correlation matrices of %s (%.2f–%.2fs after stim onset).fig', subject_number, Stimulus_type, depth_for_analysis, temporal_structure, stim_window(1), stim_window(2));
                    savefig(fig2, fullfile(pwd, fig_filename2));
                end
                
                %%% Make a second correlation matrix plot in which units are clustered for B, C, D according to their hierarchical clustering for A
                if contains(temporal_structure, 'timebinned spikes')

                    first_ori = ordered_oris(1);
                    trial_mask = trial_labels == first_ori;
                    responses_first_ori = z_trial_responses(trial_mask, :, :);
                    
                    [num_trials, num_bins, num_neurons] = size(responses_first_ori);
                    reshaped = permute(responses_first_ori, [3, 1, 2]);  % neurons x trials x timebins
                    reshaped = reshape(reshaped, num_neurons, num_trials * num_bins)';  % (trials*timebins) x neurons
                    
                    %exclude neurons with extremely low variance.
                    neuron_var = var(reshaped, 0, 1);  % 1 x neurons
                    threshold_var = 1e-30;
                    valid_neurons = neuron_var > threshold_var;
                    fprintf('Excluding %d neurons with very low variance\n', sum(~valid_neurons));
                    valid_cluster_id = cluster_id(valid_neurons);
                    % Apply filter to original 3D data
                    responses_first_ori = responses_first_ori(:, :, valid_neurons);
                    [num_trials, num_bins, num_neurons] = size(responses_first_ori);
                    reshaped = permute(responses_first_ori, [3, 1, 2]);
                    reshaped = reshape(reshaped, num_neurons, num_trials * num_bins)';

                    corrmat_first = corr(reshaped);
                    corrmat_first(1:size(corrmat_first,1)+1:end) = 1;  
                    Y_first = 1 - corrmat_first;
                    D_first = squareform(Y_first, 'tovector');
                    Z_first = linkage(D_first, 'average');
                    fixed_cluster_order = optimalleaforder(Z_first, D_first);
                    
                    % Precompute color scale limits (again) for consistency
                    all_corr_values_fixed = [];
                    
                    for i = 1:length(ordered_oris)
                        ori = ordered_oris(i);
                        trial_mask = trial_labels == ori;
                        responses_this_ori = z_trial_responses(trial_mask, :, :);
                        %use same valid_neurons as for ori A
                        responses_this_ori = responses_this_ori(:, :, valid_neurons);  % use fixed mask from ori A

                        [num_trials, num_bins, num_neurons] = size(responses_this_ori);
                        reshaped = permute(responses_this_ori, [3, 1, 2]);
                        reshaped = reshape(reshaped, num_neurons, num_trials * num_bins)';                                         

                        corrmat = corr(reshaped);
                        off_diagonal = corrmat(~eye(size(corrmat)));
                        all_corr_values_fixed = [all_corr_values_fixed; off_diagonal(:)];
                    end
                    
                    vmin_fixed = min(all_corr_values_fixed);
                    vmax_fixed = max(all_corr_values_fixed);
                    
                    % Now plot with fixed clustering using fig5
                    fig5 = figure;
                    t_fixed = tiledlayout(ceil(sqrt(length(ordered_oris))), ceil(sqrt(length(ordered_oris))), ...
                        'TileSpacing', 'compact', 'Padding', 'compact');
                    
                    for i = 1:length(ordered_oris)
                        nexttile;
                        ori = ordered_oris(i);
                        trial_mask = trial_labels == ori;
                        responses_this_ori = z_trial_responses(trial_mask, :, :);
                        
                        % Apply fixed neuron filter before reshaping
                        responses_this_ori = responses_this_ori(:, :, valid_neurons);

                        [num_trials, num_bins, num_neurons] = size(responses_this_ori);
                        reshaped = permute(responses_this_ori, [3, 1, 2]);
                        reshaped = reshape(reshaped, num_neurons, num_trials * num_bins)';
                    
                        corrmat = corr(reshaped);
                        corrmat(1:size(corrmat,1)+1:end) = NaN;  % mask diagonal just for plot
                    
                        % Apply fixed clustering
                        corrmat_reordered = corrmat(fixed_cluster_order, fixed_cluster_order);
                    
                        imagesc(corrmat_reordered, 'AlphaData', ~isnan(corrmat_reordered), [vmin_fixed vmax_fixed]);
                        axis square;
                        colormap('parula');
                        cb = colorbar;
                        ylabel(cb, 'Pearson Correlation');
                        if contains(Stimulus_type, 'TRAIN')
                            title([letters(i), ': orientation ', num2str(round(rad2deg(ori))), '°']);
                        else 
                            title([letters(i), ': orientation ', num2str(ori), '°']);
                        end
    
                        xticks(1:length(fixed_cluster_order));
                        yticks(1:length(fixed_cluster_order));
                        xticklabels(string(valid_cluster_id(fixed_cluster_order)));
                        yticklabels(string(valid_cluster_id(fixed_cluster_order)));
                        xtickangle(90);
                        set(gca, 'FontSize', 8);
                    end
                    
                    sgtitle(sprintf('%s - %s: %s cluster correlation matrices of spiking across time (%.2f–%.2fs after stim onset) Z-scored over %s [fixed clustering from ori A] Day %d', ...
                            subject_number, Stimulus_type, depth_for_analysis, stim_window(1), stim_window(2), z_score_period, experiment_info(nsession).date), ...
                            'Interpreter', 'none');
                    fig_filename5 = sprintf('%s - %s - %s cluster correlation matrices clustered Awise of %s (%.2f–%.2fs after stim onset).fig', subject_number, Stimulus_type, depth_for_analysis, temporal_structure, stim_window(1), stim_window(2));
                        savefig(fig5, fullfile(pwd, fig_filename5));

                elseif contains(temporal_structure, 'trial spikecounts')
                    first_ori = ordered_oris(1);
                    trial_mask = trial_labels == first_ori;
                    responses_first_ori = z_trial_responses(trial_mask, :);
                
                    %exclude neurons with extremely low variance.
                    neuron_var = var(responses_first_ori);  % variance across trials for each neuron
                    threshold_var = 1e-30;
                    valid_neurons = neuron_var > threshold_var;
                    fprintf('Excluding %d neurons with very low variance\n', sum(~valid_neurons));
                    valid_cluster_id = cluster_id(valid_neurons);
                    responses_first_ori = responses_first_ori(:, valid_neurons);
                    
                    corrmat_first = corr(responses_first_ori);
                    corrmat_first(1:size(corrmat_first,1)+1:end) = 1;  
                    Y_first = 1 - corrmat_first;
                    D_first = squareform(Y_first, 'tovector');
                    Z_first = linkage(D_first, 'average');
                    fixed_cluster_order = optimalleaforder(Z_first, D_first);
                
                    % Precompute color scale limits
                    all_corr_values_fixed = [];
                
                    for i = 1:length(ordered_oris)
                        ori = ordered_oris(i);
                        trial_mask = trial_labels == ori;
                        responses_this_ori = z_trial_responses(trial_mask, :);
                        
                        %use same valid_neurons as for ori A
                        responses_this_ori = responses_this_ori(:, valid_neurons);  % use fixed mask from ori A

                        corrmat = corr(responses_this_ori);
                        off_diagonal = corrmat(~eye(size(corrmat)));
                        all_corr_values_fixed = [all_corr_values_fixed; off_diagonal(:)];
                    end
                
                    vmin_fixed = min(all_corr_values_fixed);
                    vmax_fixed = max(all_corr_values_fixed);
                
                    fig5 = figure;
                    t_fixed = tiledlayout(ceil(sqrt(length(ordered_oris))), ceil(sqrt(length(ordered_oris))), ...
                        'TileSpacing', 'compact', 'Padding', 'compact');
                
                    for i = 1:length(ordered_oris)
                        nexttile;
                        ori = ordered_oris(i);
                        trial_mask = trial_labels == ori;
                        responses_this_ori = z_trial_responses(trial_mask, :);
                
                        responses_this_ori = responses_this_ori(:, valid_neurons); 

                        corrmat = corr(responses_this_ori);
                        corrmat(1:size(corrmat,1)+1:end) = NaN;
                
                        corrmat_reordered = corrmat(fixed_cluster_order, fixed_cluster_order);
                
                        imagesc(corrmat_reordered, 'AlphaData', ~isnan(corrmat_reordered), [vmin_fixed vmax_fixed]);
                        axis square;
                        colormap('parula');
                        cb = colorbar;
                        ylabel(cb, 'Pearson Correlation');
                        if contains(Stimulus_type, 'TRAIN')
                            title([letters(i), ': orientation ', num2str(round(rad2deg(ori))), '°']);
                        else 
                            title([letters(i), ': orientation ', num2str(ori), '°']);
                        end
                        xticks(1:length(fixed_cluster_order));
                        yticks(1:length(fixed_cluster_order));
                        xticklabels(string(valid_cluster_id(fixed_cluster_order)));
                        yticklabels(string(valid_cluster_id(fixed_cluster_order)));
                        xtickangle(90);
                        set(gca, 'FontSize', 8);
                    end
                
                    sgtitle(sprintf('%s - %s: %s cluster correlation matrices of trial spikecounts (%.2f–%.2fs after stim onset) Z-scored over %s [fixed clustering from ori A] Day %d', ...
                        subject_number, Stimulus_type, depth_for_analysis, stim_window(1), stim_window(2), z_score_period, experiment_info(nsession).date), ...
                        'Interpreter', 'none');
                    fig_filename5 = sprintf('%s - %s - %s cluster correlation matrices clustered Awise of %s (%.2f–%.2fs after stim onset).fig', subject_number, Stimulus_type, depth_for_analysis, temporal_structure, stim_window(1), stim_window(2));
                    savefig(fig5, fullfile(pwd, fig_filename5));
                end

                %%% SILHOUETTE SCORING
                % Test how many PCs you need before the silhouette score plateaus
                maxPCs = length(cluster_id);  % or however many you want to test
                % Preallocate silhouette score arrays
                sil_overall = zeros(maxPCs, 1);
                sil_by_type = zeros(maxPCs, length(ordered_oris));  % one column per orientation
                
                % Loop through PC counts
                for k = 1:maxPCs
                    reduced_data = score_z(:, 1:k);
                    sil_vals = silhouette(reduced_data, trial_labels);  % one value per trial
                    
                    % Overall silhouette score
                    sil_overall(k) = mean(sil_vals);
                    
                    % Trial-type-specific silhouette scores
                    for t = 1:length(ordered_oris)
                        label = ordered_oris(t);
                        idx = trial_labels == label;
                        sil_by_type(k, t) = mean(sil_vals(idx));
                    end
                end
                
                % ----- Plot -----
                fig3 = figure;
                hold on;
                
                % Plot overall silhouette score in black
                plot(1:maxPCs, sil_overall, 'k-o', 'LineWidth', 2);
                
                % Plot each trial type
                for t = 1:length(ordered_oris)
                    plot(1:maxPCs, sil_by_type(:, t), '-', ...
                         'Color', colors(t,:), ...
                         'LineWidth', 2, ...
                         'DisplayName', ori_labels{t});  % For legend
                end
                
                xlabel('Number of Principal Components');
                ylabel('Mean Silhouette Score');
                sgtitle(sprintf('%s - %s - %s %s (%.2f–%.2fs after stim onset) Silhouette Clustering Quality vs. Number of PCs (by Trial Type) Day %d', subject_number, Stimulus_type, depth_for_analysis, temporal_structure, stim_window(1), stim_window(2), experiment_info(nsession).date), 'Interpreter', 'none');
                legend([{'Overall'}, ori_labels'], 'Location', 'best');
                grid on;
                fig3_filename = sprintf('%s - %s - %s %s (%.2f–%.2fs after stim onset) Silhouette score progression with increasing PCs, coloured by trial-type.fig', subject_number, Stimulus_type, depth_for_analysis, temporal_structure, stim_window(1), stim_window(2));
                fig3_save_path = fullfile(pwd, fig3_filename);
                savefig(fig3, fig3_save_path);

                                                
                % ---- Silhouette Plot Colored by Trial Type ----
                % Compute silhouette values
                PCs_for_silhouette = 6;  % based on elbow point and silhouette analysis
                silhouette_score_z = score_z(:, 1:PCs_for_silhouette);  % retain only first X principal components
                
                % Compute silhouette values for all points
                silhouette_vals = silhouette(silhouette_score_z, trial_labels);
                
                silh_vals = silhouette(silhouette_score_z, trial_labels);  % Prevent auto-plotting
                
                % Create new figure
                fig4 = figure;
                hold on;
                
                % Sort by trial label for grouping
                % Create an index map for custom ordering
                [~, label_order] = ismember(trial_labels, ordered_oris);
                [~, sorted_idx] = sort(label_order);
                sorted_labels = trial_labels(sorted_idx);
                disp('Silhouette plot trial order:');
                disp(ordered_oris);    

                sorted_silh_vals = silh_vals(sorted_idx);
                
                % Plot each bar colored by trial type
                nTrials = length(silh_vals);
                y_pos = 1;
                for i = 1:nTrials
                    trial_label = sorted_labels(i);
                    color_idx = find(ordered_oris == trial_label);
                    barh(y_pos, sorted_silh_vals(i), 1, 'FaceColor', colors(color_idx,:), ...
                         'EdgeColor', 'none');
                    y_pos = y_pos + 1;
                end
                
                % Add color-matched legend entries
                dummy_handles = gobjects(length(ordered_oris), 1);
                for i = 1:length(ordered_oris)
                    dummy_handles(i) = barh(nTrials + 10, 0, 'FaceColor', colors(i,:), 'Visible', 'off');
                end
                legend(dummy_handles, ori_labels, 'Location', 'best');
                
                % Compute and show mean silhouette score
                mean_silhouette = mean(silh_vals);
                disp(['Mean silhouette score: ', num2str(mean_silhouette)]);
                
                % Add vertical line for mean silhouette score
                xline(mean_silhouette, '--k', ...
                    sprintf('Mean = %.2f', mean_silhouette), ...
                    'LabelHorizontalAlignment', 'left', ...
                    'LabelVerticalAlignment', 'middle', 'HandleVisibility', 'off');
                

                % Annotate
                xlabel('Silhouette Value');
                ylabel('Trial index (sorted by orientation)');
                title(sprintf('%s - %s - Silhouette plot using first %d PCs of %s Z-scored %s (%.2f–%.2fs after stim onset), coloured by trial-type', subject_number, Stimulus_type, PCs_for_silhouette, depth_for_analysis, temporal_structure, stim_window(1), stim_window(2)), 'Interpreter', 'none');
                xlim([-1 1]);
                grid on;
                
                fig4_filename = sprintf('%s - %s - Silhouette plot using first %d PCs of %s Z-scored %s (%.2f–%.2fs after stim onset), coloured by trial-type.fig', subject_number, Stimulus_type, PCs_for_silhouette, depth_for_analysis, temporal_structure, stim_window(1), stim_window(2));
                fig4_save_path = fullfile(pwd, fig4_filename);
                savefig(fig4, fig4_save_path);


                
                % for i = 1:length(ordered_oris)
                %     ori = ordered_oris(i);
                % 
                %     % Find trials where this orientation was shown
                %     trial_mask = trial_labels == ori;
                % 
                %     % Create figure for this orientation
                %     fig_pc_corr = figure('Name', sprintf('PC correlations for orientation %d°', ori));
                % 
                %     for pc = 1:4
                %         % Reconstruct trial x neuron activity using only this PC
                %         proj = score_z(trial_mask, pc) * coeff_z(:, pc)';  % (trials x neurons)
                % 
                %         % Compute correlation matrix across neurons
                %         corrmat = corr(proj);
                %         corrmat(1:size(corrmat,1)+1:end) = NaN;  % remove autocorrelations
                % 
                %         % Plot
                %         subplot(2, 2, pc);
                %         imagesc(corrmat, 'AlphaData', ~isnan(corrmat));
                %         axis square;
                %         colormap('parula');
                %         cb = colorbar;
                %         ylabel(cb, 'Pearson Correlation');
                %         title(sprintf('PC%d only', pc));
                %         xticks(1:length(cluster_id));
                %         yticks(1:length(cluster_id));
                %         xticklabels(string(cluster_id));
                %         yticklabels(string(cluster_id));
                %         xtickangle(90);
                %     end
                % 
                %     sgtitle(sprintf('%s - %s: %s cluster PC-wise correlations of spiking (%.2f–%.2fs after stim onset) Z-scored over %s, Orientation %d°, Day %d', subject_number, Stimulus_type, depth_for_analysis, stim_window(1), stim_window(2), z_score_period, ori, experiment_info(nsession).date), 'Interpreter', 'none');
                %     % Save figure if desired
                %     fig_filename = sprintf('%s - %s - PC correlations, Ori %d', subject_number, Stimulus_type, ori);
                %     savefig(fig_pc_corr, fullfile(pwd, fig_filename));
                % end

            end
        end
    end
end    




        
       
%clusters_probe = select_clusters(clusters(nprobe),params); %only look at good clusters

