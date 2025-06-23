%%% For assessing the coding of ABCD orientations in the spiking responses of 3 tuned putative single units

addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

%% setting metrics to screen good clusters
clear all
% Choose your putative single units of interest
tuned_units = [490 508 513]; % cluster IDs of three neurons whose firing rates together best define A, B, C and D distinctly, manually chosen

SUBJECTS = {'M00013'};

option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);

%%% 3/5
Stimulus_type = 'GAVNIK_ABCD'; % OMIT
%plot_choice = 'aggregate'; % auto-curated 'single_units' or in 'aggregate' or uncurated 'MUA'; MUA includes all clusters from kilosort, unfiltered
z_score_period = 'entire_session'; % z score either over 'entire_session' or 'first30secs' or 'none' (for every stimulus recording
% session from 20250205 onward, I presented grey screen to the mouse for at least 30s before starting the stimulus. 'none' may be useful 
% to try for the aggregate TRAIN case across days)
%nprobe = 1;
%base_folder='V:\Ellie\DATA\SUBJECTS';
cd('V:\Ellie\DATA\SUBJECTS\M00013\analysis\20250217\GAVNIK_ABCD') % 4/5 files will be saved here in the cd


for nsession = 11 %5/5 row number of recording date in "experiment_info" 
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
       
        psthBinSize = 0.01; % but use 1ms for raster plots
        stim_window = [0.03, 0.18];  % in seconds: [start, end]

       if (contains(Stimulus_type, 'GAVNIK_ABCD'))
            cluster_id = tuned_units;
            % For each probe (if multiple), find the spike times and IDs for only the tuned_units
            for nprobe = 1:length(clusters) % clusters has a length of 1 if there is only one probe
                
                spike_mask = ismember(clusters(nprobe).spike_id, tuned_units);
                clusters(nprobe).spike_times = clusters(nprobe).spike_times(spike_mask);
                clusters(nprobe).spike_id = clusters(nprobe).spike_id(spike_mask);
                clusters(nprobe).cluster_id = tuned_units;
            
                baseline_window = [0 30]; % in seconds - first 30s of recording is grey screen - can use for z-scoring
    
                % Define time_edges depending on z_score_period
                if contains(z_score_period, 'entire_session')
                    time_edges = 0:psthBinSize:max(clusters(nprobe).spike_times); %time bin edges from 0 to the max time in steps of ptshBinSize
                elseif contains(z_score_period, 'first30secs')
                    time_edges = baseline_window(1):psthBinSize:baseline_window(2);
                end
                
                              
                % Count total number of stimulus presentations
                total_presentations = length(Task_info.stim_onset);
                                
                z_trial_responses = zeros(total_presentations, length(cluster_id));
                trial_labels = zeros(total_presentations, 1);
                trial_counter = 0;
                                       
                for ori = 1:length(ordered_oris) % loop through each unique stimulus orientation shown to the mouse during this task
                    onset_times = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));
                    n_trials = length(onset_times);
                     
                    for trial_idx = 1:n_trials
                        trial_counter = trial_counter + 1;
                        trial_labels(trial_counter) = ordered_oris(ori);   
                        trial_vector = zeros(1, length(cluster_id)); % one value per neuron per trial (x, y and z dimensions)   
    
                        for nCluster = 1:length(cluster_id) % loop through each chosen cluster
                            spike_times_this_cluster = clusters(nprobe).spike_times(clusters(nprobe).spike_id == cluster_id(nCluster)); % extract the spike times for this cluster    
                            % Compute appropriate histogram counts for z-scoring
                            if contains(z_score_period, 'entire_session')
                                zscore_counts = histcounts(spike_times_this_cluster, time_edges);
                            elseif contains(z_score_period, 'first30secs')
                                baseline_spikes = spike_times_this_cluster(spike_times_this_cluster >= baseline_window(1) & spike_times_this_cluster <= baseline_window(2));
                                zscore_counts = histcounts(baseline_spikes, time_edges);
                            end
                                % Get spike counts per trial (1 row per trial)
                            onset = onset_times(trial_idx);
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  onset, stim_window, psthBinSize); % for this orientation, get the binnedArray (matrix of spike counts per timebin), for 30-150ms following stimulus onset
                            
                            % Now z-score. With z-scoring, all units contribute equally to the PCA, regardless of firing rate.
                            z_binnedArray = (binnedArray - mean(zscore_counts))./std(zscore_counts);

                            trial_vector(nCluster) = mean(z_binnedArray, 2); % mean binned spike count during this trial window for this cluster
                            
                        end
                        % Fill in the current row of trial_responses (i.e., current trial) with the activity from all selected clusters (neurons).
                        z_trial_responses(trial_counter, :) = trial_vector; % trial_counter is a running index that tracks which trial we're on (regardless of orientation), : refers to all columns
                            
                    end
                end               
                
    
                % Centroid calculations to characterise separation of trial-type clusters
                centroids = zeros(length(ordered_oris), size(z_trial_responses, 2)); % n neuron dimensions for each orientation. If there are three neurons, there will be x, y, z dimensions for each orientation
                for i = 1:length(ordered_oris)
                    idx = trial_labels == ordered_oris(i);
                    centroids(i, :) = mean(z_trial_responses(idx, :), 1);
                end

                centroid_distances = pdist(centroids, 'euclidean');  % vector of pairwise distances
                mean_centroid_distance = mean(centroid_distances);  % mean separation between clusters
                disp(['Mean centroid-to-centroid distance: ', num2str(mean_centroid_distance)]);

                within_cluster_spread = zeros(length(ordered_oris), 1);
                for i = 1:length(ordered_oris)
                    idx = trial_labels == ordered_oris(i);
                    distances_to_centroid = sqrt(sum((z_trial_responses(idx, :) - centroids(i, :)).^2, 2));
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
                    spread_i = mean(sqrt(sum((z_trial_responses(idx_i, :) - centroid_i).^2, 2)));
                
                    % Mean distance to all *other* centroids
                    other_centroids = centroids(setdiff(1:end, i), :);
                    dist_to_others = sqrt(sum((other_centroids - centroid_i).^2, 2));
                    mean_dist_to_others = mean(dist_to_others);
                
                    % Separation score for this class
                    per_class_sep_scores(i) = mean_dist_to_others / spread_i;
                end
                
                 % Define unique orientations and colors
                colors = lines(length(ordered_oris));  % One distinct color per orientation in order
                
                % Generate legend labels: A 30°, B 60°, etc.
                letters = 'A':'Z';
                ori_labels = arrayfun(@(i, ori) sprintf('%s %.0f°', letters(i), ori), ...
                                      (1:length(ordered_oris))', ordered_oris, 'UniformOutput', false);

                % Display per-class separation scores
                disp('Per-class separation scores:');
                for i = 1:length(ordered_oris)
                    fprintf('%s: %.4f\n', ori_labels{i}, per_class_sep_scores(i));
                end
                
                % ---- Conditional 3D Scatter Plot Using Z-scored Spike Counts from Selected Neurons ---
                if size(z_trial_responses, 2) == 3 
                    fig1 = figure;
                    hold on;
                    for i = 1:length(ordered_oris)
                        idx = trial_labels == ordered_oris(i);
                        % use a jitter as many dots overlap
                        jitter_x = 0.1;  % For z_trial_responses(:,1)
                        jitter_y = 0.015;   % For z_trial_responses(:,2)
                        jitter_z = 0.3;  % For z_trial_responses(:,3)
    
                        scatter3(z_trial_responses(idx,1) + randn(sum(idx),1)*jitter_x, ...
                                 z_trial_responses(idx,2) + randn(sum(idx),1)*jitter_y, ...
                                 z_trial_responses(idx,3) + randn(sum(idx),1)*jitter_z, ...
                                 36, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.6);
                        
                    end
                    title(sprintf('%s - %s - Z-scored spikecount of clusters %d, %d, %d across trials', subject_number, Stimulus_type, tuned_units(1), tuned_units(2), tuned_units(3)), 'Interpreter', 'None');
                    xlabel(sprintf('unit ID %d (x-axis)', tuned_units(1)));
                    ylabel(sprintf('unit ID %d (y-axis)', tuned_units(2)));
                    zlabel(sprintf('unit ID %d (z-axis)', tuned_units(3)));
                    %legend(ori_labels, 'Location', 'best');
                    grid on;
                    view(3);
                    axis equal;
                    hold on;
                    scatter3(centroids(:,1), centroids(:,2), centroids(:,3), 100, colors, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                    legend([ori_labels', {'Centroids'}], 'Location', 'best');
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
                    fig1_filename = sprintf('%s - %s - Z-scored spikecount of clusters %d, %d, %d across trials with trial_type separation scores.fig', subject_number, Stimulus_type, tuned_units(1), tuned_units(2), tuned_units(3));
                    fig1_save_path = fullfile(pwd, fig1_filename);
                    savefig(fig1, fig1_save_path);
                end

                % ---- Silhouette Plot Colored by Trial Type ----
                % Compute silhouette values
                silh_vals = silhouette(z_trial_responses, trial_labels);  % Prevent auto-plotting
                
                % Create new figure
                fig2 = figure;
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
                title(sprintf('%s - %s - Silhouette plot of Z-scored spikecount of clusters %d, %d, %d coloured by trial-type', subject_number, Stimulus_type, tuned_units(1), tuned_units(2), tuned_units(3)), 'Interpreter', 'none');
                xlim([-1 1]);
                grid on;
                
                fig2_filename = sprintf('%s - %s - Silhouette plot of Z-scored spikecount of clusters %d, %d, %d coloured by trial-type.fig', subject_number, Stimulus_type, tuned_units(1), tuned_units(2), tuned_units(3));
                fig2_save_path = fullfile(pwd, fig2_filename);
                savefig(fig2, fig2_save_path);

            end
        end
    end
end    




        
       
%clusters_probe = select_clusters(clusters(nprobe),params); %only look at good clusters

