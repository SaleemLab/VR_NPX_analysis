%%% For analysis of unit spiking in response to visual stimuli

addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

%% setting metrics to screen good clusters
clear all
% Choose your probe depth of interest
session_specific_L4 = {4, 4440:4580; 5, 4440:4580; 6, 4460:4600; 7, 4480:4620; 8, 4500:4640}; % 1/5. um. Set for each SESSION based on CSD +/- 70um
% 4, 4440:4580; 5, 4440:4580; 6, 4460:4600; 7, 4480:4620; 8, 4500:4640; %% M00013 TRAIN
% 11, 4500:4640; 12, 4510:4650; 13, 4510:4650; 14, 4510:4650; 15, 4510:4650; %% M00013 GAVNIK
% 

SUBJECTS = {'M00013'};

params = create_cluster_selection_params('sorting_option','ellie');
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);

%%% 2/5
Stimulus_type = 'TRAIN'; 
plot_choice = 'aggregate'; % curated 'single_units' or in 'aggregate' or uncurated 'MUA'; MUA includes all clusters from kilosort, unfiltered
plot_type = 'FR'; % 'FR' firing rate or 'raster'
neuron_type = 'All'; % for GAVNIK stimuli (coded for so far) set this to 'PYR only' if you want to include only putative pyramidal neurons. distinguish between putative PYR (wide waveform, lower tau rise than SOM), PV (narrow waveform) and SOM (wide waveform, higher tau rise in the ACG i.e. probability of spiking again increases quite slowly)
z_score_period = 'entire_session'; % 'none' = no z scoring or z score either over 'entire_session' or 'first30secs' (for every stimulus recording
% session from 20250205 onward, I presented grey screen to the mouse for at least 30s before starting the stimulus).
cd('V:\Ellie\DATA\SUBJECTS\M00013\analysis') % 3/5 files will be saved here in the cd

sessions_to_plot = [4, 5, 6, 7, 8]; %4/5 row numbers of recording dates in "experiment_info" 4, 5, 6, 7, 8   11, 12, 13, 14, 15
colors = [
    1.0,  0.8,  0.0;  % 1st group - golden yellow (1,0, 1.0, 0.0 for pure yellow)
    0.4,  0.8,  0.0;  % 2nd group - yellow-green (more yellowish)
    0.2,  0.6,  0.4;  % 3rd group - pure green
    0.2,  0.4,  0.6;  % 4th group - greenish blue (not too blue)
    0.0,  0.0,  0.8];   % 5th group - deep blue


% 5/5 Define windows after time zero to calc. the peak and mean FR - Depends on stimulus type              
if contains(Stimulus_type, 'GAVNIK')
    stim_window_starts = 0:0.15:0.45;
    stim_window_ends = 0.15:0.15:0.6;
    grey_window_starts = 0.65; % look from 50ms-150ms after stim D offset
    grey_window_ends = 0.75; % look from 50ms-150ms after stim D offset
else
    stim_window_starts = 0:0.3:0.9;
    stim_window_ends = 0.15:0.3:1.05;
    grey_window_starts = [0.2 0.5 0.8 1.1]; % look from 50ms-150ms after stim offset
    grey_window_ends = [0.3 0.6 0.9 1.2]; % look from 50ms-150ms after stim offset
end

fig1 = figure; % for traces
fig2 = figure; % for bar charts
all_peak_FR_by_stimwindow = zeros(length(sessions_to_plot), length(stim_window_starts));
all_mean_FR_by_stimwindow = zeros(length(sessions_to_plot), length(stim_window_starts));
all_peak_FR_by_greywindow = zeros(length(sessions_to_plot), length(grey_window_starts));
all_mean_FR_by_greywindow = zeros(length(sessions_to_plot), length(grey_window_starts));

sem_peak_FR_by_stimwindow = zeros(length(sessions_to_plot), length(stim_window_starts));
sem_mean_FR_by_stimwindow = zeros(length(sessions_to_plot), length(stim_window_starts));
sem_peak_FR_by_greywindow = zeros(length(sessions_to_plot), length(grey_window_starts));
sem_mean_FR_by_greywindow = zeros(length(sessions_to_plot), length(grey_window_starts));

for nsession = sessions_to_plot 
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
  

    for n = 1:length(session_info) 
        options = session_info(n).probe(1);
        subject_number = session_info(n).probe(1).SUBJECT;
        DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters = clusters_ks4;
        
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));

        all_orientations = unique(Task_info.stim_orientation);
        %Task_info.stim_onset
        
       
        params = create_cluster_selection_params('sorting_option','ellie');
        %params.orientation_tuned = ...
        psthBinSize = 0.01; % but use 1ms for raster plots
        
        idx = find([session_specific_L4{:,1}] == nsession);
        L4_depth_range = session_specific_L4{idx, 2};
        depth_range = L4_depth_range;
                

        if contains(Stimulus_type, 'TRAIN') && contains(plot_choice, 'aggregate') && contains(plot_type, 'FR')        
            
            for nprobe = 1:length(clusters)
                selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters, which pass the set parameters
                
                depth_selected_clusters = selected_clusters; % initialize with same structure
                for np = 1:length(clusters)
                    sc = selected_clusters(np);
                    peak_depths = clusters(np).peak_depth(sc.cluster_id); % get depths of selected clusters

                    % Find cluster_ids within selected depth range
                    keep_mask = peak_depths >= min(depth_range) & peak_depths <= max(depth_range);
                    depth_cluster_ids = sc.cluster_id(keep_mask)

                    % Keep only spike times and IDs for clusters at selected depth
                    depth_selected_clusters(np).cluster_id = depth_cluster_ids;
                    depth_selected_clusters(np).spike_times = sc.spike_times(ismember(sc.spike_id, depth_cluster_ids));
                    depth_selected_clusters(np).spike_id = sc.spike_id(ismember(sc.spike_id, depth_cluster_ids));
                end
                
                baseline_window = [0 30]; % in seconds - first 30s of recording is grey screen - can use for z-scoring

                % Define time_edges depending on z_score_period
                if contains(z_score_period, 'entire_session')
                    time_edges = 0:psthBinSize:max(depth_selected_clusters(nprobe).spike_times); %time bin edges from 0 to the max time in steps of ptshBinSize
                elseif contains(z_score_period, 'first30secs')
                    time_edges = baseline_window(1):psthBinSize:baseline_window(2);
                end
                

                % Combine spike times for all clusters at selected depth (i.e., for curated MUA)
                all_spike_times = [];
                all_spike_ids = [];
                for np = 1:length(depth_selected_clusters)
                    all_spike_times = [all_spike_times; depth_selected_clusters(np).spike_times];
                    all_spike_ids = [all_spike_ids; depth_selected_clusters(np).spike_id];
                end
                
                % Compute appropriate histogram counts for z-scoring
                if contains(z_score_period, 'entire_session')
                    zscore_counts = histcounts(all_spike_times, time_edges);
                elseif contains(z_score_period, 'first30secs')
                    baseline_spikes = all_spike_times(all_spike_times >= baseline_window(1) & all_spike_times <= baseline_window(2));
                    zscore_counts = histcounts(baseline_spikes, time_edges); %histcounts counts the number of datapoints in specified time bins
                end

                ordered_oris = unique(Task_info.stim_orientation, 'stable');          
                ori = 1; % make plots from onset of first stim in sequence
                    
                stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));              
                [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.30 1.8], psthBinSize);    
                
                mean_trace = mean(binnedArray, 1); % average over trials
                figure(fig1);
                hold on;

                if contains(z_score_period, 'none')
                    % Plot raw mean firing rate trace
                    plot(bins, mean_trace, 'Color', colors(find(sessions_to_plot == nsession),:), 'LineWidth', 1.5, 'DisplayName', sprintf('Day %d, L4 %d–%d µm', experiment_info(nsession).date, min(depth_range), max(depth_range)));
                    ylabel('Mean firing rate (Hz)', 'FontSize', 14);
                    ylim ([0 16]);
                else 
                    baseline_mean = mean(zscore_counts);
                    baseline_std = std(zscore_counts);
                    zscored_trials = (binnedArray - baseline_mean) / baseline_std; 
                    z_trace = mean(zscored_trials, 1);  
                    sem_trace = std(zscored_trials, 0, 1) / sqrt(size(zscored_trials, 1)); % SEM
                    
                    % Plot SEM shading
                    fill([bins, fliplr(bins)], ...
                        [z_trace + sem_trace, fliplr(z_trace - sem_trace)], ...
                        colors(find(sessions_to_plot == nsession), :), ...
                        'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');

                    plot(bins, z_trace, 'Color', colors(find(sessions_to_plot == nsession),:), 'LineWidth', 1.5, 'DisplayName', sprintf('Day %d, L4 %d–%d µm', experiment_info(nsession).date, min(depth_range), max(depth_range)));

                    if contains(z_score_period, 'entire_session')
                        ylabel('Z-scored FR (z-scored over entire session)', 'FontSize', 14);
                        ylim([-1 6]);
                    elseif contains(z_score_period, 'first30secs')
                        ylabel('Z-scored FR (z-scored over first 30s baseline)', 'FontSize', 14);
                        ylim([-1 8]);
                    end                    
                end    
                               
                peak_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                mean_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                peak_FR_by_greywindow = zeros(1, length(grey_window_starts)); % preallocate
                mean_FR_by_greywindow = zeros(1, length(grey_window_starts)); % preallocate
                  
                % per trial values for SEM calculation
                peak_FR_trials_by_stimwindow = zeros(size(zscored_trials, 1), length(stim_window_starts)); % preallocate
                mean_FR_trials_by_stimwindow = zeros(size(zscored_trials, 1), length(stim_window_starts)); % preallocate
                peak_FR_trials_by_greywindow = zeros(size(zscored_trials, 1), length(grey_window_starts)); % preallocate
                mean_FR_trials_by_greywindow = zeros(size(zscored_trials, 1), length(grey_window_starts)); % preallocate

                for i = 1:length(stim_window_starts)
                    idx_in_stimwindow = bins >= stim_window_starts(i) & bins < stim_window_ends(i);
                    if contains(z_score_period, 'none')                                         
                        peak_FR_by_stimwindow(i) = max(mean_trace(idx_in_stimwindow));
                        mean_FR_by_stimwindow(i) = mean(mean_trace(idx_in_stimwindow));
                        %for SEM
                        peak_FR_trials_by_stimwindow(:, i) = max(binnedArray(:, idx_in_stimwindow), [], 2);
                        mean_FR_trials_by_stimwindow(:, i) = mean(binnedArray(:, idx_in_stimwindow), 2);
                        
                    else                
                        peak_FR_by_stimwindow(i) = max(z_trace(idx_in_stimwindow));
                        mean_FR_by_stimwindow(i) = mean(z_trace(idx_in_stimwindow));
                        %for SEM
                        mean_FR_trials_by_stimwindow(:, i) = mean(zscored_trials(:, idx_in_stimwindow), 2);
                        peak_FR_trials_by_stimwindow(:, i) = max(zscored_trials(:, idx_in_stimwindow), [], 2);
                    end
                end
                
                for i = 1:length(grey_window_starts)
                    idx_in_greywindow = bins >= grey_window_starts(i) & bins < grey_window_ends(i); 
                    if contains(z_score_period, 'none')                                          
                        peak_FR_by_greywindow(i) = max(mean_trace(idx_in_greywindow));
                        mean_FR_by_greywindow(i) = mean(mean_trace(idx_in_greywindow));
                        %for SEM
                        peak_FR_trials_by_greywindow(:, i) = max(binnedArray(:, idx_in_greywindow), [], 2);
                        mean_FR_trials_by_greywindow(:, i) = mean(binnedArray(:, idx_in_greywindow), 2); 
                    else              
                        peak_FR_by_greywindow(i) = max(z_trace(idx_in_greywindow));
                        mean_FR_by_greywindow(i) = mean(z_trace(idx_in_greywindow));
                        %for SEM
                        peak_FR_trials_by_greywindow(:, i) = max(zscored_trials(:, idx_in_greywindow), [], 2);
                        mean_FR_trials_by_greywindow(:, i) = mean(zscored_trials(:, idx_in_greywindow), 2);
                        
                    end
                end


                session_idx = find(sessions_to_plot == nsession); % Map actual session number to index in preallocated array
                
                all_peak_FR_by_stimwindow(session_idx, :) = peak_FR_by_stimwindow;
                all_mean_FR_by_stimwindow(session_idx, :) = mean_FR_by_stimwindow;
                all_peak_FR_by_greywindow(session_idx, :) = peak_FR_by_greywindow;
                all_mean_FR_by_greywindow(session_idx, :) = mean_FR_by_greywindow;

                sem_peak_FR_by_stimwindow(session_idx, :) = std(peak_FR_trials_by_stimwindow, 0, 1) / sqrt(size(peak_FR_trials_by_stimwindow,1));
                sem_mean_FR_by_stimwindow(session_idx, :) = std(mean_FR_trials_by_stimwindow, 0, 1) / sqrt(size(mean_FR_trials_by_stimwindow,1));
                sem_peak_FR_by_greywindow(session_idx, :) = std(peak_FR_trials_by_greywindow, 0, 1) / sqrt(size(peak_FR_trials_by_greywindow,1));
                sem_mean_FR_by_greywindow(session_idx, :) = std(mean_FR_trials_by_greywindow, 0, 1) / sqrt(size(mean_FR_trials_by_greywindow,1));
                
                xlim([-0.5 2.0]);
                xticks(-0.4:0.2:1.8);
                set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 14);
                xlabel('Time (s) since onset of A', 'FontSize', 14)
                
                legend(flipud(findobj(gca,'-property','DisplayName')), 'Location', 'northeast', 'FontSize', 14);
                hold on;                
            end
                    % Define grey intervals
            grey_intervals = [-0.5 0; 0.15 0.3; 0.45 0.6; 0.75 0.9; 1.05 2];
                
                          
            % Get current y-axis limits for full vertical shading
            yl = ylim;
                
            % Shade each interval
            for i = 1:size(grey_intervals, 1)
                x = [grey_intervals(i,1), grey_intervals(i,2), grey_intervals(i,2), grey_intervals(i,1)];
                y = [yl(1), yl(1), yl(2), yl(2)];
                fill(x, y, [0.7 0.7 0.7], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % grey color with transparency
            end

            xline(0, 'k', (sprintf('A %d%s onset', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
            xline(0.30, 'k', (sprintf('B %d%s onset', round(rad2deg(ordered_oris(2))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
            xline(0.60, 'k', (sprintf('C %d%s onset', round(rad2deg(ordered_oris(3))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
            xline(0.90, 'k', (sprintf('D %d%s onset', round(rad2deg(ordered_oris(4))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                
            sgtitle(sprintf('%s - %s: Aggregate L4 single unit activity across days', subject_number, Stimulus_type), 'Interpreter', 'none');
            filename = sprintf('%s - %s - Multiplot L4 FR traces across %d days.png', subject_number, Stimulus_type, length(sessions_to_plot));
            save_path = fullfile(pwd, filename);
            exportgraphics(fig1, save_path); 

            % Also save as .fig
            fig_filename = sprintf('%s - %s - Multiplot L4 FR traces across %d days.fig', subject_number, Stimulus_type, length(sessions_to_plot));
            fig_save_path = fullfile(pwd, fig_filename);
            savefig(fig1, fig_save_path);
            
            
            % Plot bar chart showing development of peak and mean firing across training days
            figure(fig2); % switch to bar chart figure
            hold on;
                        
            num_sessions = length(sessions_to_plot);
            cumulative_peak_stim = sum(all_peak_FR_by_stimwindow(:, :), 2);
            cumulative_mean_stim = sum(all_mean_FR_by_stimwindow(:, :), 2);
            cumulative_peak_grey = sum(all_peak_FR_by_greywindow(:, :), 2);
            cumulative_mean_grey = sum(all_mean_FR_by_greywindow(:, :), 2);
                      
            % Set bar width and spacing
            num_bars = 4;
            bar_width = 0.15;
            x_offsets = ((1:num_bars) - (num_bars+1)/2); % e.g., [-1.5, -0.5, 0.5, 1.5]
            x = 1:num_sessions;
            
            for i = 1:num_sessions
                xpos = x(i);
            
                % Stimulus peak
                b1 = bar(xpos + x_offsets(1)*bar_width, cumulative_peak_stim(i), bar_width, 'FaceColor', colors(i,:));
                % Stimulus mean
                b2 = bar(xpos + x_offsets(2)*bar_width, cumulative_mean_stim(i), bar_width, 'FaceColor', colors(i,:));
                % Grey peak
                b3 = bar(xpos + x_offsets(3)*bar_width, cumulative_peak_grey(i), bar_width, 'FaceColor', [0.7 0.7 0.7]);
                % Grey mean
                b4 = bar(xpos + x_offsets(4)*bar_width, cumulative_mean_grey(i), bar_width, 'FaceColor', [0.7 0.7 0.7]);
            
                % Add dividing lines for stim windows on stimulus bars
                cum_val = 0;
                for w = 1:size(all_peak_FR_by_stimwindow, 2)
                    cum_val = cum_val + all_peak_FR_by_stimwindow(i, w);
                    plot([xpos + x_offsets(1)*bar_width - bar_width/2, xpos + x_offsets(1)*bar_width + bar_width/2], ...
                         [cum_val, cum_val], 'k-', 'LineWidth', 1)
                    if ~strcmp(z_score_period, 'none')
                        % Stimulus peak error bars
                        y_offset = 0;
                        for w = 1:size(all_peak_FR_by_stimwindow, 2)
                            y_val = y_offset + all_peak_FR_by_stimwindow(i, w);
                            y_err = sem_peak_FR_by_stimwindow(i, w);
                        
                            errorbar(xpos + x_offsets(1)*bar_width, y_val, y_err, ...
                                     'k.', 'CapSize', 8, 'LineWidth', 1);
                            y_offset = y_val;
                        end
                    end
                end
            
                cum_val = 0;
                for w = 1:size(all_mean_FR_by_stimwindow, 2)
                    cum_val = cum_val + all_mean_FR_by_stimwindow(i, w);
                    plot([xpos + x_offsets(2)*bar_width - bar_width/2, xpos + x_offsets(2)*bar_width + bar_width/2], ...
                         [cum_val, cum_val], 'k-', 'LineWidth', 1)
                    if ~strcmp(z_score_period, 'none')
                        % Stimulus mean error bars
                        y_offset = 0;
                        for w = 1:size(all_mean_FR_by_stimwindow, 2)
                            y_val = y_offset + all_mean_FR_by_stimwindow(i, w);
                            y_err = sem_mean_FR_by_stimwindow(i, w);
                        
                            errorbar(xpos + x_offsets(2)*bar_width, y_val, y_err, ...
                                     'k.', 'CapSize', 8, 'LineWidth', 1);
                            y_offset = y_val;
                        end
                    end    
                end
            
                cum_val = 0;
                for w = 1:size(all_peak_FR_by_greywindow, 2)
                    cum_val = cum_val + all_peak_FR_by_greywindow(i, w);
                    plot([xpos + x_offsets(3)*bar_width - bar_width/2, xpos + x_offsets(3)*bar_width + bar_width/2], ...
                         [cum_val, cum_val], 'k-', 'LineWidth', 1)
                    if ~strcmp(z_score_period, 'none')
                        % Grey peak error bars
                        y_offset = 0;
                        for w = 1:size(all_peak_FR_by_greywindow, 2)
                            y_val = y_offset + all_peak_FR_by_greywindow(i, w);
                            y_err = sem_peak_FR_by_greywindow(i, w);
                        
                            errorbar(xpos + x_offsets(3)*bar_width, y_val, y_err, ...
                                     'k.', 'CapSize', 8, 'LineWidth', 1);
                            y_offset = y_val;
                        end
                    end 
                end
            
                cum_val = 0;
                for w = 1:size(all_mean_FR_by_greywindow, 2)
                    cum_val = cum_val + all_mean_FR_by_greywindow(i, w);
                    plot([xpos + x_offsets(4)*bar_width - bar_width/2, xpos + x_offsets(4)*bar_width + bar_width/2], ...
                         [cum_val, cum_val], 'k-', 'LineWidth', 1)
                    if ~strcmp(z_score_period, 'none')
                        % Grey mean error bars
                        y_offset = 0;
                        for w = 1:size(all_mean_FR_by_greywindow, 2)
                            y_val = y_offset + all_mean_FR_by_greywindow(i, w);
                            y_err = sem_mean_FR_by_greywindow(i, w);
                        
                            errorbar(xpos + x_offsets(4)*bar_width, y_val, y_err, ...
                                     'k.', 'CapSize', 8, 'LineWidth', 1);
                            y_offset = y_val;
                        end
                    end    
                    
                end
                hold on;
            end
            
            xticks(1:num_sessions);
            xticklabels(arrayfun(@(x) sprintf('Day %d', experiment_info(x).date), sessions_to_plot, 'UniformOutput', false));
            set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 14);
            ylabel(sprintf('Cumulative Z-scored FR (z-scored over %s)', z_score_period), 'Interpreter', 'none', 'FontSize', 14);
            sgtitle(sprintf('%s - %s: L4 Cumulative Peak and Mean Firing Rates Across Days', subject_number, Stimulus_type), 'Interpreter', 'none');
            legend({'Stim Peak', 'Stim Mean', 'Grey Peak', 'Grey Mean'}, 'Location', 'northeast', 'FontSize', 14);
            grid on;

            % Save bar chart figure
            saveas(fig2, fullfile(pwd, sprintf('%s %s - Multiplot cumulative peak and mean FRs across %d days.fig', subject_number, Stimulus_type, length(sessions_to_plot))));
            exportgraphics(fig2, fullfile(pwd, sprintf('%s %s - Multiplot cumulative peak and mean FRs across %d days.png', subject_number, Stimulus_type, length(sessions_to_plot))));


        end

       
        


        if contains(Stimulus_type, 'GAVNIK_ABCD') && contains(plot_choice, 'aggregate') && contains(plot_type, 'FR')
            for nprobe = 1:length(clusters)
                selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters, which pass the set parameters
                
                depth_selected_clusters = selected_clusters; % initialize with same structure
                for np = 1:length(clusters)
                    sc = selected_clusters(np);
                    peak_depths = clusters(np).peak_depth(sc.cluster_id); % get depths of selected clusters

                    % Find cluster_ids within selected depth range
                    keep_mask = peak_depths >= min(depth_range) & peak_depths <= max(depth_range);
                    depth_cluster_ids = sc.cluster_id(keep_mask);

                    % If filtering for pyramidal neurons
                    if contains(neuron_type, 'PYR only')
                        % Identify putative pyramidal neurons
                        is_pyramidal = (clusters(np).cell_type == 1); % 1 is PYR, 2 is PV, 3 is SOM. putatively
                        pyramidal_cluster_ids = find(is_pyramidal);  % These are indices in clusters(np)
                
                        % Only keep depth_cluster_ids that are also pyramidal
                        depth_cluster_ids = intersect(depth_cluster_ids, pyramidal_cluster_ids);
                    end
                
                    % Keep only spike times and IDs for clusters at selected depth (and type, if filtered)
                    depth_selected_clusters(np).cluster_id = depth_cluster_ids;
                    depth_selected_clusters(np).spike_times = sc.spike_times(ismember(sc.spike_id, depth_cluster_ids));
                    depth_selected_clusters(np).spike_id = sc.spike_id(ismember(sc.spike_id, depth_cluster_ids));
                end
                
                baseline_window = [0 30]; % in seconds - first 30s of recording is grey screen - can use for z-scoring

                % Define time_edges depending on z_score_period
                if contains(z_score_period, 'entire_session')
                    time_edges = 0:psthBinSize:max(depth_selected_clusters(nprobe).spike_times); %time bin edges from 0 to the max time in steps of ptshBinSize
                elseif contains(z_score_period, 'first30secs')
                    time_edges = baseline_window(1):psthBinSize:baseline_window(2);
                end

                % Combine spike times for all clusters at selected depth (i.e., for curated MUA)
                all_spike_times = [];
                all_spike_ids = [];
                for np = 1:length(depth_selected_clusters)
                    all_spike_times = [all_spike_times; depth_selected_clusters(np).spike_times];
                    all_spike_ids = [all_spike_ids; depth_selected_clusters(np).spike_id];
                end

                % Compute appropriate histogram counts for z-scoring
                if contains(z_score_period, 'entire_session')
                    zscore_counts = histcounts(all_spike_times, time_edges);
                elseif contains(z_score_period, 'first30secs')
                    baseline_spikes = all_spike_times(all_spike_times >= baseline_window(1) & all_spike_times <= baseline_window(2));
                    zscore_counts = histcounts(baseline_spikes, time_edges); %histcounts counts the number of datapoints in specified time bins
                end

                ordered_oris = unique(Task_info.stim_orientation, 'stable');               
                ori = 1;
                stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));

                [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.3 1.35], psthBinSize);  
                mean_trace = mean(binnedArray, 1); % average over trials
                
                figure(fig1);
                hold on;

                if contains(z_score_period, 'none')
                    % Plot raw mean firing rate trace
                    plot(bins, mean_trace, 'Color', colors(find(sessions_to_plot == nsession),:), 'LineWidth', 1.5, 'DisplayName', sprintf('Day %d, L4 %d–%d µm', experiment_info(nsession).date, min(depth_range), max(depth_range)));
                    ylabel('Mean firing rate (Hz)', 'FontSize', 14);
                    ylim ([0 16]);
                else 
                    baseline_mean = mean(zscore_counts);
                    baseline_std = std(zscore_counts);
                    
                    % Z-score each trial individually so that trial to variability remains visible and SEM can be calculated
                    zscored_trials = (binnedArray - baseline_mean) / baseline_std;  % [trials x time]
                    z_trace = mean(zscored_trials, 1);                             % mean z-scored trace
                    sem_trace = std(zscored_trials, 0, 1) / sqrt(size(zscored_trials, 1));  % SEM

                    % Plot SEM shading
                    fill([bins, fliplr(bins)], ...
                         [z_trace + sem_trace, fliplr(z_trace - sem_trace)], ...
                         colors(find(sessions_to_plot == nsession), :), ...
                         'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                    
                    plot(bins, z_trace, 'Color', colors(find(sessions_to_plot == nsession),:), 'LineWidth', 1.5, 'DisplayName', sprintf('Day %d, L4 %d–%d µm', experiment_info(nsession).date, min(depth_range), max(depth_range)));

                    if contains(z_score_period, 'entire_session')
                        ylabel('Z-scored FR (z-scored over entire session)', 'FontSize', 14);
                        ylim([-1 5]);
                    elseif contains(z_score_period, 'first30secs')
                        ylabel('Z-scored FR (z-scored over first 30s baseline)', 'FontSize', 14);
                        ylim([-1 8]);
                    end
                                        
                end    
                
                peak_FR_by_stimwindow = zeros(1, length(stim_window_starts));
                mean_FR_by_stimwindow = zeros(1, length(stim_window_starts));
                peak_FR_by_greywindow = zeros(1, length(grey_window_starts));
                mean_FR_by_greywindow = zeros(1, length(grey_window_starts));

                % per trial values for SEM calculation
                peak_FR_trials_by_stimwindow = zeros(size(zscored_trials, 1), length(stim_window_starts)); % preallocate
                mean_FR_trials_by_stimwindow = zeros(size(zscored_trials, 1), length(stim_window_starts)); % preallocate
                peak_FR_trials_by_greywindow = zeros(size(zscored_trials, 1), length(grey_window_starts)); % preallocate
                mean_FR_trials_by_greywindow = zeros(size(zscored_trials, 1), length(grey_window_starts)); % preallocate
                                
                for i = 1:length(stim_window_starts)
                    idx_in_stimwindow = bins >= stim_window_starts(i) & bins < stim_window_ends(i);
                    if contains(z_score_period, 'none')
                        peak_FR_by_stimwindow(i) = max(mean_trace(idx_in_stimwindow));
                        mean_FR_by_stimwindow(i) = mean(mean_trace(idx_in_stimwindow));
                        %for SEM
                        peak_FR_trials_by_stimwindow(:, i) = max(binnedArray(:, idx_in_stimwindow), [], 2);
                        mean_FR_trials_by_stimwindow(:, i) = mean(binnedArray(:, idx_in_stimwindow), 2);
                                                
                    else
                        peak_FR_by_stimwindow(i) = max(z_trace(idx_in_stimwindow));
                        mean_FR_by_stimwindow(i) = mean(z_trace(idx_in_stimwindow));
                        %for SEM
                        mean_FR_trials_by_stimwindow(:, i) = mean(zscored_trials(:, idx_in_stimwindow), 2);
                        peak_FR_trials_by_stimwindow(:, i) = max(zscored_trials(:, idx_in_stimwindow), [], 2);
                    end
                end
                
                for i = 1:length(grey_window_starts)
                    idx_in_greywindow = bins >= grey_window_starts(i) & bins < grey_window_ends(i);
                    if contains(z_score_period, 'none')              
                        peak_FR_by_greywindow(i) = max(mean_trace(idx_in_greywindow));
                        mean_FR_by_greywindow(i) = mean(mean_trace(idx_in_greywindow));
                        %for SEM
                        peak_FR_trials_by_greywindow(:, i) = max(binnedArray(:, idx_in_greywindow), [], 2);
                        mean_FR_trials_by_greywindow(:, i) = mean(binnedArray(:, idx_in_greywindow), 2);                       
                    else                  
                        peak_FR_by_greywindow(i) = max(z_trace(idx_in_greywindow));
                        mean_FR_by_greywindow(i) = mean(z_trace(idx_in_greywindow));
                        %for SEM
                        peak_FR_trials_by_greywindow(:, i) = max(zscored_trials(:, idx_in_greywindow), [], 2);
                        mean_FR_trials_by_greywindow(:, i) = mean(zscored_trials(:, idx_in_greywindow), 2);
                        
                    end
                end
                
                session_idx = find(sessions_to_plot == nsession); % Map actual session number to index in preallocated array
                
                all_peak_FR_by_stimwindow(session_idx, :) = peak_FR_by_stimwindow;
                all_mean_FR_by_stimwindow(session_idx, :) = mean_FR_by_stimwindow;
                all_peak_FR_by_greywindow(session_idx, :) = peak_FR_by_greywindow;
                all_mean_FR_by_greywindow(session_idx, :) = mean_FR_by_greywindow;
                
                sem_peak_FR_by_stimwindow(session_idx, :) = std(peak_FR_trials_by_stimwindow, 0, 1) / sqrt(size(peak_FR_trials_by_stimwindow,1));
                sem_mean_FR_by_stimwindow(session_idx, :) = std(mean_FR_trials_by_stimwindow, 0, 1) / sqrt(size(mean_FR_trials_by_stimwindow,1));
                sem_peak_FR_by_greywindow(session_idx, :) = std(peak_FR_trials_by_greywindow, 0, 1) / sqrt(size(peak_FR_trials_by_greywindow,1));
                sem_mean_FR_by_greywindow(session_idx, :) = std(mean_FR_trials_by_greywindow, 0, 1) / sqrt(size(mean_FR_trials_by_greywindow,1));
          
                xlim([-0.5 1.5]);
                xticks(-0.4:0.2:1.4);
                set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 14);
                xlabel('Time (s) since onset of A', 'FontSize', 14)
                
                legend(flipud(findobj(gca,'-property','DisplayName')), 'Location', 'northeast', 'FontSize', 14);
                hold on;                
            end
                    % Define grey intervals
            grey_intervals = [-0.5 0; 0.6 1.5];
                           
            % Get current y-axis limits for full vertical shading
            yl = ylim;
                
            % Shade each interval
            for i = 1:size(grey_intervals, 1)
                x = [grey_intervals(i,1), grey_intervals(i,2), grey_intervals(i,2), grey_intervals(i,1)];
                y = [yl(1), yl(1), yl(2), yl(2)];
                fill(x, y, [0.7 0.7 0.7], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % grey color with transparency
            end

            xline(0, 'k', (sprintf('A %d%s onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
            xline(0.15, 'k', (sprintf('B %d%s onset', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
            xline(0.30, 'k', (sprintf('C %d%s onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
            xline(0.45, 'k', (sprintf('D %d%s onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);

            sgtitle(sprintf('%s - %s: Aggregate L4 single unit activity across days', subject_number, Stimulus_type), 'Interpreter', 'none');
            filename = sprintf('%s - %s - Multiplot L4 FR traces across %d days.png', subject_number, Stimulus_type, length(sessions_to_plot));
            save_path = fullfile(pwd, filename);
            exportgraphics(fig1, save_path); 

            % Also save as .fig
            fig_filename = sprintf('%s - %s - Multiplot L4 FR traces across %d days.fig', subject_number, Stimulus_type, length(sessions_to_plot));
            fig_save_path = fullfile(pwd, fig_filename);
            savefig(fig1, fig_save_path);
            
            
            % Plot bar chart showing development of peak and mean firing across training days
            figure(fig2); % switch to bar chart figure
            hold on;
                        
            num_sessions = length(sessions_to_plot);
            cumulative_peak_stim = sum(all_peak_FR_by_stimwindow(:, :), 2);
            cumulative_mean_stim = sum(all_mean_FR_by_stimwindow(:, :), 2);
            cumulative_peak_grey = sum(all_peak_FR_by_greywindow(:, :), 2);
            cumulative_mean_grey = sum(all_mean_FR_by_greywindow(:, :), 2);
                      
            % Set bar width and spacing
            num_bars = 4;
            bar_width = 0.15;
            x_offsets = ((1:num_bars) - (num_bars+1)/2); % e.g., [-1.5, -0.5, 0.5, 1.5]
            x = 1:num_sessions;
            
            for i = 1:num_sessions
                xpos = x(i);
            
                % Stimulus peak
                b1 = bar(xpos + x_offsets(1)*bar_width, cumulative_peak_stim(i), bar_width, 'FaceColor', colors(i,:));
                % Stimulus mean
                b2 = bar(xpos + x_offsets(2)*bar_width, cumulative_mean_stim(i), bar_width, 'FaceColor', colors(i,:));
                % Grey peak
                b3 = bar(xpos + x_offsets(3)*bar_width, cumulative_peak_grey(i), bar_width, 'FaceColor', [0.7 0.7 0.7]);
                % Grey mean
                b4 = bar(xpos + x_offsets(4)*bar_width, cumulative_mean_grey(i), bar_width, 'FaceColor', [0.7 0.7 0.7]);
            
                % Add dividing lines for stim windows on stimulus bars
                cum_val = 0;
                for w = 1:size(all_peak_FR_by_stimwindow, 2)
                    cum_val = cum_val + all_peak_FR_by_stimwindow(i, w);
                    plot([xpos + x_offsets(1)*bar_width - bar_width/2, xpos + x_offsets(1)*bar_width + bar_width/2], ...
                         [cum_val, cum_val], 'k-', 'LineWidth', 1)
                    
                    if ~strcmp(z_score_period, 'none')
                        % Stimulus peak error bars
                        y_offset = 0;
                        for w = 1:size(all_peak_FR_by_stimwindow, 2)
                            y_val = y_offset + all_peak_FR_by_stimwindow(i, w);
                            y_err = sem_peak_FR_by_stimwindow(i, w);
                        
                            errorbar(xpos + x_offsets(1)*bar_width, y_val, y_err, ...
                                     'k.', 'CapSize', 8, 'LineWidth', 1);
                            y_offset = y_val;
                        end
                    end
                end
            
                cum_val = 0;
                for w = 1:size(all_mean_FR_by_stimwindow, 2)
                    cum_val = cum_val + all_mean_FR_by_stimwindow(i, w);
                    plot([xpos + x_offsets(2)*bar_width - bar_width/2, xpos + x_offsets(2)*bar_width + bar_width/2], ...
                         [cum_val, cum_val], 'k-', 'LineWidth', 1)
                    if ~strcmp(z_score_period, 'none')
                        % Stimulus mean error bars
                        y_offset = 0;
                        for w = 1:size(all_mean_FR_by_stimwindow, 2)
                            y_val = y_offset + all_mean_FR_by_stimwindow(i, w);
                            y_err = sem_mean_FR_by_stimwindow(i, w);
                        
                            errorbar(xpos + x_offsets(2)*bar_width, y_val, y_err, ...
                                     'k.', 'CapSize', 8, 'LineWidth', 1);
                            y_offset = y_val;
                        end
                    end    
                end
            
                cum_val = 0;
                for w = 1:size(all_peak_FR_by_greywindow, 2)
                    cum_val = cum_val + all_peak_FR_by_greywindow(i, w);
                    plot([xpos + x_offsets(3)*bar_width - bar_width/2, xpos + x_offsets(3)*bar_width + bar_width/2], ...
                         [cum_val, cum_val], 'k-', 'LineWidth', 1)
                    if ~strcmp(z_score_period, 'none')
                        % Grey peak error bars
                        y_offset = 0;
                        for w = 1:size(all_peak_FR_by_greywindow, 2)
                            y_val = y_offset + all_peak_FR_by_greywindow(i, w);
                            y_err = sem_peak_FR_by_greywindow(i, w);
                        
                            errorbar(xpos + x_offsets(3)*bar_width, y_val, y_err, ...
                                     'k.', 'CapSize', 8, 'LineWidth', 1);
                            y_offset = y_val;
                        end
                    end    
                end
            
                cum_val = 0;
                for w = 1:size(all_mean_FR_by_greywindow, 2)
                    cum_val = cum_val + all_mean_FR_by_greywindow(i, w);
                    plot([xpos + x_offsets(4)*bar_width - bar_width/2, xpos + x_offsets(4)*bar_width + bar_width/2], ...
                         [cum_val, cum_val], 'k-', 'LineWidth', 1)
                    if ~strcmp(z_score_period, 'none')
                        % Grey mean error bars
                        y_offset = 0;
                        for w = 1:size(all_mean_FR_by_greywindow, 2)
                            y_val = y_offset + all_mean_FR_by_greywindow(i, w);
                            y_err = sem_mean_FR_by_greywindow(i, w);
                        
                            errorbar(xpos + x_offsets(4)*bar_width, y_val, y_err, ...
                                     'k.', 'CapSize', 8, 'LineWidth', 1);
                            y_offset = y_val;
                        end
                    end    
                end
                hold on;
            end
            
            xticks(1:num_sessions);
            ylabel(sprintf('Cumulative Z-scored FR (z-scored over %s)', z_score_period), 'Interpreter', 'none', 'FontSize', 14);
            xlabel('Training Day', 'FontSize', 14)
            set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 14);
            sgtitle(sprintf('%s - %s: L4 Cumulative Peak and Mean Firing Rates Across Days', subject_number, Stimulus_type), 'Interpreter', 'none');
            legend({'Stim Peak', 'Stim Mean', 'Grey Peak', 'Grey Mean'}, 'Location', 'northeast', 'FontSize', 14);
            grid on;

            % Save bar chart figure
            saveas(fig2, fullfile(pwd, sprintf('%s %s - Multiplot cumulative peak and mean FRs across %d days.fig', subject_number, Stimulus_type, length(sessions_to_plot))));
            exportgraphics(fig2, fullfile(pwd, sprintf('%s %s - Multiplot cumulative peak and mean FRs across %d days.png', subject_number, Stimulus_type, length(sessions_to_plot))));


        end
                    
               
    end
end
hold off;

