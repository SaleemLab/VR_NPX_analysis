%%% For analysis of unit spiking in response to visual stimuli

addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

%% setting metrics to screen good clusters
clear all
% Choose your probe depth of interest
session_specific_L4 = {4, 4440:4580; 5, 4440:4580; 6, 4460:4600; 7, 4480:4620; 8, 4500:4640;}; % 1/5. um. Set for each SESSION based on CSD +/- 70um
 
SUBJECTS = {'M00013'};

params = create_cluster_selection_params('sorting_option','ellie');
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);

%%% 2/5
Stimulus_type = 'TRAIN'; 
plot_choice = 'aggregate'; % curated 'single_units' or in 'aggregate' or uncurated 'MUA'; MUA includes all clusters from kilosort, unfiltered
plot_type = 'FR'; % 'FR' firing rate or 'raster'
z_score_period = 'first30secs'; % z score either over 'entire_session' or 'first30secs' (for every stimulus recording
% session from 20250205 onward, I presented grey screen to the mouse for at least 30s before starting the stimulus)
cd('V:\Ellie\DATA\SUBJECTS\M00013\analysis') % 3/5 files will be saved here in the cd

% 4/5 Define 150 ms windows after time zero to calc. the peak and mean FR - Depends on stimulus type              
window_starts = 0:0.15:1.35;
window_ends = 0.15:0.15:1.5;

sessions_to_plot = [4, 5, 6, 7, 8]; %5/5 row numbers of recording dates in "experiment_info" 
colors = [
    1.0,  0.8,  0.0;  % 1st group - golden yellow (1,0, 1.0, 0.0 for pure yellow)
    0.4,  0.8,  0.0;  % 2nd group - yellow-green (more yellowish)
    0.2,  0.6,  0.4;  % 3rd group - pure green
    0.2,  0.4,  0.6;  % 4th group - greenish blue (not too blue)
    0.0,  0.0,  0.8];   % 5th group - deep blue

fig1 = figure; % for traces
fig2 = figure; % for bar charts
all_peak_zFR_by_window = zeros(length(sessions_to_plot), length(window_starts));
all_mean_zFR_by_window = zeros(length(sessions_to_plot), length(window_starts));

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
                    depth_cluster_ids = sc.cluster_id(keep_mask);

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
                z_trace = (mean_trace - mean(zscore_counts)) / std(zscore_counts); %mean(zscore_counts) gives the mean spikes per timebin in the reference period; this is then deducted from the spikecount of each trial-averaged timebin
                
                figure(fig1);
                hold on;
                plot(bins, z_trace, 'Color', colors(find(sessions_to_plot == nsession),:), 'LineWidth', 1.5, 'DisplayName', sprintf('Day %d, L4 %d–%d µm', experiment_info(nsession).date, min(depth_range), max(depth_range)));
                
                peak_zFR_by_window = zeros(1, length(window_starts)); % preallocate
                mean_zFR_by_window = zeros(1, length(window_starts)); % preallocate
                ylim([-1 8]);
                                
                for i = 1:length(window_starts)
                    % Find indices of bins within current window
                    idx_in_window = bins >= window_starts(i) & bins < window_ends(i);                    
                    peak_zFR_by_window(i) = max(z_trace(idx_in_window));
                    mean_zFR_by_window(i) = mean(z_trace(idx_in_window));
                    
                end
                session_idx = find(sessions_to_plot == nsession); % Map actual session number to index in preallocated array
                all_peak_zFR_by_window(session_idx, :) = peak_zFR_by_window;
                all_mean_zFR_by_window(session_idx, :) = mean_zFR_by_window;
                
                xlim([-0.5 2.0]);
                xticks(-0.4:0.2:1.8);
                ylim([-1 8]);
                xlabel('Time (s) since onset of A')
                ylabel('Z-scored FR');
                legend(flipud(findobj(gca,'-property','DisplayName')), 'Location', 'northeast');
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

            xline(0, 'k', (sprintf('A onset; %d%s', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
            xline(0.30, 'k', (sprintf('B onset; %d%s', round(rad2deg(ordered_oris(2))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
            xline(0.60, 'k', (sprintf('C onset; %d%s', round(rad2deg(ordered_oris(3))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
            xline(0.90, 'k', (sprintf('D onset; %d%s', round(rad2deg(ordered_oris(4))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
                
            sgtitle(sprintf('%s - %s: Aggregate L4 single unit activity across days', subject_number, Stimulus_type), 'Interpreter', 'none');
            filename = sprintf('%s - Multiplot - aggregate L4 single unit FRs.png', Stimulus_type);
            save_path = fullfile(pwd, filename);
            exportgraphics(fig1, save_path); 

            % Also save as .fig
            fig_filename = sprintf('%s - %s - Multiplot aggregate L4 single unit FRs.fig', subject_number, Stimulus_type);
            fig_save_path = fullfile(pwd, fig_filename);
            savefig(fig1, fig_save_path);
            
            
            % Plot bar chart showing development of peak and mean firing across training days
            figure(fig2); % switch to bar chart figure
            hold on;
            stim_windows = [1, 3, 5, 7];
            grey_windows = [2, 4, 6, 8, 9, 10];
            
            num_sessions = length(sessions_to_plot);
            cumulative_peak_stim = sum(all_peak_zFR_by_window(:, stim_windows), 2);
            cumulative_mean_stim = sum(all_mean_zFR_by_window(:, stim_windows), 2);
            cumulative_peak_grey = sum(all_peak_zFR_by_window(:, grey_windows), 2);
            cumulative_mean_grey = sum(all_mean_zFR_by_window(:, grey_windows), 2);
                      
            % Set bar width and spacing
            bar_width = 0.15;
            x_offsets = [-1.5, -0.5, 0.5, 1.5];
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
                for w = stim_windows
                    cum_val = cum_val + all_peak_zFR_by_window(i, w);
                    plot([xpos + x_offsets(1)*bar_width - bar_width/2, xpos + x_offsets(1)*bar_width + bar_width/2], ...
                         [cum_val, cum_val], 'k-', 'LineWidth', 1)
                end
            
                cum_val = 0;
                for w = stim_windows
                    cum_val = cum_val + all_mean_zFR_by_window(i, w);
                    plot([xpos + x_offsets(2)*bar_width - bar_width/2, xpos + x_offsets(2)*bar_width + bar_width/2], ...
                         [cum_val, cum_val], 'k-', 'LineWidth', 1)
                end
            
                cum_val = 0;
                for w = grey_windows
                    cum_val = cum_val + all_peak_zFR_by_window(i, w);
                    plot([xpos + x_offsets(3)*bar_width - bar_width/2, xpos + x_offsets(3)*bar_width + bar_width/2], ...
                         [cum_val, cum_val], 'k-', 'LineWidth', 1)
                end
            
                cum_val = 0;
                for w = grey_windows
                    cum_val = cum_val + all_mean_zFR_by_window(i, w);
                    plot([xpos + x_offsets(4)*bar_width - bar_width/2, xpos + x_offsets(4)*bar_width + bar_width/2], ...
                         [cum_val, cum_val], 'k-', 'LineWidth', 1)
                end
                hold on;
            end
            
            xticks(1:num_sessions);
            xticklabels(arrayfun(@(x) sprintf('Day %d', experiment_info(x).date), sessions_to_plot, 'UniformOutput', false));
            ylabel('Cumulative z-scored FR');
            sgtitle(sprintf('%s - %s: L4 Cumulative Peak and Mean Z-scored Firing Rates Across Days', subject_number, Stimulus_type), 'Interpreter', 'none');
            legend({'Stim Peak', 'Stim Mean', 'Grey Peak', 'Grey Mean'}, 'Location', 'northwest');
            grid on;

            % Save bar chart figure
            saveas(fig2, fullfile(pwd, sprintf('%s %s - Cumulative peak and mean FRs by day.fig', subject_number, Stimulus_type)));
            exportgraphics(fig2, fullfile(pwd, sprintf('%s %s - Cumulative peak and mean FRs by day.png', subject_number, Stimulus_type)));


        end

       
        


        if (contains(Stimulus_type, 'GAVNIK_ABCD') || contains(Stimulus_type, 'GAVNIK DCBA')) && contains(plot_choice, 'aggregate') && contains(plot_type, 'FR')
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
                fig = figure;
                fig.Name = sprintf('%s: Aggregate activity: %s depth range %d–%d μm', Stimulus_type, depth_for_analysis, min(depth_range), max(depth_range));
                fig.Position = [114 90 770 650];
                tiledlayout(5,1);

                for ori = 1:length(ordered_oris)
                    stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));

                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.02 0.17], psthBinSize);
                    
                    mean_trace = mean(binnedArray, 1); % average over trials
                    z_trace = (mean_trace - mean(zscore_counts)) / std(zscore_counts); % z-score using session-wide stats

                    nexttile;
                    plot(bins, z_trace, 'k', 'LineWidth', 1.5);
                    xline(0,'r', 'LineWidth', 1);
                    xlim([-0.02 0.17]);
                    xticks([0, 0.05, 0.1, 0.15]);
                    ylabel('Z-scored FR');
                    ori_deg = round(ordered_oris(ori));
                    title(sprintf('Orientation %d%s', ori_deg, char(176)));
                    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
                end

                % Extra tile to show post-stimulus oscillations following final stim in seq
                if length(ordered_oris) >= 4
                    stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(4));

                    [~, bins_long, ~, ~, ~, binnedArray_long] = psthAndBA(all_spike_times, stim_onsets, [-0.02 0.75], psthBinSize);
    
                    mean_trace_long = mean(binnedArray_long, 1);
                    z_trace_long = (mean_trace_long - mean(zscore_counts)) / std(zscore_counts);

                    nexttile;
                    plot(bins_long, z_trace_long, 'k', 'LineWidth', 1.5);
                    xline(0,'r', 'LineWidth', 1);
                    xlim([-0.02 0.75]);
                    xticks(0:0.05:0.75);
                    ylabel('Z-scored FR');
                    title(sprintf('Orientation %d%s (Extended to show post-stimulus oscillations)', round(ordered_oris(4)), char(176)));
                    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
                end

                sgtitle(sprintf('%s - %s: Aggregate single unit activity: %s depth range %d–%d μm', subject_number, Stimulus_type, depth_for_analysis, min(depth_range), max(depth_range)), 'Interpreter', 'none');

                filename = sprintf('Aggregate single unit %s FRs.pdf', depth_for_analysis);
                save_path = fullfile(pwd, filename);
                saveas(fig, save_path);
            end
        end




        
        
        if (contains(Stimulus_type, 'GAVNIK_A_CD') || contains(Stimulus_type, 'GAVNIK_E_CD')) && contains(plot_choice, 'aggregate') && contains(plot_type, 'FR')
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
                fig = figure;
                fig.Name = sprintf('%s: Aggregate activity: %s depth range %d–%d μm', Stimulus_type, depth_for_analysis, min(depth_range), max(depth_range));
                fig.Position = [114 90 770 650];
                tiledlayout(5,1);

                for ori = 1:length(ordered_oris)
                    stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));

                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.02 0.17], psthBinSize);
                    
                    mean_trace = mean(binnedArray, 1); % average over trials
                    z_trace = (mean_trace - mean(zscore_counts)) / std(zscore_counts); % from each trial-averaged timebin of spikes, deduct the mean spikes per timebin in the reference period (baseline or entire session)

                    nexttile;
                    plot(bins, z_trace, 'k', 'LineWidth', 1.5);
                    xline(0,'r', 'LineWidth', 1);
                    xlim([-0.02 0.17]);
                    xticks([0, 0.05, 0.1, 0.15]);
                    ylabel('Z-scored FR');
                    
                    ori_deg = round(ordered_oris(ori));
                    if ori == 2
                        title('Grey screen'); % in the GAVNIK_A_CD and GAVNIK_E_CD conditions, the second stimulus is absent
                    else
                        title(sprintf('Orientation %d%s', ori_deg, char(176)));
                    end

                    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
                end

                % Extra tile to show post-stimulus oscillations following final stim in seq
                if length(ordered_oris) >= 4
                    stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(4));

                    [~, bins_long, ~, ~, ~, binnedArray_long] = psthAndBA(all_spike_times, stim_onsets, [-0.02 0.75], psthBinSize);
    
                    mean_trace_long = mean(binnedArray_long, 1);
                    z_trace_long = (mean_trace_long - mean(zscore_counts)) / std(zscore_counts);

                    nexttile;
                    plot(bins_long, z_trace_long, 'k', 'LineWidth', 1.5);
                    xline(0,'r', 'LineWidth', 1);
                    xlim([-0.02 0.75]);
                    xticks(0:0.05:0.75);
                    ylabel('Z-scored FR');
                    title(sprintf('Orientation %d%s (Extended to show post-stimulus oscillations)', round(ordered_oris(4)), char(176)));
                    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
                end

                sgtitle(sprintf('%s - %s: Aggregate single unit activity: %s depth range %d–%d μm', subject_number, Stimulus_type, depth_for_analysis, min(depth_range), max(depth_range)), 'Interpreter', 'none');
                
                filename = sprintf('Aggregate single unit %s FRs.pdf', depth_for_analysis);
                save_path = fullfile(pwd, filename);
                saveas(fig, save_path);
            end
        end
       
    end
end
hold off;

