%%% For analysis of unit spiking in response to visual stimuli

addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

%% setting metrics to screen good clusters
clear all
% Choose your probe depth of interest
session_specific_L4 = {15, 4510:4650;}; % 1/5. um. Set for each SESSION based on CSD +/- 70um
% 4, 4440:4580; 5, 4440:4580; 6, 4460:4600; 7, 4480:4620; 8, 4500:4640; 
% 11, 4500:4640; 12, 4510:4650; 13, 4510:4650; 14, 4510:4650; 15, 4510:4650;
SUBJECTS = {'M00013'};

params = create_cluster_selection_params('sorting_option','ellie');
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);

%%% 2/5
Stimulus_types = {'GAVNIK_ABCD', 'GAVNIK_A_CD', 'GAVNIK_E_CD'}; % IN ORDER with abcd first. 'GAVNIK_A_CD', 'GAVNIK_E_CD', 'GAVNIK DCBA'
plot_type = 'FR'; % 'FR' firing rate or 'raster'
z_score_period = 'entire_session'; % 'none' = no z scoring or z score either over 'entire_session' or 'first30secs' (for every stimulus recording
% session from 20250205 onward, I presented grey screen to the mouse for at least 30s before starting the stimulus).
cd('V:\Ellie\DATA\SUBJECTS\M00013\analysis\20250221') % 3/5 files will be saved here in the cd

sessions_to_plot = [15]; %4/5 row numbers of recording dates in "experiment_info" 11, 12, 13, 14, 15

if contains(Stimulus_types{2}, 'GAVNIK DCBA')
    stimulus_colors = [
        0.0,  0.0,  1.0;  % 1st group - blue - GAVNIK_ABCD
        1.0,  0.0,  0.0;];  % 2nd group - red - GAVNIK_DCBA    
elseif contains(Stimulus_types{2}, 'GAVNIK_A_CD')
    stimulus_colors = [
        0.0,  0.0,  1.0;  % 1st group - blue - GAVNIK_ABCD
        1.0,  0.0,  0.0;  % 2nd group - red - GAVNIK_A_CD
        0.0,  1.0,  0.0;];  % 3rd group - green - GAVNIK_E_CD  
end    

% 5/5 Define windows after time zero to calc. the peak and mean FR - Depends on stimulus type              
if contains(Stimulus_types{1}, 'GAVNIK')
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
%fig2 = figure; % for bar charts
all_peak_FR_by_stimwindow = zeros(length(sessions_to_plot), length(stim_window_starts));
all_mean_FR_by_stimwindow = zeros(length(sessions_to_plot), length(stim_window_starts));
all_peak_FR_by_greywindow = zeros(length(sessions_to_plot), length(grey_window_starts));
all_mean_FR_by_greywindow = zeros(length(sessions_to_plot), length(grey_window_starts));


for nsession = sessions_to_plot 
    for stim_idx = 1:length(Stimulus_types)
        Stimulus_type = Stimulus_types{stim_idx};
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
                    
    
            if contains(Stimulus_types{2}, 'GAVNIK DCBA') % only make this plot if GAVNIK DCBA is present
                if (contains(Stimulus_type, 'GAVNIK DCBA') || contains(Stimulus_type, 'GAVNIK_ABCD')) && contains(plot_type, 'FR')
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
                        ori = 1;
                        stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));
        
                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.3 1.35], psthBinSize);  
                        mean_trace = mean(binnedArray, 1); % average over trials
                        
                        figure(fig1);
                        hold on;
        
                        if contains(z_score_period, 'none')
                            % Plot raw mean firing rate trace
                            plot(bins, mean_trace, 'Color', stimulus_colors(stim_idx,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%s (200x)', Stimulus_types{stim_idx}));
                            ylabel('Mean firing rate (Hz)', 'FontSize', 14);
                            ylim ([0 16]);
                        else 
                            baseline_mean = mean(zscore_counts);
                            baseline_std = std(zscore_counts);
                            % Z-score each trial individually so that trial to variability remains visible and SEM can be calculated
                            zscored_trials = (binnedArray - baseline_mean) / baseline_std;  % [trials x time]
                            z_trace = mean(zscored_trials, 1);                             % mean z-scored trace
                            sem_trace = std(zscored_trials, 0, 1) / sqrt(size(zscored_trials, 1));  % SEM
                            
                            plot(bins, z_trace, 'Color', stimulus_colors(stim_idx,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%s (200x)', Stimulus_types{stim_idx}));
                               % Plot SEM shading
                            fill([bins, fliplr(bins)], ...
                                [z_trace + sem_trace, fliplr(z_trace - sem_trace)], ...
                                stimulus_colors(stim_idx,:), ...
                                'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');

                            if contains(z_score_period, 'entire_session')
                                ylabel('Z-scored FR (z-scored over entire session)', 'FontSize', 14);
                                ylim([-1 5]);
                            elseif contains(z_score_period, 'first30secs')
                                ylabel('Z-scored FR (z-scored over first 30s baseline)', 'FontSize', 14);
                                ylim([-1 8]);
                            end
                                                   
                        end    
                                       
                        peak_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                        mean_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                        peak_FR_by_greywindow = zeros(1, length(grey_window_starts)); % preallocate
                        mean_FR_by_greywindow = zeros(1, length(grey_window_starts)); % preallocate
                                        
                        for i = 1:length(stim_window_starts)
                            if contains(z_score_period, 'none')
                                % Find indices of bins within current window
                                idx_in_stimwindow = bins >= stim_window_starts(i) & bins < stim_window_ends(i);                    
                                peak_FR_by_stimwindow(i) = max(mean_trace(idx_in_stimwindow));
                                mean_FR_by_stimwindow(i) = mean(mean_trace(idx_in_stimwindow));
                                
                            else
                                % Find indices of bins within current window
                                idx_in_stimwindow = bins >= stim_window_starts(i) & bins < stim_window_ends(i);                    
                                peak_FR_by_stimwindow(i) = max(z_trace(idx_in_stimwindow));
                                mean_FR_by_stimwindow(i) = mean(z_trace(idx_in_stimwindow));
                                
                            end
                        end
                        
                        for i = 1:length(grey_window_starts)
                            if contains(z_score_period, 'none')
                                % Find indices of bins within current window
                                idx_in_greywindow = bins >= grey_window_starts(i) & bins < grey_window_ends(i);                    
                                peak_FR_by_greywindow(i) = max(mean_trace(idx_in_greywindow));
                                mean_FR_by_greywindow(i) = mean(mean_trace(idx_in_greywindow));
                                
                            else
                                % Find indices of bins within current window
                                idx_in_greywindow = bins >= grey_window_starts(i) & bins < grey_window_ends(i);                    
                                peak_FR_by_greywindow(i) = max(z_trace(idx_in_greywindow));
                                mean_FR_by_greywindow(i) = mean(z_trace(idx_in_greywindow));
                                
                            end
                        end
                        
                        session_idx = find(sessions_to_plot == nsession); % Map actual session number to index in preallocated array
                        all_peak_FR_by_stimwindow(session_idx, :) = peak_FR_by_stimwindow;
                        all_mean_FR_by_stimwindow(session_idx, :) = mean_FR_by_stimwindow;
                        all_peak_FR_by_greywindow(session_idx, :) = peak_FR_by_greywindow;
                        all_mean_FR_by_greywindow(session_idx, :) = mean_FR_by_greywindow;
                        
                        xlim([-0.5 1.5]);
                        xticks(-0.4:0.2:1.4);
                        xlabel('Time (s)', 'FontSize', 14)
                        set(gca, 'FontSize', 14);  % Tick labels font size
                        legend(flipud(findobj(gca,'-property','DisplayName')), 'Location', 'northeast', 'Interpreter', 'none');
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
                    if (contains(Stimulus_types{stim_idx}, 'GAVNIK DCBA')) % to label the plot only once
                        xline(0, 'k', (sprintf('A %d%s or D onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                        xline(0.15, 'k', (sprintf('B %d%s or C onset', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                        xline(0.30, 'k', (sprintf('C %d%s or B onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                        xline(0.45, 'k', (sprintf('D %d%s or A onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                        
                    end
                    sgtitle(sprintf('%s Aggregate L4 single unit activity %s vs %s', subject_number, Stimulus_types{1}, Stimulus_types{2}), 'Interpreter', 'none', 'FontSize', 16);
                    filename = sprintf('%s L4 FRs %s vs %s.png', subject_number, Stimulus_types{1}, Stimulus_types{2});
                    save_path = fullfile(pwd, filename);
                    exportgraphics(fig1, save_path); 
    
                    % Also save as .fig
                    fig_filename = sprintf('%s L4 FRs %s vs %s.fig', subject_number, Stimulus_types{1}, Stimulus_types{2});
                    fig_save_path = fullfile(pwd, fig_filename);
                    savefig(fig1, fig_save_path);
                end
            end    
        end            
        

        
        if contains(Stimulus_types{2}, 'GAVNIK_A_CD') % only make this plot if GAVNIK_A_CD is present
            if (contains(Stimulus_type, 'GAVNIK_ABCD') || contains(Stimulus_type, 'GAVNIK_A_CD') || contains(Stimulus_type, 'GAVNIK_E_CD')) && contains(plot_type, 'FR')
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
                    
                    ori = 1; %plot from onset of first stimulus in sequence
                    stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));
    
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.3 1.35], psthBinSize);   
                    mean_trace = mean(binnedArray, 1); % average over trials
    
                    figure(fig1);
                    hold on;
                    
                    if contains(z_score_period, 'none')
                        % Plot raw mean firing rate trace
                        plot(bins, mean_trace, 'Color', stimulus_colors(stim_idx,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%s (200x)', Stimulus_types{stim_idx}));
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
                             stimulus_colors(stim_idx,:), ...
                             'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                        
                        plot(bins, z_trace, 'Color', stimulus_colors(stim_idx,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%s (200x)', Stimulus_types{stim_idx}));
                        if contains(z_score_period, 'entire_session')
                            ylabel('Z-scored FR (z-scored over entire session)', 'FontSize', 14);
                            ylim([-1 5]);
                        elseif contains(z_score_period, 'first30secs')
                            ylabel('Z-scored FR (z-scored over first 30s baseline)', 'FontSize', 14);
                            ylim([-1 8]);
                        end                   
                    end    
                                       
                    peak_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                    mean_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                    peak_FR_by_greywindow = zeros(1, length(grey_window_starts)); % preallocate
                    mean_FR_by_greywindow = zeros(1, length(grey_window_starts)); % preallocate
                                    
                    for i = 1:length(stim_window_starts)
                        if contains(z_score_period, 'none')
                            % Find indices of bins within current window
                            idx_in_stimwindow = bins >= stim_window_starts(i) & bins < stim_window_ends(i);                    
                            peak_FR_by_stimwindow(i) = max(mean_trace(idx_in_stimwindow));
                            mean_FR_by_stimwindow(i) = mean(mean_trace(idx_in_stimwindow));
                            
                        else
                            % Find indices of bins within current window
                            idx_in_stimwindow = bins >= stim_window_starts(i) & bins < stim_window_ends(i);                    
                            peak_FR_by_stimwindow(i) = max(z_trace(idx_in_stimwindow));
                            mean_FR_by_stimwindow(i) = mean(z_trace(idx_in_stimwindow));
                            
                        end
                    end
                    
                    for i = 1:length(grey_window_starts)
                        if contains(z_score_period, 'none')
                            % Find indices of bins within current window
                            idx_in_greywindow = bins >= grey_window_starts(i) & bins < grey_window_ends(i);                    
                            peak_FR_by_greywindow(i) = max(mean_trace(idx_in_greywindow));
                            mean_FR_by_greywindow(i) = mean(mean_trace(idx_in_greywindow));
                            
                        else
                            % Find indices of bins within current window
                            idx_in_greywindow = bins >= grey_window_starts(i) & bins < grey_window_ends(i);                    
                            peak_FR_by_greywindow(i) = max(z_trace(idx_in_greywindow));
                            mean_FR_by_greywindow(i) = mean(z_trace(idx_in_greywindow));
                            
                        end
                    end
                    
                    session_idx = find(sessions_to_plot == nsession); % Map actual session number to index in preallocated array
                    all_peak_FR_by_stimwindow(session_idx, :) = peak_FR_by_stimwindow;
                    all_mean_FR_by_stimwindow(session_idx, :) = mean_FR_by_stimwindow;
                    all_peak_FR_by_greywindow(session_idx, :) = peak_FR_by_greywindow;
                    all_mean_FR_by_greywindow(session_idx, :) = mean_FR_by_greywindow;
                    
                    xlim([-0.5 1.5]);
                    ylim([-1 5]);
                    xticks(-0.4:0.2:1.4);
                    xlabel('Time (s)', 'FontSize', 14)
                    set(gca, 'FontSize', 14);  % Tick labels font size
                    legend(flipud(findobj(gca,'-property','DisplayName')), 'Location', 'northeast', 'Interpreter', 'none');
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
                if (contains(Stimulus_types{stim_idx}, 'GAVNIK_A_CD')) % to label the plot only once
                    xline(0, 'k', (sprintf('A %d%s or E novel onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                    xline(0.15, 'k', (sprintf('B %d%s or grey onset', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                    xline(0.30, 'k', (sprintf('C %d%s onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                    xline(0.45, 'k', (sprintf('D %d%s onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                    
                end
                sgtitle(sprintf('%s Aggregate L4 single unit activity %s vs %s vs %s', subject_number, Stimulus_types{1}, Stimulus_types{2}, Stimulus_types{3}), 'Interpreter', 'none', 'FontSize', 16);
                filename = sprintf('%s L4 FRs %s vs %s vs %s.png', subject_number, Stimulus_types{1}, Stimulus_types{2}, Stimulus_types{3});
                save_path = fullfile(pwd, filename);
                exportgraphics(fig1, save_path); 

                % Also save as .fig
                fig_filename = sprintf('%s L4 FRs %s vs %s vs %s.fig', subject_number, Stimulus_types{1}, Stimulus_types{2}, Stimulus_types{3});
                fig_save_path = fullfile(pwd, fig_filename);
                savefig(fig1, fig_save_path);

            end
        end
    end
end
hold off;

