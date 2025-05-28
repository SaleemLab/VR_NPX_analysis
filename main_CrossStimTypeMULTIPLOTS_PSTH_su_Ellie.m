%%% For analysis of unit spiking in response to visual stimuli

addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

%% setting metrics to screen good clusters
clear all
% Choose your probe depth of interest
session_specific_L4 = {15, 4510:4650;}; % 1/4. um. Set for each SESSION based on CSD +/- 70um
% 4, 4440:4580; 5, 4440:4580; 6, 4460:4600; 7, 4480:4620; 8, 4500:4640; 
% 11, 4500:4640; 12, 4510:4650; 13, 4510:4650; 14, 4510:4650; 15, 4510:4650;
depth_for_analysis = 'L5_6'; % choose 'L4' or 'L2_3' (max(L4_depth_range) + 500) or 'L5_6' (min(L4_depth_range) - 400)
SUBJECTS = {'M00013'};

params = create_cluster_selection_params('sorting_option','ellie');
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);

%%% 2/4
Stimulus_types = {'GAVNIK_ABCD', 'GAVNIK_A_CD', 'GAVNIK_E_CD'}; % IN ORDER 'GAVNIK_A_CD', 'GAVNIK_E_CD', 'GAVNIK DCBA'
plot_type = 'FR'; % 'FR' firing rate or 'raster'
z_score_period = 'entire_session'; % 'none' = no z scoring or z score either over 'entire_session' or 'first30secs' (for every stimulus recording
% session from 20250205 onward, I presented grey screen to the mouse for at least 30s before starting the stimulus).
cd('V:\Ellie\DATA\SUBJECTS\M00013\analysis\20250221') % 3/4 files will be saved here in the cd

sessions_to_plot = [15]; %4/4 row numbers of recording dates in "experiment_info" 11, 12, 13, 14, 15

stimulus_colors = [
        0.0,  0.0,  1.0;  % 1st group - blue - GAVNIK_ABCD
        1.0,  0.0,  0.0;  % 2nd group - red - GAVNIK_A_CD
        0.0,  1.0,  0.0;];  % 3rd group - green - GAVNIK_E_CD  


for nsession = sessions_to_plot 
    %Create a figure handle for each cluster outside the stimulus-type loop
    cluster_figures = containers.Map('KeyType', 'int32', 'ValueType', 'any'); %% to store figure handles

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
            
            if strcmp(depth_for_analysis, 'L4') % choose 'L4' or 'L2/3' (max(L4_depth_range) + 500) or 'L5/6' (min(L4_depth_range) - 400)
                depth_range = L4_depth_range;
            elseif strcmp(depth_for_analysis, 'L2_3')
                depth_range = (max(L4_depth_range) : max(L4_depth_range) + 500);
            elseif strcmp(depth_for_analysis, 'L5_6')
                depth_range = (min(L4_depth_range) - 400 : min(L4_depth_range));
            end
            
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
                %all_spike_times = [];
                %all_spike_ids = [];
                %for np = 1:length(depth_selected_clusters)
                  %  all_spike_times = [all_spike_times; depth_selected_clusters(np).spike_times];
                  %  all_spike_ids = [all_spike_ids; depth_selected_clusters(np).spike_id];
                %end
                for cluster_idx = 1:length(depth_cluster_ids)
                    cluster_id = depth_cluster_ids(cluster_idx);
                    
                    % Extract spike times for this cluster only
                    this_cluster_mask = depth_selected_clusters(nprobe).spike_id == cluster_id;
                    cluster_spike_times = depth_selected_clusters(nprobe).spike_times(this_cluster_mask);
                    
                    % Create figure for this cluster if not already created
                    if ~isKey(cluster_figures, cluster_id)
                        cluster_figures(cluster_id) = figure('Name', sprintf('Cluster %d', cluster_id), 'NumberTitle', 'off'); %% <<< MODIFIED >>>
                        figure(cluster_figures(cluster_id)); % activates the figure
                        hold on;
                    end
                    fig = cluster_figures(cluster_id);
                    figure(fig); % Activate figure

                
                
                    % Compute appropriate histogram counts for z-scoring
                    if contains(z_score_period, 'entire_session')
                        zscore_counts = histcounts(cluster_spike_times, time_edges);
                    elseif contains(z_score_period, 'first30secs')
                        baseline_spikes = cluster_spike_times(cluster_spike_times >= baseline_window(1) & cluster_spike_times <= baseline_window(2));
                        zscore_counts = histcounts(baseline_spikes, time_edges); %histcounts counts the number of datapoints in specified time bins
                    end
    
    
                    ordered_oris = unique(Task_info.stim_orientation, 'stable');               
                    
                    ori = 1; %plot from onset of first stimulus in sequence
                    stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));
    
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(cluster_spike_times, stim_onsets, [-0.3 1.35], psthBinSize);   
                    mean_trace = mean(binnedArray, 1); % average over trials
    
                                        
                    if contains(z_score_period, 'none')
                        % Plot raw mean firing rate trace
                        plot(bins, mean_trace, 'Color', stimulus_colors(stim_idx,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%s (200x)', Stimulus_types{stim_idx}));
                        ylabel('Mean firing rate (Hz)', 'FontSize', 14);
                        ylim ([0 16]);
                    else % Compute appropriate histogram counts for z-scoring
                        if contains(z_score_period, 'entire_session')
                            zscore_counts = histcounts(cluster_spike_times, time_edges);
                            ylabel('Z-scored FR (z-scored over entire session)', 'FontSize', 14);
                            ylim([-1 5]);
                        elseif contains(z_score_period, 'first30secs')
                            baseline_spikes = cluster_spike_times(cluster_spike_times >= baseline_window(1) & cluster_spike_times <= baseline_window(2));
                            zscore_counts = histcounts(baseline_spikes, time_edges); %histcounts counts the number of datapoints in specified time bins
                            ylabel('Z-scored FR (z-scored over first 30s baseline)', 'FontSize', 14);
                            ylim([-1 8]);
                        end
                        z_trace = (mean_trace - mean(zscore_counts)) / std(zscore_counts); %mean(zscore_counts) gives the mean spikes per timebin in the reference period; this is then deducted from the spikecount of each trial-averaged timebin
                        plot(bins, z_trace, 'Color', stimulus_colors(stim_idx,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%s (200x)', Stimulus_types{stim_idx}));         
                    end   
    
                    xlim([-0.5 1.5]);
                    xticks(-0.4:0.2:1.4);
                    xlabel('Time (s)', 'FontSize', 14)
                    set(gca, 'FontSize', 14);  % Tick labels font size
                    legend(flipud(findobj(gca,'-property','DisplayName')), 'Location', 'northeast', 'Interpreter', 'none');
                    hold on;                
                
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
                    
                    if stim_idx == 1 %only label (with standard ABCD orientations!) once
                        xline(0, 'k', (sprintf('A %d%s or E novel onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                        xline(0.15, 'k', (sprintf('B %d%s or grey onset', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                        xline(0.30, 'k', (sprintf('C %d%s onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                        xline(0.45, 'k', (sprintf('D %d%s onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                    
                        sgtitle(sprintf('%s %s cluster %d FRs %s vs %s vs %s', subject_number, depth_for_analysis, cluster_id, Stimulus_types{1}, Stimulus_types{2}, Stimulus_types{3}), 'Interpreter', 'none', 'FontSize', 16);
                    end    
                    if stim_idx == 3 %only save once    
                        filename = sprintf('%s %s cluster %d FRs %s vs %s vs %s.png', subject_number, depth_for_analysis, cluster_id, Stimulus_types{1}, Stimulus_types{2}, Stimulus_types{3});
                        save_path = fullfile(pwd, filename);
                        exportgraphics(fig, save_path); 
        
                        % Also save as .fig
                        fig_filename = sprintf('%s %s cluster %d FRs %s vs %s vs %s.fig', subject_number, depth_for_analysis, cluster_id, Stimulus_types{1}, Stimulus_types{2}, Stimulus_types{3});
                        fig_save_path = fullfile(pwd, fig_filename);
                        savefig(fig, fig_save_path);
                    end    
                end
            end
        end
    end
end
