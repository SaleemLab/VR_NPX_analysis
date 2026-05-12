%%% ERB 2025 For analysis of unit spiking in response to visual stimuli
% Depths of putative L5 and putative CA1 are based on best regional power per PSD analysis (for which, use PSD_analysis_UnProcessedLFP_ellie).
% L4 depth is identified via CSD analysis (for which, run CSD_Gratings_afterFILT_ellie or CSD_GAVNIKstims_afterFILT_ellie)

addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

%% setting metrics to screen good clusters
clear all
% 1/5 Choose your probe depth of interest and mouse
depth_for_analysis = 'V1'; % choose 'L4' or 'L2_3' or 'L5_6' or 'V1' or 'CA1' or 'Sub_CA1' or 'Sub_HPC'
SUBJECTS = {'M00069'};

params = create_cluster_selection_params('sorting_option','ellie');
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);

%%% 2/5
Stimulus_types = {'GAVNIK_ABCD', 'GAVNIK_A_CD', 'GAVNIK_E_CD'}; % IN ORDER with abcd first 'GAVNIK_A_CD', 'GAVNIK_E_CD', 'GAVNIK DCBA'
plot_type = 'FR'; % 'FR' firing rate or 'raster'
z_score_period = 'stim_session'; % 'none' = no z scoring or z score either over 'stim_session' (excludes variable greyscreen periods before and after stim paradigm), 'entire_session' or 'first30secs' (for every stimulus recording
% session from 20250205 onward, I presented grey screen to the mouse for at least 30s before starting the stimulus).   

cd('V:\Ellie\DATA\SUBJECTS\M00069\analysis\20251003') % 3/5 files will be saved here in the cd

%% SET THIS 4/5***** For NPX2.0 you will use a different L4 channel for each shank. Use CSD to estimate the best channel to use in L4
probe_type = 1; % NPX1.0 is type 0, NPX2.0 is type 1.

%5/5 row numbers of recording dates in "experiment_info" 15, 9
sessions_to_plot = [8]; 

if contains(Stimulus_types{2}, 'GAVNIK DCBA')
    stimulus_colors = [
        0.0,  0.0,  1.0;  % 1st group - blue - GAVNIK_ABCD
        1.0,  0.0,  0.0;];  % 2nd group - red - GAVNIK_DCBA
     stimulus_colors_raster = [
            0.0,  0.0,  0.6;  % 1st group - darker blue for raster - GAVNIK_ABCD
            0.6,  0.0,  0.0;];  % 2nd group - darker red for raster - GAVNIK_DCBA

elseif contains(Stimulus_types{2}, 'GAVNIK_A_CD')
    stimulus_colors = [
            0.0,  0.0,  1.0;  % 1st group - blue - GAVNIK_ABCD
            1.0,  0.0,  0.0;  % 2nd group - red - GAVNIK_A_CD
            0.0,  1.0,  0.0;];  % 3rd group - green - GAVNIK_E_CD  
    stimulus_colors_raster = [
            0.0,  0.0,  0.6;  % 1st group - darker blue for raster - GAVNIK_ABCD
            0.6,  0.0,  0.0;  % 2nd group - darker red for raster - GAVNIK_A_CD
            0.0,  0.5,  0.0;];  % 3rd group - darker green for raster - GAVNIK_E_CD
end

for nsession = sessions_to_plot 
    %Create a figure handle for each cluster outside the stimulus-type loop
    cluster_figures = containers.Map('KeyType', 'int32', 'ValueType', 'any'); %% to store figure handles
    cluster_raster_figures = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    
    all_z_traces_by_stim = cell(1, length(Stimulus_types)); % for heatmaps
    resp_A_ABCD = []; % to sort by strenght of response to A
    resp_B_ABCD = [];
    bins_ref = [];

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
            load(fullfile(options.ANALYSIS_DATAPATH, '..', 'earliest_V1sink_CSD.mat'))
            load(fullfile(options.ANALYSIS_DATAPATH, '..', 'depths_from_PSD.mat'))
            %% --- Normalise PSD struct for NPX1.0 to look like NPX2.0 ---
            if probe_type == 0
                % If PSD fields are not shank-wrapped, wrap them into shank1
                if ~isfield(depths_from_PSD, 'shank1')
                    tmp = depths_from_PSD;
                    depths_from_PSD = struct();
                    depths_from_PSD.shank1 = tmp;
                end
            end
            files = dir(fullfile(options.EPHYS_DATAPATH, '*ChanMap*.mat')); %the channel map has the y coordinate of each channel
            file_to_load = fullfile(options.EPHYS_DATAPATH, files(1).name); %load() does not accept wildcards like *
            load(file_to_load);

            %% NPX2.0 you will use a different L4 channel for each shank. Use CSD to estimate the best channel to use in L4
            if probe_type == 0
                layerfour_channels = earliest_V1sink_CSD.overall_best_halfmax_channel;
                layerfour_channels = layerfour_channels(:);
                shank_ids = ones(size(layerfour_channels)); % dummy shank ID = 1
    
                assert(isfield(earliest_V1sink_CSD, 'overall_best_halfmax_depth'), ...
                    'Missing overall_best_halfmax_depth in earliest_V1sink_CSD');
            
                earliest_V1sink_CSD = struct( ...
                    'shank_id', 1, ...
                    'best_channel_this_shank', layerfour_channels(1), ...
                    'best_depth_this_shank', earliest_V1sink_CSD(1).overall_best_halfmax_depth ...
                );
               
            elseif probe_type == 1
                % NPX2.0 → multiple shanks
                % earliest_V1sink_CSD is a struct array, extract fields safely
                layerfour_channels = [earliest_V1sink_CSD.best_channel_this_shank]; 
                shank_ids = [earliest_V1sink_CSD.shank_id]; 
                % ensure column vectors
                layerfour_channels = layerfour_channels(:); 
                shank_ids = shank_ids(:);
            end

            % --- Identify shanks and their layer 4 channels ---
            [unique_shanks, ia, ~] = unique(shank_ids, 'stable');
            best_channels_per_shank = layerfour_channels(ia); % one per shank
            
            % --- Preallocate shank-specific depth ranges ---
            L4_depth_range = cell(numel(unique_shanks), 1);
            V1_depth_range = cell(numel(unique_shanks), 1);
            supra_L4_depth_range = cell(numel(unique_shanks), 1);
            sub_L4_depth_range = cell(numel(unique_shanks), 1);
    
            CA1_depth_range = cell(numel(unique_shanks), 1);
            Sub_CA1_depth_range = cell(numel(unique_shanks), 1);
                   
            % --- Compute per-shank depth ranges using CSD and PSD data ---
            for iShank = 1:numel(unique_shanks)
                this_shank   = unique_shanks(iShank);
                this_channel = best_channels_per_shank(iShank);
                
                this_shank_name = ['shank', num2str(this_shank)];  % creates e.g. 'shank2'
                Brain_surface_depth = depths_from_PSD.(this_shank_name).surface_depth_PSD;
                L4_channel_depth    = earliest_V1sink_CSD(find([earliest_V1sink_CSD.shank_id] == this_shank, 1)).best_depth_this_shank;
                L5_depth            = depths_from_PSD.(this_shank_name).L5_depth_PSD;
                CA1_depth           = depths_from_PSD.(this_shank_name).CA1_depth_PSD;
            
                L4_depth_range{iShank}   = [L4_channel_depth - 60, L4_channel_depth + 60]; % giving electrodes the full extent of 120um inclusive 
                % (hence measuring spiking over depth range greater than 120um. As 60 is divisible by 15 and 10, this will give the same effective range for NPX1.0 (which has staggered electrodes every 10um down the shank) and NPX2.0 (which has electrode rows every 15um down each shank)
                V1_depth_range{iShank}   = [L5_depth - 330, L5_depth + 700];
                supra_L4_depth_range{iShank} = [L4_channel_depth + 65, L5_depth + 700]; % begins 5um above top of inferred L4 so no channels are double-counted on either NPX1.0 or NPX2.0
                sub_L4_depth_range{iShank} = [L5_depth - 330, L4_channel_depth - 65]; % begins 5um below bottom of inferred L4 so no channels are double-counted on either NPX1.0 or NPX2.0
    
                CA1_depth_range{iShank}  = [CA1_depth - 150, CA1_depth + 150];
                Sub_CA1_depth_range{iShank} = [min(CA1_depth_range{iShank}) - 1000, min(CA1_depth_range{iShank})];  
            end
            
            all_orientations = unique(Task_info.stim_orientation);
            %Task_info.stim_onset
            
           
            params = create_cluster_selection_params('sorting_option','ellie');
            %params.orientation_tuned = ...
            psthBinSize = 0.01; % but use 1ms for raster plots
            
            switch depth_for_analysis
                case 'L4' 
                    depth_ranges = L4_depth_range;
                case 'V1'
                    depth_ranges = V1_depth_range;
                case 'supraL4' 
                    depth_ranges = supra_L4_depth_range;
                case 'subL4' 
                    depth_ranges = sub_L4_depth_range;
    
                case 'CA1' 
                    depth_ranges = CA1_depth_range;
                case 'Sub_CA1'
                    depth_ranges = Sub_CA1_depth_range;
                case 'Sub_HPC'
                    depth_ranges = Sub_HPC_depth_range;
            end
                        
            if contains(z_score_period, 'stim_session') % excludes variable grey screen period before and after the stimulus paradigm ran
                time_edges = (min(Task_info.stim_onset) - 2):psthBinSize:(max(Task_info.stim_onset) + 2);
            end


            if contains(Stimulus_types{2}, 'GAVNIK DCBA') % only make this plot if GAVNIK DCBA is present
                if (contains(Stimulus_type, 'GAVNIK DCBA') || contains(Stimulus_type, 'GAVNIK_ABCD')) && contains(plot_type, 'FR')
                    for nprobe = 1:length(clusters)
                        selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters, which pass the set parameters
                        
                        depth_selected_clusters = selected_clusters; % initialize with same structure
                        for np = 1:length(clusters)
                            sc = selected_clusters(np);
                            peak_depths = clusters(np).peak_depth(sc.cluster_id); % get depths of selected clusters
        
                            % Find cluster_ids within selected depth range
                            keep_mask = peak_depths >= min(depth_ranges) & peak_depths <= max(depth_ranges);
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
        
                        for cluster_idx = 1:length(depth_cluster_ids)
                            cluster_id = depth_cluster_ids(cluster_idx);
                            
                            % Extract spike times for this cluster only
                            this_cluster_mask = depth_selected_clusters(nprobe).spike_id == cluster_id;
                            cluster_spike_times = depth_selected_clusters(nprobe).spike_times(this_cluster_mask);
                            
                            % Create traces figure for this cluster if not already created
                            if ~isKey(cluster_figures, cluster_id)
                                cluster_figures(cluster_id) = figure('Name', sprintf('Cluster %d', cluster_id), 'NumberTitle', 'off'); %% <<< MODIFIED >>>
                                figure(cluster_figures(cluster_id)); % activates the figure
                                hold on;
                            end
        
                            % Create raster figure for this cluster if not already created
                            if ~isKey(cluster_raster_figures, cluster_id)
                                raster_fig = figure('Name', sprintf('Cluster %d Raster', cluster_id), 'NumberTitle', 'off');
                                cluster_raster_figures(cluster_id) = raster_fig;
                                tiledlayout(raster_fig, 2, 1); % 2 vertically stacked panels
                            end
        
                            fig = cluster_figures(cluster_id);
                            figure(fig); % Activate traces figure
                        
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
                                                                    
                            if contains(z_score_period, 'none')
                                mean_trace = mean(binnedArray, 1) / psthBinSize; % average over trials
                                sem_trace = std(binnedArray, 0, 1) / sqrt(size(binnedArray, 1))/ psthBinSize;  % [1 x time]
                                % Plot SEM shading
                                fill([bins, fliplr(bins)], ...
                                     [mean_trace + sem_trace, fliplr(mean_trace - sem_trace)], ...
                                     stimulus_colors(stim_idx,:), ...
                                     'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                                hold on;
        
                                % Plot raw mean firing rate trace
                                plot(bins, mean_trace, 'Color', stimulus_colors(stim_idx,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%s (200x)', Stimulus_types{stim_idx}));
                                ylabel('Mean firing rate across trials (Hz)', 'FontSize', 14);
                                %ylim ([0 16]);
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
                                xline(0, 'k', (sprintf('A %d%s or D onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                                xline(0.15, 'k', (sprintf('B %d%s or C onset', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                                xline(0.30, 'k', (sprintf('C %d%s or B onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                                xline(0.45, 'k', (sprintf('D %d%s or A onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                                sgtitle(sprintf('%s %s cluster %d FRs %s vs %s', subject_number, depth_for_analysis, cluster_id, Stimulus_types{1}, Stimulus_types{2}), 'Interpreter', 'none', 'FontSize', 16); 
                            end    
                            if stim_idx == 2 %only save once as .fig    
                                fig_filename = sprintf('%s %s cluster %d FRs %s vs %s (zscore %s).fig', subject_number, depth_for_analysis, cluster_id, Stimulus_types{1}, Stimulus_types{2}, z_score_period);
                                fig_save_path = fullfile(pwd, fig_filename);
                                savefig(fig, fig_save_path);
                            end    
                            
                            % Activate raster figure and tile
                            raster_fig = cluster_raster_figures(cluster_id);
                            figure(raster_fig); % activate figure
                            t = findobj(raster_fig, 'Type', 'tiledlayout');
                            nexttile(t, stim_idx); % Use tile based on stim_idx
                            ax = gca;
        
                            % Plot the raster
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(cluster_spike_times, stim_onsets, [-0.3 1.35], psthBinSize/10); 
                            plot(ax, rasterX, rasterY, 'Color', stimulus_colors_raster(stim_idx,:), 'LineWidth', 3, 'HandleVisibility', 'off'); %  
                            hold on;
        
                            
                            
                            if stim_idx == 1 %only label (with standard ABCD orientations!) once
                                xline(ax, 0, 'k', (sprintf('A %d%s onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                                xline(ax, 0.15, 'k', (sprintf('B %d%s onset', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                                xline(ax, 0.30, 'k', (sprintf('C %d%s onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                                xline(ax, 0.45, 'k', (sprintf('D %d%s onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                            end    
                            
                            if stim_idx == 2 %only label (with standard ABCD orientations!) once
                                xline(ax, 0, 'k', (sprintf('D %d%s onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                                xline(ax, 0.15, 'k', (sprintf('C %d%s onset', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                                xline(ax, 0.30, 'k', (sprintf('B %d%s onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                                xline(ax, 0.45, 'k', (sprintf('A %d%s onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                            end  
        
                            sgtitle(sprintf('%s %s cluster %d rasters %s vs %s', subject_number, depth_for_analysis, cluster_id, Stimulus_types{1}, Stimulus_types{2}), 'Interpreter', 'none', 'FontSize', 16);
                            ylabel('Trial', 'FontSize', 14);
                            title(Stimulus_types{stim_idx}, 'Interpreter', 'none', 'FontSize', 14);
                            set(gca, 'FontSize', 14, 'TickDir', 'out', 'box', 'off');
        
                            if stim_idx == 2 %only save once as .fig    
                                fig_filename = sprintf('%s %s cluster %d rasters %s vs %s.fig', subject_number, depth_for_analysis, cluster_id, Stimulus_types{1}, Stimulus_types{2});
                                fig_save_path = fullfile(pwd, fig_filename);
                                savefig(raster_fig, fig_save_path);
                            end  
                        end
                    end
                end    
            end







            if contains(Stimulus_types{2}, 'GAVNIK_A_CD') % only make this plot if GAVNIK_A_CD is present
                %if (contains(Stimulus_type, 'GAVNIK_ABCD') || contains(Stimulus_type, 'GAVNIK_A_CD') || contains(Stimulus_type, 'GAVNIK_E_CD')) && contains(plot_type, 'FR')
                for nprobe = 1:length(clusters)
                    selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters, which pass the set parameters
                    
                    depth_selected_clusters = selected_clusters; % initialize with same structure
                    for np = 1:length(clusters)
                        sc = selected_clusters(np);
                        
                        cluster_channels = clusters(np).peak_channel(sc.cluster_id);
                        cluster_depths = ycoords(cluster_channels); % get depths of selected clusters from ycoords (peak_depths from SI are 15 microns different..)
        
                        % Map each cluster to its shank
                        cluster_shanks = kcoords(cluster_channels);
        
                        % Initialize container for depth-filtered cluster IDs
                        depth_cluster_ids = [];
                
                        % Loop over each shank
                        for iShank = 1:numel(unique_shanks)
                            this_shank = unique_shanks(iShank);
                            this_range = depth_ranges{iShank};
                
                            % Logical mask for clusters on this shank and within this shank's depth range
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
                    
                    baseline_window = [0 30]; % in seconds - first 30s of recording is grey screen - can use for z-scoring
    
                    % Define time_edges depending on z_score_period
                    if contains(z_score_period, 'entire_session')
                        time_edges = 0:psthBinSize:max(depth_selected_clusters(nprobe).spike_times); %time bin edges from 0 to the max time in steps of ptshBinSize
                    elseif contains(z_score_period, 'first30secs')
                        time_edges = baseline_window(1):psthBinSize:baseline_window(2);
                    end
    

                    % Initialise storage of per neuron activity for response heatmap
                    all_z_traces = [];   % [neurons x time] 

                    for cluster_idx = 1:length(depth_cluster_ids)
                        cluster_id = depth_cluster_ids(cluster_idx);
                        
                        % Extract spike times for this cluster only
                        this_cluster_mask = depth_selected_clusters(nprobe).spike_id == cluster_id;
                        cluster_spike_times = depth_selected_clusters(nprobe).spike_times(this_cluster_mask);
                        
                        % Create traces figure for this cluster if not already created
                        if ~isKey(cluster_figures, cluster_id)
                            cluster_figures(cluster_id) = figure('Name', sprintf('Cluster %d', cluster_id), 'NumberTitle', 'off'); %% <<< MODIFIED >>>
                            figure(cluster_figures(cluster_id)); % activates the figure
                            hold on;
                        end
    
                        % Create raster figure for this cluster if not already created
                        if ~isKey(cluster_raster_figures, cluster_id)
                            raster_fig = figure('Name', sprintf('Cluster %d Raster', cluster_id), 'NumberTitle', 'off');
                            cluster_raster_figures(cluster_id) = raster_fig;
                            tiledlayout(raster_fig, 3, 1); % 3 vertically stacked panels
                        end
    
                        fig = cluster_figures(cluster_id);
                        figure(fig); % Activate traces figure
                    
                        % Compute appropriate histogram counts for z-scoring
                        if contains(z_score_period, 'entire_session') || contains(z_score_period, 'stim_session')
                            zscore_counts = histcounts(cluster_spike_times, time_edges);
                        elseif contains(z_score_period, 'first30secs')
                            baseline_spikes = cluster_spike_times(cluster_spike_times >= baseline_window(1) & cluster_spike_times <= baseline_window(2));
                            zscore_counts = histcounts(baseline_spikes, time_edges); %histcounts counts the number of datapoints in specified time bins
                        end
        
                        ordered_oris = unique(Task_info.stim_orientation, 'stable');               
                        
                        ori = 1; %plot from onset of first stimulus in sequence
                        stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));
        
                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(cluster_spike_times, stim_onsets, [-0.3 1.35], psthBinSize);   
                                                                
                        if contains(z_score_period, 'none')
                            mean_trace = mean(binnedArray, 1) / psthBinSize; % average over trials
                            sem_trace = std(binnedArray, 0, 1) / sqrt(size(binnedArray, 1))/ psthBinSize;  % [1 x time]
                            % Plot SEM shading
                            fill([bins, fliplr(bins)], ...
                                 [mean_trace + sem_trace, fliplr(mean_trace - sem_trace)], ...
                                 stimulus_colors(stim_idx,:), ...
                                 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                            hold on;
    
                            % Plot raw mean firing rate trace
                            plot(bins, mean_trace, 'Color', stimulus_colors(stim_idx,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%s (200x)', Stimulus_types{stim_idx}));
                            ylabel('Mean firing rate across trials (Hz)', 'FontSize', 14);
                            %ylim ([0 16]);
                        else 
                            baseline_mean = mean(zscore_counts);
                            baseline_std = std(zscore_counts);
                            % Z-score each trial individually so that trial to variability remains visible and SEM can be calculated
                            zscored_trials = (binnedArray - baseline_mean) / baseline_std;  % [trials x time]
                            z_trace = mean(zscored_trials, 1);                             % mean z-scored trace
                            
                            A_ABCD_window = bins >= 0.03 & bins < 0.08; % window for A onset responses sorting
                            if contains(Stimulus_types{1}, 'GAVNIK_ABCD')
                                B_ABCD_window = bins >= 0.18 & bins < 0.23; % window for B onset responses sorting
                            elseif contains(Stimulus_types{1}, 'GAVNIK250_ABCD')
                                B_ABCD_window = bins >= 0.28 & bins < 0.33; % window for B onset responses sorting
                            end

                            if contains(Stimulus_type, '_ABCD')
                                resp_A_ABCD(end+1) = mean(z_trace(A_ABCD_window));
                                resp_B_ABCD(end+1) = mean(z_trace(B_ABCD_window));
                            end
                                                        
                            if isempty(bins_ref)
                                bins_ref = bins;
                            end
                            all_z_traces(end+1, :) = z_trace; % store for responses heatmap
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
                            elseif contains(z_score_period, 'stim_session')
                                ylabel('Z-scored FR (z-scored over stim session)', 'FontSize', 14);
                                ylim([-1 5]);
                            elseif contains(z_score_period, 'first30secs')
                                ylabel('Z-scored FR (z-scored over first 30s baseline)', 'FontSize', 14);
                                ylim([-1 8]);
                            end
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
                        if stim_idx == 3 %only save once as .fig    
                            fig_filename = sprintf('%s %s cluster %d FRs %s vs %s vs %s (zscore %s).fig', subject_number, depth_for_analysis, cluster_id, Stimulus_types{1}, Stimulus_types{2}, Stimulus_types{3}, z_score_period);
                            fig_save_path = fullfile(pwd, fig_filename);
                            %savefig(fig, fig_save_path);
                        end    
                        
                        % Activate raster figure and tile
                        raster_fig = cluster_raster_figures(cluster_id);
                        figure(raster_fig); % activate figure
                        t = findobj(raster_fig, 'Type', 'tiledlayout');
                        nexttile(t, stim_idx); % Use tile based on stim_idx
                        ax = gca;
    
                        % Plot the raster
                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(cluster_spike_times, stim_onsets, [-0.3 1.35], psthBinSize/10); 
                        plot(ax, rasterX, rasterY, 'Color', stimulus_colors_raster(stim_idx,:), 'LineWidth', 3, 'HandleVisibility', 'off'); %  
                        hold on;
    
                        %look at the distribution of spikes in the omission window i.e. 30-150ms after onset of element 2 
                        analysis_window = [0.18 0.3];  % seconds after onset of element 1 (which comes on 0.15s before element 2)
                        trial_means = [];
                        
                        % Loop over each stimulus onset
                        for i = 1:length(stim_onsets)
                            onset = stim_onsets(i);
                            % Get spikes within the window for this trial
                            spikes_in_trial = cluster_spike_times(cluster_spike_times >= onset + analysis_window(1) & ...
                                                                  cluster_spike_times <= onset + analysis_window(2));
                            spikes_in_trial = spikes_in_trial - onset;
    
                            if ~isempty(spikes_in_trial)
                                trial_means(end+1) = mean(spikes_in_trial);
                            end
                        end
                        
                        % Only compute if spikes were found
                        if ~isempty(trial_means)
                            mean_spike_time = mean(trial_means);
                            sem_spike_time = std(trial_means) / sqrt(length(trial_means));  % SEM across trials; n = # of trials with spikes
                            ci99 = 2.576 * sem_spike_time;  % 99% CI
                            % Add a marker for the mean spike time 
                            %y_val = max(rasterY)/2;  % Place it in the middle of the raster
                           
                            hold on;
                            
                            %errorbar(ax, mean_spike_time, y_val, ci99, 'horizontal', ...
                             %   'o', ...
                             %   'Color', 'k', ...
                             %   'CapSize', 10, ...
                             %   'LineWidth', 1.0, ...
                             %   'MarkerSize', 5, ...
                             %   'MarkerFaceColor', stimulus_colors_raster(stim_idx,:), ...
                             %   'MarkerEdgeColor', stimulus_colors_raster(stim_idx,:), ...
                             %   'DisplayName', 'mean trial spike time [0.18–0.3s] 99% CI');
                        end
                        %legend show
                        %legend('Location', 'southwest');
                        xlim(ax, [-0.2 1.2]); % to highlight post stim oscillations
                        xticks(ax, 0:0.2:1.2); % to highlight post stim oscillations
                        %xlim(ax, [0.14 0.31]); % to highlight response to second element
                       % xticks(ax, 0.15:0.01:0.3); % to highlight response to second element
                        
                        xline(ax, 0.30, 'k', (sprintf('C %d%s onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                        xline(ax, 0.45, 'k', (sprintf('D %d%s onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                        
                        if stim_idx == 1 %only label (with standard ABCD orientations!) once
                            xline(ax, 0, 'k', (sprintf('A %d%s onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                            xline(ax, 0.15, 'k', (sprintf('B %d%s onset', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                        end    
                        
                        if stim_idx == 2 %only label (with standard ABCD orientations!) once
                            xline(ax, 0, 'k', (sprintf('A %d%s onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                            xline(ax, 0.15, 'k', (sprintf('grey onset')), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                        end  
    
                        if stim_idx == 3 %only label (with standard ABCD orientations!) once
                            xline(ax, 0, 'k', (sprintf('E %d%s onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                            xline(ax, 0.15, 'k', (sprintf('grey onset')), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 14);
                            xlabel('Time (s) since onset of first element', 'FontSize', 14);
                        end  
                        sgtitle(sprintf('%s %s cluster %d rasters %s vs %s vs %s', subject_number, depth_for_analysis, cluster_id, Stimulus_types{1}, Stimulus_types{2}, Stimulus_types{3}), 'Interpreter', 'none', 'FontSize', 16);
                        ylabel('Trial', 'FontSize', 14);
                        title(Stimulus_types{stim_idx}, 'Interpreter', 'none', 'FontSize', 14);
                        set(gca, 'FontSize', 14, 'TickDir', 'out', 'box', 'off');
    
                        if stim_idx == 3 %only save once as .fig    
                            fig_filename = sprintf('%s %s cluster %d rasters %s vs %s vs %s.fig', subject_number, depth_for_analysis, cluster_id, Stimulus_types{1}, Stimulus_types{2}, Stimulus_types{3});
                            fig_save_path = fullfile(pwd, fig_filename);
                            %savefig(raster_fig, fig_save_path);
                        end  
                    end
                    
                    all_z_traces_by_stim{stim_idx} = all_z_traces; % store traces for this stim type
                
                end 
            end
        end
    end
    [~, sort_idx_A] = sort(resp_A_ABCD, 'ascend');
    [~, sort_idx_B] = sort(resp_B_ABCD, 'ascend');

    figure;
    sgtitle(sprintf('%s %s cluster responses, sorted by response to A in GAVNIK_ABCD', subject_number, depth_for_analysis), 'Interpreter', 'none', 'FontSize', 16);
    
    for stim_idx = 1:length(Stimulus_types)
    
        subplot(1, length(Stimulus_types), stim_idx);
    
        traces = all_z_traces_by_stim{stim_idx};
        sorted_traces = traces(sort_idx_A, :);
    
        imagesc(bins_ref, 1:size(sorted_traces,1), sorted_traces);
        axis xy;
    
        title(Stimulus_types{stim_idx}, 'Interpreter','none');
        xlabel('Time (s)');
        ylabel('V1 Clusters');
    
        xlim([-0.3 1.35]);
    
        % Stimulus timing
        xline(0, '--k');
        xline(0.15, '--k');
        xline(0.30, '--k');
        xline(0.45, '--k');

        % Set color limits asymmetrically around zero
        minVal = -1; %z-scored FR doesn't typically go below -1
        maxVal = 2; % z-scored FR s rarely go above 4 in this dataset
        n = 256; % number of colours in colormap - standard resolution for a smooth gradient
        
        % Position of zero in colormap
        zero_pos = round(n * (0 - minVal) / (maxVal - minVal));
        
        % Blue → white (NEGATIVE values ONLY)
        blue_part = interp1([0 1], [0 0 1; 1 1 1], linspace(0,1,zero_pos));
        
        % White → red (POSITIVE values ONLY)
        red_part  = interp1([0 1], [1 1 1; 1 0 0], linspace(0,1,n-zero_pos+1));
        
        % Remove duplicate white row
        red_part = red_part(2:end,:);
        
        cmi = [blue_part; red_part];
        caxis([minVal maxVal]);
        colormap(cmi);

        % Colorbar
        cb = colorbar;
        ylabel(cb, 'Z-scored firing rate');
    
    end    
       
   

    figure;
    sgtitle(sprintf('%s %s cluster responses, sorted by response to B in GAVNIK_ABCD', subject_number, depth_for_analysis), 'Interpreter', 'none', 'FontSize', 16);
    
    for stim_idx = 1:length(Stimulus_types)
    
        subplot(1, length(Stimulus_types), stim_idx);
    
        traces = all_z_traces_by_stim{stim_idx};
        sorted_traces = traces(sort_idx_B, :);
    
        imagesc(bins_ref, 1:size(sorted_traces,1), sorted_traces);
        axis xy;
    
        title(Stimulus_types{stim_idx}, 'Interpreter','none');
        xlabel('Time (s)');
        ylabel('V1 Clusters');
    
        xlim([-0.3 1.35]);
    
        % Stimulus timing
        xline(0, '--k');
        xline(0.15, '--k');
        xline(0.30, '--k');
        xline(0.45, '--k');

        % Set color limits asymmetrically around zero
        minVal = -1; %z-scored FR doesn't typically go below -1
        maxVal = 2; % z-scored FR s rarely go above 4 in this dataset
        n = 256; % number of colours in colormap - standard resolution for a smooth gradient
        
        % Position of zero in colormap
        zero_pos = round(n * (0 - minVal) / (maxVal - minVal));
        
        % Blue → white (NEGATIVE values ONLY)
        blue_part = interp1([0 1], [0 0 1; 1 1 1], linspace(0,1,zero_pos));
        
        % White → red (POSITIVE values ONLY)
        red_part  = interp1([0 1], [1 1 1; 1 0 0], linspace(0,1,n-zero_pos+1));
        
        % Remove duplicate white row
        red_part = red_part(2:end,:);
        
        cmi = [blue_part; red_part];
        caxis([minVal maxVal]);
        colormap(cmi);

        % Colorbar
        cb = colorbar;
        ylabel(cb, 'Z-scored firing rate');
    end       
    
end
