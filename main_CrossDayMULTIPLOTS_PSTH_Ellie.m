%%% ERB 2025 For analysis of unit spiking in response to visual stimuli
% Depths of putative L5 and putative CA1 are based on best regional power per PSD analysis (for which, use PSD_analysis_UnProcessedLFP_ellie).
% L4 depth is identified via CSD analysis (for which, run CSD_Gratings_afterFILT_ellie or CSD_GAVNIKstims_afterFILT_ellie)

addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

%% setting metrics to screen good clusters
clear all
% Choose your probe depth of interest
depth_for_analysis = 'sub_L4'; % choose 'L4' or 'V1' or 'supra_L4' or 'sub_L4', 'CA1' or 'Sub_CA1' or 'Sub_HPC'

SUBJECTS = {'M00069'};

params = create_cluster_selection_params('sorting_option','ellie');
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);

%%% 2/6
Stimulus_type = 'GAVNIK_ABCD'; %%%% BUT IF THERE ARE MULTIPLE TYPES E.G. GAVNIK_ABCD_1, also NEED TO SPECIFY THIS IN 5/6
plot_choice = 'aggregate'; % curated 'single_units' or in 'aggregate' or uncurated 'MUA'; MUA includes all clusters from kilosort, unfiltered
plot_type = 'FR'; % 'FR' firing rate or 'raster'
neuron_type = 'All'; % for GAVNIK stimuli (coded for so far) set this to 'PYR only' if you want to include only putative pyramidal neurons. distinguish between putative PYR (wide waveform, lower tau rise than SOM), PV (narrow waveform) and SOM (wide waveform, higher tau rise in the ACG i.e. probability of spiking again increases quite slowly)
z_score_period = 'entire_session'; % 'none' = no z scoring or z score either over 'entire_session' or 'first30secs' (for every stimulus recording
% session from 20250205 onward, I presented grey screen to the mouse for at least 30s before starting the stimulus).

% 3/6 files will be saved here in the cd
cd('V:\Ellie\DATA\SUBJECTS\M00069\analysis') % 3/5 files will be saved here in the cd

%% SET THIS 4/6***** For NPX2.0 you will use a different L4 channel for each shank. Use CSD to estimate the best channel to use in L4
probe_type = 1; % NPX1.0 is type 0, NPX2.0 is type 1.

%%% 5/6
sessions_to_plot = [4, 5, 6, 7, 8]; %5/6 row numbers of recording dates in "experiment_info" 4, 5, 6, 7, 8   11, 12, 13, 14, 15
stimulus_choice = containers.Map( ...
    [4 5 6 7 8], ...
    {'GAVNIK_ABCD_1', 'GAVNIK_ABCD_1', 'GAVNIK_ABCD_1', 'GAVNIK_ABCD_1', 'GAVNIK_ABCD'} ...
);

colors = [
    1.0,  0.8,  0.0;  % 1st group - golden yellow (1,0, 1.0, 0.0 for pure yellow)
    0.4,  0.8,  0.0;  % 2nd group - yellow-green (more yellowish)
    0.2,  0.6,  0.4;  % 3rd group - pure green
    0.2,  0.4,  0.6;  % 4th group - greenish blue (not too blue)
    0.0,  0.0,  0.8;   % 5th group - deep blue
    %0.4,  0.2,  0.8;  % 6th group - blueish purple
    %0.8,  0.2,  0.6;  % 7th group - reddish purple
];

% 6/6 Define windows after time zero to calc. the peak and mean FR - Depends on stimulus type              
if contains(Stimulus_type, 'GAVNIK_')
    stim_window_starts = 0:0.15:0.45;
    stim_window_ends = 0.15:0.15:0.6;
    grey_window_starts = 0.65; % look from 50ms-150ms after stim D offset
    grey_window_ends = 0.75; % look from 50ms-150ms after stim D offset
elseif contains(Stimulus_type, 'GAVNIK250_ABCD')    
    stim_window_starts = 0:0.25:0.75;
    stim_window_ends = 0.25:0.25:1.0;
    grey_window_starts = 1.05; % look from 50ms-250ms after stim D offset
    grey_window_ends = 1.25; % look from 50ms-250ms after stim D offset
elseif contains(Stimulus_type, 'E_1000ms') || contains(Stimulus_type, 'A_1000ms')
    stim_window_starts = 0;
    stim_window_ends = 1;
    grey_window_starts = 1; 
    grey_window_ends = 2;   
elseif contains(Stimulus_type, 'TRAIN250') 
    stim_window_starts = 0:0.5:1.5;
    stim_window_ends = 0.25:0.5:1.75;
    grey_window_starts = [0.3 0.8 1.3 1.8]; % look from 50ms-250ms after stim offset
    grey_window_ends = [0.5 1.0 1.5 2.0]; % look from 50ms-250ms after stim offset
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
    %session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    %stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stim_names = experiment_info(nsession).StimulusName;

    desired_stim = stimulus_choice(nsession);

    match_idx = strcmp(stim_names, desired_stim);

    assert(nnz(match_idx) == 1, ...
        'Session %d: expected exactly one match for %s', ...
        nsession, desired_stim);

    session_info  = experiment_info(nsession).session(match_idx);
    stimulus_name = stim_names(match_idx);

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
        
            L4_depth_range{iShank}   = [L4_channel_depth - 70, L4_channel_depth + 70];
            V1_depth_range{iShank}   = [L5_depth - 330, L5_depth + 700];
            supra_L4_depth_range{iShank} = [L4_channel_depth + 70, L5_depth + 700];
            sub_L4_depth_range{iShank} = [L5_depth - 330, L4_channel_depth - 70];

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
            case 'supra_L4' 
                depth_ranges = supra_L4_depth_range;
            case 'sub_L4' 
                depth_ranges = sub_L4_depth_range;

            case 'CA1' 
                depth_ranges = CA1_depth_range;
            case 'Sub_CA1'
                depth_ranges = Sub_CA1_depth_range;
            case 'Sub_HPC'
                depth_ranges = Sub_HPC_depth_range;
        end
                

        if contains(Stimulus_type, 'TRAIN') || contains(Stimulus_type, 'E_1000ms') || contains(Stimulus_type, 'A_1000ms')...
                && contains(plot_choice, 'aggregate') && contains(plot_type, 'FR')        
            
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
                if(contains(Stimulus_type, 'E_1000ms') || contains(Stimulus_type, 'A_1000ms'))
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.30 2], psthBinSize);    
                elseif contains(Stimulus_type, 'TRAIN250')
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.30 2.8], psthBinSize); 
                else
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.30 1.8], psthBinSize);    
                end    
                mean_trace = mean(binnedArray, 1); % average over trials
                figure(fig1);
                hold on;

                if contains(z_score_period, 'none')
                    % Plot raw mean firing rate trace
                    plot(bins, mean_trace, 'Color', colors(find(sessions_to_plot == nsession),:), 'LineWidth', 1.5, 'DisplayName', sprintf('Day %d, %s', experiment_info(nsession).date, depth_for_analysis));
                    ylabel('Mean firing rate (Hz)', 'FontSize', 24);
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

                    plot(bins, z_trace, 'Color', colors(find(sessions_to_plot == nsession),:), 'LineWidth', 1.5, 'DisplayName', sprintf('Day %d, %s', experiment_info(nsession).date, depth_for_analysis));

                    if contains(z_score_period, 'entire_session')
                        ylabel('Z-scored FR (z-scored over entire session)', 'FontSize', 24);
                        ylim([-1 4]);
                    elseif contains(z_score_period, 'first30secs')
                        ylabel('Z-scored FR (z-scored over first 30s baseline)', 'FontSize', 24);
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
                
                if contains(Stimulus_type, 'TRAIN250')
                    xlim([-0.3 2.8]);
                    xticks(-0.4:0.2:2.8);
                else    
                    xlim([-0.3 2.0]);
                    xticks(-0.4:0.2:2.0);
                end    
                set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);
                xlabel('Time (s) since onset of A', 'FontSize', 24)
                
                legend(flipud(findobj(gca,'-property','DisplayName')), 'Location', 'northeast', 'FontSize', 24);
                hold on;                
            end
                    % Define grey intervals
            if (contains(Stimulus_type, '_1000ms'))        
                grey_intervals = [-0.5 0; 1 2];
            elseif (contains(Stimulus_type, 'TRAIN250'))    
                grey_intervals = [-0.5 0; 0.25 0.5; 0.75 1.0; 1.25 1.5; 1.75 2.8];
            else
                grey_intervals = [-0.5 0; 0.15 0.3; 0.45 0.6; 0.75 0.9; 1.05 2];
            end    
                          
            % Get current y-axis limits for full vertical shading
            yl = ylim;
                
            % Shade each interval
            for i = 1:size(grey_intervals, 1)
                x = [grey_intervals(i,1), grey_intervals(i,2), grey_intervals(i,2), grey_intervals(i,1)];
                y = [yl(1), yl(1), yl(2), yl(2)];
                fill(x, y, [0.7 0.7 0.7], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % grey color with transparency
            end
            
            if (contains(Stimulus_type, 'E_1000ms'))
                xline(0, 'k', (sprintf('E %d%s onset', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
            elseif (contains(Stimulus_type, 'A_1000ms'))
                xline(0, 'k', (sprintf('A %d%s onset', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
            elseif (contains(Stimulus_type, 'TRAIN250'))
                xline(0, 'k', (sprintf('A %d%s onset', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(0.50, 'k', (sprintf('B %d%s onset', round(rad2deg(ordered_oris(2))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(1.0, 'k', (sprintf('C %d%s onset', round(rad2deg(ordered_oris(3))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(1.5, 'k', (sprintf('D %d%s onset', round(rad2deg(ordered_oris(4))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
            else    
                xline(0, 'k', (sprintf('A %d%s onset', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(0.30, 'k', (sprintf('B %d%s onset', round(rad2deg(ordered_oris(2))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(0.60, 'k', (sprintf('C %d%s onset', round(rad2deg(ordered_oris(3))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(0.90, 'k', (sprintf('D %d%s onset', round(rad2deg(ordered_oris(4))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
            end

            sgtitle(sprintf('%s - %s: Aggregate %s single unit activity across days', subject_number, Stimulus_type, depth_for_analysis), 'Interpreter', 'none');
            filename = sprintf('%s - %s - Multiplot %s FR traces across %d days.png', subject_number, Stimulus_type, depth_for_analysis, length(sessions_to_plot));
            save_path = fullfile(pwd, filename);
            exportgraphics(fig1, save_path); 

            % Also save as .fig
            fig_filename = sprintf('%s - %s - Multiplot %s FR traces across %d days.fig', subject_number, Stimulus_type, depth_for_analysis, length(sessions_to_plot));
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
            set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);
            ylabel(sprintf('Cumulative Z-scored FR (z-scored over %s)', z_score_period), 'Interpreter', 'none', 'FontSize', 24);
            sgtitle(sprintf('%s - %s: %s Cumulative Peak and Mean Firing Rates Across Days', subject_number, Stimulus_type, depth_for_analysis), 'Interpreter', 'none');
            legend({'Stim Peak', 'Stim Mean', 'Grey Peak', 'Grey Mean'}, 'Location', 'northeast', 'FontSize', 24);
            grid on;

            % Save bar chart figure
            saveas(fig2, fullfile(pwd, sprintf('%s %s - Multiplot %s cumulative peak and mean FRs across %d days.fig', subject_number, Stimulus_type, depth_for_analysis, length(sessions_to_plot))));
            exportgraphics(fig2, fullfile(pwd, sprintf('%s %s - Multiplot %s cumulative peak and mean FRs across %d days.png', subject_number, Stimulus_type, depth_for_analysis, length(sessions_to_plot))));


        end

       
        


        if contains(Stimulus_type, 'GAVNIK_ABCD') || contains(Stimulus_type, 'GAVNIK250_ABCD') || contains(Stimulus_type, 'GAVNIK_A___')... 
                && contains(plot_choice, 'aggregate') && contains(plot_type, 'FR')
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
                
                if contains(Stimulus_type, 'GAVNIK250_ABCD')
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.3 1.75], psthBinSize);  
                else    
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.3 1.35], psthBinSize);  
                end
                
                mean_trace = mean(binnedArray, 1); % average over trials
                
                figure(fig1);
                hold on;

                if contains(z_score_period, 'none')
                    ylim ([0 16]);
                     % Define grey intervals
                    if contains(Stimulus_type, 'GAVNIK250_ABCD')
                        grey_intervals = [-0.5 0; 1.0 1.9];
                    else    
                        grey_intervals = [-0.5 0; 0.6 1.5];
                    end

                    % Get current y-axis limits for full vertical shading
                    yl = ylim;
                        
                    % Shade each interval
                    for i = 1:size(grey_intervals, 1)
                        x = [grey_intervals(i,1), grey_intervals(i,2), grey_intervals(i,2), grey_intervals(i,1)];
                        y = [yl(1), yl(1), yl(2), yl(2)];
                        fill(x, y, [0.7 0.7 0.7], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % grey color with transparency
                    end
                    
                    % Plot raw mean firing rate trace
                    plot(bins, mean_trace, 'Color', colors(find(sessions_to_plot == nsession),:), 'LineWidth', 1.5, 'DisplayName', sprintf('Day %d, %s', experiment_info(nsession).date, depth_for_analysis));
                    ylabel('Mean firing rate (Hz)', 'FontSize', 24);
                    
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
                    
                    if contains(z_score_period, 'entire_session')
                        ylabel('Z-scored FR (z-scored over entire session)', 'FontSize', 24);
                        ylim([-1 5]);
                    elseif contains(z_score_period, 'first30secs')
                        ylabel('Z-scored FR (z-scored over first 30s baseline)', 'FontSize', 24);
                        ylim([-1 8]);
                    end

                    % Define grey intervals
                    if (contains(Stimulus_type, 'GAVNIK_A___'))
                        grey_intervals = [-0.5 0; 0.15 1.5];
                    elseif contains(Stimulus_type, 'GAVNIK250_ABCD')  
                        grey_intervals = [-0.5 0; 1.0 1.9];
                    else
                        grey_intervals = [-0.5 0; 0.6 1.5];
                    end    
                                   
                    % Get current y-axis limits for full vertical shading
                    yl = ylim;
                        
                    % Shade each interval
                    for i = 1:size(grey_intervals, 1)
                        x = [grey_intervals(i,1), grey_intervals(i,2), grey_intervals(i,2), grey_intervals(i,1)];
                        y = [yl(1), yl(1), yl(2), yl(2)];
                        fill(x, y, [0.7 0.7 0.7], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % grey color with transparency
                    end

                    plot(bins, z_trace, 'Color', colors(find(sessions_to_plot == nsession),:), 'LineWidth', 1.5, 'DisplayName', sprintf('Day %d, %s', experiment_info(nsession).date, depth_for_analysis));
                                        
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
                
                if contains(Stimulus_type, 'GAVNIK250_ABCD')
                    xlim([-0.3 1.75]);
                    xticks(-0.2:0.2:1.6);
                else    
                    xlim([-0.3 1.35]);
                    xticks(-0.4:0.2:1.4);
                end
                set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);
                xlabel('Time (s) since onset of A', 'FontSize', 24)
                
                legend(flipud(findobj(gca,'-property','DisplayName')), 'Location', 'northeast', 'FontSize', 24);
                hold on;                
            end
            
            if (contains(Stimulus_type, 'GAVNIK_A___'))
                xline(0, 'k', (sprintf('A %d%s onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
            elseif (contains(Stimulus_type, 'GAVNIK250_ABCD'))
                xline(0, 'k', (sprintf('A %d%s onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(0.25, 'k', (sprintf('B %d%s onset', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(0.50, 'k', (sprintf('C %d%s onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(0.75, 'k', (sprintf('D %d%s onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
            else
                xline(0, 'k', (sprintf('A %d%s onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(0.15, 'k', (sprintf('B %d%s onset', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(0.30, 'k', (sprintf('C %d%s onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(0.45, 'k', (sprintf('D %d%s onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
            end

            sgtitle(sprintf('%s - %s: Aggregate %s single unit activity across days', subject_number, Stimulus_type, depth_for_analysis), 'Interpreter', 'none');
            filename = sprintf('%s - %s - Multiplot %s FR traces across %d days.png', subject_number, Stimulus_type, depth_for_analysis, length(sessions_to_plot));
            save_path = fullfile(pwd, filename);
            exportgraphics(fig1, save_path); 

            % Also save as .fig
            fig_filename = sprintf('%s - %s - Multiplot %s FR traces across %d days.fig', subject_number, Stimulus_type, depth_for_analysis, length(sessions_to_plot));
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
            ylabel(sprintf('Cumulative Z-scored FR (z-scored over %s)', z_score_period), 'Interpreter', 'none', 'FontSize', 24);
            xlabel('Training Day', 'FontSize', 24)
            set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);
            sgtitle(sprintf('%s - %s: %s Cumulative Peak and Mean Firing Rates Across Days', subject_number, Stimulus_type, depth_for_analysis), 'Interpreter', 'none');
            legend({'Stim Peak', 'Stim Mean', 'Grey Peak', 'Grey Mean'}, 'Location', 'northeast', 'FontSize', 24);
            grid on;

            % Save bar chart figure
            saveas(fig2, fullfile(pwd, sprintf('%s %s - Multiplot %s cumulative peak and mean FRs across %d days.fig', subject_number, Stimulus_type, depth_for_analysis, length(sessions_to_plot))));
            exportgraphics(fig2, fullfile(pwd, sprintf('%s %s - Multiplot %s cumulative peak and mean FRs across %d days.png', subject_number, Stimulus_type, depth_for_analysis, length(sessions_to_plot))));


        end
                    
               
    end
end
hold off;

