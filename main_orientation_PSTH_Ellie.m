%%% For analysis of unit spiking in response to visual stimuli. ERB 2025
% Depths of V1 L5 and CA1 are both based on best regional power per PSD analysis (for which, use PSD_analysis_UnProcessedLFP_ellie).
% V1 L4 depth is identified via CSD analysis (for which, run CSD_Gratings_afterFILT_ellie or CSD_GAVNIKstims_afterFILT_ellie

addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

%% setting metrics to screen good clusters
clear all
% 1/5 Choose your probe depth of interest 
depth_for_analysis = 'V1'; % choose 'L4' or 'V1' or 'CA1' or 'Sub_CA1' or 
SUBJECTS = {'M00087'};
params = create_cluster_selection_params('sorting_option','ellie');
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);

%%% 2/5
Stimulus_type = 'F_1000ms'; % OMIT 'GAVNIK_ABCD_1 DCBA ADCD E_CD
% plot_choice 'struct' gives output of tuning metrics.
plot_choice = 'aggregate'; % curated 'single_units' or in 'aggregate' or uncurated 'MUA'; MUA includes all clusters from kilosort, unfiltered
plot_type = 'FR'; % 'FR' firing rate or 'raster' or 'struct' (for no plotting but output of metrics).
sliced_plot_option = 'no'; % 'yes' if you want to plot traces by groups of 40 trials to look for changes during the session
z_method = 'per_neuron'; % 'per_neuron' or 'in_aggregate'
z_score_period = 'stim_session'; % z score either over 'stim_session' (excludes variable greyscreen periods before and after stim paradigm), 'entire_session' or 'first30secs' or 'none' (for every stimulus recording
% session from 20250205 onward, I presented grey screen to the mouse for at least 30s before starting the stimulus.
%nprobe = 1;
%base_folder='V:\Ellie\DATA\SUBJECTS';

% 3/5 files will be saved here in the cd
cd('V:\Ellie\DATA\SUBJECTS\M00087\analysis\20260227\F_1000ms') 

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

% 5/5 
for nsession = 20 % row number of recording date in "experiment_info" 
    session_info = experiment_info(nsession).session(strcmp(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(strcmp(experiment_info(nsession).StimulusName,Stimulus_type));
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
        load(fullfile(options.ANALYSIS_DATAPATH, '..', 'earliest_V1sink_CSD.mat'))
        load(fullfile(options.ANALYSIS_DATAPATH, '..', 'depths_from_PSD.mat'))
        files = dir(fullfile(options.EPHYS_DATAPATH, '*ChanMap*.mat')); %the channel map has the y coordinate of each channel
        file_to_load = fullfile(options.EPHYS_DATAPATH, files(1).name); %load() does not accept wildcards like *
        load(file_to_load);
        
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
            
            this_shank_name = ['shank', num2str(this_shank)];  % creates e.g. 'shank2'
            Brain_surface_depth = depths_from_PSD.(this_shank_name).surface_depth_PSD;
            L4_channel_depth    = earliest_V1sink_CSD(find([earliest_V1sink_CSD.shank_id] == this_shank, 1)).best_depth_this_shank;
            L5_depth            = depths_from_PSD.(this_shank_name).L5_depth_PSD;
            CA1_depth           = depths_from_PSD.(this_shank_name).CA1_depth_PSD;
        
            L4_depth_range{iShank}   = [L4_channel_depth - 60, L4_channel_depth + 60]; % giving electrodes the full extent of 120um inclusive 
            % (hence measuring spiking over depth range greater than 120um. As 60 is divisible by 15 and 10, this will give the same effective range for NPX1.0 (which has staggered electrodes every 10um down the shank) and NPX2.0 (which has electrode rows every 15um down each shank)
            V1_depth_range{iShank}   = [L5_depth - 330, L5_depth + 700];
            CA1_depth_range{iShank}  = [CA1_depth - 150, CA1_depth + 150];
            Sub_CA1_depth_range{iShank} = [min(CA1_depth_range{iShank}) - 1000, min(CA1_depth_range{iShank})];  
        end
             
        all_orientations = unique(Task_info.stim_orientation); % uniqe values sorted in ascending order
        %Task_info.stim_onset
               
        params = create_cluster_selection_params('sorting_option','ellie');
        %params.orientation_tuned = ...
        psthBinSize = 0.01; % but use 1ms for raster plots
        
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

        if contains(z_score_period, 'stim_session') % excludes variable grey screen period before and after the stimulus paradigm ran
            time_edges = (min(Task_info.stim_onset) - 2):psthBinSize:(max(Task_info.stim_onset) + 2);
        end

        if contains(Stimulus_type, 'OP_Tuning') && contains(plot_choice, 'single_units') 
            % For M00013, the OP_Tuning stimulus was 150ms, with 100-200ms grey screen between stimuli, and took 24 different directions, moving at 2Hz
            % For M00014 (and beyond), the OP_Tuning stimulus was 67ms (but Bonsai typically keeps it on 84ms, 1 frame longer), with 17ms grey screen 
            % between stimuli (but Bonsai typically keeps it grey 34ms, 1 frame longer), and took 12 different static orientations
            
            for nprobe = 1:length(clusters)
                selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters, which pass the set parameters

                depth_selected_clusters = selected_clusters; % initialize with same structure
                np = nprobe;
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
                
                cluster_id = depth_selected_clusters(nprobe).cluster_id; % cluster_ids of units which pass the set parameters [NB these are one count higher than per zero-based pythonic SI output cluster IDs...]
                
                baseline_window = [0 30]; % in seconds - first 30s of recording is grey screen - can use for z-scoring

                % Define time_edges depending on z_score_period
                if contains(z_score_period, 'entire_session')
                    time_edges = 0:psthBinSize:max(depth_selected_clusters(nprobe).spike_times); %time bin edges from 0 to the max time in steps of ptshBinSize
                elseif contains(z_score_period, 'first30secs')
                    time_edges = baseline_window(1):psthBinSize:baseline_window(2);
                end
                
                ori_response = [];
                z_binnedArray = [];

                for nCluster = 1:length(cluster_id) % loop through each good cluster
                
                    spike_times_this_cluster = depth_selected_clusters(nprobe).spike_times(depth_selected_clusters(nprobe).spike_id == cluster_id(nCluster)); % extract the spike times for this cluster
                    %[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset, [-0.150 0.150], psthBinSize/10); % gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
                    
                    % Compute appropriate histogram counts for z-scoring
                    if contains(z_score_period, 'entire_session') || contains(z_score_period, 'stim_session')
                        zscore_counts = histcounts(spike_times_this_cluster, time_edges);
                    elseif contains(z_score_period, 'first30secs')
                        baseline_spikes = spike_times_this_cluster(spike_times_this_cluster >= baseline_window(1) & spike_times_this_cluster <= baseline_window(2));
                        zscore_counts = histcounts(baseline_spikes, time_edges);
                    end

                    fig(nCluster)=figure; % open a figure window
                    fig(nCluster).Name=sprintf('Orientation response %s Cluster %i', depth_for_analysis, cluster_id(nCluster)); %overall figure title includes the cluster_id (one count higher than zero-based SI output)
                    fig(nCluster).Position = [114 90 770 650]; % sets the size of the figure window
                    

                    for ori = 1:length(all_orientations) % loop through each unique stimulus orientation shown to the mouse during this task
                        if isequal(session_info(n).probe(1).SUBJECT, 'M00013')
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==all_orientations(ori)), [-0.07 0.150], psthBinSize/10); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window around stimulus onset
                        else
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==all_orientations(ori)), [-0.05 0.120], psthBinSize/10); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window around stimulus onset
                        end
                        nexttile % opens a subplot tile
                        plot(rasterX,rasterY,'k','LineWidth',1) % plot the spiking raster for the current orientation
                        xline(0,'r',LineWidth=1) % put a vertical red line at time zero (stimulus onset)
                        if isequal(session_info(n).probe(1).SUBJECT, 'M00013')
                            xlim([-0.15 0.150])
                        else
                            xlim([-0.05 0.120])
                        end

                        ylim([0 sum(Task_info.stim_orientation==all_orientations(ori))]) % set the max y coord to be the total # trials with this orientation
                        title(sprintf('Orientation %d%s',round(rad2deg(all_orientations(ori))), char(176)))
                        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',24)
                        
                        onset = Task_info.stim_onset(Task_info.stim_orientation == all_orientations(ori));
                        if isequal(session_info(n).probe(1).SUBJECT, 'M00013')
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  onset, [0.02 0.150], psthBinSize);
                        else
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  onset, [0.02 0.10], psthBinSize);
                        end

                        z_binnedArray{nCluster}{ori} = (binnedArray - mean(zscore_counts))./std(zscore_counts);
                        
                                            
                    end
                    nexttile % add a new subplot tile for a summary plot showing how the z-scored firing for the cluster varies by orientation
                    % Each orientation is a different colour
                    clear PLOT
                    for ori = 1:length(all_orientations)
                        PLOT(ori) = plot(bins, mean(z_binnedArray{nCluster}{ori})); %plot the z-scored PSTH (mean across trials) for each orientation
                        hold on
                    end
                    legend(PLOT(1:end),{num2str(round(rad2deg(all_orientations)))},'box','off')
                    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',24)


                    nexttile % generate a summary tuning curve for the cluster, showing z-scored response to each orientation, with error bars 
                    mean_response=[]; %initialise array to hold the mean response across trials for each orientation
                    std_response = []; %initialise array to hold the sd of the cluster's response for each orienation
                    for ori = 1:length(all_orientations)
                        ori_response{nCluster}{ori} = mean(z_binnedArray{nCluster}{ori}(:,bins>0.02&bins<0.11), 2); % extract the mean z-scored firing 20-110ms after stim onset, for each trial
                        mean_response(ori) = mean(ori_response{nCluster}{ori}); % calc the mean of mean z-scored firing after stim-onset, across trials
                        se_response(ori) = std(ori_response{nCluster}{ori})./sqrt(length(ori_response{nCluster}{ori})); %calc the SE of the mean z-scored firing after stim-onset
                    end


                    clear PLOT
                    % plot(round(rad2deg(all_orientations)),mean_response)
                    % hold on;
                    errorbar(round(rad2deg(all_orientations)),mean_response,se_response,se_response) %plot the mean of z-scored firing (y-axis) in response to each orientation (x-axis)
                    % with error bars showing SE above and below mean z-scored response for each orientation
                    % Note that z-scoring against entire session leads to quite low z-scores as a stimulus is on the screen half the time
                    if isequal(session_info(n).probe(1).SUBJECT, 'M00013')
                        xlim([-20 370]) % x-axis runs from -20 across 360 degrees
                    else
                        xlim([-20 190]) % x-axis runs from -20 across 180 degrees
                    end
                        % legend(PLOT(1:end),{num2str(round(rad2deg(all_orientations)))},'box','off')
                    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',24)
                    cluster_depth = clusters(nprobe).peak_depth(cluster_id(nCluster));
                    sgtitle(sprintf('%s - %s - Response of %s Cluster %i (%.0f µm)', subject_number, Stimulus_type, depth_for_analysis, cluster_id(nCluster), cluster_depth), 'Interpreter', 'none');
                    fig_filename1 = sprintf('%s - %s - Response of %s Cluster %i.fig', subject_number, Stimulus_type, depth_for_analysis, cluster_id(nCluster));
                    savefig(fig(nCluster), fullfile(pwd, fig_filename1));

                end
            end
        end


            

        if contains(Stimulus_type, 'OP_Tuning') && contains(plot_choice, 'struct') 
            % For M00013, the OP_Tuning stimulus was 150ms, with 100-200ms grey screen between stimuli, and took 24 different directions, moving at 2Hz
            % For M00014 (and beyond), the OP_Tuning stimulus was 67ms (but Bonsai typically keeps it on 84ms, 1 frame longer), with 17ms grey screen 
            % between stimuli (but Bonsai typically keeps it grey 34ms, 1 frame longer), and took 12 different static orientations
            
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
                
                % Initialize matrix to hold summed spike counts: [nOrientations x nClusters]
                summedNetSpikeCounts = nan(length(all_orientations), length(cluster_id)); 
                summedSpontActivity = nan(length(all_orientations), length(cluster_id));
                
                OP_OSI = nan(1, length(cluster_id));
                preferred_oris = nan(1, length(cluster_id));
                OP_DSI = nan(1, length(cluster_id));
                preferred_dirs = nan(1, length(cluster_id));

                OP_tuning = struct( ...
                    'cluster_id', []);   

                % For (1 - CV) orientation and direction preference function:    
                th = all_orientations'; % transpose to a row vector to suit the function format
                % Preallocate output vectors
                gDSI     = zeros(1, length(cluster_id));
                DirAngle  = zeros(1, length(cluster_id));
                gOSI     = zeros(1, length(cluster_id));
                OriAngle  = zeros(1, length(cluster_id));

                for nCluster = 1:length(cluster_id) % loop through each good cluster
                
                    spike_times_this_cluster = depth_selected_clusters(nprobe).spike_times(depth_selected_clusters(nprobe).spike_id == cluster_id(nCluster)); % extract the spike times for this cluster
                                        
                    % Compute appropriate histogram counts for z-scoring
                    if contains(z_score_period, 'entire_session') || contains(z_score_period, 'stim_session')
                        zscore_counts = histcounts(spike_times_this_cluster, time_edges);
                    elseif contains(z_score_period, 'first30secs')
                        baseline_spikes = spike_times_this_cluster(spike_times_this_cluster >= baseline_window(1) & spike_times_this_cluster <= baseline_window(2));
                        zscore_counts = histcounts(baseline_spikes, time_edges);
                    end

                    
                    for ori = 1:length(all_orientations) % loop through each unique stimulus orientation shown to the mouse during this task
                        onset = Task_info.stim_onset(Task_info.stim_orientation == all_orientations(ori));
                        if isequal(session_info(n).probe(1).SUBJECT, 'M00013')
                            [~, ~, ~, ~, spikeCount, ~] = psthAndBA(spike_times_this_cluster,  onset, [0.02 0.11], psthBinSize); % for this orientation, get the spike count during 20-110ms following stimulus onset
                            [~, ~, ~, ~, spontActivity, ~] = psthAndBA(spike_times_this_cluster,  onset, [-0.07 0.02], psthBinSize); % for this orientation, get the spike count during the 90ms prior to stimulus response
                        else
                            [~, ~, ~, ~, spikeCount, ~] = psthAndBA(spike_times_this_cluster,  onset, [0.02 0.09], psthBinSize); % for this orientation, get the spike count during stimulus presentation
                        end   
                        % Sum spike counts across all trials of this orientation
                        summedSpontActivity(ori,nCluster) = sum(spontActivity);
                        summedNetSpikeCounts(ori, nCluster) = sum(spikeCount) - sum(spontActivity); % each row of summedSpikeCounts is an orientation; each column is a cluster, each value is the spikecount during all trials of that orientation            
                        
                    end

                    
                    % ===== Manual Orientation & Direction Preference Calculation =====
                    % Extract this cluster's summed spike counts across all orientations
                    thisClusterCounts = summedNetSpikeCounts(:, nCluster);  % length = number of orientations
                    dirs_rad = all_orientations(:); % in radians across the full 2 pi of a circle
                    nDirs = length(dirs_rad);
                    nOris = nDirs / 2;  % Assumes directions evenly spaced

                    % --- Direction Preference (0–2π)
                    [~, max_dir_idx] = max(thisClusterCounts);
                    preferred_dir_rad = dirs_rad(max_dir_idx);
                    preferred_dirs(nCluster) = rad2deg(mod(preferred_dir_rad, 2*pi));

                    % --- Orientation Preference (0–π): collapse across direction pairs
                    % Collapse responses from direction pairs (e.g., 0 & π, π/4 & 5π/4)
                    ori_angles = mod(dirs_rad, pi);  % 0 to π
                    bin_width = pi / nOris;
                    % Round to nearest bin center to reduce floating-point error
                    rounded_ori_angles = round(ori_angles / bin_width) * bin_width;

                    [unique_oris, ~, ori_group_idx] = unique(rounded_ori_angles);  % group indices
                    
                    ori_responses = accumarray(ori_group_idx, thisClusterCounts);  % sum both directions
                    [~, max_ori_idx] = max(ori_responses);
                    preferred_ori_rad = unique_oris(max_ori_idx);
                    preferred_oris(nCluster) = rad2deg(mod(preferred_ori_rad, pi));  % 0–180°

                    % --- Manual OSI
                    orth_ori = mod(preferred_ori_rad + pi/2, pi); % calc the orghogonal orientation (90 degrees different) from the preferred ori, in radians, within the 0 to pi range.
                    % Find closest matching index in `unique_oris`
                    [~, orth_idx] = min(abs(unique_oris - orth_ori));
                    osi_num = ori_responses(max_ori_idx) - ori_responses(orth_idx); % numerator for OSI calc.
                    osi_denom = ori_responses(max_ori_idx) + ori_responses(orth_idx); % denominator for OSI calc.
                    OP_OSI(nCluster) = osi_num / (osi_denom + eps);  % add eps to avoid division by zero - eps is the smallest positive number that can be added to 1.0 and result in a value different from 1.0 in double precision

                    % --- Manual DSI
                    opposite_dir = mod(preferred_dir_rad + pi, 2*pi);
                    [~, opp_idx] = min(abs(dirs_rad - opposite_dir));
                    dsi_num = thisClusterCounts(max_dir_idx) - thisClusterCounts(opp_idx);
                    dsi_denom = thisClusterCounts(max_dir_idx) + thisClusterCounts(opp_idx);
                    OP_DSI(nCluster) = dsi_num / (dsi_denom + eps);  % add eps to avoid division by zero

                    %%%% Using complex VECTORs to calculate orientation and direction preferences:                
                    r = summedNetSpikeCounts(:, nCluster)'; % transpose to a row vector
                    [gDSI(nCluster),DirAngle(nCluster),gOSI(nCluster),OriAngle(nCluster)] = calcCV(th,r);
                    DirAngle = mod(DirAngle, 360);  % Wraps negative angles and angles ≥360 back into [0, 360)
                    OriAngle = mod(OriAngle, 180);  % Orientation is periodic over 180°, so wrap accordingly
                end

                % Package VECTOR results into a struct array, "OP_tuning" for reuse
                OP_tuning(nprobe).cluster_id = cluster_id(:)'; % include cluster ids for reference
                OP_tuning(nprobe).DirAngle = DirAngle(:)'; % preferred direction
                OP_tuning(nprobe).OriAngle = OriAngle(:)'; % preferred orientation
                OP_tuning(nprobe).gOSI = gOSI(:)'; % (Circular Variance; 0 is untuned, 1 highly tuned)
                OP_tuning(nprobe).gDSI = gDSI(:)';
                OP_tuning(nprobe).summedNetSpikeCounts_per_ori = summedNetSpikeCounts;
                OP_tuning(nprobe).all_orientations = all_orientations(:);
                save('OP_tuning.mat', 'OP_tuning'); % save to current directory

                % Plot gOSI vs gDSI to check the data
                plot_gOSI = OP_tuning(nprobe).gOSI;
                plot_gDSI = OP_tuning(nprobe).gDSI;
                cluster_ids = OP_tuning(nprobe).cluster_id;
                fig=figure;
                scatter(plot_gOSI, plot_gDSI, 40, 'filled');
                xlabel('Orientation Selectivity (gOSI)');
                ylabel('Direction Selectivity (gDSI)');
                grid on;
                xlim([0 1]);
                ylim([0 1]);
                axis square;
                % Annotate each point with cluster ID
                hold on;
                for i = 1:length(cluster_ids)
                    text(plot_gOSI(i) + 0.01, plot_gDSI(i), num2str(cluster_ids(i)), ...
                        'FontSize', 8, 'Color', [0.3 0.3 0.3]); % small offset and light gray
                end
                sgtitle(sprintf('%s_%s: gOSI vs gDSI Day %d', subject_number, Stimulus_type, experiment_info(nsession).date), 'Interpreter', 'None');
                fig_filename = sprintf('%s_%s_%s_gOSI vs gDSI Day %d.fig', ...
                    subject_number, Stimulus_type, depth_for_analysis, experiment_info(nsession).date);
                savefig(fig, fullfile(pwd, fig_filename));
                hold off;

                % Plot preferred ori vs direction to check the data
                plot_DirAngle = OP_tuning(nprobe).DirAngle;
                plot_OriAngle = OP_tuning(nprobe).OriAngle;
                cluster_ids = OP_tuning(nprobe).cluster_id;
                fig=figure;
                scatter(plot_OriAngle, plot_DirAngle, 40, 'filled');
                xlabel('Preferred Orientation');
                ylabel('Preferred Direction');
                grid on;
                xlim([0 180]);
                ylim([0 360]);
                axis square;
                % Annotate each point with cluster ID
                hold on;
                for i = 1:length(cluster_ids)
                    text(plot_OriAngle(i) + 1, plot_DirAngle(i), num2str(cluster_ids(i)), ...
                        'FontSize', 8, 'Color', [0.3 0.3 0.3]); % small offset and light gray
                end
                sgtitle(sprintf('%s_%s: pref Ori vs pref Dir, Day %d', subject_number, Stimulus_type, experiment_info(nsession).date), 'Interpreter', 'None');
                fig_filename = sprintf('%s_%s_%s_pref Ori vs pref Dir, Day %d.fig', ...
                subject_number, Stimulus_type, depth_for_analysis, experiment_info(nsession).date);
                savefig(fig, fullfile(pwd, fig_filename));
                hold off;
              
            end
        end
        



        if contains(Stimulus_type, 'TRAIN') && contains(plot_choice, 'single_units') 
            for nprobe = 1:length(clusters)
                selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters, which pass the set parameters
                
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

                    % Keep only spike times and IDs for clusters at selected depths
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
                
                ori_response = [];
                z_binnedArray = [];

                for nCluster = 1:length(cluster_id) % loop through each good cluster in selected depth
                
                    spike_times_this_cluster = depth_selected_clusters(nprobe).spike_times(depth_selected_clusters(nprobe).spike_id == cluster_id(nCluster)); % extract the spike times for this cluster
                    %[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset, [-0.150 0.30], 0.001); % gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
                    
                    if contains(z_score_period, 'entire_session') || contains(z_score_period, 'stim_session')
                        zscore_counts = histcounts(spike_times_this_cluster, time_edges); % use all spike times
                    elseif contains(z_score_period, 'first30secs')
                        baseline_spikes = spike_times_this_cluster(spike_times_this_cluster >= baseline_window(1) & spike_times_this_cluster <= baseline_window(2));
                        zscore_counts = histcounts(baseline_spikes, time_edges); % use just baseline spikes
                    end

                    fig(nCluster)=figure; % open a figure window
                    fig(nCluster).Name=sprintf('%s Grating responses Cluster %i', Stimulus_type, cluster_id(nCluster)); %overall figure title includes the cluster_id (one count higher than zero-based SI output)
                    fig(nCluster).Position = [114 90 770 650]; % sets the size of the figure window
                
                    % Extract ordered orientations based on presentation sequence
                    ordered_oris = unique(Task_info.stim_orientation, 'stable'); % radians for TRAIN but degrees for GAVNIK
                    
                    ori = 1;  

                    if contains(plot_type, 'raster')
                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.3 1.8], psthBinSize/10); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
                        plot(rasterX,rasterY,'k','LineWidth',1) % plot the spiking raster for the current orientation
                        xline(0,'r',LineWidth=1) % put a vertical red line at time zero (stimulus onset)
                        xlim([-0.5 2])
                        ylim([0 sum(Task_info.stim_orientation==ordered_oris(ori))]) % set the max y coord to be the total # trials with this orientation
                        ylabel('Trial');
                        % imagesc(binnedArray)
                        % xticks([1.5 10.5 15.5 20.5 30.5])
                        % xticklabels([-0.150 -0.050 0 0.050 0.150])
                        % xline(15.5,'r',LineWidth=1)
                        % colorbar
                        % colormap(flipud(gray))
                        sgtitle(sprintf('%s %s cluster %d rasters %s', subject_number, depth_for_analysis, cluster_id, Stimulus_type, 'Interpreter', 'none', 'FontSize', 24));
                        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',24)
                        fig_filename1 = sprintf('%s - %s - %s Cluster_%i_raster.fig', subject_number, Stimulus_type, depth_for_analysis, cluster_id(nCluster));
                        savefig(fig(nCluster), fullfile(pwd, fig_filename1));

                    end
                    
                    
                    if contains(plot_type, 'FR')
                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.3 1.8], psthBinSize); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
                        if contains(z_score_period, 'none')
                            mean_trace = mean(binnedArray, 1) / psthBinSize; % average over trials
                            sem_trace = std(binnedArray, 0, 1) / sqrt(size(binnedArray, 1))/ psthBinSize;  % [1 x time]
                            % Plot SEM shading
                            fill([bins, fliplr(bins)], ...
                                 [mean_trace + sem_trace, fliplr(mean_trace - sem_trace)], ...
                                 'b', ...
                                 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                            hold on;
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
                            % Plot raw mean firing rate trace
                            plot(bins, mean_trace, 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', sprintf('%s (200x)', Stimulus_type));
                            ylabel('Mean firing rate across trials (Hz)', 'FontSize', 24);
                            %ylim ([0 16]);
                        else 
                            baseline_mean = mean(zscore_counts);
                            baseline_std = std(zscore_counts);
                            % Z-score each trial individually so that trial to variability remains visible and SEM can be calculated
                            zscored_trials = (binnedArray - baseline_mean) / baseline_std;  % [trials x time]
                            z_trace = mean(zscored_trials, 1);                             % mean z-scored trace
                            sem_trace = std(zscored_trials, 0, 1) / sqrt(size(zscored_trials, 1));  % SEM

                            if contains(z_score_period, 'entire_session')
                                ylabel('Z-scored FR (z-scored over entire session)', 'FontSize', 24);
                                ylim([-1 5]);
                            elseif contains(z_score_period, 'stim_session')    
                                ylabel('FR Z-scored over stim session', 'FontSize', 24);
                            elseif contains(z_score_period, 'first30secs')
                                ylabel('Z-scored FR (z-scored over first 30s baseline)', 'FontSize', 24);
                                ylim([-1 8]);
                            end
                            
                            % Define grey intervals
                            grey_intervals = [-0.5 0; 0.15 0.3; 0.45 0.6; 0.75 0.9; 1.05 1.5];
                                         
                            % Get current y-axis limits for full vertical shading
                            yl = ylim;
                                
                            % Shade each interval
                            for i = 1:size(grey_intervals, 1)
                                x = [grey_intervals(i,1), grey_intervals(i,2), grey_intervals(i,2), grey_intervals(i,1)];
                                y = [yl(1), yl(1), yl(2), yl(2)];
                                fill(x, y, [0.7 0.7 0.7], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % grey color with transparency
                            end

                            plot(bins, z_trace, 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', sprintf('%s (200x)', Stimulus_type));
                                hold on; 
                                % Plot SEM shading
                                fill([bins, fliplr(bins)], ...
                                    [z_trace + sem_trace, fliplr(z_trace - sem_trace)], ...
                                    'b', ...
                                    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');                                    
                        end   
        
                        xlim([-0.3 1.8]);
                        xticks(-0.2:0.2:1.8);
                        xlabel('Time (s)', 'FontSize', 24)
                        set(gca, 'FontSize', 24);  % Tick labels font size
                        legend(flipud(findobj(gca,'-property','DisplayName')), 'Location', 'northeast', 'Interpreter', 'none');
                        hold on;                
                        
                        xline(0, 'k', (sprintf('A %d%s onset', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                        xline(0.30, 'k', (sprintf('B %d%s onset', round(rad2deg(ordered_oris(2))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                        xline(0.60, 'k', (sprintf('C %d%s onset', round(rad2deg(ordered_oris(3))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                        xline(0.90, 'k', (sprintf('D %d%s onset', round(rad2deg(ordered_oris(4))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);
                        
                        sgtitle(sprintf('%s %s cluster %d FRs %s', subject_number, depth_for_analysis, cluster_id(nCluster), Stimulus_type), 'Interpreter', 'none', 'FontSize', 24);
                        fig_filename = sprintf('%s %s cluster %d FRs %s (zscore %s).fig', subject_number, depth_for_analysis, cluster_id(nCluster), Stimulus_type, z_score_period);
                        fig_save_path = fullfile(pwd, fig_filename);
                        savefig(fig(nCluster), fig_save_path);
                    end
                        
                end
            end
        end

        
               

        if contains(Stimulus_type, 'TRAIN')|| contains(Stimulus_type, '_1000ms') || contains(Stimulus_type, 'F_150ms')...
                || contains(Stimulus_type, 'A___') || contains(Stimulus_type, 'A_50ms') || contains(Stimulus_type, 'A_500ms')... 
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
                % all_spike_times = [];
                % all_spike_ids = [];
                % for np = 1:length(depth_selected_clusters)
                %     all_spike_times = [all_spike_times; depth_selected_clusters(np).spike_times];
                %     all_spike_ids = [all_spike_ids; depth_selected_clusters(np).spike_id];
                % end                  

                % Combine spike times separately for each shank
                for iShank = 1:numel(unique_shanks)
                    all_spike_times{iShank} = [];
                    all_spike_ids{iShank} = [];
                end
                
                % Storage for per-neuron spike times (needed if z_method is 'per_neuron')
                for iShank = 1:numel(unique_shanks)
                    spike_times_by_cluster{iShank} = {};
                end


                for np = 1:length(depth_selected_clusters)
                    sc = depth_selected_clusters(np);
                    
                    cluster_channels = clusters(np).peak_channel(sc.cluster_id);
                    cluster_shanks = kcoords(cluster_channels);
                
                    for iShank = 1:numel(unique_shanks)
                        this_shank = unique_shanks(iShank);
                
                        shank_mask = cluster_shanks == this_shank;
                        shank_cluster_ids = sc.cluster_id(shank_mask);

                        for cid = shank_cluster_ids'
                            neuron_mask = sc.spike_id == cid;
                            neuron_spikes = sc.spike_times(neuron_mask);                       
                            spike_times_by_cluster{iShank}{end+1} = neuron_spikes;                        
                        end
                
                        spike_mask = ismember(sc.spike_id, shank_cluster_ids);
                
                        all_spike_times{iShank} = [all_spike_times{iShank}; sc.spike_times(spike_mask)];
                        all_spike_ids{iShank} = [all_spike_ids{iShank}; sc.spike_id(spike_mask)];
                    end
                end


                ordered_oris = unique(Task_info.stim_orientation, 'stable');
                                              
                fig = figure;
                fig.Name = sprintf('Aggregate activity: %s', depth_for_analysis);
                fig.Position = [114 90 770 650];
                %tiledlayout(5,1);

                %for ori = 1:length(ordered_oris)
                ori = 1; % make plots from onset of first stim in sequence
                    
                stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));

                    %[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.150 0.30], psthBinSize);
                %[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.30 1.8], psthBinSize);                      
                %mean_trace = mean(binnedArray, 1); % average binned firing over trials

                for iShank = 1:numel(unique_shanks)
                    if strcmp(z_method,'in_aggregate')
                        [psth{iShank}, bins{iShank}, rasterX{iShank}, rasterY{iShank}, spikeCounts{iShank}, binnedArray{iShank}] = ...
                            psthAndBA(all_spike_times{iShank}, stim_onsets, [-0.30 1.8], psthBinSize);                    
                        mean_trace{iShank} = mean(binnedArray{iShank}, 1); % average binned firing over trials              
                    elseif strcmp(z_method,'per_neuron')
                        % ----- PSTH PER NEURON -----
                        neuron_traces{iShank} = [];
                
                        for num = 1:length(spike_times_by_cluster{iShank})
                
                            neuron_spikes = spike_times_by_cluster{iShank}{num};
                
                            [~, bins{iShank}, ~, ~, ~, binnedArray] = ...
                                psthAndBA(neuron_spikes, stim_onsets, [-0.30 1.8], psthBinSize);
                
                            neuron_trace = mean(binnedArray,1);                
                            neuron_traces{iShank}(num,:) = neuron_trace;                
                        end                
                        % population mean (used later if no z-scoring)
                        mean_trace{iShank} = mean(neuron_traces{iShank},1);
                    end                                        
                end


                for iShank = 1:numel(unique_shanks)
                    subplot(numel(unique_shanks),1,iShank)
                    title(sprintf('Shank %d', unique_shanks(iShank)))
                    hold on;

                    if contains(z_score_period, 'none')
                        % Plot raw mean firing rate trace
                        plot(bins{iShank}, mean_trace{iShank}, 'b', 'LineWidth', 1.5);
                        ylabel('Mean firing rate (Hz)');
                    %    ylim ([0 16]);
                   
                    else
                        if strcmp(z_method,'in_aggregate')
                        % Compute appropriate histogram counts for z-scoring
                            if contains(z_score_period, 'entire_session')
                                zscore_counts = histcounts(all_spike_times{iShank}, time_edges);
                                hold on;
                                ylim([-1 6]);
                                ylabel('FR Z-scored in aggregate over entire session)', 'FontSize', 24);
                            elseif contains(z_score_period, 'stim_session')  
                                zscore_counts = histcounts(all_spike_times{iShank}, time_edges);
                                hold on;
                                ylabel('FR Z-scored in aggregate over stim session)', 'FontSize', 24);
                            elseif contains(z_score_period, 'first30secs')
                                baseline_spikes = all_spike_times{iShank}(all_spike_times{iShank} >= baseline_window(1) & all_spike_times{iShank} <= baseline_window(2));
                                zscore_counts = histcounts(baseline_spikes, time_edges); %histcounts counts the number of datapoints in specified time bins
                                hold on;
                                ylim([-1 8]);
                                ylabel('FR Z-scored in aggregate over first 30s baseline)');
                            end
        
                            z_trace = (mean_trace{iShank} - mean(zscore_counts)) / std(zscore_counts); %mean(zscore_counts) gives the mean spikes per timebin in the reference period; this is then deducted from the spikecount of each trial-averaged timebin
                        
                        elseif strcmp(z_method,'per_neuron')

                            % ----- Z-SCORE EACH NEURON -----                  
                            nNeurons = size(neuron_traces{iShank},1);
                            z_neuron_traces = nan(size(neuron_traces{iShank})); % allow invalid neurons to remain NaN instead of corrupting the mean 
                            if contains(z_score_period,'entire_session')
                                ylabel('FR Z-scored per-neuron over entire session','FontSize',24);
                            elseif contains(z_score_period,'stim_session')
                                ylabel('FR Z-scored per-neuron over stim session','FontSize',24);     
                            elseif contains(z_score_period,'first30secs')
                                ylabel('FR Z-scored per-neuron over first 30s baseline','FontSize',24);
                            end
                            

                            for nneur = 1:nNeurons                    
                                neuron_spikes = spike_times_by_cluster{iShank}{nneur};                   
                                if contains(z_score_period,'entire_session') || contains(z_score_period,'stim_session')                    
                                    counts = histcounts(neuron_spikes,time_edges);
                                elseif contains(z_score_period,'first30secs')                   
                                    baseline_spikes = neuron_spikes( ...
                                        neuron_spikes >= baseline_window(1) & ...
                                        neuron_spikes <= baseline_window(2));                    
                                    counts = histcounts(baseline_spikes,time_edges); 
                                end
                    
                                mu = mean(counts);
                                sigma = std(counts);  
                                if sigma > 0 % Neurons with sigma = 0 will stay NaN and be ignored.
                                    z_neuron_traces(nneur,:) = (neuron_traces{iShank}(nneur,:) - mu) / sigma; 
                                end
                            end                  
                            % Average z-scored neurons
                            z_trace = mean(z_neuron_traces,1,'omitnan');                    
                            %ylim([-1 6]);                                             
                        end

                        plot(bins{iShank}, z_trace, 'b', 'LineWidth', 1.5);
                        
                    end
    
                    % Define windows after time zero and calc. the peak and mean FR [for exploration]
                    if contains(Stimulus_type, 'A_50ms')
                        stim_window_starts = 0;
                        stim_window_ends = 0.05;
                        grey_window_starts = 0.05; 
                        grey_window_ends = [1.8]; 
                        peak_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                        mean_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                        peak_FR_by_greywindow = zeros(1, length(stim_window_starts)); % preallocate
                        mean_FR_by_greywindow = zeros(1, length(stim_window_starts)); % preallocate
                    elseif contains(Stimulus_type, 'F_150ms')
                        stim_window_starts = 0;
                        stim_window_ends = 0.15;
                        grey_window_starts = 0.15; 
                        %window_ends = 0.15:0.15:1.5;
                        grey_window_ends = [1.8]; % window ending at 1.35s is where a stimulus would end if there was a fifth stimulus in the sequence)
                        peak_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                        mean_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                        peak_FR_by_greywindow = zeros(1, length(stim_window_starts)); % preallocate
                        mean_FR_by_greywindow = zeros(1, length(stim_window_starts)); % preallocate 
                    elseif contains(Stimulus_type, '_1000ms')
                        stim_window_starts = 0;
                        stim_window_ends = 1.0;
                        grey_window_starts = 1.0; 
                        %window_ends = 0.15:0.15:1.5;
                        grey_window_ends = [1.8]; % window ending at 1.35s is where a stimulus would end if there was a fifth stimulus in the sequence)
                        peak_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                        mean_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                        peak_FR_by_greywindow = zeros(1, length(stim_window_starts)); % preallocate
                        mean_FR_by_greywindow = zeros(1, length(stim_window_starts)); % preallocate      
                    elseif contains(Stimulus_type, 'A_500ms')
                        stim_window_starts = 0;
                        stim_window_ends = 0.5;
                        grey_window_starts = 0.5; 
                        %window_ends = 0.15:0.15:1.5;
                        grey_window_ends = [1.8]; % window ending at 1.35s is where a stimulus would end if there was a fifth stimulus in the sequence)
                        peak_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                        mean_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                        peak_FR_by_greywindow = zeros(1, length(stim_window_starts)); % preallocate
                        mean_FR_by_greywindow = zeros(1, length(stim_window_starts)); % preallocate
                    else
                        stim_window_starts = 0:0.3:0.9;
                        stim_window_ends = 0.15:0.3:1.05;
                        grey_window_starts = [0.15 0.45 0.75 1.05 1.2 1.35]; % look from 0ms after stim offset 
                        grey_window_ends = [0.3 0.6 0.9 1.2 1.35 1.5]; % window ending at 1.35s is where a stimulus would end if there was a fifth stimulus in the sequence)
                        peak_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                        mean_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                        peak_FR_by_greywindow = zeros(1, length(stim_window_starts)); % preallocate
                        mean_FR_by_greywindow = zeros(1, length(stim_window_starts)); % preallocate
                    end    
                                    
                    yl = ylim; % Get y-axis limits for text placement
                        
                    for i = 1:length(stim_window_starts)
                        if contains(z_score_period, 'none')
                            % Find indices of bins within current window
                            idx_in_stimwindow = bins{iShank} >= stim_window_starts(i) & bins{iShank} < stim_window_ends(i);                    
                            peak_FR_by_stimwindow(i) = max(mean_trace{iShank}(idx_in_stimwindow));
                            mean_FR_by_stimwindow(i) = mean(mean_trace{iShank}(idx_in_stimwindow));
                            
                        else
                            % Find indices of bins within current window
                            idx_in_stimwindow = bins{iShank} >= stim_window_starts(i) & bins{iShank} < stim_window_ends(i);                    
                            peak_FR_by_stimwindow(i) = max(z_trace(idx_in_stimwindow));
                            mean_FR_by_stimwindow(i) = mean(z_trace(idx_in_stimwindow));
                            
                        end
                    end
                    
                    for i = 1:length(grey_window_starts)
                        if contains(z_score_period, 'none')
                            % Find indices of bins within current window
                            idx_in_greywindow = bins{iShank} >= grey_window_starts(i) & bins{iShank} < grey_window_ends(i);                    
                            peak_FR_by_greywindow(i) = max(mean_trace{iShank}(idx_in_greywindow));
                            mean_FR_by_greywindow(i) = mean(mean_trace{iShank}(idx_in_greywindow));
                            
                        else
                            % Find indices of bins within current window
                            idx_in_greywindow = bins{iShank} >= grey_window_starts(i) & bins{iShank} < grey_window_ends(i);                    
                            peak_FR_by_greywindow(i) = max(z_trace(idx_in_greywindow));
                            mean_FR_by_greywindow(i) = mean(z_trace(idx_in_greywindow));
                            
                        end
                    end
                    
                    if contains(Stimulus_type, 'A_50ms')
                        % Define grey intervals
                        grey_intervals = [-0.5 0; 0.05 1.8];
                    elseif contains(Stimulus_type, 'A_500ms') 
                        grey_intervals = [-0.5 0; 0.5 1.8];
                    elseif contains(Stimulus_type, '_150ms') 
                        grey_intervals = [-0.5 0; 0.15 1.8];
                    elseif contains(Stimulus_type, '_1000ms') 
                        grey_intervals = [-0.5 0; 1.0 1.8];    
                    else
                        grey_intervals = [-0.5 0; 0.15 0.3; 0.45 0.6; 0.75 0.9; 1.05 2];
                    end    
                    
                    hold on; % Make sure current plot stays visible
                    
                    % Get current y-axis limits for full vertical shading
                    yl = ylim;
                    
                    % Shade each interval
                    for i = 1:size(grey_intervals, 1)
                        x = [grey_intervals(i,1), grey_intervals(i,2), grey_intervals(i,2), grey_intervals(i,1)];
                        y = [yl(1), yl(1), yl(2), yl(2)];
                        fill(x, y, [0.7 0.7 0.7], 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % grey color with transparency
                    end
                    
                    if contains(Stimulus_type, 'A_50ms') || contains(Stimulus_type, 'A_500ms')
                        xline(0, 'k', (sprintf('A onset; %d%s', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                    elseif contains(Stimulus_type, 'F_150ms') || contains(Stimulus_type, 'F_1000ms')
                        xline(0, 'k', (sprintf('F onset; %d%s', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                    elseif contains(Stimulus_type, 'E_1000ms')
                        xline(0, 'k', (sprintf('E onset; %d%s', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);    
                    else
                        xline(0, 'k', (sprintf('A onset; %d%s', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                        xline(0.30, 'k', (sprintf('B onset; %d%s', round(rad2deg(ordered_oris(2))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                        xline(0.60, 'k', (sprintf('C onset; %d%s', round(rad2deg(ordered_oris(3))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                        xline(0.90, 'k', (sprintf('D onset; %d%s', round(rad2deg(ordered_oris(4))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                    end
    
                    xlim([-0.3 1.8]);
                    xticks(-0.2:0.2:1.8);
                    if contains(Stimulus_type, 'F_150ms') || contains(Stimulus_type, 'F_1000ms')
                        xlabel('Time (s) since onset of F')
                    elseif contains(Stimulus_type, 'E_1000ms')
                        xlabel('Time (s) since onset of E')    
                    else
                        xlabel('Time (s) since onset of A')
                    end    
                    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);
    
                end                    
                sgtitle(sprintf('%s day %d - %s - Aggregate single unit activity: %s', subject_number, experiment_info(nsession).date, Stimulus_type, depth_for_analysis), 'Interpreter', 'none', 'FontSize', 24);
                %filename = sprintf('%s - Aggregate single unit %s FRs.png', Stimulus_type, depth_for_analysis);
                %save_path = fullfile(pwd, filename);
                %exportgraphics(fig, save_path);  % Add this line

                % Also save as .fig
                fig_filename = sprintf('%s - Aggregate single unit %s FRs.fig', Stimulus_type, depth_for_analysis);
                fig_save_path = fullfile(pwd, fig_filename);
                savefig(fig, fig_save_path);
                    
            end
        end

        


        if (contains(Stimulus_type, 'OMIT') || strcmp(Stimulus_type, 'E_CD')... 
                && contains(plot_choice, 'aggregate') && contains(plot_type, 'FR'))
            for nprobe = 1:length(clusters)
                selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters, which pass the set parameters
                
                depth_selected_clusters = selected_clusters; % initialize with same structure
                np = nprobe;
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
                if contains(z_score_period, 'entire_session') || contains(z_score_period, 'stim_session') 
                    zscore_counts = histcounts(all_spike_times, time_edges);
                elseif contains(z_score_period, 'first30secs')
                    baseline_spikes = all_spike_times(all_spike_times >= baseline_window(1) & all_spike_times <= baseline_window(2));
                    zscore_counts = histcounts(baseline_spikes, time_edges); %histcounts counts the number of datapoints in specified time bins
                end

                ordered_oris = unique(Task_info.stim_orientation, 'stable');
                
                % Each set of 4 stimuli makes a "trial" (i.e., sequence repeats every 4 stimuli)
                sequence_length = 4;
                stim_oris = Task_info.stim_orientation;
                stim_contrasts = Task_info.stim_contrast;
                stim_onsets = Task_info.stim_onset;

                fig = figure;
                fig.Name = sprintf('Aggregate activity: %s', depth_for_analysis);
                fig.Position = [114 90 770 650];
                %tiledlayout(4,1);

                % baseline or Session-wide stats for z-scoring
                
                mu = mean(zscore_counts);
                sigma = std(zscore_counts);

                stim_pos = 1;  % plot trace from onset of first stim in sequence
                blue_onsets = [];
                red_onsets = [];

                    % Go through all the stimuli of this type (i.e., every 4th starting from stim_pos)
                for i = stim_pos:sequence_length:length(stim_contrasts)
                    % Check if we have enough buffer for checking previous, next, two-back
                    curr = i;
                    prev = i - 1;
                    next = i + 1;
                    twoback = i - 2;

                        % Check if any surrounding stimulus (not current) has zero contrast
                    has_zero_surround = false;
                    if prev >= 1 && stim_contrasts(prev) == 0
                        has_zero_surround = true;
                    end
                    if next <= length(stim_contrasts) && stim_contrasts(next) == 0
                        has_zero_surround = true;
                    end
                    if twoback >= 1 && stim_contrasts(twoback) == 0
                        has_zero_surround = true;
                    end

                        % Check if current has zero contrast too 
                    if stim_contrasts(curr) == 0
                        has_zero_surround = true;
                    end

                        % Assign onset to red or blue
                    if has_zero_surround
                        red_onsets(end+1) = stim_onsets(curr);
                    else
                        blue_onsets(end+1) = stim_onsets(curr);
                    end
                end

                    % Compute PSTHs for red and blue
                [~, bins, ~, ~, ~, binned_red] = psthAndBA(all_spike_times, red_onsets, [-0.3 1.80], psthBinSize);
                [~, ~, ~, ~, ~, binned_blue] = psthAndBA(all_spike_times, blue_onsets, [-0.3 1.80], psthBinSize);

                % Z-score trial-by-trial, then average and compute SEM
                zscored_red = (binned_red - mu) / sigma;    % [trials x time]
                zscored_blue = (binned_blue - mu) / sigma;
                
                z_red = mean(zscored_red, 1);
                z_blue = mean(zscored_blue, 1);
                
                sem_red = std(zscored_red, 0, 1) / sqrt(size(zscored_red, 1));
                sem_blue = std(zscored_blue, 0, 1) / sqrt(size(zscored_blue, 1));

                    % Plot both traces
                    % Determine which stimuli correspond to this position
                curr_idx_list = stim_pos:sequence_length:length(stim_oris);

                    % Identify condition for each stimulus (red/green or blue)
                red_idxs = false(size(curr_idx_list));
                blue_idxs = false(size(curr_idx_list));

                for idx = 1:length(curr_idx_list)
                    curr = curr_idx_list(idx);
                    prev = curr - 1;
                    next = curr + 1;
                    twoback = curr - 2;

                    has_zero_surround = false;
                    if prev >= 1 && stim_contrasts(prev) == 0, has_zero_surround = true; end
                    if next <= length(stim_contrasts) && stim_contrasts(next) == 0, has_zero_surround = true; end
                    if twoback >= 1 && stim_contrasts(twoback) == 0, has_zero_surround = true; end
                    if stim_contrasts(curr) == 0, has_zero_surround = true; end

                    if has_zero_surround
                        red_idxs(idx) = true;
                    else
                        blue_idxs(idx) = true;
                    end
                end

                % Plot SEM shading
                fill([bins, fliplr(bins)], ...
                     [z_blue + sem_blue, fliplr(z_blue - sem_blue)], ...
                     [0.2 0.4 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % light blue
                     hold on;
                p1 = plot(bins, z_blue, 'b', 'LineWidth', 1.5); hold on;
                
                % Plot red or green depending on stimulus type
                if contains(Stimulus_type, 'OMIT')
                    fill([bins, fliplr(bins)], ...
                         [z_red + sem_red, fliplr(z_red - sem_red)], ...
                         [1 0.4 0.4], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % light red
                    p2 = plot(bins, z_red, 'r', 'LineWidth', 1.5);
                elseif contains(Stimulus_type, 'E_CD')
                    fill([bins, fliplr(bins)], ...
                         [z_red + sem_red, fliplr(z_red - sem_red)], ...
                         [0.5 1 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % light green
                    p2 = plot(bins, z_red, 'g', 'LineWidth', 1.5);
                end


                % Define grey intervals
                if strcmp(Stimulus_type, 'OMIT50grey')
                    grey_intervals = [-0.5 0; 0.15 0.2; 0.35 0.4; 0.55 0.6; 0.75 2];
                else
                    grey_intervals = [-0.5 0; 0.15 0.3; 0.45 0.6; 0.75 0.9; 1.05 2];
                end    
                
                hold on; % Make sure current plot stays visible
                
                % Get current y-axis limits for full vertical shading
                ylim([-1 6]);
                yl = ylim;
                
                % Shade each interval
                for i = 1:size(grey_intervals, 1)
                    x = [grey_intervals(i,1), grey_intervals(i,2), grey_intervals(i,2), grey_intervals(i,1)];
                    y = [yl(1), yl(1), yl(2), yl(2)];
                    fill(x, y, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % grey color with transparency
                end
                
                if strcmp(Stimulus_type, 'OMIT')
                    xline(0, 'k', (sprintf('A %d%s onset', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                end
                if strcmp(Stimulus_type, 'E_CD')
                    xline(0, 'k', sprintf('A %d%s or novel E %d%s onset', round(rad2deg(ordered_oris(1))), char(176), round(rad2deg(ordered_oris(5))), char(176)), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24); % E is always 5th orientation presented as E_CD is never the first sequence
                end    
                if strcmp(Stimulus_type, 'OMIT50grey')
                    xline(0, 'k', (sprintf('A %d%s onset', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                    xline(0.20, 'k', (sprintf('B %d%s or grey onset', round(rad2deg(ordered_oris(2))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                    xline(0.40, 'k', (sprintf('C %d%s onset', round(rad2deg(ordered_oris(3))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                    xline(0.60, 'k', (sprintf('D %d%s onset', round(rad2deg(ordered_oris(4))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                else
                    xline(0.30, 'k', (sprintf('B %d%s or grey onset', round(rad2deg(ordered_oris(2))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                    xline(0.60, 'k', (sprintf('C %d%s onset', round(rad2deg(ordered_oris(3))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                    xline(0.90, 'k', (sprintf('D %d%s onset', round(rad2deg(ordered_oris(4))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                end          

                xlim([-0.3 1.8]);
                xticks(-0.2:0.2:1.8);
                xlabel('Time (s)', 'FontSize', 24)
                ylabel(sprintf('Z-scored FR (z-scored over %s)', z_score_period), 'Interpreter', 'none');
                set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);

                    % Make legend labels
                if strcmp(Stimulus_type, 'OMIT') || strcmp(Stimulus_type, 'OMIT50grey') 
                    legend(sprintf('ABCD (%dx)', size(blue_onsets, 2)), sprintf('A_CD (%dx)', size(red_onsets, 2)), 'Location', 'northeast', 'FontSize', 24, 'Interpreter', 'none');
                end    
                if strcmp(Stimulus_type, 'E_CD')
                    legend(sprintf('ABCD (%dx)', size(blue_onsets, 2)), sprintf('%s (%dx)', Stimulus_type, size(red_onsets, 2)), 'Location', 'northeast', 'FontSize', 24, 'Interpreter', 'none');              
                end 
                
                sgtitle(sprintf('%s - %s - aggregate single unit activity: %s', subject_number, Stimulus_type, depth_for_analysis), 'Interpreter', 'none');
                %filename = sprintf('%s %s %s Aggregate single unit FRs.png', subject_number, Stimulus_type, depth_for_analysis);
                %save_path = fullfile(pwd, filename);
                %exportgraphics(fig, save_path);
                % Also save as .fig
                fig_filename = sprintf('%s %s - Aggregate single unit %s FRs.fig', subject_number, Stimulus_type, depth_for_analysis);
                fig_save_path = fullfile(pwd, fig_filename);
                savefig(fig, fig_save_path);
            end
        end



        
        if (contains(Stimulus_type, 'OMIT') || strcmp(Stimulus_type, 'E_CD')) && contains(plot_choice, 'single_units') 
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
                
                % Each set of 4 stimuli makes a "trial" (i.e., sequence repeats every 4 stimuli)
                sequence_length = 4;
                stim_oris = Task_info.stim_orientation;
                stim_contrasts = Task_info.stim_contrast;
                stim_onsets = Task_info.stim_onset;

                cluster_id = depth_selected_clusters(nprobe).cluster_id;

                baseline_window = [0 30]; % in seconds - first 30s of recording is grey screen - can use for z-scoring

                % Define time_edges depending on z_score_period
                if contains(z_score_period, 'entire_session')
                    time_edges = 0:psthBinSize:max(depth_selected_clusters(nprobe).spike_times); %time bin edges from 0 to the max time in steps of ptshBinSize
                elseif contains(z_score_period, 'first30secs')
                    time_edges = baseline_window(1):psthBinSize:baseline_window(2);
                end

                for nCluster = 1:length(cluster_id)
                    
                    spike_times_this_cluster = depth_selected_clusters(nprobe).spike_times(depth_selected_clusters(nprobe).spike_id == cluster_id(nCluster));
                    
                    % Compute appropriate histogram counts for z-scoring
                    if contains(z_score_period, 'entire_session') || contains(z_score_period, 'stim_session')
                        zscore_counts = histcounts(spike_times_this_cluster, time_edges);
                    elseif contains(z_score_period, 'first30secs')
                        baseline_spikes = spike_times_this_cluster(spike_times_this_cluster >= baseline_window(1) & spike_times_this_cluster <= baseline_window(2));
                        zscore_counts = histcounts(baseline_spikes, time_edges);
                    end

                    for stim_pos = 1 % make charts from onset of first stimulus in sequence
                        blue_onsets = [];
                        red_onsets = [];

                        % Go through all the stimuli of this type (i.e., every 4th starting from stim_pos)
                        for i = stim_pos:sequence_length:length(stim_contrasts)
                            % Check if we have enough buffer for checking previous, next, two-back
                            curr = i;
                            prev = i - 1;
                            next = i + 1;
                            twoback = i - 2;

                            % Check if any surrounding stimulus (not current) has zero contrast
                            has_zero_surround = false;
                            if prev >= 1 && stim_contrasts(prev) == 0
                                has_zero_surround = true;
                            end
                            if next <= length(stim_contrasts) && stim_contrasts(next) == 0
                                has_zero_surround = true;
                            end
                            if twoback >= 1 && stim_contrasts(twoback) == 0
                                has_zero_surround = true;
                            end

                            % Check if current has zero contrast too 
                            if stim_contrasts(curr) == 0
                                has_zero_surround = true;
                            end

                            % Assign onset to red or blue
                            if has_zero_surround
                                red_onsets(end+1) = stim_onsets(curr);
                            else
                                blue_onsets(end+1) = stim_onsets(curr);
                            end
                        end

                        % Compute PSTHs for red and blue
                        [~, bins, ~, ~, ~, binned_red] = psthAndBA(spike_times_this_cluster, red_onsets, [-0.3 1.8], psthBinSize);
                        [~, ~, ~, ~, ~, binned_blue] = psthAndBA(spike_times_this_cluster, blue_onsets, [-0.3 1.8], psthBinSize);
                      
                        fig1 = figure; %fig1 is for raster plots
                        fig1.Name = sprintf('%s - Spiking raster %s Cluster %i', Stimulus_type, depth_for_analysis, cluster_id(nCluster));
                        fig1.Position = [114 90 770 650];
                        
                        tiledlayout(100, 1); %Use tiledlayout with weighted row heights so the plot for the blue condition (160 trials) can be 4x taller than that for the red condition (40 trials)
                        
                        % Set raster window
                        rasterWindow = [-0.3 1.80];
                        
                        % Get raster for blue
                        [~, ~, rasterX_blue, rasterY_blue, ~, ~] = psthAndBA(spike_times_this_cluster, blue_onsets, rasterWindow, psthBinSize/10);
                        nexttile([74 1]); 
                        
                        plot(rasterX_blue, rasterY_blue, 'Color', [0, 0, 0.6]); %dark blue
                        xlim([-0.3 1.8]);
                        xticks(-0.2:0.2:1.8);
                        xlabel('Time (s)');
                        ylim([0 length(blue_onsets)]);
                        ylabel('ABCD Trial');
                        set(gca, 'TickDir', 'out', 'box', 'off', 'Color', 'none', 'FontSize', 24);

                        ordered_oris = unique(Task_info.stim_orientation, 'stable'); % orientations in order of first appearance (NB an ABCD trial always starts the block) radians for TRAIN but degrees for GAVNIK
                        xline(0, 'k', (sprintf('A %d%s onset', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                        xline(0.30, 'k', (sprintf('B %d%s onset', round(rad2deg(ordered_oris(2))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                        xline(0.60, 'k', (sprintf('C %d%s onset', round(rad2deg(ordered_oris(3))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                        xline(0.90, 'k', (sprintf('D %d%s onset', round(rad2deg(ordered_oris(4))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                                    
                        % Get raster for red (or green if E_CD)
                        [~, ~, rasterX_red, rasterY_red, ~, ~] = psthAndBA(spike_times_this_cluster, red_onsets, rasterWindow, psthBinSize/10);
                        nexttile([26 1]); 
                        
                        color_raster = [0.6, 0, 0]; %dark red
                        if contains(Stimulus_type, 'E_CD'), color_raster = [0, 0.5, 0]; end %dark green
                        plot(rasterX_red, rasterY_red, 'Color', color_raster);
                        xlim([-0.3 1.8]);
                        xticks(-0.2:0.2:1.8);
                        xlabel('Time (s)');
                        ylim([0 length(red_onsets)]);
                        ylabel('Test Trial');
                        set(gca, 'TickDir', 'out', 'box', 'off', 'Color', 'none', 'FontSize', 24);
                                              
                        if contains(Stimulus_type, 'E_CD')
                            firststim = 'E';
                            firstori = 5;
                        elseif contains(Stimulus_type, 'OMIT')
                            firststim = 'A';
                            firstori = 1;
                        end

                        xline(0, 'k', (sprintf('%s %d%s onset', firststim, round(rad2deg(ordered_oris(firstori))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                        xline(0.30, 'k', (sprintf('grey onset')), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                        xline(0.60, 'k', (sprintf('C %d%s onset', round(rad2deg(ordered_oris(3))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                        xline(0.90, 'k', (sprintf('D %d%s onset', round(rad2deg(ordered_oris(4))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                             
                        sgtitle(sprintf('%s - ABCD vs Omission condition: Raster plots %s Cluster %i', subject_number, depth_for_analysis, cluster_id(nCluster)), 'Interpreter', 'none');
                        
                        fig_filename1 = sprintf('%s - %s - %s Cluster_%i_raster.fig', subject_number, Stimulus_type, depth_for_analysis, cluster_id(nCluster));
                        savefig(fig1, fullfile(pwd, fig_filename1));

                        % Plot both traces
                        fig2 = figure; %fig2 is for FR traces
                        fig2.Name = sprintf('%s FR response %s Cluster %i', Stimulus_type, depth_for_analysis, cluster_id(nCluster));
                        fig2.Position = [114 90 770 650];
                        
                        % Get orientation/contrast for legend
                        curr_idx_list = stim_pos:sequence_length:length(stim_oris);
                        red_idxs = false(size(curr_idx_list));
                        blue_idxs = false(size(curr_idx_list));
                        for idx = 1:length(curr_idx_list)
                            curr = curr_idx_list(idx);
                            prev = curr - 1;
                            next = curr + 1;
                            twoback = curr - 2;

                            has_zero_surround = false;
                            if prev >= 1 && stim_contrasts(prev) == 0, has_zero_surround = true; end
                            if next <= length(stim_contrasts) && stim_contrasts(next) == 0, has_zero_surround = true; end
                            if twoback >= 1 && stim_contrasts(twoback) == 0, has_zero_surround = true; end
                            if stim_contrasts(curr) == 0, has_zero_surround = true; end

                            if has_zero_surround
                                red_idxs(idx) = true;
                            else
                                blue_idxs(idx) = true;
                            end
                        end

                        ori_blue = unique(stim_oris(curr_idx_list(blue_idxs)));
                        contrast_blue = unique(stim_contrasts(curr_idx_list(blue_idxs)));
                        ori_red = unique(stim_oris(curr_idx_list(red_idxs)));
                        contrast_red = unique(stim_contrasts(curr_idx_list(red_idxs)));

                        % Format legend labels
                        blue_label = sprintf('ABCD %dx', length(blue_onsets));
                        red_label = sprintf('%s %dx', Stimulus_type, length(red_onsets));
                        % Plot
                        nexttile;
                        if contains(z_score_period, 'none')
                            % Raw mean and SEM (Hz)
                            mean_blue = mean(binned_blue, 1) / psthBinSize;
                            sem_blue = std(binned_blue, 0, 1) / sqrt(size(binned_blue, 1)) / psthBinSize;
                        
                            mean_red = mean(binned_red, 1) / psthBinSize;
                            sem_red = std(binned_red, 0, 1) / sqrt(size(binned_red, 1)) / psthBinSize;
                        
                            % SEM shading
                            fill([bins, fliplr(bins)], [mean_blue + sem_blue, fliplr(mean_blue - sem_blue)], ...
                                 [0, 0, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off'); hold on;
                            p1 = plot(bins, mean_blue, 'b', 'LineWidth', 1.5);
                        
                            if contains(Stimulus_type, 'OMIT')
                                fill([bins, fliplr(bins)], [mean_red + sem_red, fliplr(mean_red - sem_red)], ...
                                     [1, 0, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                                p2 = plot(bins, mean_red, 'r', 'LineWidth', 1.5);
                            elseif contains(Stimulus_type, 'E_CD')
                                fill([bins, fliplr(bins)], [mean_red + sem_red, fliplr(mean_red - sem_red)], ...
                                     [0, 0.5, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                                p2 = plot(bins, mean_red, 'g', 'LineWidth', 1.5);
                            end
                        
                            ylabel('Mean FR across trials (Hz)', 'FontSize', 24);
                        
                        else
                            % Session-wide stats for z-scoring
                            mu = mean(zscore_counts);
                            sigma = std(zscore_counts);

                            zscored_blue_trials = (binned_blue - mu) / sigma;
                            zscored_red_trials  = (binned_red  - mu) / sigma;
                            
                            z_blue = mean(zscored_blue_trials, 1);
                            sem_blue = std(zscored_blue_trials, 0, 1) / sqrt(size(zscored_blue_trials, 1));
                            
                            z_red = mean(zscored_red_trials, 1);
                            sem_red = std(zscored_red_trials, 0, 1) / sqrt(size(zscored_red_trials, 1));
                        
                            % SEM shading
                            fill([bins, fliplr(bins)], [z_blue + sem_blue, fliplr(z_blue - sem_blue)], ...
                                 [0, 0, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off'); hold on;
                            p1 = plot(bins, z_blue, 'b', 'LineWidth', 1.5);
                        
                            if contains(Stimulus_type, 'OMIT')
                                fill([bins, fliplr(bins)], [z_red + sem_red, fliplr(z_red - sem_red)], ...
                                     [1, 0, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                                p2 = plot(bins, z_red, 'r', 'LineWidth', 1.5);
                            elseif contains(Stimulus_type, 'E_CD')
                                fill([bins, fliplr(bins)], [z_red + sem_red, fliplr(z_red - sem_red)], ...
                                     [0, 0.5, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                                p2 = plot(bins, z_red, 'g', 'LineWidth', 1.5);
                            end
                        
                            if contains(z_score_period, 'entire_session')
                                ylabel('Z-scored FR (z-scored over entire session)', 'FontSize', 24);
                            elseif contains(z_score_period, 'stim_session')
                                ylabel('FR Z-scored over stim session)', 'FontSize', 24);     
                            elseif contains(z_score_period, 'first30secs')
                                ylabel('Z-scored FR (z-scored over first 30s baseline)', 'FontSize', 24);
                            end
                        end

                        xline(0,'k','LineWidth',1);
                        xlim([-0.3 1.8]);
                        xticks(-0.2:0.2:1.8);
                        xlabel('Time (s)');
                        
                        legend([p1 p2], {blue_label, red_label}, 'Location', 'northeast', 'Interpreter', 'none');  
                        legend boxoff;
                        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);

                                % Define grey intervals
                        grey_intervals = [-0.5 0; 0.15 0.3; 0.45 0.6; 0.75 0.9; 1.05 2];
                            
                                      
                        % Get current y-axis limits for full vertical shading
                        yl = ylim;
                            
                        % Shade each interval
                        for i = 1:size(grey_intervals, 1)
                            x = [grey_intervals(i,1), grey_intervals(i,2), grey_intervals(i,2), grey_intervals(i,1)];
                            y = [yl(1), yl(1), yl(2), yl(2)];
                            fill(x, y, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % grey color with transparency
                        end
                        
                        if contains(Stimulus_type, 'E_CD')
                            firststim = sprintf('or E %d%s onset', round(rad2deg(ordered_oris(5))), char(176));
                        elseif contains(Stimulus_type, 'OMIT')
                            firststim = 'onset';
                        end

                        xline(0, 'k', (sprintf('A %d%s %s', round(rad2deg(ordered_oris(1))), char(176), firststim)), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                        xline(0.30, 'k', (sprintf('B %d%s or grey onset', round(rad2deg(ordered_oris(2))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                        xline(0.60, 'k', (sprintf('C %d%s onset', round(rad2deg(ordered_oris(3))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                        xline(0.90, 'k', (sprintf('D %d%s onset', round(rad2deg(ordered_oris(4))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                             
                        sgtitle(sprintf('%s - %s: FR %s Cluster %i', subject_number, Stimulus_type, depth_for_analysis, cluster_id(nCluster)), 'Interpreter', 'none');
                        fig_filename2 = sprintf('%s - %s - %s Cluster_%i_FR.fig', subject_number, Stimulus_type, depth_for_analysis, cluster_id(nCluster));
                        savefig(fig2, fullfile(pwd, fig_filename2));

                    end                               
                end
            end
        end



        
        if strcmp(Stimulus_type, 'DCBA')... 
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
                if contains(z_score_period, 'entire_session') || contains(z_score_period, 'stim_session')
                    zscore_counts = histcounts(all_spike_times, time_edges);
                elseif contains(z_score_period, 'first30secs')
                    baseline_spikes = all_spike_times(all_spike_times >= baseline_window(1) & all_spike_times <= baseline_window(2));
                    zscore_counts = histcounts(baseline_spikes, time_edges); %histcounts counts the number of datapoints in specified time bins
                end

                % Define expected normal and reversed orientation sequences
                normal_sequence = round(rad2deg(unique(Task_info.stim_orientation, 'stable')))'; % test blocks per my protocol always begin with the trained sequence
                reversed_sequence = fliplr(normal_sequence);

                sequence_length = 4;
                num_sequences = floor(length(Task_info.stim_orientation) / sequence_length);
                normal_onsets_by_pos = cell(1, sequence_length);
                reversed_onsets_by_pos = cell(1, sequence_length);

                for seq = 1:num_sequences
                    idx_range = (seq - 1) * sequence_length + (1:sequence_length);
                    oris = Task_info.stim_orientation(idx_range);
                    onsets = Task_info.stim_onset(idx_range);

                    % Round to degrees to avoid floating-point errors
                    rounded_oris = round(rad2deg(oris(:)))';

                    if isequal(rounded_oris, normal_sequence)
                        for pos = 1:sequence_length
                            normal_onsets_by_pos{pos}(end+1) = onsets(pos);
                        end
                    elseif isequal(rounded_oris, reversed_sequence)
                        for pos = 1:sequence_length
                            reversed_onsets_by_pos{pos}(end+1) = onsets(pos);
                        end
                    end
                end

                fig = figure;
                fig.Name = sprintf('Aggregate activity: %s depth range %d–%d μm', depth_for_analysis, min(depth_range), max(depth_range));
                fig.Position = [114 90 770 650];
                
                % baseline or Session-wide stats
                mu = mean(zscore_counts);
                sigma = std(zscore_counts);

                pos = 1;
                % Normal trials (blue)
                [~, bins, ~, ~, ~, ba_normal] = psthAndBA(all_spike_times, normal_onsets_by_pos{pos}, [-0.3 1.80], psthBinSize);
                %z_normal = (mean(ba_normal,1) - mu) / sigma;
                % Reversed trials (red)
                [~, ~, ~, ~, ~, ba_reversed] = psthAndBA(all_spike_times, reversed_onsets_by_pos{pos}, [-0.3 1.80], psthBinSize);
                %z_reversed = (mean(ba_reversed,1) - mu) / sigma;
                
                % Z-score each trial for normal and reversed
                zscored_normal = (ba_normal - mu) / sigma;
                zscored_reversed = (ba_reversed - mu) / sigma;
                
                % Mean and SEM across trials
                z_normal = mean(zscored_normal, 1);
                sem_normal = std(zscored_normal, 0, 1) / sqrt(size(zscored_normal, 1));
                
                z_reversed = mean(zscored_reversed, 1);
                sem_reversed = std(zscored_reversed, 0, 1) / sqrt(size(zscored_reversed, 1));
                
                % Plot SEM shading
                fill([bins, fliplr(bins)], ...
                     [z_normal + sem_normal, fliplr(z_normal - sem_normal)], ...
                     [0.6 0.8 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % lighter blue
                hold on;     
                plot(bins, z_normal, 'b', 'LineWidth', 1.5); hold on;
                
                fill([bins, fliplr(bins)], ...
                     [z_reversed + sem_reversed, fliplr(z_reversed - sem_reversed)], ...
                     [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % lighter red
                plot(bins, z_reversed, 'r', 'LineWidth', 1.5); hold on;

                % Define grey intervals
                grey_intervals = [-0.5 0; 0.15 0.3; 0.45 0.6; 0.75 0.9; 1.05 2];
                
                hold on; % Make sure current plot stays visible
                
                % Get current y-axis limits for full vertical shading
                ylim([-1 6]);
                yl = ylim;
                
                % Shade each interval
                for i = 1:size(grey_intervals, 1)
                    x = [grey_intervals(i,1), grey_intervals(i,2), grey_intervals(i,2), grey_intervals(i,1)];
                    y = [yl(1), yl(1), yl(2), yl(2)];
                    fill(x, y, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % grey color with transparency
                end

                xline(0, 'k', (sprintf('A %d%s or D onset', normal_sequence(1), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                xline(0.30, 'k', (sprintf('B %d%s or C onset', normal_sequence(2), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                xline(0.60, 'k', (sprintf('C %d%s or B onset', normal_sequence(3), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);
                xline(0.90, 'k', (sprintf('D %d%s or A onset', normal_sequence(4), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'FontSize', 24);

                xlim([-0.3 1.8]);
                xticks(-0.2:0.2:1.8);
                xlabel('Time (s)', 'Fontsize', 24)
                ylabel(sprintf('Z-scored FR (z-scored over %s)', z_score_period), 'Interpreter', 'none');
                    
                legend(sprintf('ABCD (%dx)', size(normal_onsets_by_pos{1}, 2)), sprintf('DCBA (%dx)', size(reversed_onsets_by_pos{1}, 2)), 'Location', 'northeast', 'Box', 'on');
                set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);
                
                sgtitle(sprintf('%s - %s - aggregate single unit activity, %s depth range %d–%d μm', subject_number, Stimulus_type, depth_for_analysis, min(depth_range), max(depth_range)), 'Interpreter', 'none', 'FontSize', 24);

                %filename = sprintf('%s %s %s Aggregate single unit FRs.png', subject_number, Stimulus_type, depth_for_analysis);
                %save_path = fullfile(pwd, filename);
                %exportgraphics(fig, save_path);
                % Also save as .fig
                fig_filename = sprintf('%s %s - Aggregate single unit %s FRs.fig', subject_number, Stimulus_type, depth_for_analysis);
                fig_save_path = fullfile(pwd, fig_filename);
                savefig(fig, fig_save_path);
            end
        end

        


        if (contains(Stimulus_type, 'GAVNIK_ABCD') || contains(Stimulus_type, 'GAVNIK200_ABCD') || contains(Stimulus_type, 'GAVNIK250_ABCD') || contains(Stimulus_type, 'GAVNIK DCBA')) && contains(plot_choice, 'aggregate') && contains(plot_type, 'FR')
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
                if contains(z_score_period, 'entire_session') || contains(z_score_period, 'stim_session')
                    zscore_counts = histcounts(all_spike_times, time_edges);
                elseif contains(z_score_period, 'first30secs')
                    baseline_spikes = all_spike_times(all_spike_times >= baseline_window(1) & all_spike_times <= baseline_window(2));
                    zscore_counts = histcounts(baseline_spikes, time_edges); %histcounts counts the number of datapoints in specified time bins
                end

                ordered_oris = unique(Task_info.stim_orientation, 'stable');               
                fig = figure;
                fig.Name = sprintf('%s: Aggregate activity: %s', Stimulus_type, depth_for_analysis);
                fig.Position = [114 90 770 650];
                
                ori = 1;
                stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));

                if (contains(Stimulus_type, 'GAVNIK250_ABCD'))
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.3 2.0], psthBinSize);
                elseif (contains(Stimulus_type, 'GAVNIK200_ABCD'))
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.3 1.8], psthBinSize);
                else
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.3 1.35], psthBinSize);
                end

                mean_trace = mean(binnedArray, 1); % average over trials
                z_trace = (mean_trace - mean(zscore_counts)) / std(zscore_counts); % z-score using session-wide stats

                plot(bins, z_trace, 'b', 'LineWidth', 1.5);

                % Define windows after time zero and calc. the peak and mean FR [for exploration]             
                if (contains(Stimulus_type, 'GAVNIK250_ABCD'))
                    stim_window_starts = 0:0.25:0.75;
                    stim_window_ends = 0.25:0.25:1.0;
                    grey_window_starts = 1.0:0.25:1.75; % look in each 250ms window
                    grey_window_ends = 1.25:0.25:2.0; % look in each 250ms window
                elseif (contains(Stimulus_type, 'GAVNIK200_ABCD'))
                    stim_window_starts = 0:0.2:0.6;
                    stim_window_ends = 0.2:0.2:0.8;
                    grey_window_starts = 0.8:0.2:1.6; % look in each 200ms window
                    grey_window_ends = 1.0:0.2:1.8; % look in each 200ms window
                else
                    stim_window_starts = 0:0.15:0.45;
                    stim_window_ends = 0.15:0.15:0.6;
                    grey_window_starts = 0.6:0.15:1.2; % look in each 150ms window
                    grey_window_ends = 0.75:0.15:1.35; % look in each 150ms window
                end  
                peak_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                mean_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                peak_FR_by_greywindow = zeros(1, length(stim_window_starts)); % preallocate
                mean_FR_by_greywindow = zeros(1, length(stim_window_starts)); % preallocate

                                                  
                for i = 1:length(stim_window_starts)
                    % Find indices of bins within current window
                    idx_in_stimwindow = bins >= stim_window_starts(i) & bins < stim_window_ends(i);                    
                    peak_FR_by_stimwindow(i) = max(z_trace(idx_in_stimwindow));
                    mean_FR_by_stimwindow(i) = mean(z_trace(idx_in_stimwindow));                    
                end
                
                for i = 1:length(grey_window_starts)
                    % Find indices of bins within current window
                    idx_in_greywindow = bins >= grey_window_starts(i) & bins < grey_window_ends(i);                    
                    peak_FR_by_greywindow(i) = max(z_trace(idx_in_greywindow));
                    mean_FR_by_greywindow(i) = mean(z_trace(idx_in_greywindow));
                end
                
                % Define grey intervals
                if (contains(Stimulus_type, 'GAVNIK250_ABCD'))
                    grey_intervals = [-0.5 0; 1.0 2.0];
                elseif (contains(Stimulus_type, 'GAVNIK200_ABCD'))
                    grey_intervals = [-0.5 0; 0.8 1.8];
                else
                    grey_intervals = [-0.5 0; 0.6 1.5];
                end    

                hold on; % Make sure current plot stays visible
                
                % Get current y-axis limits for full vertical shading
                ylim([-1 5]);
                yl = ylim;
                
                % Shade each interval
                for i = 1:size(grey_intervals, 1)
                    x = [grey_intervals(i,1), grey_intervals(i,2), grey_intervals(i,2), grey_intervals(i,1)];
                    y = [yl(1), yl(1), yl(2), yl(2)];
                    fill(x, y, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % grey color with transparency
                end

                xline(stim_window_starts(1), 'k', (sprintf('A onset; %d%s', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(stim_window_starts(2), 'k', (sprintf('B onset; %d%s', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(stim_window_starts(3), 'k', (sprintf('C onset; %d%s', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(stim_window_starts(4), 'k', (sprintf('D onset; %d%s', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                
                if (contains(Stimulus_type, 'GAVNIK250_ABCD'))
                    xlim([-0.3 2.0]);
                    xticks(-0.2:0.2:2.0);
                elseif (contains(Stimulus_type, 'GAVNIK200_ABCD'))
                    xlim([-0.3 1.8]);
                    xticks(-0.2:0.2:1.8);    
                else
                    xlim([-0.3 1.35]);
                    xticks(-0.2:0.2:1.2);
                end

                ylabel('Z-scored FR');
                xlabel('Time since onset of A (s)')
                set(gca, "TickDir", "out", 'Color', 'none', 'FontSize', 24);
                
                sgtitle(sprintf('%s - %s day %d - Aggregate single unit activity: %s', subject_number, Stimulus_type, experiment_info(nsession).date, depth_for_analysis), 'Interpreter', 'none');

                %filename = sprintf('%s - Aggregate single unit %s FRs.png', Stimulus_type, depth_for_analysis);
                %save_path = fullfile(pwd, filename);
                %exportgraphics(fig, save_path);  

                % Also save as .fig
                fig_filename = sprintf('%s - Aggregate single unit %s FRs.fig', Stimulus_type, depth_for_analysis);
                fig_save_path = fullfile(pwd, fig_filename);
                savefig(fig, fig_save_path);

            end
        end




        if (contains(Stimulus_type, 'GAVNIK_ABCD') || contains(Stimulus_type, 'GAVNIK200_ABCD') || contains(Stimulus_type, 'GAVNIK250_ABCD') || contains(Stimulus_type, 'GAVNIK_DCBA')) &&...
           contains(plot_choice, 'single_units') 
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
                    time_edges = 0:psthBinSize:max(depth_selected_clusters(nprobe).spike_times);
                elseif contains(z_score_period, 'first30secs')
                    time_edges = baseline_window(1):psthBinSize:baseline_window(2);
                end
                
                %ori_response = [];
                %z_binnedArray = [];

                for nCluster = 1:length(cluster_id) % loop through each good cluster at selected depth
                
                    spike_times_this_cluster = depth_selected_clusters(nprobe).spike_times(depth_selected_clusters(nprobe).spike_id == cluster_id(nCluster)); % extract the spike times for this cluster
                    
                    % Compute appropriate histogram counts for z-scoring
                    if contains(z_score_period, 'entire_session') || contains(z_score_period, 'stim_session')
                        zscore_counts = histcounts(spike_times_this_cluster, time_edges);
                    elseif contains(z_score_period, 'first30secs')
                        baseline_spikes = spike_times_this_cluster(spike_times_this_cluster >= baseline_window(1) & spike_times_this_cluster <= baseline_window(2));
                        zscore_counts = histcounts(baseline_spikes, time_edges);
                    end    

                    fig(nCluster)=figure; % open a figure window
                    fig(nCluster).Name=sprintf('%s Grating responses Cluster %i', Stimulus_type, cluster_id(nCluster)); %overall figure title includes the cluster_id (one count higher than zero-based SI output)
                    fig(nCluster).Position = [114 90 770 650]; % sets the size of the figure window
                
                    % Extract ordered orientations based on presentation sequence
                    ordered_oris = unique(Task_info.stim_orientation, 'stable'); % radians for TRAIN but degrees for GAVNIK
                    
                    for ori = 1 % plot from onset of first stimulus in sequence 

                        if contains(plot_type, 'raster')
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.3 1.35], psthBinSize/10); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
                            plot(rasterX,rasterY,'k','LineWidth',1) % plot the spiking raster for the current orientation
                            xline(0,'r',LineWidth=1) % put a vertical red line at time zero (stimulus onset)
                            xlim([-0.2 1.2])
                            xticks(0:0.2:1.2);
                            ylim([0 sum(Task_info.stim_orientation==ordered_oris(ori))]) % set the max y coord to be the total # trials with this orientation
                            ylabel('Trial');
                                                        
                            sgtitle(sprintf('%s %s cluster %d rasters %s', subject_number, depth_for_analysis, cluster_id, Stimulus_type, 'Interpreter', 'none', 'FontSize', 24));
                            set(gca,"TickDir","out", 'Color','none','FontSize',24)
                        end
                        
                        
                        if contains(plot_type, 'FR')
                            if contains(Stimulus_type, 'GAVNIK250_ABCD')
                                [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.3 2.0], psthBinSize); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
                            elseif contains(Stimulus_type, 'GAVNIK200_ABCD')
                                [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.3 1.8], psthBinSize); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
                            else
                                [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.3 1.35], psthBinSize); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
                            end

                            if contains(z_score_period, 'none')
                                mean_trace = mean(binnedArray, 1) / psthBinSize; % average over trials
                                sem_trace = std(binnedArray, 0, 1) / sqrt(size(binnedArray, 1))/ psthBinSize;  % [1 x time]
                                % Plot SEM shading
                                fill([bins, fliplr(bins)], ...
                                     [mean_trace + sem_trace, fliplr(mean_trace - sem_trace)], ...
                                     'b', ...
                                     'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                                hold on;
                                    % Define grey intervals
                                if contains(Stimulus_type, 'GAVNIK250_ABCD')
                                    grey_intervals = [-0.5 0; 1.0 2.0];
                                elseif contains(Stimulus_type, 'GAVNIK200_ABCD')
                                    grey_intervals = [-0.5 0; 0.8 1.8];    
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
                                plot(bins, mean_trace, 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', sprintf('%s (200x)', Stimulus_type));
                                ylabel('Mean firing rate across trials (Hz)', 'FontSize', 24);
                                %ylim ([0 16]);
                            else 
                                baseline_mean = mean(zscore_counts);
                                baseline_std = std(zscore_counts);
                                % Z-score each trial individually so that trial to variability remains visible and SEM can be calculated
                                zscored_trials = (binnedArray - baseline_mean) / baseline_std;  % [trials x time]
                                z_trace = mean(zscored_trials, 1);                             % mean z-scored trace
                                sem_trace = std(zscored_trials, 0, 1) / sqrt(size(zscored_trials, 1));  % SEM

                                if contains(z_score_period, 'entire_session')
                                    ylabel('Z-scored FR (z-scored over entire session)', 'FontSize', 24);
                                    ylim([-1 5]);
                                elseif contains(z_score_period, 'stim_session')
                                    ylabel('FR Z-scored over stim session', 'FontSize', 24);     
                                elseif contains(z_score_period, 'first30secs')
                                    ylabel('Z-scored FR (z-scored over first 30s baseline)', 'FontSize', 24);
                                    ylim([-1 8]);
                                end
                                     % Define grey intervals
                                if (contains(Stimulus_type, 'GAVNIK250_ABCD'))     
                                    grey_intervals = [-0.5 0; 1.0 2.0];
                                elseif (contains(Stimulus_type, 'GAVNIK200_ABCD'))     
                                    grey_intervals = [-0.5 0; 0.8 1.8];    
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
                                
                                plot(bins, z_trace, 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', sprintf('%s (200x)', Stimulus_type));
                                hold on; 
                                % Plot SEM shading
                                fill([bins, fliplr(bins)], ...
                                    [z_trace + sem_trace, fliplr(z_trace - sem_trace)], ...
                                    'b', ...
                                    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
                                
                            end   
                            if (contains(Stimulus_type, 'GAVNIK250_ABCD'))
                                xlim([-0.3 2.0]);
                                xticks(-0.2:0.2:2.0);
                            elseif (contains(Stimulus_type, 'GAVNIK200_ABCD'))
                                xlim([-0.3 1.8]);
                                xticks(-0.2:0.2:1.8);    
                            else
                                xlim([-0.3 1.5]);
                                xticks(-0.2:0.2:1.4);
                            end    
                            xlabel('Time (s)', 'FontSize', 24)
                            set(gca, 'FontSize', 24);  % Tick labels font size
                            legend(flipud(findobj(gca,'-property','DisplayName')), 'Location', 'northeast', 'Interpreter', 'none');
                            hold on;                
                        
                                                       
                            xline(0, 'k', (sprintf('A %d%s onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                            if (contains(Stimulus_type, 'GAVNIK250_ABCD'))
                                xline(0.25, 'k', (sprintf('B %d%s onset', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                                xline(0.50, 'k', (sprintf('C %d%s onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                                xline(0.75, 'k', (sprintf('D %d%s onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                            elseif (contains(Stimulus_type, 'GAVNIK200_ABCD'))
                                xline(0.2, 'k', (sprintf('B %d%s onset', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                                xline(0.4, 'k', (sprintf('C %d%s onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                                xline(0.6, 'k', (sprintf('D %d%s onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);    
                            else
                                xline(0.15, 'k', (sprintf('B %d%s onset', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                                xline(0.30, 'k', (sprintf('C %d%s onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                                xline(0.45, 'k', (sprintf('D %d%s onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                            end    
                            sgtitle(sprintf('%s %s cluster %d FRs %s', subject_number, depth_for_analysis, cluster_id(nCluster), Stimulus_type), 'Interpreter', 'none', 'FontSize', 24);
                            fig_filename = sprintf('%s %s cluster %d FRs %s (zscore %s).fig', subject_number, depth_for_analysis, cluster_id(nCluster), Stimulus_type, z_score_period);
                            fig_save_path = fullfile(pwd, fig_filename);
                            savefig(fig(nCluster), fig_save_path);
                            
                        end
                                          
                    end
                end
            end
        end



        
        if (contains(Stimulus_type, 'GAVNIK_A_CD') || contains(Stimulus_type, 'GAVNIK_E_CD')) && contains(plot_choice, 'aggregate') && contains(plot_type, 'FR')
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
                if contains(z_score_period, 'entire_session') || contains(z_score_period, 'stim_session')
                    zscore_counts = histcounts(all_spike_times, time_edges);
                elseif contains(z_score_period, 'first30secs')
                    baseline_spikes = all_spike_times(all_spike_times >= baseline_window(1) & all_spike_times <= baseline_window(2));
                    zscore_counts = histcounts(baseline_spikes, time_edges); %histcounts counts the number of datapoints in specified time bins
                end


                ordered_oris = unique(Task_info.stim_orientation, 'stable');               
                
                fig = figure;
                fig.Name = sprintf('%s: Aggregate activity: %s depth range %d–%d μm', Stimulus_type, depth_for_analysis, min(depth_range), max(depth_range));
                fig.Position = [114 90 770 650];
                
                ori = 1; % plot from onset of first stimulus in sequence
                stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));

                [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.3 1.35], psthBinSize);
                
                mean_trace = mean(binnedArray, 1); % average over trials

                if contains(z_score_period, 'none')
                    % Plot raw mean firing rate trace
                    plot(bins, mean_trace, 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', sprintf('%s (200x)', Stimulus_type));
                    ylabel('Mean firing rate (Hz)', 'FontSize', 24);
                    ylim ([0 16]);
                else 
                    baseline_mean = mean(zscore_counts);
                    baseline_std = std(zscore_counts);
                    % Z-score each trial individually so that trial to variability remains visible and SEM can be calculated
                    zscored_trials = (binnedArray - baseline_mean) / baseline_std;  % [trials x time]

                    if strcmp(sliced_plot_option, 'yes')
                        num_trials = size(binnedArray, 1);
                        num_groups = 5;
                        group_size = num_trials/num_groups;
                        % Define manually chosen colors for each group
                        colors = [
                            1.0,  0.8,  0.0;  % 1st group - golden yellow (1.0,  1.0,  0.0; for pure yellow)
                            0.4,  0.8,  0.0;  % 2nd group - yellow-green (more yellowish)
                            0.2,  0.6,  0.4;  % 3rd group - pure green
                            0.2,  0.4,  0.6;  % 4th group - greenish blue (not too blue)
                            0.0,  0.0,  0.8;  % 5th group - deep blue
                        ];
                        for g = 1:num_groups
                            idx_start = (g - 1) * group_size + 1;
                            idx_end = g * group_size;
                        
                            group_trials = binnedArray(idx_start:idx_end, :);
                            
                            % Z-score this group of trials (use the same baseline if specified globally)
                            zscored_group = (group_trials - baseline_mean) / baseline_std;
                            z_trace = mean(zscored_group, 1);
                            hold on;
                            plot(bins, z_trace, 'Color', colors(g, :), 'LineWidth', 1.5, ...
                                 'DisplayName', sprintf('%s trials %d–%d', Stimulus_type, idx_start, idx_end));

                        end
                    else


                        z_trace = mean(zscored_trials, 1);                             % mean z-scored trace
                        sem_trace = std(zscored_trials, 0, 1) / sqrt(size(zscored_trials, 1));  % SEM
                        
                         % Plot SEM shading
                        hold on;
                        fill([bins, fliplr(bins)], ...
                             [z_trace + sem_trace, fliplr(z_trace - sem_trace)], ...
                             [1 0.4 0.4], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % light red
                        hold on;
                        plot(bins, z_trace, 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', sprintf('%s (200x)', Stimulus_type));
                    end    
                    if contains(z_score_period, 'entire_session')
                        ylabel('Z-scored FR (z-scored over entire session)', 'FontSize', 24);
                        ylim([-1 5]);
                    elseif contains(z_score_period, 'stim_session')
                        ylabel('FR Z-scored over stim session)', 'FontSize', 24);     
                    elseif contains(z_score_period, 'first30secs')
                        ylabel('Z-scored FR (z-scored over first 30s baseline)', 'FontSize', 24);
                        ylim([-1 8]);
                    end                   
                end  


                xlim([-0.3 1.5]);
                ylim([-1 5]);
                xticks(-0.2:0.2:1.4);
                xlabel('Time (s)', 'FontSize', 24)
                set(gca, 'FontSize', 24);  % Tick labels font size
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
                
                xline(0, 'k', (sprintf('A %d%s onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(0.15, 'k', (sprintf('grey onset')), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(0.30, 'k', (sprintf('C %d%s onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                xline(0.45, 'k', (sprintf('D %d%s onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                
                sgtitle(sprintf('%s - %s: Aggregate single unit activity: %s depth range %d–%d μm', subject_number, Stimulus_type, depth_for_analysis, min(depth_range), max(depth_range)), 'Interpreter', 'none');
                %filename = sprintf('%s L4 FRs %s.png', subject_number, Stimulus_type);
                %save_path = fullfile(pwd, filename);
                %exportgraphics(fig, save_path); 

                % Also save as .fig
                fig_filename = sprintf('%s L4 FRs %s.fig', subject_number, Stimulus_type);
                fig_save_path = fullfile(pwd, fig_filename);
                savefig(fig, fig_save_path);
                           

                
            end
        end



        
        if (contains(Stimulus_type, 'GAVNIK_A_CD') || contains(Stimulus_type, 'GAVNIK_E_CD')) &&...
           contains(plot_choice, 'single_units') 
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
                    time_edges = 0:psthBinSize:max(depth_selected_clusters(nprobe).spike_times);
                elseif contains(z_score_period, 'first30secs')
                    time_edges = baseline_window(1):psthBinSize:baseline_window(2);
                end
                
                ori_response = [];
                z_binnedArray = [];

                for nCluster = 1:length(cluster_id) % loop through each good cluster at selected depth
                
                    spike_times_this_cluster = depth_selected_clusters(nprobe).spike_times(depth_selected_clusters(nprobe).spike_id == cluster_id(nCluster)); % extract the spike times for this cluster
                    %[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset, [-0.150 0.30], 0.001); % gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
             
                    % Compute appropriate histogram counts for z-scoring
                    if contains(z_score_period, 'entire_session') || contains(z_score_period, 'stim_session')
                        zscore_counts = histcounts(spike_times_this_cluster, time_edges);
                    elseif contains(z_score_period, 'first30secs')
                        baseline_spikes = spike_times_this_cluster(spike_times_this_cluster >= baseline_window(1) & spike_times_this_cluster <= baseline_window(2));
                        zscore_counts = histcounts(baseline_spikes, time_edges);
                    end

                    fig(nCluster)=figure; % open a figure window
                    fig(nCluster).Name=sprintf('%s Grating responses Cluster %i', Stimulus_type, cluster_id(nCluster)); %overall figure title includes the cluster_id (one count higher than zero-based SI output)
                    fig(nCluster).Position = [114 90 770 650]; % sets the size of the figure window
                
                    tiledlayout(5,1); % vertical layout (5 rows × 1 column)
                    % Extract ordered orientations based on presentation sequence
                    ordered_oris = unique(Task_info.stim_orientation, 'stable'); % radians for TRAIN but degrees for GAVNIK
                    
                    for ori = 1:length(ordered_oris)  

                        if contains(plot_type, 'raster')
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.02 0.17], psthBinSize/10); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
                            nexttile % opens a subplot tile
                            plot(rasterX,rasterY,'k','LineWidth',1) % plot the spiking raster for the current orientation
                            xline(0,'r',LineWidth=1) % put a vertical red line at time zero (stimulus onset)
                            xlim([-0.02 0.17])
                            xticks([0, 0.05, 0.1, 0.15]);
                            ylim([0 sum(Task_info.stim_orientation==ordered_oris(ori))]) % set the max y coord to be the total # trials with this orientation
                            ylabel('Trial');
                            % imagesc(binnedArray)
                            % xticks([1.5 10.5 15.5 20.5 30.5])
                            % xticklabels([-0.150 -0.050 0 0.050 0.150])
                            % xline(15.5,'r',LineWidth=1)
                            % colorbar
                            % colormap(flipud(gray))
                            
                            if ori == 2
                                title('Grey screen'); % in the GAVNIK_A_CD and GAVNIK_E_CD conditions, the second stimulus is absent
                            else
                                title(sprintf('Orientation %d%s', round(ordered_oris(ori)), char(176)));
                            end
                            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',24)
                        end
                        
                        
                        if contains(plot_type, 'FR')
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.02 0.17], psthBinSize); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
                            
                            mean_trace = mean(binnedArray, 1);
                            z_trace = (mean_trace - mean(zscore_counts)) / std(zscore_counts);
                            nexttile;
                            plot(bins, z_trace, 'k', 'LineWidth', 1.5);
                            xline(0, 'r', 'LineWidth', 1);
                            xlim([-0.02 0.17]);
                            xticks([0, 0.05, 0.1, 0.15]);
                            ylabel('Z-scored FR');
                            if ori == 2
                                title('Grey screen'); % in the GAVNIK_A_CD and GAVNIK_E_CD conditions, the second stimulus is absent
                            else
                                title(sprintf('Orientation %d%s', round(ordered_oris(ori)), char(176)));
                            end
                            
                            set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);
                        end

                        %[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.150 0.30], psthBinSize);
                        %z_binnedArray{nCluster}{ori} = (binnedArray - mean(spikecounts_cluster))./std(spikecounts_cluster); %z-score normalisation of the binned Array for this orientation: subtracts the mean and divides by the s.d. of the spike count histogram for the whole session
                        % hold on
                        % plot(bins,mean(z_binnedArray));
             
                        % time_selected 
                    
                    end
                    
                    % === Extra tile for post-stimulus activity for 4th orientation ===
                    if length(ordered_oris) >= 4 && contains(plot_type, 'FR')
                        stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(4));
                        [~, bins_long, ~, ~, ~, binnedArray_long] = psthAndBA(spike_times_this_cluster, stim_onsets, [-0.02 0.75], psthBinSize);
                        mean_trace_long = mean(binnedArray_long, 1);
                        z_trace_long = (mean_trace_long - mean(zscore_counts)) / std(zscore_counts);

                        nexttile;
                        plot(bins_long, z_trace_long, 'k', 'LineWidth', 1.5);
                        xline(0, 'r', 'LineWidth', 1);
                        xlim([-0.02 0.75]);
                        xticks(0:0.05:0.75);
                        ylabel('Z-scored FR');
                        
                        title(sprintf('Orientation %d%s (Extended to show post-stimulus oscillations)', round(ordered_oris(4)), char(176)));
                        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);
                    end
                    
                    if length(ordered_oris) >= 4 && contains(plot_type, 'raster')
                        stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(4));
    
                        [~, ~, rasterX_long, rasterY_long, ~, ~] = psthAndBA(spike_times_this_cluster, stim_onsets, [-0.02 0.75], psthBinSize/10);
    
                        nexttile;
                        plot(rasterX_long, rasterY_long, 'k', 'LineWidth', 1);
                        xline(0, 'r', 'LineWidth', 1);
                        xlim([-0.02 0.75]);
                        ylim([0 sum(Task_info.stim_orientation == ordered_oris(4))]);
                        xticks(0:0.05:0.75);
                        
                        title(sprintf('Orientation %d%s (Extended to show post-stimulus oscillations)', round(ordered_oris(4)), char(176)));
                        ylabel('Trial');
                        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);
                    end


                    % Sanitize Stimulus_type for filenames
                    safeStimulusType = regexprep(Stimulus_type, '[:\\/*?"<>| ]', '_');
                    cluster_depth = clusters(nprobe).peak_depth(cluster_id(nCluster));
                    sgtitle(sprintf('%s - %s: Response of %s Cluster %i (%.0f µm)', subject_number, Stimulus_type, depth_for_analysis, cluster_id(nCluster), cluster_depth), 'Interpreter', 'None');

                    if contains(plot_type, 'raster')
                        exportgraphics(fig(nCluster), sprintf('%s_%s_Cluster_%i_raster.pdf', safeStimulusType, depth_for_analysis, cluster_id(nCluster)));
                    end

                    if contains(plot_type, 'FR')
                        exportgraphics(fig(nCluster), sprintf('%s_%s_Cluster_%i_FR.pdf', safeStimulusType, depth_for_analysis, cluster_id(nCluster)));
                    end
                end

                %save_all_figures(save_path,[],'ContentType',ContentType)
            end
        end


    end
end

%clusters_probe = select_clusters(clusters(nprobe),params); %only look at good clusters
%[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikeTimes, eventTimes, window, psthBinSize);
