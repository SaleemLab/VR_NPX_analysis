%%% ERB 2025. This code performs sleep-staging and calls functions to
% detect slow waves, ripples and spindles, and saves the results. It plots V1 peri-ripplepeak spiking activity. 
% For LFP analyses this script uses one channel from putative L5 and one from putative
% CA1, both based on best regional power per PSD analysis (for which, use PSD_analysis_UnProcessedLFP_ellie).
% L4 depth is identified via CSD analysis (for which, run CSD_Gratings_afterFILT_ellie or CSD_GAVNIKstims_afterFILT_ellie)

clear all
addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

SUBJECTS = {'M00087'};  %% set this - 1/4
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);
params = create_cluster_selection_params('sorting_option','ellie');
psthBinSize = 0.01; % this script divides this by 10 (to 0.001s) for raster plots

%% 2/4
Stimulus_type = 'Sleep_Box_2'; % 'Sleep_Box' 'Sleep_Box_1' 'Sleep_Box_2' 'Sleep_Box_3'
cd('V:\Ellie\DATA\SUBJECTS\M00087\analysis\20260211\Sleep_Box_2')
V1_sleepstaging_shank = 1; %*************************visually inspect the PSD plots and choose whichever shank best captures V1
HPC_sleepstaging_shank = 4;   %*****************************visually inspect the PSD plots and choose whichever shank best captures CA1

%% 3/4***** For NPX2.0 you will use a different L4 channel for each shank. Use CSD to estimate the best channel to use in L4
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

%% 4/4
for nsession = 6 %5/5 row number of recording date in "experiment_info" 
    session_info = experiment_info(nsession).session(strcmp(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(strcmp(experiment_info(nsession).StimulusName,Stimulus_type));

    for n = 1:length(session_info) % How many recording sessions 
        options = session_info(n).probe(1);
        subject_number = session_info(n).probe(1).SUBJECT;
        
        %DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters = clusters_ks4;
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
        best_V1_channel = cell(numel(unique_shanks), 1);
        best_HPC_channel = cell(numel(unique_shanks), 1);

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
            V1_depth_range{iShank}   = [L5_depth - 330, L5_depth + 700]; % Senzai 2019 - distance from mid L5 to lower L6 appears to be ~260um. Senzai 2019 - mid L5 appears to fall ~610um below the brain surface                         
            CA1_depth_range{iShank}  = [CA1_depth - 150, CA1_depth + 150];
            Sub_CA1_depth_range{iShank} = [min(CA1_depth_range{iShank}) - 1000, min(CA1_depth_range{iShank})];  
            
            best_V1_channel{iShank} = find((kcoords == this_shank) & (ycoords == L5_depth), 1); 
            best_HPC_channel{iShank} = find((kcoords == this_shank) & (ycoords == CA1_depth), 1);                     
        end
              
        
        DIR = dir(fullfile(options.EPHYS_DATAPATH,'*lf.bin')); % Locate the file containing LF data
        % Load the channel map here 
        files = dir(fullfile(options.EPHYS_DATAPATH, '*ChanMap*.mat')); %the channel map has the y coordinate of each channel. The dir function lists the contents of a folder or provides information about files and directories matching a specified pattern
        % Construct the full path to the channel map file
        channel_map_filename = fullfile(files(1).folder, files(1).name); % files(1) refers to the first file in the list of matches
        load(channel_map_filename);
        
        % Select all channels across V1 and CA1, by shank
        selected_channels = [];
        for iShank = 1:numel(unique_shanks)
            this_shank = unique_shanks(iShank);
            idx = (kcoords == this_shank) & (ycoords >= CA1_depth_range{iShank}(1)) & (ycoords <= V1_depth_range{iShank}(2));
            selected_channels = [selected_channels; find(idx)];
        end
        
        selected_channels = unique(selected_channels);

        best_HPC_channel_idx = cell(numel(unique_shanks), 1);
        best_V1_channel_idx  = cell(numel(unique_shanks), 1);
        
        for iShank = 1:numel(unique_shanks)
            best_HPC_channel_idx{iShank} = ...
                find(selected_channels == best_HPC_channel{iShank}, 1);
        
            best_V1_channel_idx{iShank} = ...
                find(selected_channels == best_V1_channel{iShank}, 1);
        end

        % Extract LFP
        [raw_LFP,tvec,SR,chan_config,~] = load_LFP_NPX(options,[],'selected_channels', selected_channels);
        rms_values = rms(raw_LFP, 2); % for detection of noisy/bad channels
        channel_std = std(raw_LFP, 0, 2);

        %%%%% Find good CA1 channels
        good_HPC_channels = [];

        for iShank = 1:numel(unique_shanks)
            this_shank = unique_shanks(iShank);
                    
            % Absolute channel indices % First select channels within 150um of the mid-CA1 depth identified by PSD analysis
            HPC_channels = find((kcoords == this_shank) & (ycoords >= CA1_depth_range{iShank}(1)) & (ycoords <= CA1_depth_range{iShank}(2)));
        
            % Find the indices of HPC_channels *within* selected_channels (which raw_LFP lines up with)
            HPC_relative_idx = find(ismember(selected_channels, HPC_channels));        
    
            % Preallocate bad channel mask
            bad_channels = false(length(HPC_relative_idx), 1);
            
            % Loop through just the relevant LFP rows
            for i = 1:numel(HPC_relative_idx)
                idx = HPC_relative_idx(i);  % row in raw_LFP
                bad_channels(i) = rms_values(idx) > mean(rms_values) + 3 * std(rms_values) || ...
                                  channel_std(idx) > 5 * median(channel_std);
            end

             % Select good channels
            good_HPC_channels = [good_HPC_channels; HPC_channels(~bad_channels)];
        end
        
        % Find best CA1 channel with highest ripple power
        %[~, best_channel] = max(power{nprobe}(good_HPC_channels, 6)); %running extract_PSD_profile caused matlab to become unresponsive
        %best_HPC_channel = good_HPC_channels(best_channel);

        
        %%%%% Find good V1 LFP channels
        good_V1_channels = [];
        for iShank = 1:numel(unique_shanks)
            this_shank = unique_shanks(iShank);
            V1_channels = find((kcoords == this_shank) & (ycoords >= V1_depth_range{iShank}(1)) & (ycoords <= V1_depth_range{iShank}(2)) );
           % Find the indices of V1_channels *within* selected_channels (which raw_LFP lines up with)
            V1_relative_idx = find(ismember(selected_channels, V1_channels));
        
            % Preallocate bad channel mask
            bad_channels = false(length(V1_relative_idx), 1);
        
            % Loop through just the relevant LFP rows
            for i = 1:numel(V1_relative_idx)
                idx = V1_relative_idx(i);  % row in raw_LFP
                % ERB -I have relaxed the bad channel 3 * std(rms_values) criterion to 4* as 3* was excluding the best spiking L5 channel..
                bad_channels(i) = rms_values(idx) > mean(rms_values) + 4 * std(rms_values) || ...
                                  channel_std(idx) > 5 * median(channel_std);
            end
            % Select good channels
            good_V1_channels = [good_V1_channels; V1_channels(~bad_channels)];
        end
        
        % Find best V1 channel with highest high freq power per PSD
        %[~, best_channel] = max(power{nprobe}(good_V1_channels, 7)); %running extract_PSD_profile caused matlab to become unresponsive
        %best_V1_channel = good_V1_channels(best_channel);    
        

        %%%%%  Sleep detection
        speedTreshold = 1;

        %detect moments when the smoothed frame-by-frame pixel change (mobility) exceeds a threshold, i.e. when the change between consecutive time points is large (above 2000).
        mobility_thresholded = abs([0 diff(movmean(Behaviour.mobility,1/mean(diff(Behaviour.tvec))))])>2000;%1 second movemean windows

        %smooth out brief pauses; treat them as part of the movement        
        mob_index = find(mobility_thresholded==1);
        for nindex = 1:length(mob_index)
            if nindex<length(mob_index)
                if Behaviour.tvec(mob_index(nindex+1))-Behaviour.tvec(mob_index(nindex))<1 % if immobility less than 1 second, treats that as movement period
                    mobility_thresholded(mob_index(nindex):mob_index(nindex+1))=1;
                end
            end
        end

        %             figure; plot(Behaviour.tvec,Behaviour.mobility);hold on;
        %             figure; plot(Behaviour.tvec,abs([0 diff(movmean(Behaviour.mobility,1/mean(diff(Behaviour.tvec))))]))
        Behaviour.mobility_thresholded = mobility_thresholded;

        % mobility_thresholded = interp1(Behaviour.tvec,double(mobility_thresholded),LFP(probe_no).tvec,'linear');
        
        Behaviour.mobility_zscore = abs([0 diff(movmean(Behaviour.mobility,1/mean(diff(Behaviour.tvec))))]);% Diff of pixel change
        Behaviour.mobility_zscore(isnan(Behaviour.mobility_zscore))=mean(Behaviour.mobility_zscore,'omitnan');
        Behaviour.mobility_zscore=zscore(Behaviour.mobility_zscore);
        
        %             figure; plot(Behaviour.mobility_zscore)

        speed = double(mobility_thresholded); % Convert from logical to numeric
        %speed = double(session_clusters.mobility_thresholded{1});
        %Behaviour.mobility_zscore=session_clusters.mobility_zscore{1};
        speed = interp1(Behaviour.tvec,speed,tvec,'linear'); %  interpolate the speed signal (movement binary vector) to the tvec timeline

        %%%% Classify the sleep state of the mouse
        iV1  = find(unique_shanks == V1_sleepstaging_shank, 1);
        iHPC = find(unique_shanks == HPC_sleepstaging_shank, 1);
        
        [freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
        [tvec' raw_LFP(best_V1_channel_idx{iV1},:)'],[tvec' raw_LFP(best_HPC_channel_idx{iHPC},:)'],...
        [tvec' speed'],speedTreshold);

        % Plot the behavioural state data to sense check it
        % Initialize state vectors
        movement_trace    = zeros(size(tvec));
        SWS_trace         = zeros(size(tvec));
        REM_trace         = zeros(size(tvec));
        quietWake_trace   = zeros(size(tvec));
        
        % Fill state traces manually
        for i = 1:size(movement,1)
            movement_trace(tvec >= movement(i,1) & tvec <= movement(i,2)) = 1;
        end
        for i = 1:size(SWS,1)
            SWS_trace(tvec >= SWS(i,1) & tvec <= SWS(i,2)) = 1;
        end
        for i = 1:size(REM,1)
            REM_trace(tvec >= REM(i,1) & tvec <= REM(i,2)) = 1;
        end
        for i = 1:size(quietWake,1)
            quietWake_trace(tvec >= quietWake(i,1) & tvec <= quietWake(i,2)) = 1;
        end
        
        % Plot
        figure; hold on;
        
        % Movement
        for i = 1:size(movement,1)
            patch([movement(i,1) movement(i,2) movement(i,2) movement(i,1)], ...
                  [3.5 3.5 4.5 4.5], 'k', 'EdgeColor','none','FaceAlpha',0.6);
        end
        
        % Quiet Wake
        for i = 1:size(quietWake,1)
            patch([quietWake(i,1) quietWake(i,2) quietWake(i,2) quietWake(i,1)], ...
                  [2.5 2.5 3.5 3.5], 'g', 'EdgeColor','none','FaceAlpha',0.6);
        end


        % SWS
        for i = 1:size(SWS,1)
            patch([SWS(i,1) SWS(i,2) SWS(i,2) SWS(i,1)], ...
                  [1.5 1.5 2.5 2.5], 'b', 'EdgeColor','none','FaceAlpha',0.6);
        end
        
        % REM
        for i = 1:size(REM,1)
            patch([REM(i,1) REM(i,2) REM(i,2) REM(i,1)], ...
                  [0.5 0.5 1.5 1.5], 'r', 'EdgeColor','none','FaceAlpha',0.6);
        end
        
        
        yticks([1 2 3 4])
        ylim([0 5]);    
        % Remove default labels
        yticklabels({})
        
        % Add custom colored labels to the left of the axis
        text(min(xlim)-0.01*range(xlim), 3, 'Quiet Wake', 'Color','g', ...
            'FontSize',14, 'HorizontalAlignment','right', 'VerticalAlignment','middle')
        text(min(xlim)-0.01*range(xlim), 1, 'REM', 'Color','r', ...
            'FontSize',14, 'HorizontalAlignment','right', 'VerticalAlignment','middle')
        text(min(xlim)-0.01*range(xlim), 2, 'SWS', 'Color','b', ...
            'FontSize',14, 'HorizontalAlignment','right', 'VerticalAlignment','middle')
        text(min(xlim)-0.01*range(xlim), 4, 'Movement', 'Color','k', ...
            'FontSize',14, 'HorizontalAlignment','right', 'VerticalAlignment','middle')
        xlabel('Time (s)');
        set(gca,"TickDir","out",'box', 'off', 'FontSize', 14)
        sgtitle(sprintf('%s Sleep/Wake States over time Day %s %s ', subject_number, session_info.probe.SESSION, Stimulus_type), 'Interpreter', 'none', 'FontSize', 16);
        fig_filename = sprintf('%s SleepWake States over time Day %s %s ', subject_number, session_info.probe.SESSION, Stimulus_type);
        fig_save_path = fullfile(pwd, fig_filename);
        savefig(fig_save_path);   

        %%%% Save the behavioural states to a struct
        behavioural_state = struct();
        behavioural_state.quietWake = quietWake;
        behavioural_state.SWS = SWS;
        behavioural_state.REM = REM;
        behavioural_state.movement = movement;
        save(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state.mat'),'behavioural_state')


        %%%% Detect Slow Waves and Spindles and Ripples
        % Extract V1 MUA - use strict filtering if this yields >50 cells. Otherwise drop the filtering but cap the firing rates at 99%
        V1_clusters = struct();  % Predefine        
        for nprobe = 1:length(clusters)
            % filter the clusters
            selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters, which pass the set parameters
            % get the depths of the clusters which pass the filtering
            sc = selected_clusters(nprobe);
            % Get depth and shank from channel map
            cluster_channels = clusters(nprobe).peak_channel(sc.cluster_id);
            cluster_depths = ycoords(cluster_channels);
            cluster_shanks = kcoords(cluster_channels);        
            % Preallocate keep mask
            keep_mask = false(size(sc.cluster_id));
            
            % Apply shank-specific depth selection
            for iShank = 1:numel(unique_shanks)
                this_shank = unique_shanks(iShank);               
        
                keep_mask = keep_mask | ...
                    (cluster_shanks == this_shank & ...
                     cluster_depths >= V1_depth_range{iShank}(1) & ...
                     cluster_depths <= V1_depth_range{iShank}(2));
            end
            
            V1_cluster_ids = sc.cluster_id(keep_mask);

            % Keep only spike times and IDs for clusters at selected depth
            V1_clusters(nprobe).cluster_id = V1_cluster_ids;
            V1_clusters(nprobe).spike_times = sc.spike_times(ismember(sc.spike_id, V1_cluster_ids));
            V1_clusters(nprobe).spike_id = sc.spike_id(ismember(sc.spike_id, V1_cluster_ids));            
        end
        % Detecting slow waves using one shank here- could do separately for each shank?
        % detects neocortical slow waves using a combination of a positive deflection in the LFP (delta wave) and a dip in gamma power.
        temp = DetectSlowWaves_masa('time',tvec,'lfp',raw_LFP(best_V1_channel_idx{iV1},:),'spikes',V1_clusters(nprobe),'NREMInts',behavioural_state.SWS);
        
        timevec = tvec';
        % Finding spindles using one shank here- could do separately for each shank?
        [spindles(nprobe)] = FindSpindles_masa(raw_LFP(best_V1_channel_idx{iV1},:), timevec,'behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(timevec)),...
            'noise',[],'passband',[9 17],'thresholds',[1 3],'show','on');
        % Finding ripples using one shank here - could do separately for each shank?
        [ripples(nprobe)] = FindRipples_masa(raw_LFP(best_HPC_channel_idx{iHPC},:), timevec,'behaviour',Behaviour,'minDuration',30,'durations',[30 200],'frequency',mean(1./diff(timevec)),...
        'noise',[],'passband',[125 300],'thresholds',[2 5],'show','on'); % ,'best_channel',best_HPC_channel
        

        % Extract Ripple peak times to assign to a behavioural state
        ripple_times = ripples(nprobe).peaktimes;
        
        % Preallocate state labels
        ripple_state = strings(size(ripple_times));
        
        %% --- Classify SWS ripples ---
        for iInt = 1:size(SWS,1)
            in_state = ripple_times >= SWS(iInt,1) & ...
                       ripple_times <= SWS(iInt,2);
        
            ripple_state(in_state) = "SWS";
        end
        
        %% --- Classify REM ripples ---
        for iInt = 1:size(REM,1)
            in_state = ripple_times >= REM(iInt,1) & ...
                       ripple_times <= REM(iInt,2);
        
            ripple_state(in_state) = "REM";
        end
        
        %% --- Classify quiet wake ripples ---
        for iInt = 1:size(quietWake,1)
            in_state = ripple_times >= quietWake(iInt,1) & ...
                       ripple_times <= quietWake(iInt,2);
        
            ripple_state(in_state) = "quietWake";
        end
        
        %% --- Classify movement ripples ---
        for iInt = 1:size(movement,1)
            in_state = ripple_times >= movement(iInt,1) & ...
                       ripple_times <= movement(iInt,2);
        
            ripple_state(in_state) = "movement";
        end
        
        ripples(nprobe).state = ripple_state;

        disp(['SWS ripples: ', num2str(sum(ripple_state=="SWS"))]);
        disp(['REM ripples: ', num2str(sum(ripple_state=="REM"))]);
        disp(['quietWake ripples: ', num2str(sum(ripple_state=="quietWake"))]);
        disp(['movement ripples: ', num2str(sum(ripple_state=="movement"))]);
        disp(['Unclassified ripples: ', num2str(sum(ripple_state==""))]);

        slowwave_file = fullfile(options.ANALYSIS_DATAPATH, 'slowwaves.mat');
        save(slowwave_file, 'temp');
        spindle_file = fullfile(options.ANALYSIS_DATAPATH, 'spindles.mat');
        save(spindle_file, 'spindles');
        ripple_file  = fullfile(options.ANALYSIS_DATAPATH, 'ripples.mat');
        save(ripple_file, 'ripples');       

        %%%% Look at V1 spiking around SWS ripple peaks

        % ----- Summary plot for all V1 clusters -----
        % Combine all V1 cluster spike times
        all_spike_times = [];
        for nCluster = 1:length(V1_cluster_ids)
            this_spike_times = V1_clusters.spike_times(V1_clusters.spike_id == V1_cluster_ids(nCluster));
            all_spike_times = [all_spike_times; this_spike_times];
        end
        
        % Call psthAndBA on all spike times
        [psth, bins, ~, ~, ~, binnedArray] = psthAndBA(all_spike_times, ripples(nprobe).peaktimes(ripples(nprobe).state == "SWS"), [-1 1], psthBinSize);
        mean_FR = mean(binnedArray, 1) / psthBinSize;
        
        % Plot
        summaryFig = figure;
        summaryFig.Name = 'Aggregate FR during SWS from all V1 clusters';

        % plot(bins, mean_FR, 'k', 'LineWidth', 2);
        % xline(0, 'r', 'LineWidth', 1);
        % xlabel('Time from ripple peak (s)');
        % ylabel('Aggregate firing rate (Hz)');
        % title(sprintf('%s: Aggregate V1 peri-ripplepeak FR (%i clusters) %s Day %s', subject_number, length(V1_cluster_ids), Stimulus_type, session_info.probe.SESSION), 'Interpreter', 'none');
        % 
        %%% CHANGED: explicit axes handle
        axSummary = axes('Parent', summaryFig);                        %%% CHANGED
        plot(axSummary, bins, mean_FR, 'k', 'LineWidth', 2);           %%% CHANGED
        xline(axSummary, 0, 'r', 'LineWidth', 1);                      %%% CHANGED
        xlabel(axSummary, 'Time from ripple peak (s)');                %%% CHANGED
        ylabel(axSummary, 'Aggregate firing rate (Hz)');               %%% CHANGED
        title(axSummary, sprintf( ...
            '%s: Aggregate V1 SWS peri-ripplepeak FR (%i clusters) %s Day %s', subject_number, length(V1_cluster_ids), Stimulus_type, session_info.probe.SESSION), ...
            'Interpreter','none');                                     %%% CHANGED

        % Save
        savefig(summaryFig, fullfile(pwd, sprintf('%s Aggregate V1 SWS peri-ripplepeak FR %s Day %s.fig', subject_number, Stimulus_type, session_info.probe.SESSION)));
        
        
        allCluster_binnedArrays = {};
        for nCluster = 1:length(V1_cluster_ids) % loop through each V1 cluster
            spike_times_this_cluster = V1_clusters.spike_times(V1_clusters.spike_id == V1_cluster_ids(nCluster));
            
            fig(nCluster)=figure; % open a figure window
            fig(nCluster).Name=sprintf('SWS peri-ripplepeak response V1 Cluster %i', V1_cluster_ids(nCluster)); %overall figure title includes the cluster_id (one count higher than zero-based SI output)
            fig(nCluster).Position = [114 90 770 650]; % sets the size of the figure window
            
            %clear psthAndBA
            axRaster = subplot(2,1,1,'Parent',fig(nCluster));         %%% CHANGED
            axFR     = subplot(2,1,2,'Parent',fig(nCluster));        %%% CHANGED
            
            [~, ~, rasterX, rasterY, ~, ~] = psthAndBA(spike_times_this_cluster,  ripples(nprobe).peaktimes(ripples(nprobe).state == "SWS"), [-1 1], psthBinSize/10); % get the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window around ripple peak
            
            plot(axRaster, rasterX, rasterY, 'k', 'LineWidth', 2);    %%% CHANGED
            xline(axRaster, 0, 'r', LineWidth=1);                     %%% CHANGED
            xlim(axRaster, [-1 1]);                                   %%% CHANGED
            xlabel(axRaster, 'Time from ripple peak (s)');            %%% CHANGED
            ylabel(axRaster, 'Ripple number');                         %%% CHANGED
    
            % Plot raster
            % figure(fig(nCluster));
            % sgtitle(sprintf('%s V1 cluster %i Peri-ripplepeak spiking (%ix CA1 Ripples) %s Day %s', subject_number, V1_cluster_ids(nCluster), length(ripples.peaktimes), Stimulus_type, session_info.probe.SESSION), 'Interpreter', 'none', 'FontSize', 16);
            % subplot(2,1,1);
            % plot(rasterX,rasterY,'k','LineWidth',2) % plot the spiking raster 
            % %line(rasterX,rasterY,'k','LineWidth',2) % plot the spiking raster
            % xline(0,'r',LineWidth=1) % put a vertical red line at time zero (peak ripple power)
            % xlim([-1 1])
            % xlabel('Time from ripple peak (s)');
            % ylabel('Ripple number');
            %set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',14)
            
            %%% CHANGED: sgtitle explicitly tied to figure
            sgtitle(fig(nCluster), sprintf( ...
                '%s V1 cluster %i SWS peri-ripplepeak spiking (%ix CA1 Ripples) %s Day %s', subject_number, V1_cluster_ids(nCluster), ...
                length(ripples(nprobe).peaktimes(ripples(nprobe).state == "SWS")), Stimulus_type, session_info.probe.SESSION), 'Interpreter','none','FontSize',16);                   %%% CHANGED

            % Plot FR
            [psth, bins, ~, ~, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster, ripples(nprobe).peaktimes(ripples(nprobe).state == "SWS"), [-1 1], psthBinSize); % get the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window around ripple peak
            allCluster_binnedArrays{end+1} = binnedArray; % each is [numRipples x numBins]
            mean_trace = mean(binnedArray, 1) / psthBinSize; % average over ripple events
            
            % subplot(2,1,2);
            % hold on;
            % plot(bins, mean_trace, 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', sprintf('Ripples (%ix)', length(ripples.peaktimes)));
            % xlim([-1 1])
            % xline(0,'r',LineWidth=1) % put a vertical red line at time zero (peak ripple power)
            % xlabel('Time from ripple peak (s)');
            % ylabel('Mean firing rate across ripples (Hz)', 'FontSize', 14);
            
            %%% CHANGED: target FR plot to explicit axes
            plot(axFR, bins, mean_trace, 'b', 'LineWidth', 1.5, ...   %%% CHANGED
                'DisplayName', sprintf('Ripples (%ix)', length(ripples(nprobe).peaktimes(ripples(nprobe).state == "SWS"))));
            hold(axFR, 'on');                                         %%% CHANGED
            xline(axFR, 0, 'r', LineWidth=1);                          %%% CHANGED
            xlim(axFR, [-1 1]);                                        %%% CHANGED
            xlabel(axFR, 'Time from ripple peak (s)');                 %%% CHANGED
            ylabel(axFR, 'Mean firing rate across SWS ripples (Hz)');      %%% CHANGED


            fig_filename = sprintf('%s Cluster %i SWS Peri-Ripplepeak spiking (%ix CA1 Ripples) %s Day %s.fig', subject_number, V1_cluster_ids(nCluster), length(ripples(nprobe).peaktimes(ripples(nprobe).state == "SWS")), Stimulus_type, session_info.probe.SESSION);
            fig_save_path = fullfile(pwd, fig_filename);
            savefig(fig(nCluster), fig_save_path);

        end

        
        
        %% Extract and save sleep data for projection into pls space
        SleepData = struct();
        SleepData.subject = subject_number;
        SleepData.session = session_info.probe.SESSION;
        SleepData.StimulusType = Stimulus_type;
        SleepData.behavioural_state = behavioural_state;

        SleepData.V1 = struct();  % initialize per probe

        for nprobe = 1:length(V1_clusters)

            % Initialize
            nNeurons = length(V1_clusters(nprobe).cluster_id);
            nRipples = length(ripples(nprobe).peaktimes(ripples(nprobe).state == "SWS"));
            nTimeBins = size(allCluster_binnedArrays{1},2); % bins from psthAndBA
            
            X_sleep = nan(nRipples, nTimeBins, nNeurons);
            
            for nCluster = 1:nNeurons
                binnedArray = allCluster_binnedArrays{nCluster};  % [nRipples x nTimeBins]
                all_counts = binnedArray(:);  % vector of all counts
                meanCounts = mean(all_counts);
                stdCounts  = std(all_counts);
                if stdCounts == 0, stdCounts = 1; % avoid divide by zero if neuron silent
                end
                
                % Apply z-scoring to each ripple separately
                z_binned = (binnedArray - meanCounts) ./ stdCounts;
                X_sleep(:,:,nCluster) = z_binned;
            end
        
            SleepData.V1(nprobe).cluster_id = V1_clusters(nprobe).cluster_id;
            SleepData.V1(nprobe).X_sleep = X_sleep; % [ripples x timebins x neurons]
            SleepData.V1(nprobe).X_sleep_reshaped = reshape(X_sleep, nRipples, nTimeBins*nNeurons);
            SleepData.V1(nprobe).X_sleep_mean = squeeze(mean(X_sleep, 2)); % [nRipples x nNeurons]  
            
            SleepData.V1(nprobe).binnedArrays = allCluster_binnedArrays; % {cluster}{event x timebin}
            SleepData.V1(nprobe).bins = bins;
            SleepData.V1(nprobe).psthBinSize = psthBinSize;
            SleepData.V1(nprobe).ripple_peaktimes = ripples(nprobe).peaktimes(ripples(nprobe).state == "SWS");
        end    
        
        SleepData.ripples = ripples;
        SleepData.spindles = spindles;
        SleepData.slowwaves = temp;
    
        save(fullfile(options.ANALYSIS_DATAPATH,'SleepData.mat'),'SleepData','-v7.3');

    end    
end    
       
