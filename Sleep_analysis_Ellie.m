%%% ERB 2025. This code performs sleep-staging and calls functions to
% detect slow waves, ripples and spindles, and saves the results. It plots V1 peri-ripplepeak spiking activity. 
% For LFP analyses this script uses one channel from putative L5 and one from putative
% CA1, both based on best regional power per PSD analysis (for which, use PSD_analysis_UnProcessedLFP_ellie).
% L4 channel pair is identified via CSD analysis (for which, run CSD_Gratings_afterFILT_ellie or CSD_GAVNIKstims_afterFILT_ellie)

clear all
addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

SUBJECTS = {'M00013'};  %% set this - 1/4
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);
params = create_cluster_selection_params('sorting_option','ellie');
%% 2/4
Stimulus_type = 'Sleep_Box_2'; % 'Sleep_Box' 'Sleep_Box_1' 'Sleep_Box_2' 'Sleep_Box_3'
%% 3/4
cd('V:\Ellie\DATA\SUBJECTS\M00013\analysis\20250204\Sleep_Box_2')
%% 4/4
for nsession = 5 %5/5 row number of recording date in "experiment_info" 
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

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

        Brain_surface_depth = depths_from_PSD.surface_depth_PSD;
        L4_channel_pair = earliest_V1sink_CSD(1).overall_median_channel_pair; % from CSD analysis
        L4_channel_pair_depth = median(earliest_V1sink_CSD(1).overall_median_pair_depth); % take the median of the depths of the two channels
        %L4_channel_pair_depth = 4820; % enter manually if CSD analysis is not avaiable.
        L5_depth = depths_from_PSD.L5_depth_PSD; % strongest spiking depth in V1 region, from PSD analysis
        best_V1_channel = find(ycoords == L5_depth, 1); % for simplicity, return the first channel at this depth (i.e. either xcoord 11 or 27)
        CA1_depth = depths_from_PSD.CA1_depth_PSD;
        best_HPC_channel = find(ycoords == CA1_depth, 1); % for simplicity, return the first channel at this depth (i.e. either xcoord 11 or 27)
        
        V1_depth_range = [L5_depth - 400 , L4_channel_pair_depth + 570]; % low to high
                
        DIR = dir(fullfile(options.EPHYS_DATAPATH,'*lf.bin')); % Locate the file containing LF data
        % Load the channel map here 
        files = dir(fullfile(options.EPHYS_DATAPATH, '*ChanMap*.mat')); %the channel map has the y coordinate of each channel. The dir function lists the contents of a folder or provides information about files and directories matching a specified pattern
        % Construct the full path to the channel map file
        channel_map_filename = fullfile(files(1).folder, files(1).name); % files(1) refers to the first file in the list of matches
        load(channel_map_filename);
        
        % Select all channels across V1 and CA1
        selected_channels = find(ycoords >= (CA1_depth - 150) & ycoords <= max(V1_depth_range));
        best_HPC_channel_idx = find(selected_channels == best_HPC_channel);
        best_V1_channel_idx = find(selected_channels == best_V1_channel);

        % Extract LFP
        [raw_LFP,tvec,SR,chan_config,~] = load_LFP_NPX(options,[],'selected_channels', selected_channels);
        rms_values = rms(raw_LFP, 2); % for detection of noisy/bad channels
        channel_std = std(raw_LFP, 0, 2);

        %%%%% Find good CA1 channels
        % First select channels within 150um of the mid-CA1 depth identified by PSD analysis
        HPC_channels = find(ycoords >= (CA1_depth - 150) & ycoords <= (CA1_depth + 150));

        % Find the indices of HPC_channels *within* selected_channels (which raw_LFP lines up with)
        HPC_relative_idx = find(ismember(selected_channels, HPC_channels));
        
        % Preallocate bad channel mask
        bad_channels = false(length(HPC_relative_idx), 1);
        
        % Loop through just the relevant LFP rows
        for i = 1:length(HPC_relative_idx)
            idx = HPC_relative_idx(i);  % row in raw_LFP
            bad_channels(i) = rms_values(idx) > mean(rms_values) + 3 * std(rms_values) || ...
                              channel_std(idx) > 5 * median(channel_std);
        end


        % Select good channels
        good_HPC_channels = HPC_channels(~bad_channels);
        
        % Find best CA1 channel with highest ripple power
        %[~, best_channel] = max(power{nprobe}(good_HPC_channels, 6)); %running extract_PSD_profile caused matlab to become unresponsive
        %best_HPC_channel = good_HPC_channels(best_channel);

        
        %%%%% Find good V1 LFP channels
        % First select channels within the putative range of V1
        V1_channels = find(ycoords >= min(V1_depth_range) & ycoords <= max(V1_depth_range));

        % Find the indices of V1_channels *within* selected_channels (which raw_LFP lines up with)
        V1_relative_idx = find(ismember(selected_channels, V1_channels));
        
        % Preallocate bad channel mask
        bad_channels = false(length(V1_relative_idx), 1);
        
        % Loop through just the relevant LFP rows
        for i = 1:length(V1_relative_idx)
            idx = V1_relative_idx(i);  % row in raw_LFP
            % ERB -I have relaxed the bad channel 3 * std(rms_values) criterion to 4* as 3* was excluding the best spiking L5 channel..
            bad_channels(i) = rms_values(idx) > mean(rms_values) + 4 * std(rms_values) || ...
                              channel_std(idx) > 5 * median(channel_std);
        end


        % Select good channels
        good_V1_channels = V1_channels(~bad_channels);
        
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
        [freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
        [tvec' raw_LFP(best_V1_channel_idx,:)'],[tvec' raw_LFP(best_HPC_channel_idx,:)'],...
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
        area(tvec, movement_trace * 4, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        area(tvec, SWS_trace * 3,      'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        area(tvec, REM_trace * 2,      'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        area(tvec, quietWake_trace,    'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.6);

        yticks([1 2 3 4])
        yticklabels({'Quiet Wake', 'REM', 'SWS', 'Movement'})
        xlabel('Time (s)');
        ylabel('State');
        set(gca,"TickDir","out",'box', 'off', 'FontSize', 14)
        sgtitle(sprintf('%s Sleep/Wake States over time Day %s %s ', subject_number, session_info.probe.SESSION, Stimulus_type), 'Interpreter', 'none', 'FontSize', 16);
        ylim([0 5]);       

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
            peak_depths = clusters(nprobe).peak_depth(sc.cluster_id); % get depths of selected clusters

            % Find cluster_ids within selected depth range
            keep_mask = peak_depths >= min(V1_depth_range) & peak_depths <= max(V1_depth_range);
            V1_cluster_ids = sc.cluster_id(keep_mask);

            % Keep only spike times and IDs for clusters at selected depth
            V1_clusters(nprobe).cluster_id = V1_cluster_ids;
            V1_clusters(nprobe).spike_times = sc.spike_times(ismember(sc.spike_id, V1_cluster_ids));
            V1_clusters(nprobe).spike_id = sc.spike_id(ismember(sc.spike_id, V1_cluster_ids));
            
        end
        
        temp = DetectSlowWaves_masa('time',tvec,'lfp',raw_LFP(best_V1_channel_idx,:),'spikes',V1_clusters(nprobe),'NREMInts',behavioural_state.SWS);
        
        timevec = tvec';
        
        [spindles(nprobe)] = FindSpindles_masa(raw_LFP(best_V1_channel_idx,:), timevec,'behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(timevec)),...
            'noise',[],'passband',[9 17],'thresholds',[1 3],'show','on');

        [ripples(nprobe)] = FindRipples_masa(raw_LFP(best_HPC_channel_idx,:), timevec,'behaviour',Behaviour,'minDuration',30,'durations',[30 200],'frequency',mean(1./diff(timevec)),...
        'noise',[],'passband',[125 300],'thresholds',[2 5],'show','on','best_channel',best_HPC_channel);
        

        slowwave_file = fullfile(options.ANALYSIS_DATAPATH, 'slowwaves.mat');
        save(slowwave_file, 'temp');
        spindle_file = fullfile(options.ANALYSIS_DATAPATH, 'spindles.mat');
        save(spindle_file, 'spindles');
        ripple_file  = fullfile(options.ANALYSIS_DATAPATH, 'ripples.mat');
        save(ripple_file, 'ripples');
        

        %%%% Look at V1 spiking around ripple peaks
        psthBinSize = 0.01; % but convert to 1ms for raster plots
        allCluster_binnedArrays = {};
        for nCluster = 1:length(V1_cluster_ids) % loop through each V1 cluster
            spike_times_this_cluster = V1_clusters.spike_times(V1_clusters.spike_id == V1_cluster_ids(nCluster));
            
            fig(nCluster)=figure; % open a figure window
            fig(nCluster).Name=sprintf('Peri-ripplepeak response V1 Cluster %i', V1_cluster_ids(nCluster)); %overall figure title includes the cluster_id (one count higher than zero-based SI output)
            fig(nCluster).Position = [114 90 770 650]; % sets the size of the figure window
            
            [~, ~, rasterX, rasterY, ~, ~] = psthAndBA(spike_times_this_cluster,  ripples.peaktimes, [-1 1], psthBinSize/10); % get the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window around ripple peak
            
            % Plot raster
            figure(fig(nCluster));
            sgtitle(sprintf('%s V1 cluster %i Peri-ripplepeak spiking (%ix CA1 Ripples) %s Day %s', subject_number, V1_cluster_ids(nCluster), length(ripples.peaktimes), Stimulus_type, session_info.probe.SESSION), 'Interpreter', 'none', 'FontSize', 16);
            subplot(2,1,1);
            plot(rasterX,rasterY,'k','LineWidth',2) % plot the spiking raster 
            xline(0,'r',LineWidth=1) % put a vertical red line at time zero (peak ripple power)
            xlim([-1 1])
            xlabel('Time from ripple peak (s)');
            ylabel('Ripple number');
            %set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',14)
                        
            % Plot FR
            [psth, bins, ~, ~, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster, ripples.peaktimes, [-1 1], psthBinSize); % get the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window around ripple peak
            allCluster_binnedArrays{end+1} = binnedArray; % each is [numRipples x numBins]
            mean_trace = mean(binnedArray, 1) / psthBinSize; % average over ripple events
            subplot(2,1,2);
            hold on;
            plot(bins, mean_trace, 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', sprintf('Ripples (%ix)', length(ripples.peaktimes)));
            xlim([-1 1])
            xline(0,'r',LineWidth=1) % put a vertical red line at time zero (peak ripple power)
            xlabel('Time from ripple peak (s)');
            ylabel('Mean firing rate across ripples (Hz)', 'FontSize', 14);
            
            fig_filename = sprintf('%s Cluster %i Peri-Ripplepeak spiking (%ix CA1 Ripples) %s Day %s.fig', subject_number, V1_cluster_ids(nCluster), length(ripples.peaktimes), Stimulus_type, session_info.probe.SESSION);
            fig_save_path = fullfile(pwd, fig_filename);
            savefig(fig(nCluster), fig_save_path);

        end

        % ----- Summary plot for all V1 clusters -----
        % Combine all V1 cluster spike times
        all_spike_times = [];
        for nCluster = 1:length(V1_cluster_ids)
            this_spike_times = V1_clusters.spike_times(V1_clusters.spike_id == V1_cluster_ids(nCluster));
            all_spike_times = [all_spike_times; this_spike_times];
        end
        
        % Call psthAndBA on all spike times
        [psth, bins, ~, ~, ~, binnedArray] = psthAndBA(all_spike_times, ripples.peaktimes, [-1 1], psthBinSize);
        mean_FR = mean(binnedArray, 1) / psthBinSize;
        
        % Plot
        summaryFig = figure;
        summaryFig.Name = 'Aggregate FR from all V1 clusters';
        plot(bins, mean_FR, 'k', 'LineWidth', 2);
        xline(0, 'r', 'LineWidth', 1);
        xlabel('Time from ripple peak (s)');
        ylabel('Aggregate firing rate (Hz)');
        title(sprintf('%s: Aggregate V1 peri-ripplepeak FR (%i clusters) %s Day %s', subject_number, length(V1_cluster_ids), Stimulus_type, session_info.probe.SESSION), 'Interpreter', 'none');
        
        % Save
        savefig(summaryFig, fullfile(pwd, sprintf('%s Aggregate V1 peri-ripplepeak FR %s Day %s.fig', subject_number, Stimulus_type, session_info.probe.SESSION)));
        
    end    
end    
       