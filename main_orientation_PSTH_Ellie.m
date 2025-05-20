%%% For analysis of unit spiking in response to visual stimuli

addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

%% setting metrics to screen good clusters
clear all
% Choose your probe depth of interest
L4_depth_range = 4510:4650; % 1/5. um. Set for each SESSION based on CSD +/- 70um
V1_depth_range = (min(L4_depth_range) - 400) : (max(L4_depth_range) + 500); 
CA1_depth_range = 3640:3940; % 2/5. um. Set for each SESSION based on PSD; ~300um around Ripple power "bump"
Sub_CA1_depth_range = 1550:(min(CA1_depth_range));
depth_for_analysis = 'L4'; % choose 'L4' or 'V1' or 'CA1' or 'Sub_CA1'

SUBJECTS = {'M00013'};

params = create_cluster_selection_params('sorting_option','ellie');
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);

%%% 3/5
Stimulus_type = 'GAVNIK_ABCD'; 
plot_choice = 'aggregate'; % curated 'single_units' or in 'aggregate' or uncurated 'MUA'; MUA includes all clusters from kilosort, unfiltered
plot_type = 'FR'; % 'FR' firing rate or 'raster'
z_score_period = 'entire_session'; % z score either over 'entire_session' or 'first30secs' or 'none' (for every stimulus recording
% session from 20250205 onward, I presented grey screen to the mouse for at least 30s before starting the stimulus. 'none' may be useful 
% to try for the aggregate TRAIN case across days)
%nprobe = 1;
%base_folder='V:\Ellie\DATA\SUBJECTS';
cd('V:\Ellie\DATA\SUBJECTS\M00013\analysis\20250221\GAVNIK_ABCD') % 4/5 files will be saved here in the cd


for nsession = 15 %5/5 row number of recording date in "experiment_info" 
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

        all_orientations = unique(Task_info.stim_orientation);
        %Task_info.stim_onset
        
       
        params = create_cluster_selection_params('sorting_option','ellie');
        %params.orientation_tuned = ...
        psthBinSize = 0.01; % but use 1ms for raster plots
        
        switch depth_for_analysis
            case 'L4' 
                depth_range = L4_depth_range;
            case 'V1'
                depth_range = V1_depth_range;
            case 'CA1' 
                depth_range = CA1_depth_range;
            case 'Sub_CA1'
                depth_range = Sub_CA1_depth_range;
        end

        if contains(Stimulus_type, 'OP_Tuning') 
            % For M00013, the OP_Tuning stimulus was 150ms, with 150ms grey screen between stimuli, and took 24 different directions, moving at 2Hz
            % For M00014 (and beyond), the OP_Tuning stimulus was 67ms (but Bonsai typically keeps it on 84ms, 1 frame longer), with 17ms grey screen 
            % between stimuli (but Bonsai typically keeps it grey 34ms, 1 frame longer), and took 12 different static orientations
            
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
                    if contains(z_score_period, 'entire_session')
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
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==all_orientations(ori)), [-0.15 0.150], psthBinSize/10); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window around stimulus onset
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
                        % imagesc(binnedArray)
                        % xticks([1.5 10.5 15.5 20.5 30.5])
                        % xticklabels([-0.150 -0.050 0 0.050 0.150])
                        % xline(15.5,'r',LineWidth=1)
                        % colorbar
                        % colormap(flipud(gray))
                        title(sprintf('Orientation %d%s',round(rad2deg(all_orientations(ori))), char(176)))
                        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
                       
                        if isequal(session_info(n).probe(1).SUBJECT, 'M00013')
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==all_orientations(ori)), [-0.150 0.150], psthBinSize);
                        else
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==all_orientations(ori)), [-0.05 0.120], psthBinSize);
                        end

                        z_binnedArray{nCluster}{ori} = (binnedArray - mean(zscore_counts))./std(zscore_counts);
                        
                        % hold on
                        % plot(bins,mean(z_binnedArray));
             
                        % time_selected 
                    
                    end
                    nexttile % add a new subplot tile for a summary plot showing how the z-scored firing for the cluster varies by orientation 
                    clear PLOT
                    for ori = 1:length(all_orientations)
                        PLOT(ori) = plot(bins, mean(z_binnedArray{nCluster}{ori})); %plot the z-scored PSTH (mean across trials) for each orientation
                        hold on
                    end
                    legend(PLOT(1:end),{num2str(round(rad2deg(all_orientations)))},'box','off')
                    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


                    nexttile % generate a summary tuning curve for the cluster, showing z-scored response to each orientation, with error bars 
                    mean_response=[]; %initialise array to hold the mean response across trials for each orientation
                    std_response = []; %initialise array to hold the sd of the cluster's response for each orienation
                    for ori = 1:length(all_orientations)
                        ori_response{nCluster}{ori} = max(z_binnedArray{nCluster}{ori}(:,bins>0.02&bins<0.1)'); % extract the max z-score of firing 20-100ms after stim onset, for each trial
                        mean_response(ori) = mean(ori_response{nCluster}{ori}); % calc the mean of max z-scored firing after stim-onset, across trials
                        se_response(ori) = std(ori_response{nCluster}{ori})./sqrt(length(ori_response{nCluster}{ori})); %calc the SE of the max z-scored firing after stim-onset
                    end


                    clear PLOT
                    % plot(round(rad2deg(all_orientations)),mean_response)
                    % hold on;
                    errorbar(round(rad2deg(all_orientations)),mean_response,se_response,se_response) %plot the mean of max z-scored firing (y-axis) in response to each orientation (x-axis)
                    % with error bars showing SE above and below mean z-scored max response for each orientation
                    if isequal(session_info(n).probe(1).SUBJECT, 'M00013')
                        xlim([-20 370]) % x-axis runs from -20 across 360 degrees
                    else
                        xlim([-20 190]) % x-axis runs from -20 across 180 degrees
                    end
                        % legend(PLOT(1:end),{num2str(round(rad2deg(all_orientations)))},'box','off')
                    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
                    cluster_depth = clusters(nprobe).peak_depth(cluster_id(nCluster));
                    sgtitle(sprintf('%s - %s: Response of %s Cluster %i (%.0f µm)', subject_number, Stimulus_type, depth_for_analysis, cluster_id(nCluster), cluster_depth), 'Interpreter', 'none');
                    exportgraphics(fig(nCluster), sprintf('%s_Cluster_%i.pdf', depth_for_analysis, cluster_id(nCluster))); %save as a pdf in cd 
                end

                %save_all_figures(save_path,[],'ContentType',ContentType)

            end
        end




        if contains(Stimulus_type, 'TRAIN') && contains(plot_choice, 'single_units') 
            for nprobe = 1:length(clusters)
                selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters, which pass the set parameters
                
                depth_selected_clusters = selected_clusters; % initialize with same structure
                for np = 1:length(clusters)
                    sc = selected_clusters(np);
                    peak_depths = clusters(np).peak_depth(sc.cluster_id); % get depths of selected clusters

                    % Find cluster_ids within selected depth range
                    keep_mask = peak_depths >= min(depth_range) & peak_depths <= max(depth_range);
                    depth_cluster_ids = sc.cluster_id(keep_mask);

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
                    
                    if contains(z_score_period, 'entire_session')
                        zscore_counts = histcounts(spike_times_this_cluster, time_edges); % use all spike times
                    elseif contains(z_score_period, 'first30secs')
                        baseline_spikes = spike_times_this_cluster(spike_times_this_cluster >= baseline_window(1) & spike_times_this_cluster <= baseline_window(2));
                        zscore_counts = histcounts(baseline_spikes, time_edges); % use just baseline spikes
                    end

                    fig(nCluster)=figure; % open a figure window
                    fig(nCluster).Name=sprintf('%s Grating responses Cluster %i', Stimulus_type, cluster_id(nCluster)); %overall figure title includes the cluster_id (one count higher than zero-based SI output)
                    fig(nCluster).Position = [114 90 770 650]; % sets the size of the figure window
                
                    tiledlayout(5,1); % vertical layout (5 rows × 1 column)
                    % Extract ordered orientations based on presentation sequence
                    ordered_oris = unique(Task_info.stim_orientation, 'stable'); % radians for TRAIN but degrees for GAVNIK
                    
                    for ori = 1:length(ordered_oris)  

                        if contains(plot_type, 'raster')
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.150 0.30], psthBinSize/10); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
                            nexttile % opens a subplot tile
                            plot(rasterX,rasterY,'k','LineWidth',1) % plot the spiking raster for the current orientation
                            xline(0,'r',LineWidth=1) % put a vertical red line at time zero (stimulus onset)
                            xlim([-0.15 0.30])
                            ylim([0 sum(Task_info.stim_orientation==ordered_oris(ori))]) % set the max y coord to be the total # trials with this orientation
                            ylabel('Trial');
                            % imagesc(binnedArray)
                            % xticks([1.5 10.5 15.5 20.5 30.5])
                            % xticklabels([-0.150 -0.050 0 0.050 0.150])
                            % xline(15.5,'r',LineWidth=1)
                            % colorbar
                            % colormap(flipud(gray))
                            title(sprintf('Direction %d%s', round(rad2deg(ordered_oris(ori))), char(176)))
                            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
                        end
                        
                        
                        if contains(plot_type, 'FR')
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.150 0.30], psthBinSize); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
                            
                            mean_trace = mean(binnedArray, 1); %here binned array is just for this orientation
                            z_trace = (mean_trace - mean(zscore_counts)) / std(zscore_counts); %dynamically z-scores to baseline period or entire session depending on definition of z_score_period
                            nexttile;
                            plot(bins, z_trace, 'k', 'LineWidth', 1.5);
                            xline(0, 'r', 'LineWidth', 1);
                            xlim([-0.15 0.30]);
                            ylabel('Z-scored FR');
                            
                            title(sprintf('Direction %d%s', round(rad2deg(ordered_oris(ori))), char(176)));
                            set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
                        end

                        %[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.150 0.30], psthBinSize);
                        z_binnedArray{nCluster}{ori} = (binnedArray - mean(zscore_counts)) ./ std(zscore_counts); %z-score normalisation of the binned Array for this orientation: subtracts the mean and divides by the s.d. of the spike count histogram for the baseline period or whole session, depending on definition of z_score_period
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
                        
                        title(sprintf('Direction %d%s (Extended to show post-stimulus oscillations)', round(rad2deg(ordered_oris(4))), char(176)));
                        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
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
                        
                        title(sprintf('Direction %d%s (Extended to show post-stimulus oscillations)', round(rad2deg(ordered_oris(4))), char(176)));
                        ylabel('Trial');
                        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
                    end


                    % Sanitize Stimulus_type for filenames
                    safeStimulusType = regexprep(Stimulus_type, '[:\\/*?"<>| ]', '_');
                    cluster_depth = clusters(nprobe).peak_depth(cluster_id(nCluster));
                    sgtitle(sprintf('%s - %s: Response of %s Cluster %i (%.0f µm)', subject_number, Stimulus_type, depth_for_analysis, cluster_id(nCluster), cluster_depth), 'Interpreter', 'none');

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

                ordered_oris = unique(Task_info.stim_orientation, 'stable');
                                              
                fig = figure;
                fig.Name = sprintf('Aggregate activity: %s depth range %d–%d μm', depth_for_analysis, min(depth_range), max(depth_range));
                fig.Position = [114 90 770 650];
                %tiledlayout(5,1);

                %for ori = 1:length(ordered_oris)
                ori = 1; % make plots from onset of first stim in sequence
                    
                stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));

                    %[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.150 0.30], psthBinSize);
                [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.30 1.8], psthBinSize);    
                    
                mean_trace = mean(binnedArray, 1); % average binned firing over trials

                if contains(z_score_period, 'none')
                    % Plot raw mean firing rate trace
                    plot(bins, mean_trace, 'b', 'LineWidth', 1.5);
                    ylabel('Mean firing rate (Hz)');
                    ylim ([0 16]);
                else
                    % Compute appropriate histogram counts for z-scoring
                    if contains(z_score_period, 'entire_session')
                        zscore_counts = histcounts(all_spike_times, time_edges);
                        hold on;
                        ylim([-1 6]);
                        ylabel('Z-scored FR (z-scored over entire session)');
                    elseif contains(z_score_period, 'first30secs')
                        baseline_spikes = all_spike_times(all_spike_times >= baseline_window(1) & all_spike_times <= baseline_window(2));
                        zscore_counts = histcounts(baseline_spikes, time_edges); %histcounts counts the number of datapoints in specified time bins
                        hold on;
                        ylim([-1 8]);
                        ylabel('Z-scored FR (z-scored over first 30s baseline)');
                    end

                    z_trace = (mean_trace - mean(zscore_counts)) / std(zscore_counts); %mean(zscore_counts) gives the mean spikes per timebin in the reference period; this is then deducted from the spikecount of each trial-averaged timebin

                    plot(bins, z_trace, 'b', 'LineWidth', 1.5);
                    
                end

                % Define windows after time zero and calc. the peak and mean FR              
                %window_starts = 0:0.15:1.35;
                stim_window_starts = 0:0.3:0.9;
                stim_window_ends = 0.15:0.3:1.05;
                grey_window_starts = [0.2 0.5 0.8 1.1 1.2 1.4]; % look from 50ms after stim offset except at 1.2s (when there would be a stimulus if there was a fifth stimulus)
                %window_ends = 0.15:0.15:1.5;
                grey_window_ends = [0.3 0.6 0.9 1.2 1.35 1.5]; % window ending at 1.35s is where a stimulus would end if there was a fifth stimulus in the sequence)
                peak_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                mean_FR_by_stimwindow = zeros(1, length(stim_window_starts)); % preallocate
                peak_FR_by_greywindow = zeros(1, length(stim_window_starts)); % preallocate
                mean_FR_by_greywindow = zeros(1, length(stim_window_starts)); % preallocate
                                
                yl = ylim; % Get y-axis limits for text placement
                    
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
                                % Define grey intervals
                grey_intervals = [-0.5 0; 0.15 0.3; 0.45 0.6; 0.75 0.9; 1.05 2];
                
                hold on; % Make sure current plot stays visible
                
                % Get current y-axis limits for full vertical shading
                yl = ylim;
                
                % Shade each interval
                for i = 1:size(grey_intervals, 1)
                    x = [grey_intervals(i,1), grey_intervals(i,2), grey_intervals(i,2), grey_intervals(i,1)];
                    y = [yl(1), yl(1), yl(2), yl(2)];
                    fill(x, y, [0.7 0.7 0.7], 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % grey color with transparency
                end

                xline(0, 'k', (sprintf('A onset; %d%s', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(0.30, 'k', (sprintf('B onset; %d%s', round(rad2deg(ordered_oris(2))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(0.60, 'k', (sprintf('C onset; %d%s', round(rad2deg(ordered_oris(3))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(0.90, 'k', (sprintf('D onset; %d%s', round(rad2deg(ordered_oris(4))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');

                xlim([-0.5 2.0]);
                xticks(-0.4:0.2:1.8);
                xlabel('Time (s) since onset of A')
                                
                sgtitle(sprintf('%s day %d - %s - Aggregate single unit activity: %s depth range %d–%d μm', subject_number, experiment_info(nsession).date, Stimulus_type, depth_for_analysis, min(depth_range), max(depth_range)), 'Interpreter', 'none', 'FontSize', 12);
                filename = sprintf('%s - Aggregate single unit %s FRs.png', Stimulus_type, depth_for_analysis);
                save_path = fullfile(pwd, filename);
                exportgraphics(fig, save_path);  % Add this line

                % Also save as .fig
                fig_filename = sprintf('%s - Aggregate single unit %s FRs.fig', Stimulus_type, depth_for_analysis);
                fig_save_path = fullfile(pwd, fig_filename);
                savefig(fig, fig_save_path);
            end
        end

        


        if (contains(Stimulus_type, 'OMIT') || strcmp(Stimulus_type, 'E_CD')) && contains(plot_choice, 'aggregate') && contains(plot_type, 'FR')
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
                
                % Each set of 4 stimuli makes a "trial" (i.e., sequence repeats every 4 stimuli)
                sequence_length = 4;
                stim_oris = Task_info.stim_orientation;
                stim_contrasts = Task_info.stim_contrast;
                stim_onsets = Task_info.stim_onset;

                fig = figure;
                fig.Name = sprintf('Aggregate activity: %s depth range %d–%d μm', depth_for_analysis, min(depth_range), max(depth_range));
                fig.Position = [114 90 770 650];
                %tiledlayout(4,1);

                % Session-wide stats for z-scoring
                
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

                    % Z-score
                z_red = (mean(binned_red,1) - mu) / sigma;
                z_blue = (mean(binned_blue,1) - mu) / sigma;

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

                % Plot
                p1 = plot(bins, z_blue, 'b', 'LineWidth', 1.5); hold on;
                if contains(Stimulus_type, 'OMIT')
                    p2 = plot(bins, z_red, 'r', 'LineWidth', 1.5);
                elseif contains(Stimulus_type, 'E_CD')
                    p2 = plot(bins, z_red, 'g', 'LineWidth', 1.5);
                end

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

                xline(0, 'k', (sprintf('A %d%s onset', round(rad2deg(ordered_oris(1))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(0.30, 'k', (sprintf('B %d%s or grey onset', round(rad2deg(ordered_oris(2))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(0.60, 'k', (sprintf('C %d%s onset', round(rad2deg(ordered_oris(3))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(0.90, 'k', (sprintf('D %d%s onset', round(rad2deg(ordered_oris(4))), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');

                xlim([-0.5 2.0]);
                xticks(-0.4:0.2:1.8);
                xlabel('Time (s)')
                ylabel('Z-scored FR');
                set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);

                    % Make legend labels
                legend(sprintf('ABCD (160x)'), sprintf('%s (40x)', Stimulus_type), 'Location', 'northeast');              
                          
                sgtitle(sprintf('%s - %s - aggregate single unit activity: %s depth range %d–%d μm', subject_number, Stimulus_type, depth_for_analysis, min(depth_range), max(depth_range)), 'Interpreter', 'none');
                filename = sprintf('%s %s %s Aggregate single unit FRs.png', subject_number, Stimulus_type, depth_for_analysis);
                save_path = fullfile(pwd, filename);
                exportgraphics(fig, save_path);
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
                    peak_depths = clusters(np).peak_depth(sc.cluster_id); % get depths of selected clusters

                    % Find cluster_ids within selected depth range
                    keep_mask = peak_depths >= min(depth_range) & peak_depths <= max(depth_range);
                    depth_cluster_ids = sc.cluster_id(keep_mask);

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
                    if contains(z_score_period, 'entire_session')
                        zscore_counts = histcounts(spike_times_this_cluster, time_edges);
                    elseif contains(z_score_period, 'first30secs')
                        baseline_spikes = spike_times_this_cluster(spike_times_this_cluster >= baseline_window(1) & spike_times_this_cluster <= baseline_window(2));
                        zscore_counts = histcounts(baseline_spikes, time_edges);
                    end

                    fig = figure;
                    fig.Name = sprintf('Response %s Cluster %i', depth_for_analysis, cluster_id(nCluster));
                    fig.Position = [114 90 770 650];
                    tiledlayout(5,3);

                    % Session-wide stats for z-scoring
                    mu = mean(zscore_counts);
                    sigma = std(zscore_counts);

                    for stim_pos = 1:sequence_length  % go through each of the 4 stimulus positions
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
                        [~, bins, ~, ~, ~, binned_red] = psthAndBA(spike_times_this_cluster, red_onsets, [-0.15 0.30], psthBinSize);
                        [~, ~, ~, ~, ~, binned_blue] = psthAndBA(spike_times_this_cluster, blue_onsets, [-0.15 0.30], psthBinSize);

                        % Z-score
                        z_red = (mean(binned_red,1) - mu) / sigma;
                        z_blue = (mean(binned_blue,1) - mu) / sigma;
                        
                        % Set raster window
                        rasterWindow = [-0.15 0.30];
                        
                        % Get raster for blue
                        [~, ~, rasterX_blue, rasterY_blue, ~, ~] = psthAndBA(spike_times_this_cluster, blue_onsets, rasterWindow, psthBinSize/10);
                        nexttile;
                        plot(rasterX_blue, rasterY_blue, 'b', 'LineWidth', 1);
                        xline(0, 'k', 'LineWidth', 1);
                        xlim(rasterWindow);
                        ylim([0 length(blue_onsets)]);
                        ylabel('Trial');
                        title(sprintf('Stim %d: raster', stim_pos));
                        set(gca, 'TickDir', 'out', 'box', 'off', 'Color', 'none', 'FontSize', 12);
                        
                        % Get raster for red (or green if E_CD)
                        [~, ~, rasterX_red, rasterY_red, ~, ~] = psthAndBA(spike_times_this_cluster, red_onsets, rasterWindow, psthBinSize/10);
                        nexttile;
                        color_raster = 'r';
                        if contains(Stimulus_type, 'E_CD'), color_raster = 'g'; end
                        plot(rasterX_red, rasterY_red, color_raster, 'LineWidth', 1);
                        xline(0, 'k', 'LineWidth', 1);
                        xlim(rasterWindow);
                        ylim([0 length(red_onsets)]);
                        ylabel('Trial');
                        title(sprintf('Stim %d: raster', stim_pos));
                        set(gca, 'TickDir', 'out', 'box', 'off', 'Color', 'none', 'FontSize', 12);


                        % Plot both traces
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
                        blue_label = sprintf('%d° (contrast %d)', round(rad2deg(ori_blue)), round(contrast_blue));
                        red_label = sprintf('%d° (contrast %d)', round(rad2deg(ori_red)), round(contrast_red));
                        % Plot
                        nexttile;
                        p1 = plot(bins, z_blue, 'b', 'LineWidth', 1.5); hold on;
                        if contains(Stimulus_type, 'OMIT')
                            p2 = plot(bins, z_red, 'r', 'LineWidth', 1.5);
                        elseif contains(Stimulus_type, 'E_CD')
                            p2 = plot(bins, z_red, 'g', 'LineWidth', 1.5);
                        end
                        xline(0,'k','LineWidth',1);
                        xlim([-0.15 0.30]);
                        ylabel('Z-scored FR');
                        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
                        legend([p1 p2], {blue_label, red_label}, 'Location', 'northeast');  
                        legend boxoff;
                        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
                    end
                    
                    % === 5th subplot: Extended response to final stimulus (stim_pos = 4) ===
                    stim_pos = 4;  % last stimulus position
                    blue_onsets = [];
                    red_onsets = [];

                    for i = stim_pos:sequence_length:length(stim_contrasts)
                        curr = i;
                        prev = i - 1;
                        next = i + 1;
                        twoback = i - 2;
                    
                        has_zero_surround = false;
                        if prev >= 1 && stim_contrasts(prev) == 0, has_zero_surround = true; end
                        if next <= length(stim_contrasts) && stim_contrasts(next) == 0, has_zero_surround = true; end
                        if twoback >= 1 && stim_contrasts(twoback) == 0, has_zero_surround = true; end
                        if stim_contrasts(curr) == 0, has_zero_surround = true; end
                    
                        if has_zero_surround
                            red_onsets(end+1) = stim_onsets(curr);
                        else
                            blue_onsets(end+1) = stim_onsets(curr);
                        end
                    end
                    
                    nexttile([1 3]);
                    
                    % Longer PSTH window for extended post-stimulus response
                    long_window = [-0.02 0.75]; 
                    [~, bins_long, ~, ~, ~, binned_red_long] = psthAndBA(spike_times_this_cluster, red_onsets, long_window, psthBinSize);
                    [~, ~, ~, ~, ~, binned_blue_long] = psthAndBA(spike_times_this_cluster, blue_onsets, long_window, psthBinSize);
                    
                    z_red_long = (mean(binned_red_long,1) - mu) / sigma;
                    z_blue_long = (mean(binned_blue_long,1) - mu) / sigma;
                    
                    p1 = plot(bins_long, z_blue_long, 'b', 'LineWidth', 1.5); hold on;
                    if contains(Stimulus_type, 'OMIT')
                        p2 = plot(bins_long, z_red_long, 'r', 'LineWidth', 1.5);
                    elseif contains(Stimulus_type, 'E_CD')
                        p2 = plot(bins_long, z_red_long, 'g', 'LineWidth', 1.5);
                    end
                    xline(0,'k','LineWidth',1);
                    xlim(long_window);
                    xticks(0:0.05:0.75);
                    ylabel('Z-scored FR');
                    title(sprintf('Stimulus %d (Extended to show post-stimulus oscillations)', stim_pos));
                    legend([p1 p2], {'Blue', 'Red'}, 'Location', 'northeast'); legend boxoff;
                    
                   
                    
                    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
                    cluster_depth = clusters(nprobe).peak_depth(cluster_id(nCluster));
                    sgtitle(sprintf('%s - %s: Response of %s Cluster %i (%.0f µm)', subject_number, Stimulus_type, depth_for_analysis, cluster_id(nCluster), cluster_depth), 'Interpreter', 'none');
                    saveas(fig, sprintf('%s_%s_Cluster_%i_FR.pdf', Stimulus_type, depth_for_analysis, cluster_id(nCluster))); %saveas saves to cd
                end
            end
        end



        
        if strcmp(Stimulus_type, 'DCBA') && contains(plot_choice, 'aggregate') && contains(plot_type, 'FR')
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
                z_normal = (mean(ba_normal,1) - mu) / sigma;

                    % Reversed trials (red)
                [~, ~, ~, ~, ~, ba_reversed] = psthAndBA(all_spike_times, reversed_onsets_by_pos{pos}, [-0.3 1.80], psthBinSize);
                z_reversed = (mean(ba_reversed,1) - mu) / sigma;

                    % Plot
                plot(bins, z_normal, 'b', 'LineWidth', 1.5); hold on;
                plot(bins, z_reversed, 'r', 'LineWidth', 1.5);
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

                xline(0, 'k', (sprintf('A %d%s or D onset', normal_sequence(1), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(0.30, 'k', (sprintf('B %d%s or C onset', normal_sequence(2), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(0.60, 'k', (sprintf('C %d%s or B onset', normal_sequence(3), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(0.90, 'k', (sprintf('D %d%s or A onset', normal_sequence(4), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');

                xlim([-0.5 2.0]);
                xticks(-0.4:0.2:1.8);
                xlabel('Time (s)')
                ylabel('Z-scored FR');
                    
                legend(sprintf('ABCD (160x)'), sprintf('DCBA (40x)'), 'Location', 'northeast');
                set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
                
                sgtitle(sprintf('%s - %s - aggregate single unit activity, %s depth range %d–%d μm', subject_number, Stimulus_type, depth_for_analysis, min(depth_range), max(depth_range)), 'Interpreter', 'none', 'FontSize', 12);

                filename = sprintf('%s %s %s Aggregate single unit FRs.png', subject_number, Stimulus_type, depth_for_analysis);
                save_path = fullfile(pwd, filename);
                exportgraphics(fig, save_path);
                % Also save as .fig
                fig_filename = sprintf('%s %s - Aggregate single unit %s FRs.fig', subject_number, Stimulus_type, depth_for_analysis);
                fig_save_path = fullfile(pwd, fig_filename);
                savefig(fig, fig_save_path);
            end
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
                
                ori = 1;
                stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));

                [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.3 1.35], psthBinSize);
                    
                mean_trace = mean(binnedArray, 1); % average over trials
                z_trace = (mean_trace - mean(zscore_counts)) / std(zscore_counts); % z-score using session-wide stats

                plot(bins, z_trace, 'b', 'LineWidth', 1.5);

                 % Define windows after time zero and calc. the peak and mean FR              
                stim_window_starts = 0:0.15:0.45;
                stim_window_ends = 0.15:0.15:0.6;
                grey_window_starts = 0.6:0.15:1.2; % look in each 150ms window
                grey_window_ends = 0.75:0.15:1.35; % look in each 150ms window
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
                grey_intervals = [-0.5 0; 0.6 1.5];
                
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

                xline(0, 'k', (sprintf('A onset; %d%s', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(0.15, 'k', (sprintf('B onset; %d%s', round(ordered_oris(2)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(0.30, 'k', (sprintf('C onset; %d%s', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');
                xline(0.45, 'k', (sprintf('D onset; %d%s', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left');

                xlim([-0.5 1.5]);
                xticks(-0.4:0.2:1.4);
                ylabel('Z-scored FR');
                xlabel('Time since onset of A (s)')
                set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
                
                sgtitle(sprintf('%s - %s day %d - Aggregate single unit activity: %s depth range %d–%d μm', subject_number, Stimulus_type, experiment_info(nsession).date, depth_for_analysis, min(depth_range), max(depth_range)), 'Interpreter', 'none');

                filename = sprintf('%s - Aggregate single unit %s FRs.png', Stimulus_type, depth_for_analysis);
                save_path = fullfile(pwd, filename);
                exportgraphics(fig, save_path);  

                % Also save as .fig
                fig_filename = sprintf('%s - Aggregate single unit %s FRs.fig', Stimulus_type, depth_for_analysis);
                fig_save_path = fullfile(pwd, fig_filename);
                savefig(fig, fig_save_path);

            end
        end




        if (contains(Stimulus_type, 'GAVNIK_ABCD') || contains(Stimulus_type, 'GAVNIK_DCBA')) &&...
           contains(plot_choice, 'single_units') 
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
                    if contains(z_score_period, 'entire_session')
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
                            
                            title(sprintf('Orientation %d%s', round(ordered_oris(ori)), char(176)));
                            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
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
                            
                            title(sprintf('Orientation %d%s', round(ordered_oris(ori)), char(176)));                            
                            set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
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
                        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
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
                        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
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



        
        if (contains(Stimulus_type, 'GAVNIK_A_CD') || contains(Stimulus_type, 'GAVNIK_E_CD')) &&...
           contains(plot_choice, 'single_units') 
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
                    if contains(z_score_period, 'entire_session')
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
                            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
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
                            
                            set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
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
                        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
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
                        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
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
