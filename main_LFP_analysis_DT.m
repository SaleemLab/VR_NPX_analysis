%% Main pipeline for extracting and saving and analysing LFP including ripple event and theta cycle detection
% This is a higher-order multi-purpose pipeline that calls dependent functions for
% different analysis piepline

% For mapping of visual receptive field, Please refer to
% Sparse_noise_RF_mapping_masa.mat (subject to change)


%% Set path
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))


%% spatial cell

clear all
SUBJECTS = {'M23087'};
option = 'bilateral';

SUBJECTS = {'M23034'};
option = 'V1-MEC';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'Track'; % analysing 3rd session of M23034


for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'))

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'))
        clusters = clusters_ks3;
        sorting_option = 'spikeinterface';

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'))
        %         %plot speed of each lap
        %         no_lap = size(Task_info.start_time_all,1);
        %         no_bin = 70;
        %         psth_speed = zeros(no_lap,no_bin);
        %         for iLap = 1:no_lap
        %             lap_index = Behaviour.tvec >= Task_info.start_time_all(iLap) & Behaviour.tvec <= Task_info.end_time_all(iLap);
        %             timevec_lap = Behaviour.tvec(lap_index);
        %             position_lap = Behaviour.position(lap_index);
        %             speed_lap = Behaviour.speed(lap_index);
        %             [N,edges,bin] = histcounts(position_lap,no_bin);
        %             for iBin = 1:no_bin
        %                 psth_speed(iLap,iBin) = median(speed_lap(bin == iBin));
        %             end
        %         end
        %         figure;hold on;
        %         for iBlock = 1:5
        %
        %             plot(1:2:139,mean(psth_speed((iBlock-1)*40+1:iBlock*40,:),'omitnan'))
        %         end
        options.importMode = 'KS';
        [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,1);

        for nprobe = 1:length(clusters)
            clusters(nprobe).region = strings(length(clusters(nprobe).cluster_id),1);
            clusters(nprobe).region(:) = 'n.a';
            options = session_info(n).probe(nprobe);

            if options.probe_V1 == options.probe_id
                V1_channels = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by probe');
                HPC_channels = determine_region_channels(best_channels{nprobe},options,'region','HPC','group','by probe');
                clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,V1_channels))) = 'V1';
                clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,HPC_channels))) = 'HPC';

            elseif options.probe_MEC  == options.probe_id
                %             V1_channels = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by probe');
                MEC_channels = determine_region_channels(best_channels{nprobe},options,'region','MEC_entry','group','by probe');
                %             clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,V1_channels))) = 'V1';
                clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,MEC_channels))) = 'MEC';
            end
        end

        %%%%%%%%%%%%%%%%%% place holder for merging units
        %%%%%%%%%%%%%%%%%%%% load match id

        metric_param = create_cluster_selection_params('sorting_option',sorting_option);
        metric_param.unstable_ids = @(x) x==0;

        place_field = struct();

        for nprobe = 1:length(clusters)
            options = session_info(n).probe(nprobe);

            load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('probe%ium_merge_suggestion.mat',options.probe_id)))

            %             [~,cluster_id] = select_clusters(clusters_combined,metric_param);
            merged_cluster_temp  = merge_cluster(clusters(nprobe),match_ids);
            merged_clusters(nprobe) = select_clusters(merged_cluster_temp,metric_param);
        end
        
        % save merged cluster variables (useful more reactivation activity detection)
        if contains(stimulus_name{n},'Masa')
            save(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'merged_clusters');
        else
            save(fullfile(options.ANALYSIS_DATAPATH,'merged_clusters.mat'))
        end
        %

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end
        
        place_fields = calculate_place_fields_masa_NPX_against_shuffle(clusters_combined,Task_info,Behaviour,[0 140],2,[]);
        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_place_fields.mat'),'place_fields')

        % Plotting
        metric_param = create_cluster_selection_params('sorting_option',sorting_option);
        metric_param.unstable_ids = @(x) x==0;
        for nprobe = 1:length(merged_clusters)

            [C,ia,ic] = unique(merged_clusters(nprobe).merged_cluster_id);

            plot_raster_both_track(merged_clusters(nprobe).spike_times,merged_clusters(nprobe).merged_spike_id,Task_info,Behaviour,[5 1],[0 140],2,...
                'unit_depth',merged_clusters(nprobe).peak_depth(ia),'unit_region',merged_clusters(nprobe).region(ia),'unit_id',C);

            plot_raster_both_track_extended(merged_clusters(nprobe).spike_times,merged_clusters(nprobe).merged_spike_id,Task_info,Behaviour,[5 1],[0 140],2,...
                'unit_depth',merged_clusters(nprobe).peak_depth(ia),'unit_region',merged_clusters(nprobe).region(ia),'unit_id',C);

            plot_raster_end_of_track(merged_clusters(nprobe).spike_times,merged_clusters(nprobe).merged_spike_id,Task_info,Behaviour,[5 1],[0 3],2,...
                'unit_depth',merged_clusters(nprobe).peak_depth(ia),'unit_region',merged_clusters(nprobe).region(ia),'unit_id',C);

% plot_raster_both_track_extended(HPC_clusters.spike_times,HPC_clusters.merged_spike_id,Task_info,Behaviour,[3 1],[0 140],2,...
%                 'unit_depth',HPC_clusters.peak_depth(ia),'unit_region',HPC_clusters.region(ia),'unit_id',C,'place_fields',place_fields);
        end


        % plot spatial raster plot
        if ~isempty(V1_clusters_L)
            %             plot_raster_both_track(V1_clusters_L,Task_info,Behaviour,[5 1],[0 140],5);
            %             plot_spatial_CCG(V1_clusters_L,Task_info,Behaviour,[3 1],[0 140],2)
        end

        if ~isempty(V1_clusters_R)
            %             plot_raster_both_track(V1_clusters_R,Task_info,Behaviour,[5 1],[0 140],5);
            %             plot_spatial_CCG(V1_clusters_R,Task_info,Behaviour,[3 1],[0 140],2)
        end

        % HPC
        if ~isempty(HPC_clusters_L)
            %             plot_raster_both_track(HPC_clusters_L,Task_info,Behaviour,[5 1],[0 140],5);
            %             plot_spatial_CCG(HPC_clusters_L,Task_info,Behaviour,[3 1],[0 140],2)
        end

        if ~isempty(HPC_clusters_R)
            %             plot_raster_both_track(HPC_clusters_R,Task_info,Behaviour,[5 1],[0 140],5);
            %             plot_spatial_CCG(HPC_clusters_R,Task_info,Behaviour,[3 1],[0 140],2)
        end

        %         plot_raster_single_track(HPC_clusters_L,Task_info,Behaviour,[5 1],[0 140],5);

        if ~isempty(HPC_clusters_combined)

        end

        if ~isempty(V1_clusters_combined)
            %             plot_raster_both_track_extended(HPC_clusters_combined,Task_info,Behaviour,[3 1],[0 140],2);
            plot_raster_both_track(HPC_clusters_combined,Task_info,Behaviour,[3 1],[0 140],2);
            %             plot_spatial_CCG(V1_clusters_combined,Task_info,Behaviour,[3 1],[0 140],2)
        end




        % Spatial modulation
        x_bin_size = mean(diff(place_fields_V1_L(1).x_bin_centres));
        SMI = calculate_spatial_modulation_index(V1_clusters_L,Task_info,Behaviour,[0 140],x_bin_size,'place_fields',place_fields_V1_L,'subplot_xy',[3 1],'plot_option',1)
        SMI = calculate_spatial_modulation_index(V1_clusters_R,Task_info,Behaviour,[0 140],x_bin_size,'place_fields',place_fields_V1_R,'subplot_xy',[3 1],'plot_option',1)


        metric_param.region = @(x) contains(x,'HPC_L');
        [HPC_clusters_L,cluster_id] = select_clusters(clusters(1),metric_param);

        SMI = calculate_spatial_modulation_index(HPC_clusters_L,Task_info,Behaviour,[0 140],x_bin_size,'place_fields',place_fields_HPC_L,'subplot_xy',[3 1],'plot_option',1)
        SMI = calculate_spatial_modulation_index(HPC_clusters_R,Task_info,Behaviour,[0 140],x_bin_size,'place_fields',place_fields_HPC_R,'subplot_xy',[3 1],'plot_option',1)
        %         calculate_spatial_modulation_index(place_fields_V1_L);

        % plot population vector correlation
        options.region = 'HPC';
        [normalised_raw_matrix,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = ...
            plot_place_cell_map_correlation(HPC_clusters_combined,place_fields_HPC_combined,...
            Task_info,Behaviour,options);

        save_all_figures((fullfile(options.ANALYSIS_DATAPATH,'PV correlation')),[])
        options = rmfield(options,'probe_combined');
        close all


        % Bayesian decoding 10 fold cross validated
        DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'RUN Bayesian decoding'));
        if isempty(DIR)
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'RUN Bayesian decoding'))
        end

        options.region = 'HPC Left';
        [probability_ratio_RUN_lap_HPC_L,estimated_position_lap_CV_HPC_L.track] = bayesian_decoding_RUN_lap_cross_validation(HPC_clusters_L,place_fields_HPC_L,Behaviour,Task_info,options);
        %
        save_all_figures((fullfile(options.ANALYSIS_DATAPATH,'RUN Bayesian decoding')),[])

        options.region = 'HPC Right';
        [probability_ratio_RUN_lap_HPC_R,estimated_position_lap_CV_HPC_R.track] = bayesian_decoding_RUN_lap_cross_validation(HPC_clusters_R,place_fields_HPC_R,Behaviour,Task_info,options);
        save_all_figures((fullfile(options.ANALYSIS_DATAPATH,'RUN Bayesian decoding')),[])

        if ~isempty(HPC_clusters_combined)
            options.region = 'HPC Combined';
            [probability_ratio_RUN_lap_HPC_combined,estimated_position_lap_CV_HPC_combined.track] = bayesian_decoding_RUN_lap_cross_validation(HPC_clusters_combined,place_fields_HPC_combined,Behaviour,Task_info,options);
            save_all_figures((fullfile(options.ANALYSIS_DATAPATH,'RUN Bayesian decoding')),[])

        else
            probability_ratio_RUN_lap_HPC_combined = [];
            estimated_position_lap_CV_HPC_combined = [];
        end

        save('estimated_position_lap_CV_HPC.mat','estimated_position_lap_CV_HPC_combined','estimated_position_lap_CV_HPC_L','estimated_position_lap_CV_HPC_R')
        save('probability_ratio_RUN_lap_HPC.mat','probability_ratio_RUN_lap_HPC_combined',"probability_ratio_RUN_lap_HPC_L","probability_ratio_RUN_lap_HPC_R")

        if ~isempty(V1_clusters_combined)
            options.region = 'V1 combined';
            V1_clusters_combined = combine_clusters_from_multiple_probes(V1_clusters_L,V1_clusters_R);
            [probability_ratio_RUN_lap_V1_combined,estimated_position_lap_CV_V1_combined.track] = bayesian_decoding_RUN_lap_cross_validation(V1_clusters_L,place_fields_V1_L,Behaviour,Task_info,options);
            save_all_figures((fullfile(options.ANALYSIS_DATAPATH,'RUN Bayesian decoding')),[])
        else
            probability_ratio_RUN_lap_V1_combined = [];
            estimated_position_lap_CV_V1_combined = [];
        end

        options.region = 'V1 Right';
        [probability_ratio_RUN_lap_V1_R,estimated_position_lap_CV_V1_R.track] = bayesian_decoding_RUN_lap_cross_validation(V1_clusters_R,place_fields_V1_R,Behaviour,Task_info,options);
        save_all_figures((fullfile(options.ANALYSIS_DATAPATH,'RUN Bayesian decoding')),[])

        options.region = 'V1 Left';
        [probability_ratio_RUN_lap_V1_combined,estimated_position_lap_CV_V1_combined.track] = bayesian_decoding_RUN_lap_cross_validation(V1_clusters_combined,place_fields_V1_combined,Behaviour,Task_info,options);
        save_all_figures((fullfile(options.ANALYSIS_DATAPATH,'RUN Bayesian decoding')),[])

        save('estimated_position_lap_CV_V1.mat','estimated_position_lap_CV_V1_L',"estimated_position_lap_CV_V1_R","estimated_position_lap_CV_V1_combined")
        save('probability_ratio_RUN_lap.mat','probability_ratio_RUN_lap_V1_L','probability_ratio_RUN_lap_V1_R',"probability_ratio_RUN_lap_V1_combined")


        %         %GLM analysis
        %         spatial_modulation_GLM_analysis(V1_clusters_L,[],Behaviour,Task_info);
        %
        %         spatial_modulation_GLM_analysis(V1_clusters_L,place_fields_V1_L(1).cluster_id(place_fields_V1_L(1).all_good_cells_LIBERAL),Behaviour,Task_info);
    end
end

%% Extract and save LFP from noise channel, cortical L5 channel

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'..','extracted_PSD.mat'),'power');
        raw_LFP = [];
        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);
            selected_channels = [];
            channel_regions = [];
            shank_id = [];

            all_fields = fieldnames(best_channels{nprobe});
            all_fields = {all_fields{contains(all_fields,'depth')}};
            all_fields = erase(all_fields,'_depth');
            for nregion = 1:length(all_fields)
                region_name = all_fields{nregion};
                [~,channels_temp] = determine_region_channels(best_channels{nprobe},options,'region',region_name,'group','by shank');
                %                 noise_channel{options.probe_hemisphere} = channels_temp;
                selected_channels = [selected_channels channels_temp];
                channel_regions = [channel_regions nregion*ones(1,length(channels_temp))];
                shank_id = [shank_id unique(ceil(best_channels{nprobe}.xcoord/250))];
            end
            channel_regions(isnan(selected_channels)) = []; % remove nan channel (Missing best channels for some shanks e.g. only 3 shanks with CA1)
            shank_id(isnan(selected_channels)) = [];
            selected_channels(isnan(selected_channels)) = [];

            %     column = 1;
            if nprobe ~= 1
                % Information about AP band probe 1 sample size to align probe
                % 2 LFP traces.
                session_info(n).probe(1).importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP
                [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(session_info(n).probe(1),[]);
                [raw_LFP{nprobe} tvec SR chan_config ~] = load_LFP_NPX(options,[],'probe_no',nprobe,'probe_1_total_sample',imecMeta.nFileSamp,'selected_channels',selected_channels);
            else
                %         session_info.probe(1).importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP
                %         [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(session_info.probe(1),[]);
                [raw_LFP{nprobe} tvec SR chan_config ~] = load_LFP_NPX(options,[],'selected_channels',selected_channels);
            end

            % Save downsampled LFP from key channels
            options.importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP
            [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(options,[]);

            if isfield(options,'probe_hemisphere')
                LFP(nprobe).probe_hemisphere = options.probe_hemisphere;
                LFP(nprobe).probe_id = options.probe_id;
            else
                LFP(nprobe).probe_id = options.probe_id;
            end

            LFP(nprobe).tvec = tvec;
            for nregion = 1:length(all_fields)
                LFP(nprobe).(all_fields{nregion}) = raw_LFP{nprobe}(channel_regions == nregion,:);
                LFP(nprobe).(sprintf('%s_shank_id',all_fields{nregion})) = shank_id(channel_regions == nregion); % only avaliable shanks
                LFP(nprobe).(sprintf('%s_channel',all_fields{nregion})) = selected_channels(channel_regions == nregion);
                LFP(nprobe).(sprintf('%s_depth',all_fields{nregion})) = chan_config.Ks_ycoord(selected_channels(channel_regions == nregion));
                LFP(nprobe).(sprintf('%s_power',all_fields{nregion})) = power{nprobe}(selected_channels(channel_regions == nregion),:);
                %                 LFP(nprobe).(sprintf('%s_depth',all_fields{nregion})) = chan_config.Ks_ycoord(selected_channels(channel_regions == nregion));
            end

        end
        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP')

        clear replay reactivations ripples raw_LFP

        if contains(stimulus_name{n},'Masa')
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        else
            load(fullfile(options.ANALYSIS_DATAPATH,'merged_clusters.mat'))
        end

        if length(merged_clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end

        %%%%%%%%%%%%%%%%%%
        % Ripple and candidate events detection
        %%%%%%%%%%%%%%%%%%
        for nprobe = 1:length(session_info(n).probe)
            
            if session_info(n).probe(1).probe_id == 0 
                if isfield(session_info(n).probe(1),'probe_V1') & session_info(n).probe(1).probe_V1 == 1
                   nprobe = 2; 
                end
            end
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
            %                 Behavioural state detection
            speed_interp = interp1(Behaviour.tvec,Behaviour.speed,LFP(nprobe).tvec','linear');
            speedTreshold = 1;
            if isfield(LFP(nprobe),'L4')
                cortex_LFP = LFP(nprobe).L4;
            elseif isfield(LFP(nprobe),'L5')
                cortex_LFP = LFP(nprobe).L5;
            elseif isfield(LFP(nprobe),'MEC')

            end

            if isfield(LFP(nprobe),'CA1')
                CA1_LFP = LFP(nprobe).CA1;
                cortex_channel = best_channels{nprobe}.L5_channel;
                CA1_channel = best_channels{nprobe}.CA1_channel;
                [freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
                    [LFP(nprobe).tvec' cortex_LFP'],[LFP(nprobe).tvec' CA1_LFP'],...
                    [LFP(nprobe).tvec' speed_interp],speedTreshold);
            end
            %             behavioural_state.freezing = freezing;
            behavioural_state.quietWake = quietWake;
            behavioural_state.SWS = SWS;
            behavioural_state.REM = REM;
            behavioural_state.movement = movement;


            % Detect CA1 populational bursting events (Candidate events)
            zscore_min = 0;
            zscore_max = 3;
            
            HPC_channels = determine_region_channels(best_channels{2},options,'region','CA1','group','by probe');
%             merged_clusters.region
            sorting_option = 'spikeinterface';
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.peak_channel = @(x) ismember(x,HPC_channels);
            CA1_clusters = select_clusters(merged_clusters(2),metric_param);            
            [replay,reactivations] = detect_candidate_events_masa(LFP(2).tvec,CA1_LFP,...
                [CA1_clusters.merged_spike_id CA1_clusters.spike_times],Behaviour,zscore_min,zscore_max,options);

            if length(reactivations.probe(probe_no).onset) < 50
                HPC_channels = determine_region_channels(best_channels{nprobe},options,'region','HPC','group','by probe');
                %             merged_clusters.region
                sorting_option = 'spikeinterface';
                metric_param = create_cluster_selection_params('sorting_option',sorting_option);
                metric_param.peak_channel = @(x) ismember(x,HPC_channels);
                CA1_clusters(nprobe) = select_clusters(merged_clusters(nprobe),metric_param);
                [replay(probe_no),reactivations(probe_no)] = detect_candidate_events_masa(LFP(nprobe).tvec,CA1_LFP,...
                    [CA1_clusters.merged_spike_id CA1_clusters.spike_times],Behaviour,zscore_min,zscore_max,options);
            end


            MEC_channels = determine_region_channels(best_channels{1},options,'region','MEC_entry','group','by probe');
            %             merged_clusters.region
            sorting_option = 'spikeinterface';
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.peak_channel = @(x) ismember(x,MEC_channels);
            metric_param.cell_type = @(x) x==1;
            MEC_non_ripple_clusters = select_clusters(merged_clusters(1),metric_param);

            MEC_channels = determine_region_channels(best_channels{1},options,'region','MEC_ripple','group','by probe');
            %             merged_clusters.region
            sorting_option = 'spikeinterface';
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.peak_channel = @(x) ismember(x,MEC_channels);
            metric_param.cell_type = @(x) x==1;
            MEC_ripple_clusters = select_clusters(merged_clusters(1),metric_param);


            V1_channels = determine_region_channels(best_channels{2},options,'region','V1','group','by probe');
            %             merged_clusters.region
            sorting_option = 'spikeinterface';
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.peak_channel = @(x) ismember(x,V1_channels);
%             metric_param.cell_type = @(x) x==1;
            V1_clusters = select_clusters(merged_clusters(2),metric_param);


            %
            %             [reactivations.probe(nprobe).awake_offset,reactivations.probe(nprobe).awake_index] = RestrictInts(reactivations.probe(nprobe).offset,behavioural_state.quietWake);
            %             reactivations.probe(nprobe).awake_onset = reactivations.probe(nprobe).onset(reactivations.probe(nprobe).awake_index);
            %             reactivations.probe(nprobe).awake_peaktimes = reactivations.probe(nprobe).peaktimes(reactivations.awake_index);

            % Detect CA1 ripple events
            [ripples] = FindRipples_masa(CA1_LFP',LFP(nprobe).tvec','behaviour',Behaviour,'minDuration',20,'durations',[30 200],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
                'noise',LFP(nprobe).surface','passband',[125 300],'thresholds',[3 5],'show','on')

            subplot(2,4,1)
            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(CA1_clusters.spike_times, reactivations.onset, [-0.5 1], 0.001);
            plot(rasterX,rasterY); hold on
            subplot(2,4,5)
            plot(bins,zscore(psth));
            title('HPC')

            subplot(2,4,2)
            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(MEC_ripple_clusters.spike_times, reactivations.onset, [-0.5 1], 0.001);
            plot(rasterX,rasterY); hold on
            subplot(2,4,6)
            plot(bins,zscore(psth));
            title('MEC ripple')

            subplot(2,4,3)
            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(MEC_non_ripple_clusters.spike_times, reactivations.onset, [-0.5 1], 0.001);
            plot(rasterX,rasterY); hold on
            subplot(2,4,7)
            plot(bins,zscore(psth));
            title('MEC non ripple')
            
            subplot(2,4,4)
          [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(V1_clusters.spike_times, reactivations.onset, [-0.5 1], 0.001);
            plot(rasterX,rasterY); hold on
            subplot(2,4,8)
            plot(bins,zscore(psth));
            title('V1')
            
            close all

        end

        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

            for event = 1:length(ripples.probe(probe_no).onset)
                ripples.probe(probe_no).speed(event) = mean(Behaviour.speed(find(Behaviour.sglxTime >= ripples.probe(probe_no).onset(event) & Behaviour.sglxTime <= ripples.probe(probe_no).offset(event))));
            end

            if contains(Stimulus_type,'RUN') % If reactivation events during lap running
                [reactivations.probe(probe_no).T1_offset,reactivations.probe(probe_no).T1_index] = RestrictInts(reactivations.probe(probe_no).offset',[lap_times(1).start' lap_times(1).end']); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                reactivations.probe(probe_no).T1_onset = reactivations.probe(probe_no).onset(reactivations.probe(probe_no).T1_index);
                reactivations.probe(probe_no).T1_midpoint = reactivations.probe(probe_no).midpoint(reactivations.probe(probe_no).T1_index);

                [reactivations.probe(probe_no).T2_offset,reactivations.probe(probe_no).T2_index] = RestrictInts(reactivations.probe(probe_no).offset',[lap_times(2).start' lap_times(2).end']); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                reactivations.probe(probe_no).T2_onset = reactivations.probe(probe_no).onset(reactivations.probe(probe_no).T2_index);
                reactivations.probe(probe_no).T2_midpoint = reactivations.probe(probe_no).midpoint(reactivations.probe(probe_no).T2_index);
            end

            if contains(Stimulus_type,'RUN') % If reactivation events during lap running
                [ripples.probe(probe_no).T1_offset,ripples.probe(probe_no).T1_index] = RestrictInts(ripples.probe(probe_no).offset',[lap_times(1).start' lap_times(1).end']); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                ripples.probe(probe_no).T1_onset = ripples.probe(probe_no).onset(ripples.probe(probe_no).T1_index);
                ripples.probe(probe_no).T1_peaktimes = ripples.probe(probe_no).peaktimes(ripples.probe(probe_no).T1_index);

                [ripples.probe(probe_no).T2_offset,ripples.probe(probe_no).T2_index] = RestrictInts(ripples.probe(probe_no).offset',[lap_times(2).start' lap_times(2).end']); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                ripples.probe(probe_no).T2_onset = ripples.probe(probe_no).onset(ripples.probe(probe_no).T2_index);
                ripples.probe(probe_no).T2_peaktimes = ripples.probe(probe_no).peaktimes(ripples.probe(probe_no).T2_index);
            end
        end

        save(sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks')),'replay','reactivations')

        delete(sprintf('extracted_ripples_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        save(sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks')),'ripples')


    end
end