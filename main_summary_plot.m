%% Script for plotting key summary figures of the experiment 


addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))


%% Cluster density plots

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','extracted_PSD.mat'),'power');

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

        if exist(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            clusters = merged_clusters;
            sorting_option = 'spikeinterface';
        elseif exist(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            clusters = clusters_ks3;
            sorting_option = 'spikeinterface';
        else
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            sorting_option = 'old';
        end

        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        
        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);
            options.importMode = 'KS';
            [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(options,[]);
            xcoord = chan_config.Ks_xcoord;
            best_channels{nprobe}.xcoord
            xcoord_avaliable = best_channels{nprobe}.xcoord(~isnan(best_channels{nprobe}.surface_depth));
            
            plot_cluster_density_profile(power{nprobe}(xcoord == xcoord_avaliable(1),:),chan_config,chan_config(xcoord == xcoord_avaliable(1),:),best_channels{nprobe},clusters,options);
        end

        if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','probe_depth_profile'))== 0
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','probe_depth_profile'))
        end

        save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','probe_depth_profile'),[])
        
        
        % Plotting raster plot
        metric_param = create_cluster_selection_params('sorting_option',sorting_option);
        metric_param.unstable_ids = @(x) x==0;
        for nprobe = 1:length(merged_clusters)

            [C,ia,ic] = unique(merged_clusters(nprobe).merged_cluster_id);

            plot_raster_both_track(merged_clusters(nprobe).spike_times,merged_clusters(nprobe).merged_spike_id,Task_info,Behaviour,[5 1],[0 140],2,...
                'unit_depth',merged_clusters(nprobe).peak_depth(ia),'unit_region',merged_clusters(nprobe).region(ia),'unit_id',C);

            plot_raster_both_track_extended(merged_clusters(nprobe).spike_times,merged_clusters(nprobe).merged_spike_id,Task_info,Behaviour,[5 1],[0 140],2,...
                'unit_depth',merged_clusters(nprobe).peak_depth(ia),'unit_region',merged_clusters(nprobe).region(ia),'unit_id',C);

            plot_perievent_spiketimes(merged_clusters(nprobe).spike_times,merged_clusters(nprobe).merged_spike_id,Task_info,Behaviour,[5 1],[-1 3],0.02,...
                'unit_depth',merged_clusters(nprobe).peak_depth(ia),'unit_region',merged_clusters(nprobe).region(ia),'unit_id',C,'event_times',Task_info.end_time_all',...
                'event_label',{'Track 1','Track 2'},'event_id',Task_info.track_ID_all);

            plot_perievent_spiketimes(merged_clusters(nprobe).spike_times,merged_clusters(nprobe).merged_spike_id,Task_info,Behaviour,[5 1],[-1 3],0.02,...
                'unit_depth',merged_clusters(nprobe).peak_depth(ia),'unit_region',merged_clusters(nprobe).region(ia),'unit_id',C,'event_times',Task_info.start_time_all,...
                'event_label',{'Track 1','Track 2'},'event_id',Task_info.track_ID_all);

            plot_perievent_spiketimes(merged_clusters(nprobe).spike_times,merged_clusters(nprobe).merged_spike_id,Task_info,Behaviour,[5 1],[-1 3],0.02,...
                'unit_depth',merged_clusters(nprobe).peak_depth(ia),'unit_region',merged_clusters(nprobe).region(ia),'unit_id',C,'event_times',Task_info.start_time_all,...
                'event_label',{'Track 1','Track 2'},'event_id',Task_info.track_ID_all);
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


%% Spatial bayesian decoding

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

        if exist(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            clusters = clusters_ks3;
            sorting_option = 'spikeinterface';
        else
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            sorting_option = 'old';
        end

        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));


        metric_param = create_cluster_selection_params('sorting_option',sorting_option);
        clusters_R= [];
        clusters_L= [];

        HPC_clusters_R = [];
        V1_clusters_R = [];

        HPC_clusters_L = [];
        V1_clusters_L = [];

        for nprobe = 1:length(clusters)

            if clusters(nprobe).probe_hemisphere == 1
                clusters_L = clusters(nprobe);
                V1_channels_L = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by probe');
                HPC_channels_L = determine_region_channels(best_channels{nprobe},options,'region','HPC','group','by probe');

                metric_param.peak_channel = @(x) ismember(x,V1_channels_L);
                [V1_clusters_L] = select_clusters(clusters_L,metric_param);

                metric_param.peak_channel = @(x) ismember(x,HPC_channels_L); % metric_param.depth_range = [] -- full range?
                [HPC_clusters_L] = select_clusters(clusters_L,metric_param);

            elseif clusters(nprobe).probe_hemisphere == 2
                clusters_R = clusters(nprobe);
                V1_channels_R = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by probe');
                HPC_channels_R = determine_region_channels(best_channels{nprobe},options,'region','HPC','group','by probe');

                metric_param.peak_channel = @(x) ismember(x,V1_channels_R); % metric_param.depth_range = [] -- full range?
                [V1_clusters_R] = select_clusters(clusters_R,metric_param);

                metric_param.peak_channel = @(x) ismember(x,HPC_channels_R); % metric_param.depth_range = [] -- full range?
                [HPC_clusters_R] = select_clusters(clusters_R,metric_param);

            end
        end


        % generate spike count for each lap


    end
end

