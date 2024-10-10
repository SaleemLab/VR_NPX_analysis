%% Main place cell and V1 spatial modulation analysis code

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))
%% Spatial raster plot and spatial tuning curves & spatial modulation analysis

clear all
% SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
SUBJECTS={'M24016','M24017','M24018'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
experiment_info=experiment_info([6 9 14 21 22 27 35 38 40]);
Stimulus_type = 'RUN1';

for nsession = 1
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);

        if contains(stimulus_name{n},'Masa2tracks')
            load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
        end

        % reconstruct cluster structure and place field structure from
        % session clusters

        clusters_combined= session_clusters;
        clusters_combined.spike_id=vertcat(session_clusters.spike_id{:});
        clusters_combined.spike_times=vertcat(session_clusters.spike_times{:});
        [clusters_combined.spike_times,index] =sort(clusters_combined.spike_times);
        clusters_combined.spike_id=clusters_combined.spike_id(index);
            
        % Cell with spatial tuning
        ia = find((clusters_combined.peak_percentile(:,1)>0.95&clusters_combined.odd_even_stability(:,1)>0.95) ...
            | (clusters_combined.peak_percentile(:,2)>0.95&clusters_combined.odd_even_stability(:,2)>0.95));

        % ib = find((clusters_combined.peak_percentile(:,1)>0.95&clusters_combined.odd_even_stability(:,1)>0.95) ...
        %     | (clusters_combined.peak_percentile(:,2)>0.95&clusters_combined.odd_even_stability(:,2)>0.95));
        % ia = find(contains(clusters_combined.region,'HPC'))
        %         [C,ia,ic] = unique(clusters_combined.cluster_id);
        C = clusters_combined.cluster_id(ia);
        
        clear place_fields
        x_bin_size =2;
        spatial_response = calculate_raw_spatial_response(clusters_combined.spike_id,clusters_combined.cluster_id,clusters_combined.spike_times,clusters_combined.tvec{1},...
            clusters_combined.position{1},clusters_combined.speed{1},clusters_combined.track_ID_all{1},clusters_combined.start_time_all{1},clusters_combined.end_time_all{1},x_bin_size);
        for track_id = 1:max(session_clusters.track_ID_all{1})
            place_fields(track_id).x_bin_edges = 0:x_bin_size:140;
            place_fields(track_id).x_bin_centres = x_bin_size/2:x_bin_size:140-x_bin_size/2;
            place_fields(track_id).raw = spatial_response(:,track_id);
        end


        Task_info.start_time_all = clusters_combined.start_time_all{1};
        Task_info.end_time_all = clusters_combined.end_time_all{1};
        Task_info.track_ID_all = clusters_combined.track_ID_all{1};

        Behaviour.tvec = clusters_combined.tvec{1};
        Behaviour.position = clusters_combined.position{1};
        Behaviour.speed = clusters_combined.speed{1};
%         place_fields=struct();
%         Behaviour.tvec =Behaviour.sglxTime_uncorrected ;% check if photodiode correction is causing this issue...
        plot_raster_both_track(clusters_combined.spike_times,clusters_combined.spike_id,Task_info,Behaviour,[5 1],[0 140],2,...
            'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C);
%    plot_raster_both_track(clusters_combined.spike_times,clusters_combined.spike_id,Task_info,Behaviour,[5 1],[0 140],2,...
%             'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'place_fields',place_fields);
        if  contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','Spatial PSTH',sprintf(erase(stimulus_name{n},'Masa2tracks_'))))
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','Spatial PSTH',sprintf(erase(stimulus_name{n},'Masa2tracks_'))),[])
        else
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','Spatial PSTH'),[])
        end


        %
        %         % Spatial modulation
        %         x_bin_size = mean(diff(place_fields_V1_L(1).x_bin_centres));
        %         SMI = calculate_spatial_modulation_index(V1_clusters_L,Task_info,Behaviour,[0 140],x_bin_size,'place_fields',place_fields_V1_L,'subplot_xy',[3 1],'plot_option',1)
        %         SMI = calculate_spatial_modulation_index(V1_clusters_R,Task_info,Behaviour,[0 140],x_bin_size,'place_fields',place_fields_V1_R,'subplot_xy',[3 1],'plot_option',1)
        %
        %
        %         metric_param.region = @(x) contains(x,'HPC_L');
        %         [HPC_clusters_L,cluster_id] = select_clusters(clusters(1),metric_param);
        %
        %         SMI = calculate_spatial_modulation_index(HPC_clusters_L,Task_info,Behaviour,[0 140],x_bin_size,'place_fields',place_fields_HPC_L,'subplot_xy',[3 1],'plot_option',1)
        %         SMI = calculate_spatial_modulation_index(HPC_clusters_R,Task_info,Behaviour,[0 140],x_bin_size,'place_fields',place_fields_HPC_R,'subplot_xy',[3 1],'plot_option',1)
        %         %         calculate_spatial_modulation_index(place_fields_V1_L);

        spatial_cell_id = find((clusters_combined.peak_percentile(:,1)>0.95&clusters_combined.odd_even_stability(:,1)>0.95) ...
            | (clusters_combined.peak_percentile(:,2)>0.95&clusters_combined.odd_even_stability(:,2)>0.95));

        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);

            % plot populational map and PV correlation
            %%%%%% HPC
            options.region = 'HPC';
            if options.probe_hemisphere==1
                cluster_id=intersect(spatial_cell_id,find(contains(clusters_combined.region,'HPC_L')));
            else
                cluster_id= intersect(spatial_cell_id,find(contains(clusters_combined.region,'HPC_R')));
            end

            [~,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = ...
                plot_place_cell_map_correlation(place_fields,cluster_id,Task_info,Behaviour,options); % Roughly 6-7 mins for shuffle and plotting

            if options.probe_hemisphere == 1
                if contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
                    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('population_vector_corr_HPC_L%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                else
                    save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_HPC_L.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                end
            elseif options.probe_hemisphere == 2
                if contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
                    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('population_vector_corr_HPC_R%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                else
                    save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_HPC_R.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                end
            end

            % plot populational map and PV correlation
            %%%%%% V1
            options.region = 'V1';

            if options.probe_hemisphere==1
                cluster_id=intersect(spatial_cell_id,find(contains(clusters_combined.region,'V1_L')));
            else
                cluster_id= intersect(spatial_cell_id,find(contains(clusters_combined.region,'V1_R')));
            end
            [~,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = ...
                plot_place_cell_map_correlation(place_fields,cluster_id,Task_info,Behaviour,options); % Roughly 6-7 mins for shuffle and plotting

            if options.probe_hemisphere == 1
                if contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
                    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('population_vector_corr_V1_L%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                else
                    save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_V1_L.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                end
            elseif options.probe_hemisphere == 2
                if contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
                    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('population_vector_corr_V1_R%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                else
                    save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_V1_R.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                end
            end

        end

        if length(session_info(n).probe) > 1
            % plot populational map and PV correlation
            %%%%%% HPC
            options.probe_combined = 1;

            options.region = 'HPC';
            cluster_id= intersect(spatial_cell_id,find(contains(clusters_combined.region,'HPC')))

            [~,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = ...
                plot_place_cell_map_correlation(place_fields,cluster_id,Task_info,Behaviour,options); % Roughly 6-7 mins for shuffle and plotting

            if contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('population_vector_corr_HPC_combined%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
            else
                save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_HPC_combined.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
            end

            % plot populational map and PV correlation
            %%%%%% V1
            options.region = 'V1';
            cluster_id= intersect(spatial_cell_id,find(contains(clusters_combined.region,'V1')))

            [~,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = ...
                plot_place_cell_map_correlation(place_fields,cluster_id,Task_info,Behaviour,options); % Roughly 6-7 mins for shuffle and plotting

            if contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('population_vector_corr_V1_combined%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
            else
                save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_V1_combined.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
            end
            options = rmfield(options,'probe_combined');
        end

        if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','populational_map'))== 0
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','populational_map'))
        end
        
        if  contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','populational_map',sprintf(erase(stimulus_name{n},'Masa2tracks_'))))
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','populational_map',sprintf(erase(stimulus_name{n},'Masa2tracks_'))),[])
        else
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','populational_map'),[])
        end

        close all

        %         %GLM analysis
        %         spatial_modulation_GLM_analysis(V1_clusters_L,[],Behaviour,Task_info);
        %
        %         spatial_modulation_GLM_analysis(V1_clusters_L,place_fields_V1_L(1).cluster_id(place_fields_V1_L(1).all_good_cells_LIBERAL),Behaviour,Task_info);
    end
end


%% Spatial bayesian decoding
%% Within session
clear all
% SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
SUBJECTS={'M24016','M24017','M24018'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
experiment_info=experiment_info([6 9 14 21 22 27 35 38 40]);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]

for nsession = 1
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);

        if contains(stimulus_name{n},'Masa2tracks')
            load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        end

        clusters_combined= session_clusters;
        clusters_combined.spike_id=vertcat(session_clusters.spike_id{:});
        clusters_combined.spike_times=vertcat(session_clusters.spike_times{:});
        [clusters_combined.spike_times,index] =sort(clusters_combined.spike_times);
        clusters_combined.spike_id=clusters_combined.spike_id(index);

        % Cell with spatial tuning
        ia = find((clusters_combined.peak_percentile(:,1)>0.95&clusters_combined.odd_even_stability(:,1)>0.95) ...
            | (clusters_combined.peak_percentile(:,2)>0.95&clusters_combined.odd_even_stability(:,2)>0.95));
        %         [C,ia,ic] = unique(clusters_combined.cluster_id);
        C = clusters_combined.cluster_id(ia);

        clear place_fields
        x_bin_size =5;
        spatial_response = calculate_raw_spatial_response(clusters_combined.spike_id,clusters_combined.cluster_id,clusters_combined.spike_times,clusters_combined.tvec{1},...
            clusters_combined.position{1},clusters_combined.speed{1},clusters_combined.track_ID_all{1},clusters_combined.start_time_all{1},clusters_combined.end_time_all{1},x_bin_size);
        for track_id = 1:max(session_clusters.track_ID_all{1})
            place_fields(track_id).x_bin_edges = 0:x_bin_size:140;
            place_fields(track_id).x_bin_centres = x_bin_size/2:x_bin_size:140-x_bin_size/2;
            place_fields(track_id).raw = spatial_response(:,track_id);
        end

        probability_ratio_RUN_lap_HPC_combined= [];
        probability_ratio_RUN_lap_V1_combined= [];

        estimated_position_lap_CV_shuffled_HPC_combined.track = [];
        estimated_position_lap_CV_shuffled_V1_combined.track = [];
        estimated_position_lap_CV_HPC_combined.track = [];
        estimated_position_lap_CV_V1_combined.track = [];

        for nprobe = 1:length(clusters)
            probability_ratio_RUN_lap_V1{nprobe}= [];
            probability_ratio_RUN_lap_HPC{nprobe}= [];
            estimated_position_lap_CV_HPC(nprobe).track = [];
            estimated_position_lap_CV_V1(nprobe).track = [];
            estimated_position_lap_CV_shuffled_HPC(nprobe).track = [];
            estimated_position_lap_CV_shuffled_V1(nprobe).track = [];

            options = session_info(n).probe(nprobe);

            % Bayesian decoding 10 fold cross validated
            DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'RUN Bayesian decoding'));
            if isempty(DIR)
                mkdir(fullfile(options.ANALYSIS_DATAPATH,'RUN Bayesian decoding'))
            end
            
            metric_param = [];
            if options.probe_hemisphere == 1
                metric_param.region = @(x) contains(x,'HPC_L');
                options.region = 'HPC Left';
            elseif options.probe_hemisphere == 2
                metric_param.region = @(x) contains(x,'HPC_R');
                options.region = 'HPC Right';
            end

            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            %             [probability_ratio_RUN_lap_HPC{nprobe},estimated_position_lap_CV_HPC(nprobe).track] = bayesian_decoding_RUN_lap_cross_validation(selected_clusters,place_fields,Behaviour,Task_info,options);
            [probability_ratio_RUN_lap_HPC{nprobe},estimated_position_lap_CV_HPC(nprobe).track,estimated_position_lap_CV_shuffled_HPC(nprobe).track] =...
                bayesian_decoding_RUN_lap_cross_validation_all(selected_clusters,place_fields,Behaviour,Task_info,options);


            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.unstable_ids = @(x) x==0;
            if clusters(nprobe).probe_hemisphere == 1
                metric_param.region = @(x) contains(x,'V1_L');
                options.region = 'V1 Left';
            elseif clusters(nprobe).probe_hemisphere == 2
                metric_param.region = @(x) contains(x,'V1_R');
                options.region = 'V1 Right';
            end

            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            [probability_ratio_RUN_lap_V1{nprobe},estimated_position_lap_CV_V1(nprobe).track,estimated_position_lap_CV_shuffled_V1(nprobe).track] = bayesian_decoding_RUN_lap_cross_validation_all(selected_clusters,place_fields,Behaviour,Task_info,options);

        end

        if length(session_info(n).probe) > 1
            options.region = 'HPC Combined';
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.unstable_ids = @(x) x==0;
            metric_param.region = @(x) contains(x,'HPC');

            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            [probability_ratio_RUN_lap_HPC_combined,estimated_position_lap_CV_HPC_combined.track,estimated_position_lap_CV_shuffled_HPC_combined.track] = bayesian_decoding_RUN_lap_cross_validation_all(selected_clusters,place_fields,Behaviour,Task_info,options);

            options.region = 'V1 Combined';
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.unstable_ids = @(x) x==0;
            metric_param.region = @(x) contains(x,'V1');

            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            [probability_ratio_RUN_lap_V1_combined,estimated_position_lap_CV_V1_combined.track,estimated_position_lap_CV_shuffled_V1_combined.track] = bayesian_decoding_RUN_lap_cross_validation_all(selected_clusters,place_fields,Behaviour,Task_info,options);

        end

        save(fullfile(options.ANALYSIS_DATAPATH,'estimated_position_lap_CV_HPC.mat'),'estimated_position_lap_CV_HPC_combined','estimated_position_lap_CV_HPC','estimated_position_lap_CV_shuffled_HPC','estimated_position_lap_CV_shuffled_HPC_combined')
        save(fullfile(options.ANALYSIS_DATAPATH,'probability_ratio_RUN_lap_HPC.mat'),'probability_ratio_RUN_lap_HPC_combined','probability_ratio_RUN_lap_HPC')

        save(fullfile(options.ANALYSIS_DATAPATH,'estimated_position_lap_CV_V1.mat'),'estimated_position_lap_CV_V1',"estimated_position_lap_CV_V1_combined",'estimated_position_lap_CV_shuffled_V1','estimated_position_lap_CV_shuffled_V1_combined')
        save(fullfile(options.ANALYSIS_DATAPATH,'probability_ratio_RUN_lap_V1.mat'),'probability_ratio_RUN_lap_V1','probability_ratio_RUN_lap_V1_combined')

    end
end


%% Decoding trajectory plotting
lap_times = [];
clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]

for nsession = [1 2 3 4 6 7 8 9 10 12 14]
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

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


        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_place_fields.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end

        load(fullfile(options.ANALYSIS_DATAPATH,'estimated_position_lap_CV_HPC.mat'))
        load(fullfile(options.ANALYSIS_DATAPATH,'probability_ratio_RUN_lap_HPC.mat'))

        load(fullfile(options.ANALYSIS_DATAPATH,'estimated_position_lap_CV_V1.mat'))
        load(fullfile(options.ANALYSIS_DATAPATH,'probability_ratio_RUN_lap_V1.mat'))
        

        % plotting decoded run trajectory
        for nprobe = 1:length(clusters)
            options = session_info(n).probe(nprobe);

            if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'))== 0
                mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'))
            end

            if clusters(nprobe).probe_hemisphere == 1
                options.region = 'HPC Left';
            elseif clusters(nprobe).probe_hemisphere == 2
                options.region = 'HPC Right';
            end

            plot_decoding_RUN_trajectory(estimated_position_lap_CV_HPC(nprobe).track,Task_info,options)
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])

            if clusters(nprobe).probe_hemisphere == 1
                metric_param.region = @(x) contains(x,'V1_L');
                options.region = 'V1 Left';
            elseif clusters(nprobe).probe_hemisphere == 2
                metric_param.region = @(x) contains(x,'V1_R');
                options.region = 'V1 Right';
            end

            plot_decoding_RUN_trajectory(estimated_position_lap_CV_V1(nprobe).track,Task_info,options)
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])
        end

        if length(session_info(n).probe) > 1
            options.region = 'HPC Combined';
            plot_decoding_RUN_trajectory(estimated_position_lap_CV_HPC_combined.track,Task_info,options)
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])

            options.region = 'V1 Combined';
            plot_decoding_RUN_trajectory(estimated_position_lap_CV_V1_combined.track,Task_info,options)
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])
        end
    end
end


%% Decoding error and log odds for V1 and HPC
lap_times = [];
clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]

for nsession = [1 2 3 4 6 7 8 9 10 12 14]
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

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


        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_place_fields.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end

        load(fullfile(options.ANALYSIS_DATAPATH,'estimated_position_lap_CV_HPC.mat'))
        load(fullfile(options.ANALYSIS_DATAPATH,'probability_ratio_RUN_lap_HPC.mat'))

        load(fullfile(options.ANALYSIS_DATAPATH,'estimated_position_lap_CV_V1.mat'))
        load(fullfile(options.ANALYSIS_DATAPATH,'probability_ratio_RUN_lap_V1.mat'))
        

        % plotting decoded run trajectory
        for nprobe = 1:length(clusters)
            options = session_info(n).probe(nprobe);

            if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'))== 0
                mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'))
            end

            if clusters(nprobe).probe_hemisphere == 1
                options.region = 'HPC Left';
            elseif clusters(nprobe).probe_hemisphere == 2
                options.region = 'HPC Right';
            end

            plot_decoding_RUN_trajectory(estimated_position_lap_CV_HPC(nprobe).track,Task_info,options)
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])

            if clusters(nprobe).probe_hemisphere == 1
                metric_param.region = @(x) contains(x,'V1_L');
                options.region = 'V1 Left';
            elseif clusters(nprobe).probe_hemisphere == 2
                metric_param.region = @(x) contains(x,'V1_R');
                options.region = 'V1 Right';
            end

            plot_decoding_RUN_trajectory(estimated_position_lap_CV_V1(nprobe).track,Task_info,options)
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])
        end

        if length(session_info(n).probe) > 1
            options.region = 'HPC Combined';
            plot_decoding_RUN_trajectory(estimated_position_lap_CV_HPC_combined.track,Task_info,options)
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])

            options.region = 'V1 Combined';
            plot_decoding_RUN_trajectory(estimated_position_lap_CV_V1_combined.track,Task_info,options)
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])
        end


        % initialise all variables
        for track_id = 1:length(place_fields)
            for temp_track = 1:length(place_fields)
                actual_position{nsession}{track_id} = [];
                actual_speed{nsession}{track_id} = [];
                VR_speed{nsession}{track_id} = [];
                decoded_position_lap_id{nsession}{track_id}  = [];


                decoded_position_HPC_combined{nsession}{track_id} = [];
                decoded_error_HPC_combined{nsession}{track_id}{temp_track}  = [];

                decoded_position_HPC_combined_shuffled{nsession}{track_id} = [];
                decoded_error_HPC_combined_shuffled{nsession}{track_id}{temp_track}  = [];

                decoded_position_V1_combined{nsession}{track_id} = [];
                decoded_error_V1_combined{nsession}{track_id}{temp_track}  = [];

                for nprobe = 1:length(session_info(n).probe)
                    probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere;
                    decoded_position_V1{probe_hemisphere}{nsession}{track_id} = [];
                    decoded_error_V1{probe_hemisphere}{nsession}{track_id}{temp_track} = [];
                    decoded_position_HPC{probe_hemisphere}{nsession}{track_id} = [];
                    decoded_error_HPC{probe_hemisphere}{nsession}{track_id}{temp_track} = [];

                    decoded_position_V1_shuffled{probe_hemisphere}{nsession}{track_id} = [];
                    decoded_error_V1_shuffled{probe_hemisphere}{nsession}{track_id}{temp_track} = [];
                    decoded_position_HPC_shuffled{probe_hemisphere}{nsession}{track_id} = [];
                    decoded_error_HPC_shuffled{probe_hemisphere}{nsession}{track_id}{temp_track} = [];
                end
            end
        end

        % Position bins from both tracks
        position_bin_across_tracks = estimated_position_lap_CV_V1(1).track(1).lap(1).track(1).position_bin_centres;
        position_bin_across_tracks = [position_bin_across_tracks position_bin_across_tracks+1000];

        for track_id = 1:length(place_fields)
            for temp_track = 1:length(place_fields)
                for lap_id = 1:length(estimated_position_lap_CV_V1(nprobe).track(track_id).lap)

                    if temp_track == 1 %Do not need to loop twice
                        decoded_position_lap_id{nsession}{track_id} = [decoded_position_lap_id{nsession}{track_id} ...
                            lap_id*ones(1,length(estimated_position_lap_CV_V1(1).track(track_id).lap(lap_id).track(track_id).run_actual_position))];
                        actual_position{nsession}{track_id} = [actual_position{nsession}{track_id} ...
                            estimated_position_lap_CV_V1(1).track(track_id).lap(lap_id).track(track_id).run_actual_position];
                        actual_speed{nsession}{track_id} = [actual_speed{nsession}{track_id} ...
                            estimated_position_lap_CV_V1(1).track(track_id).lap(lap_id).track(track_id).actual_run_speed];
                        VR_speed{nsession}{track_id} = [VR_speed{nsession}{track_id} ...
                            estimated_position_lap_CV_V1(1).track(track_id).lap(lap_id).track(track_id).run_speed];
                    end

%                     if length(session_info(n).probe)==1
%                         decoded_error_HPC{nsession}{track_id}{temp_track} = [decoded_error_HPC{nsession}{track_id}{temp_track} ...
%                             estimated_position_lap_CV_HPC.track(track_id).lap(lap_id).track(temp_track).peak_position...
%                             - estimated_position_lap_CV_HPC.track(track_id).lap(lap_id).track(track_id).run_actual_position];
%                         if temp_track == 1 %Do not need to loop twice
%                             [~,index] = max([estimated_position_lap_CV_HPC.track(track_id).lap(lap_id).track(1).run; ...
%                                 estimated_position_lap_CV_HPC.track(track_id).lap(lap_id).track(2).run]);
%                             decoded_position_HPC{nsession}{track_id} = [decoded_position_HPC{nsession}{track_id} index];
%                         end
%                     end

                    for nprobe = 1:length(session_info(n).probe)
                        probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere;
%                         decoded_position_V1_shuffled{probe_hemisphere}{nsession}{track_id} = [];
%                         decoded_error_V1_shuffled{probe_hemisphere}{nsession}{track_id}{temp_track} = [];
%                         decoded_position_HPC_shuffled{probe_hemisphere}{nsession}{track_id} = [];
%                         decoded_error_HPC_shuffled{probe_hemisphere}{nsession}{track_id}{temp_track} = [];

                        decoded_error_V1{probe_hemisphere}{nsession}{track_id}{temp_track} = [decoded_error_V1{probe_hemisphere}{nsession}{track_id}{temp_track} ...
                            estimated_position_lap_CV_V1(nprobe).track(track_id).lap(lap_id).track(temp_track).peak_position...
                            - estimated_position_lap_CV_V1(nprobe).track(track_id).lap(lap_id).track(track_id).run_actual_position];

                        decoded_error_HPC{probe_hemisphere}{nsession}{track_id}{temp_track} = [decoded_error_HPC{probe_hemisphere}{nsession}{track_id}{temp_track} ...
                            estimated_position_lap_CV_HPC(nprobe).track(track_id).lap(lap_id).track(temp_track).peak_position...
                            - estimated_position_lap_CV_HPC(nprobe).track(track_id).lap(lap_id).track(track_id).run_actual_position];

                        for nshuffle = 1:100
                            [~,index] = max(estimated_position_lap_CV_shuffled_V1(nprobe).track(track_id).lap(lap_id).track(temp_track).run);
                            decoded_error_V1_shuffled{probe_hemisphere}{nsession}{track_id}{temp_track} = [decoded_error_V1_shuffled{probe_hemisphere}{nsession}{track_id}{temp_track} ...
                                estimated_position_lap_CV_shuffled_V1(1).track(1).lap(1).track(1).position_bin_centres(index)...
                                - estimated_position_lap_CV_shuffled_V1(nprobe).track(track_id).lap(lap_id).track(track_id).run_actual_position];

                            [~,index] = max(estimated_position_lap_CV_shuffled_HPC(nprobe).track(track_id).lap(lap_id).track(temp_track).run);
                            decoded_error_HPC_shuffled{probe_hemisphere}{nsession}{track_id}{temp_track} = [decoded_error_HPC_shuffled{probe_hemisphere}{nsession}{track_id}{temp_track} ...
                                estimated_position_lap_CV_shuffled_HPC(1).track(1).lap(1).track(1).position_bin_centres(index)...
                                - estimated_position_lap_CV_shuffled_HPC(nprobe).track(track_id).lap(lap_id).track(track_id).run_actual_position];
                        end
                        


                        if temp_track == 1 %Do not need to loop twice
                            [~,index] = max([estimated_position_lap_CV_V1(nprobe).track(track_id).lap(lap_id).track(1).run;...
                                estimated_position_lap_CV_V1(nprobe).track(track_id).lap(lap_id).track(2).run]);
                            decoded_position_V1{probe_hemisphere}{nsession}{track_id} = [decoded_position_V1{probe_hemisphere}{nsession}{track_id} index];

                            [~,index] = max([estimated_position_lap_CV_HPC(nprobe).track(track_id).lap(lap_id).track(1).run; ...
                                estimated_position_lap_CV_HPC(nprobe).track(track_id).lap(lap_id).track(2).run]);
                            decoded_position_HPC{probe_hemisphere}{nsession}{track_id} = [decoded_position_HPC{probe_hemisphere}{nsession}{track_id} index];

                            for nshuffle = 1:100
                                [~,index] = max([estimated_position_lap_CV_shuffled_V1(nprobe).track{nshuffle}(track_id).lap(lap_id).track(1).run;...
                                    estimated_position_lap_CV_shuffled_V1(nprobe).track(track_id).lap(lap_id).track(2).run]);
                                decoded_position_V1_shuffled{probe_hemisphere}{nsession}{track_id} = [decoded_position_V1_shuffled{probe_hemisphere}{nsession}{track_id} index];

                                [~,index] = max([estimated_position_lap_CV_shuffled_HPC(nprobe).track(track_id).lap(lap_id).track(1).run; ...
                                    estimated_position_lap_CV_shuffled_HPC(nprobe).track(track_id).lap(lap_id).track(2).run]);
                                decoded_position_HPC_shuffled{probe_hemisphere}{nsession}{track_id} = [decoded_position_HPC_shuffled{probe_hemisphere}{nsession}{track_id} index];
                            end

                        end
                    end

                    if length(session_info(n).probe)>1
                        decoded_error_HPC_combined{nsession}{track_id}{temp_track} = [decoded_error_HPC_combined{nsession}{track_id}{temp_track} ...
                            estimated_position_lap_CV_HPC_combined.track(track_id).lap(lap_id).track(temp_track).peak_position...
                            - estimated_position_lap_CV_HPC_combined.track(track_id).lap(lap_id).track(track_id).run_actual_position];

                        decoded_error_V1_combined{nsession}{track_id}{temp_track} = [decoded_error_V1_combined{nsession}{track_id}{temp_track} ...
                            estimated_position_lap_CV_V1_combined.track(track_id).lap(lap_id).track(temp_track).peak_position...
                            - estimated_position_lap_CV_V1_combined.track(track_id).lap(lap_id).track(track_id).run_actual_position];

                        if temp_track == 1 %Do not need to loop twice
                            [~,index] = max([estimated_position_lap_CV_HPC_combined.track(track_id).lap(lap_id).track(1).run;...
                                estimated_position_lap_CV_HPC_combined.track(track_id).lap(lap_id).track(2).run]);
                            decoded_position_HPC_combined{nsession}{track_id} = [decoded_position_HPC_combined{nsession}{track_id} index];

                            [~,index] = max([estimated_position_lap_CV_V1_combined.track(track_id).lap(lap_id).track(1).run;...
                                estimated_position_lap_CV_V1_combined.track(track_id).lap(lap_id).track(2).run]);
                            decoded_position_V1_combined{nsession}{track_id} = [decoded_position_V1_combined{nsession}{track_id} index];
                        end

                    end
                end
            end
        end


        % Plotting decoding performance
        if length(session_info(n).probe)>1
            decoding_performance = plot_within_session_decoded_error_HPC_V1(decoded_position_lap_id,VR_speed,actual_position,estimated_position_lap_CV_V1,decoded_position_V1,decoded_position_HPC,decoded_position_HPC_combined,decoded_position_V1_combined,...
                decoded_error_V1,decoded_error_HPC,decoded_error_HPC_combined,decoded_error_V1_combined,place_fields,session_info(n),nsession)
        else
            decoding_performance = plot_within_session_decoded_error_HPC_V1(VR_speed,actual_position,estimated_position_lap_CV_V1,decoded_position_V1,decoded_position_HPC,[],[],...
                decoded_error_V1,decoded_error_HPC,[],[],place_fields,session_info(n),nsession)
        end
       
        save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells'),[])

    end
end




    decoded_position_V1{probe_hemisphere}{nsession}{track_id};

    scatter(decoded_position_HPC{nsession}{2}(VR_speed{nsession}{2}>5)...
        ,decoded_position_V1{1}{nsession}{2}(VR_speed{nsession}{2}>5),'blue','filled','MarkerFaceAlpha',0.05)

    [N,Xedges,Yedges,binX,binY] = histcounts2(decoded_position_HPC{nsession}{2}(VR_speed{nsession}{2}>5),decoded_position_V1{1}{nsession}{2}(VR_speed{nsession}{2}>5),14);
    imagesc(flip(N)/max(max(N)))
   
