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
experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
Stimulus_type = 'RUN1';

for nsession = 1:length(experiment_info)
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
%         ia = find((clusters_combined.peak_percentile(:,1)>0.95&clusters_combined.odd_even_stability(:,1)>0.95) ...
%             | (clusters_combined.peak_percentile(:,2)>0.95&clusters_combined.odd_even_stability(:,2)>0.95));
        ia = find(clusters_combined.odd_even_stability(:,1)>0.95 ...
            | clusters_combined.odd_even_stability(:,2)>0.95);
        % ib = find((clusters_combined.peak_percentile(:,1)>0.95&clusters_combined.odd_even_stability(:,1)>0.95) ...
        %     | (clusters_combined.peak_percentile(:,2)>0.95&clusters_combined.odd_even_stability(:,2)>0.95));
        % ia = find(contains(clusters_combined.region,'HPC'))
        %         [C,ia,ic] = unique(clusters_combined.cluster_id);
        C = clusters_combined.cluster_id(ia);
        
        clear place_fields
        speed = clusters_combined.speed{1};
        speed(isnan(speed))=0;
        w = gausswin(9);
        w = w / sum(w);
        speed = filtfilt(w,1,speed')';
        x_window = [0 140];
        x_bin_width =2;
        place_fields = calculate_spatial_cells(clusters_combined,clusters_combined.tvec{1},...
            clusters_combined.position{1},speed,clusters_combined.track_ID_all{1},clusters_combined.start_time_all{1},clusters_combined.end_time_all{1},x_window,x_bin_width);


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
experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]

for nsession = 1:length(experiment_info)
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
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            clusters = clusters_ks3;
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        end

        clusters_combined= session_clusters;
        clusters_combined.spike_id=vertcat(session_clusters.spike_id{:});
        clusters_combined.spike_times=vertcat(session_clusters.spike_times{:});
        [clusters_combined.spike_times,index] =sort(clusters_combined.spike_times);
        clusters_combined.spike_id=clusters_combined.spike_id(index);

        % Cell with spatial tuning
        ia = find(clusters_combined.odd_even_stability(:,1)>0.95 ...
            | clusters_combined.odd_even_stability(:,2)>0.95);
%                 ia = find((clusters_combined.peak_percentile(:,1)>0.95&clusters_combined.odd_even_stability(:,1)>0.95) ...
%             | (clusters_combined.peak_percentile(:,2)>0.95&clusters_combined.odd_even_stability(:,2)>0.95));
        %         [C,ia,ic] = unique(clusters_combined.cluster_id);
        C = clusters_combined.cluster_id(ia);
        
        clear place_fields_BAYESIAN
        

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

            speed = clusters_combined.speed{1};
            speed(isnan(speed))=0;
            w = gausswin(9);
            w = w / sum(w);
            speed = filtfilt(w,1,speed')';
            x_window = [0 140];
            x_bin_width =5;
            place_fields_BAYESIAN = calculate_spatial_cells(selected_clusters,selected_clusters.tvec{1},...
                selected_clusters.position{1},speed,selected_clusters.track_ID_all{1},selected_clusters.start_time_all{1},selected_clusters.end_time_all{1},x_window,x_bin_width);

            %             [probability_ratio_RUN_lap_HPC{nprobe},estimated_position_lap_CV_HPC(nprobe).track] = bayesian_decoding_RUN_lap_cross_validation(selected_clusters,place_fields,Behaviour,Task_info,options);
            [probability_ratio_RUN_lap_HPC{nprobe},estimated_position_lap_CV_HPC(nprobe).track,estimated_position_lap_CV_shuffled_HPC(nprobe).track] =...
                bayesian_decoding_RUN_lap_cross_validation_all(selected_clusters,place_fields_BAYESIAN,Behaviour,Task_info,options);
            clear place_fields_BAYESIAN

            metric_param = [];
            if clusters(nprobe).probe_hemisphere == 1
                metric_param.region = @(x) contains(x,'V1_L');
                options.region = 'V1 Left';
            elseif clusters(nprobe).probe_hemisphere == 2
                metric_param.region = @(x) contains(x,'V1_R');
                options.region = 'V1 Right';
            end

            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            place_fields_BAYESIAN = calculate_spatial_cells(selected_clusters,selected_clusters.tvec{1},...
                selected_clusters.position{1},speed,selected_clusters.track_ID_all{1},selected_clusters.start_time_all{1},selected_clusters.end_time_all{1},x_window,x_bin_width);

            [probability_ratio_RUN_lap_V1{nprobe},estimated_position_lap_CV_V1(nprobe).track,estimated_position_lap_CV_shuffled_V1(nprobe).track] =...
                bayesian_decoding_RUN_lap_cross_validation_all(selected_clusters,place_fields_BAYESIAN,Behaviour,Task_info,options);
            clear place_fields_BAYESIAN
        end

        if length(session_info(n).probe) > 1
            options.region = 'HPC Combined';
            metric_param=[]
            metric_param.region = @(x) contains(x,'HPC');

            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            place_fields_BAYESIAN = calculate_spatial_cells(selected_clusters,selected_clusters.tvec{1},...
                selected_clusters.position{1},speed,selected_clusters.track_ID_all{1},selected_clusters.start_time_all{1},selected_clusters.end_time_all{1},x_window,x_bin_width);
            [probability_ratio_RUN_lap_HPC_combined,estimated_position_lap_CV_HPC_combined.track,estimated_position_lap_CV_shuffled_HPC_combined.track] = ...
                bayesian_decoding_RUN_lap_cross_validation_all(selected_clusters,place_fields_BAYESIAN,Behaviour,Task_info,options);
            clear place_fields_BAYESIAN

            options.region = 'V1 Combined';
            metric_param=[];
            metric_param.region = @(x) contains(x,'V1');

            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            place_fields_BAYESIAN = calculate_spatial_cells(selected_clusters,selected_clusters.tvec{1},...
                selected_clusters.position{1},speed,selected_clusters.track_ID_all{1},selected_clusters.start_time_all{1},selected_clusters.end_time_all{1},x_window,x_bin_width);
            [probability_ratio_RUN_lap_V1_combined,estimated_position_lap_CV_V1_combined.track,estimated_position_lap_CV_shuffled_V1_combined.track] =...
                bayesian_decoding_RUN_lap_cross_validation_all(selected_clusters,place_fields_BAYESIAN,Behaviour,Task_info,options);
            clear place_fields_BAYESIAN

        end

        save(fullfile(options.ANALYSIS_DATAPATH,'estimated_position_lap_CV_HPC.mat'),'estimated_position_lap_CV_HPC_combined','estimated_position_lap_CV_HPC','estimated_position_lap_CV_shuffled_HPC','estimated_position_lap_CV_shuffled_HPC_combined')
        save(fullfile(options.ANALYSIS_DATAPATH,'probability_ratio_RUN_lap_HPC.mat'),'probability_ratio_RUN_lap_HPC_combined','probability_ratio_RUN_lap_HPC')

        save(fullfile(options.ANALYSIS_DATAPATH,'estimated_position_lap_CV_V1.mat'),'estimated_position_lap_CV_V1',"estimated_position_lap_CV_V1_combined",'estimated_position_lap_CV_shuffled_V1','estimated_position_lap_CV_shuffled_V1_combined')
        save(fullfile(options.ANALYSIS_DATAPATH,'probability_ratio_RUN_lap_V1.mat'),'probability_ratio_RUN_lap_V1','probability_ratio_RUN_lap_V1_combined')

    end
end
%% Log odds


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

        if contains(stimulus_name{n},'Masa2tracks')
            load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            clusters = clusters_ks3;
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        end

        clusters_combined= session_clusters;
        clusters_combined.spike_id=vertcat(session_clusters.spike_id{:});
        clusters_combined.spike_times=vertcat(session_clusters.spike_times{:});
        [clusters_combined.spike_times,index] =sort(clusters_combined.spike_times);
        clusters_combined.spike_id=clusters_combined.spike_id(index);

        % Cell with spatial tuning
        ia = find(clusters_combined.odd_even_stability(:,1)>0.95 ...
            | clusters_combined.odd_even_stability(:,2)>0.95);
%                 ia = find((clusters_combined.peak_percentile(:,1)>0.95&clusters_combined.odd_even_stability(:,1)>0.95) ...
%             | (clusters_combined.peak_percentile(:,2)>0.95&clusters_combined.odd_even_stability(:,2)>0.95));
        %         [C,ia,ic] = unique(clusters_combined.cluster_id);
        C = clusters_combined.cluster_id(ia);

        speed = clusters_combined.speed{1};
        speed(isnan(speed))=0;
        w = gausswin(9);
        w = w / sum(w);
        speed = filtfilt(w,1,speed')';
        x_window = [0 140];
        x_bin_width =5;
        place_fields_BAYESIAN = calculate_spatial_cells(clusters_combined,clusters_combined.tvec{1},...
            clusters_combined.position{1},speed,clusters_combined.track_ID_all{1},clusters_combined.start_time_all{1},clusters_combined.end_time_all{1},x_window,x_bin_width);


        load(fullfile(options.ANALYSIS_DATAPATH,'estimated_position_lap_CV_HPC.mat'))
        load(fullfile(options.ANALYSIS_DATAPATH,'probability_ratio_RUN_lap_HPC.mat'))

        load(fullfile(options.ANALYSIS_DATAPATH,'estimated_position_lap_CV_V1.mat'))
        load(fullfile(options.ANALYSIS_DATAPATH,'probability_ratio_RUN_lap_V1.mat'))

        % Plotting Log odds

        if length(session_info(n).probe)>1
            
            estimated_position_lap_CV_HPC = estimated_position_lap_CV_HPC_combined;
        else
            load estimated_position_lap_CV_HPC
        end

        

        z_log_odds = [];
        track_label = [];
        HPC_bayesian_bias = [];
        V1_z_log_odds= [];
        V1_track_label= [];
        V1_bayesian_bias = [];
        T1_T2_ratio_shuffled=[];

        % RUN lap log odds (for track 2 use 1/track 2 bias )
        % Because decoding now speed thresholded (nan for wrong track)
        probability_ratio_RUN_lap=probability_ratio_RUN_lap_HPC_combined;
        for track_id = 1:length(probability_ratio_RUN_lap{1})
            for nlap = 1:length(probability_ratio_RUN_lap{1}{track_id}{1})
                data = log(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap));

                if track_id == 1
                    data = log(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap));
                else
                    data = log(1/probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap));
                end

                for nshuffle = 1:1000
                    if track_id == 1
                        T1_T2_ratio_shuffled(nshuffle) = cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{1}(nlap));
                    else
                        T1_T2_ratio_shuffled(nshuffle) = 1/cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{1}(nlap));
                    end
                end

                shuffled_data = log(T1_T2_ratio_shuffled);
                z_log_odds{track_id}(nlap) = (data-mean(shuffled_data))/std(shuffled_data);
                track_label{track_id}(nlap) = track_id;

                HPC_bayesian_bias{track_id}(nlap) = nansum(estimated_position_lap_CV_HPC.track(track_id).lap(nlap).track(1).run_bias)/...
                    (nansum(estimated_position_lap_CV_HPC.track(track_id).lap(nlap).track(1).run_bias)+nansum(estimated_position_lap_CV_HPC.track(track_id).lap(nlap).track(2).run_bias)) ;

            end
        end


        % probability_ratio_RUN_lap{1}{track_id}
        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);
            options.ROOTPATH = ROOTPATH;
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere ;
            options.probe_no = probe_no;

            probability_ratio_RUN_lap = probability_ratio_RUN_lap_V1{probe_no};

            V1_z_log_odds{probe_hemisphere} = [];
            V1_track_label{probe_hemisphere} = [];
            V1_bayesian_bias{probe_hemisphere} = [];

            for track_id = 1:length(probability_ratio_RUN_lap{1})
                for nlap = 1:length(probability_ratio_RUN_lap{1}{track_id}{1})
                    data = log(cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));

                    if track_id == 1
                        data = log(cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
                    else
                        data = log(1/cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
                    end

                    for nshuffle = 1:1000
                        if track_id == 1
                            T1_T2_ratio_shuffled(nshuffle) = cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{track_id}(nlap));
                        else
                            T1_T2_ratio_shuffled(nshuffle) = 1/cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{track_id}(nlap));
                        end
                    end

                    shuffled_data = log(T1_T2_ratio_shuffled);
                    V1_z_log_odds{probe_hemisphere}{track_id}(nlap) = (data-mean(shuffled_data))/std(shuffled_data);
                    V1_track_label{probe_hemisphere}{track_id}(nlap) = track_id;


                    V1_bayesian_bias{probe_hemisphere}{track_id}(nlap) = nansum(estimated_position_lap_CV_V1(probe_no).track(track_id).lap(nlap).track(1).run_bias)/...
                        (nansum(estimated_position_lap_CV_V1(probe_no).track(track_id).lap(nlap).track(1).run_bias)+nansum(estimated_position_lap_CV_V1(probe_no).track(track_id).lap(nlap).track(2).run_bias)) ;

                end
            end

        end


        z_log_odds = [z_log_odds{1}, z_log_odds{2}];
        HPC_bayesian_bias = [HPC_bayesian_bias{1} HPC_bayesian_bias{2}];
        if ~isempty(V1_z_log_odds{1})
            V1_z_log_odds{1} = [V1_z_log_odds{1}{1}, V1_z_log_odds{1}{2}];
            V1_bayesian_bias{1} = [V1_bayesian_bias{1}{1}, V1_bayesian_bias{1}{2}];
        end

        if ~isempty(V1_z_log_odds{2})
            V1_z_log_odds{2} = [V1_z_log_odds{2}{1}, V1_z_log_odds{2}{2}];
            V1_bayesian_bias{2} = [V1_bayesian_bias{2}{1}, V1_bayesian_bias{2}{2}];
        end





        %         subplot(2,3,2)
        %         bar(V1_bayesian_bias{2}(sorted_id),'b','FaceAlpha',0.3,'EdgeColor','none');hold on;
        %
        %         bar(HPC_bayesian_bias(sorted_id),'k','FaceAlpha',0.3,'EdgeColor','none');hold on;
        %         for nlap = 1:length(track_orders)
        %             if track_orders(nlap) == 1
        %                 scatter(nlap,track_orders(nlap) ,'r')
        %             elseif track_orders(nlap) == 2
        %                 scatter(nlap,track_orders(nlap) -2.2,'b')
        %             end
        %         end

        [~,sorted_id] = sort([lap_times(1).start  lap_times(2).start]);
        track_orders = [ones(1,length(lap_times(1).start))  2*ones(1,length(lap_times(2).start))];
        track_orders = track_orders(sorted_id);

        fig(nsession) = figure(nsession);
        fig(nsession).Position = [500 100 1200 900];
        fig(nsession).Name = sprintf('%s %s RUN lap bayesian bias and log odds',options.SUBJECT,options.SESSION)
        colour_lines = {'b','r'};
        clear h s
        sgtitle(sprintf('%s %s RUN lap bayesian bias and log odds',options.SUBJECT,options.SESSION))
        subplot(2,2,1)
        for nprobe = 1:length(session_info(n).probe)
            probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere

            h(nprobe) = plot(V1_bayesian_bias{probe_hemisphere}(sorted_id),colour_lines{probe_hemisphere});hold on;
        end
        %         plot(HPC_bayesian_bias(sorted_id),'k');hold on;

        h(3) = bar(HPC_bayesian_bias(sorted_id),'k','FaceAlpha',0.3,'EdgeColor','none');hold on;
        for nlap = 1:length(track_orders)
            if track_orders(nlap) == 1
                s(1) = scatter(nlap,track_orders(nlap) ,3,'r','filled','MarkerFaceAlpha',1)
            elseif track_orders(nlap) == 2
                s(2) = scatter(nlap,track_orders(nlap) - 2.1,3,'b','filled','MarkerFaceAlpha',1)
            end
        end
        ylabel('Bayesian Bias')
        xlabel('lap id')
        yline(0.5,'--')
        set(gca,"TickDir","out",'box', 'off','Color','none')

        %         legend([h(1:3),s(1),s(2)],{'V1 Left','V1 Right','HPC','Track 1','Track 2'})



        subplot(2,2,2)
        for nprobe = 1:length(session_info(n).probe)
            probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere

            h(probe_hemisphere) = plot(V1_z_log_odds{probe_hemisphere}(sorted_id),colour_lines{probe_hemisphere});hold on;
        end
        %         plot(HPC_bayesian_bias(sorted_id),'k');hold on;

        h(3) = bar(z_log_odds(sorted_id),'k','FaceAlpha',0.3,'EdgeColor','none');hold on;
        for nlap = 1:length(track_orders)
            if track_orders(nlap) == 1
                s(1) = scatter(nlap,track_orders(nlap) +2.5,3,'r','filled','MarkerFaceAlpha',1)
            elseif track_orders(nlap) == 2
                s(2) = scatter(nlap,track_orders(nlap) -5.5,3,'b','filled','MarkerFaceAlpha',1)
            end
        end
        ylabel('Log odds (z)')
        xlabel('lap id')
        yline(0,'--')
        set(gca,"TickDir","out",'box', 'off','Color','none')

        if isempty(V1_z_log_odds{2})
            legend([h(1),h(3),s(1),s(2)],{'V1 Left','HPC','Track 1','Track 2'})
        elseif isempty(V1_z_log_odds{1})
            legend([h(2:3),s(1),s(2)],{'V1 Right','HPC','Track 1','Track 2'})
        else
            legend([h(1:3),s(1),s(2)],{'V1 Left','V1 Right','HPC','Track 1','Track 2'})
        end


        track_orders = [ones(1,length(lap_times(1).start))  2*ones(1,length(lap_times(2).start))];

        subplot(2,4,5)
        if ~isempty(V1_bayesian_bias{1})
            scatter(HPC_bayesian_bias(track_orders == 1), V1_bayesian_bias{1}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.2)
            hold on
            scatter(HPC_bayesian_bias(track_orders == 2), V1_bayesian_bias{1}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.2)

            mdl = fitlm(HPC_bayesian_bias',V1_bayesian_bias{1});
            [pval,~,~] = coefTest(mdl);
            x =[min(HPC_bayesian_bias') max(HPC_bayesian_bias')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'r:')
                %         xlim([-30 120])

                %     title(sprintf('Session %i',s),'Color','red')
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
            else
                plot(x,y_est,'k:')
                %         xlim([-30 120])
                %     title(sprintf('Session %i',s))
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end
            ylabel('V1 Left Bayesian bias')
            xlabel('HPC Bayesian bias')
            yline(0.5,'--')
            xline(0.5,'--')
        end
        set(gca,"TickDir","out",'box', 'off','Color','none')

        subplot(2,4,6)
        if ~isempty(V1_bayesian_bias{2})
            scatter(HPC_bayesian_bias(track_orders == 1), V1_bayesian_bias{2}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.2)
            hold on
            scatter(HPC_bayesian_bias(track_orders == 2), V1_bayesian_bias{2}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.2)

            mdl = fitlm(HPC_bayesian_bias',V1_bayesian_bias{2});
            [pval,~,~] = coefTest(mdl);
            x =[min(HPC_bayesian_bias') max(HPC_bayesian_bias')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'r:')
                %         xlim([-30 120])

                %     title(sprintf('Session %i',s),'Color','red')
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
            else
                plot(x,y_est,'k:')
                %         xlim([-30 120])
                %     title(sprintf('Session %i',s))
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end
            ylabel('V1 Right Bayesian bias')
            xlabel('HPC Bayesian bias')
            yline(0.5,'--')
            xline(0.5,'--')
        end
        set(gca,"TickDir","out",'box', 'off','Color','none')

        subplot(2,4,7)
        if ~isempty(V1_z_log_odds{1})
            scatter(z_log_odds(track_orders == 1), V1_z_log_odds{1}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.2)
            hold on
            scatter(z_log_odds(track_orders == 2), V1_z_log_odds{1}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.2)

            mdl = fitlm(z_log_odds',V1_z_log_odds{1});
            [pval,~,~] = coefTest(mdl);
            x =[min(z_log_odds') max(z_log_odds')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'r:')
                %         xlim([-30 120])

                %     title(sprintf('Session %i',s),'Color','red')
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
            else
                plot(x,y_est,'k:')
                %         xlim([-30 120])
                %     title(sprintf('Session %i',s))
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end
            ylabel('V1 Left Log odds (z)')
            xlabel('HPC log odds (z)')
            yline(0,'--')
            xline(0,'--')
        end
        set(gca,"TickDir","out",'box', 'off','Color','none')

        subplot(2,4,8)
        if ~isempty(V1_z_log_odds{2})
            scatter(z_log_odds(track_orders == 1), V1_z_log_odds{2}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.2)
            hold on
            scatter(z_log_odds(track_orders == 2), V1_z_log_odds{2}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.2)

            mdl = fitlm(z_log_odds',V1_z_log_odds{2});
            [pval,~,~] = coefTest(mdl);
            x =[min(z_log_odds') max(z_log_odds')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'r:')
                %         xlim([-30 120])

                %     title(sprintf('Session %i',s),'Color','red')
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
            else
                plot(x,y_est,'k:')
                %         xlim([-30 120])
                %     title(sprintf('Session %i',s))
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end
            ylabel('V1 Right Log odds (z)')
            xlabel('HPC log odds (z)')
            yline(0,'--')
            xline(0,'--')
            %             track_orders = [ones(1,length(lap_times(1).start))  2*ones(1,length(lap_times(2).start))];
            %             scatter3(z_log_odds(track_orders == 1), V1_z_log_odds{1}(track_orders == 1),V1_z_log_odds{2}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.5)
            %             hold on
            %             scatter3(z_log_odds(track_orders == 2), V1_z_log_odds{1}(track_orders == 2),V1_z_log_odds{2}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.5)
            %             xlabel('HPC log odds')
            %             ylabel('left V1 log odds')
            %             zlabel('right V1 log odds')
            %             set(gca,"TickDir","out",'box', 'off','Color','none')
            %             legend([h(1:3)],{'V1 Left','V1 Right','HPC'})

        end
        set(gca,"TickDir","out",'box', 'off','Color','none')
        fontsize(gcf,14,"points")


        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);
            options.ROOTPATH = ROOTPATH;
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no;
            probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere

            data = V1_bayesian_bias{probe_hemisphere};


            FPR = [];
            TPR = [];
            AUC = [];
            data_resampled = [];
            track_label_resampled = [];

            for nboot = 1:1000
                s = RandStream('mrg32k3a','Seed',nboot); % Set random seed for resampling

                index = datasample(s,1:length(data),length(data));
                data_resampled(nboot,:) = data(index);
                track_label_resampled(nboot,:) = track_orders(index)-1;
                [X,Y,T,A] = perfcurve(track_orders(index),data(index),1,'XVals',0:0.05:1,'NBoot',1);
                %             [X,Y,T,A] = perfcurve(track_label(index),z_log_odds(index),1,'XVals',0:0.05:1);

                FPR = X;
                TPR(nboot,:) = Y(:,1);
                AUC(nboot) = A(1);
            end

            if session_info(n).probe(nprobe).probe_hemisphere == 1
                fig(11) = figure(11);
                fig(11).Position = [500 100 1200 900];
                fig(11).Name = 'lap Bayesian Bias ROC two track discrimination in V1 for left probe';
                % ROC#
                subplot(2,5,nsession)
                hold on
                x = FPR';
                CI_shuffle = prctile(TPR,[2.5 97.5]);
                plot(x, CI_shuffle(2,:), 'r--', 'LineWidth', 1);
                plot(x, CI_shuffle(1,:), 'r--', 'LineWidth', 1);
                x2 = [x, fliplr(x)];
                inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
                h(2) = fill(x2, inBetween, 'r','FaceAlpha',0.2);

                h(1) = plot([0 1],[0 1],'k--')
                set(gca,"TickDir","out",'box', 'off','Color','none')
                legend([h(2) h(1)],{'Real','chance'})
                title(sprintf('Session %i AUC %.2f',nsession,mean(AUC)))
                sgtitle('lap Bayesian Bias ROC two track discrimination in V1 for left probe')
                fontsize(gcf,14,"points")
            elseif session_info(n).probe(nprobe).probe_hemisphere == 2
                fig(12) = figure(12);
                fig(12).Position = [500 100 1200 900];
                fig(12).Name = 'lap Bayesian Bias ROC two track discrimination in V1 for right probe';
                % ROC#
                subplot(2,5,nsession)
                hold on
                x = FPR';
                CI_shuffle = prctile(TPR,[2.5 97.5]);
                plot(x, CI_shuffle(2,:), 'r--', 'LineWidth', 1);
                plot(x, CI_shuffle(1,:), 'r--', 'LineWidth', 1);
                x2 = [x, fliplr(x)];
                inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
                h(2) = fill(x2, inBetween, 'r','FaceAlpha',0.2);

                h(1) = plot([0 1],[0 1],'k--')
                set(gca,"TickDir","out",'box', 'off','Color','none')
                legend([h(2) h(1)],{'Real','chance'})
                title(sprintf('Session %i AUC %.2f',nsession,mean(AUC)))
                sgtitle('lap Bayesian Bias ROC two track discrimination in V1 for right probe')
                fontsize(gcf,14,"points")
            end
        end

        data = HPC_bayesian_bias;

        FPR = [];
        TPR = [];
        AUC = [];
        data_resampled = [];
        track_label_resampled = [];

        for nboot = 1:1000
            s = RandStream('mrg32k3a','Seed',nboot); % Set random seed for resampling

            index = datasample(s,1:length(data),length(data));
            data_resampled(nboot,:) = data(index);
            track_label_resampled(nboot,:) = track_orders(index)-1;
            [X,Y,T,A] = perfcurve(track_orders(index),data(index),1,'XVals',0:0.05:1,'NBoot',1);
            %             [X,Y,T,A] = perfcurve(track_label(index),z_log_odds(index),1,'XVals',0:0.05:1);

            FPR = X;
            TPR(nboot,:) = Y(:,1);
            AUC(nboot) = A(1);
        end

        fig(13) = figure(13);
        fig(13).Position = [500 100 1200 900];
        fig(13).Name = 'lap Bayesian Bias ROC two track discrimination in HPC';

        % ROC#
        subplot(2,5,nsession)
        hold on
        x = FPR';
        CI_shuffle = prctile(TPR,[2.5 97.5]);
        plot(x, CI_shuffle(2,:), 'r--', 'LineWidth', 1);
        plot(x, CI_shuffle(1,:), 'r--', 'LineWidth', 1);
        x2 = [x, fliplr(x)];
        inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
        h(2) = fill(x2, inBetween, 'r','FaceAlpha',0.2);

        h(1) = plot([0 1],[0 1],'k--')
        set(gca,"TickDir","out",'box', 'off','Color','none')
        legend([h(2) h(1)],{'Real','chance'})
        title(sprintf('Session %i AUC %.2f',nsession,mean(AUC)))
        sgtitle('lap Bayesian Bias ROC two track discrimination in HPC')
        fontsize(gcf,14,"points")

    


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
   
