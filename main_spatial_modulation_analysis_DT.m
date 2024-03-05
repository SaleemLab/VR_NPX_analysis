%% Main place cell and V1 spatial modulation analysis code

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))
%% Spatial raster plot and spatial tuning curves & spatial modulation analysis

clear all
SUBJECTS = {'M23032','M23034','M23037','M23038'};
option = 'V1-MEC';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'Track';

for nsession = 3:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));


        options = session_info.probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'));
            clusters = clusters_ks3;
            sorting_option = 'spikeinterface';


        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));

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
        %
        for nprobe = 1:length(clusters)
            clusters(nprobe).region = strings(length(clusters(nprobe).cluster_id),1);
            clusters(nprobe).region(:) = 'n.a';
            options = session_info.probe(nprobe);
            if options.probe_id == options.probe_V1
                V1_channels = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by probe');
                HPC_channels = determine_region_channels(best_channels{nprobe},options,'region','HPC','group','by probe');
                clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,V1_channels))) = 'V1';
                clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,HPC_channels))) = 'HPC';

            elseif options.probe_id == options.probe_MEC 
                %             V1_channels = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by probe');
                MEC_channels = determine_region_channels(best_channels{nprobe},options,'region','MEC_entry','group','by probe');
                HVA_channels = determine_region_channels(best_channels{nprobe},options,'region','HVA','group','by probe');
                %             clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,V1_channels))) = 'V1';
                clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,MEC_channels))) = 'MEC_entry';
                clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,HVA_channels))) = 'HVA';
            end
        end

        %%%%%%%%%%%%%%%%%% place holder for merging units
        %%%%%%%%%%%%%%%%%%%% load match id

        metric_param = create_cluster_selection_params('sorting_option',sorting_option);
        metric_param.unstable_ids = @(x) x==0;

        place_field = struct();
        clear merged_clusters

        for nprobe = 1:length(clusters)
            options = session_info.probe(nprobe);

            if  exist(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('probe%ium_merge_suggestion.mat',options.probe_id)))
                load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('probe%ium_merge_suggestion.mat',options.probe_id)))
            elseif exist(fullfile(options.ANALYSIS_DATAPATH,'..','..','all_unit_match.mat'))==2
                load(fullfile(options.ANALYSIS_DATAPATH,'..','..','all_unit_match.mat'))
                for i = 1:length(match_ids)
                    unit_difference = size(match_ids{i},1)-length(clusters.cluster_id);
                    if sum(ismember(clusters.cluster_id,match_ids{i}(:,1))) == length(clusters.cluster_id)
                        match_ids = match_ids{i}; % find the session with the same unit id
                        
                        break
                    end
                end
                match_ids = match_ids(:,1:3);
            end

            %             [~,cluster_id] = select_clusters(clusters_combined,metric_param);
            merged_cluster_temp  = merge_cluster(clusters(nprobe),match_ids);
            merged_clusters(nprobe) = select_clusters(merged_cluster_temp,metric_param);
        end
        
        % save merged cluster variables (useful more reactivation activity detection)

            save(fullfile(options.ANALYSIS_DATAPATH,'merged_clusters.mat'))

        %

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end
        
        if exist(fullfile(options.ANALYSIS_DATAPATH,"extracted_place_fields.mat"))==2
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_place_fields.mat'),'place_fields')
        else
            place_fields = calculate_place_fields_masa_NPX_against_shuffle(clusters_combined,Task_info,Behaviour,[0 140],2,[]);
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_place_fields.mat'),'place_fields')
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

        for nprobe = 1:length(session_info.probe)
            options = session_info.probe(nprobe);

            % plot populational map and PV correlation
            %%%%%% HPC
            
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.unstable_ids = @(x) x==0;

            if options.probe_id == options.probe_V1 
                metric_param.region = @(x) contains(x,'HPC');
                options.region = 'HPC';
                options.text = {'HPC'};
            elseif options.probe_id == options.probe_MEC 
                metric_param.region = @(x) contains(x,'MEC_entry');
                options.region = 'MEC_entry';
                options.text = {'MEC'};
            end
            [~,cluster_id] = select_clusters(clusters_combined,metric_param);
            
            % Unique merged cell id
            cluster_id = unique(clusters_combined.merged_cluster_id(cluster_id));

            [~,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = ...
                plot_place_cell_map_correlation_DT(place_fields,cluster_id,Task_info,Behaviour,options); % Roughly 6-7 mins for shuffle and plotting

            if options.probe_id == options.probe_V1 
                save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_HPC.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
            elseif options.probe_id == options.probe_MEC 
                save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_MEC.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
            end

            % plot populational map and PV correlation
            %%%%%% V1
            
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.unstable_ids = @(x) x==0;

            if options.probe_id == options.probe_V1 
                metric_param.region = @(x) contains(x,'V1');
                options.region = 'V1';
                options.text = {'V1'};
            elseif options.probe_id == options.probe_MEC 
                metric_param.region = @(x) contains(x,'HVA');
                options.region = 'HVA';
                options.text = {'HVA'};
            end

            [~,cluster_id] = select_clusters(clusters_combined,metric_param);
            if sum(cluster_id) ~= 0
            % Unique merged cell id
            cluster_id = unique(clusters_combined.merged_cluster_id(cluster_id));

            [~,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = ...
                plot_place_cell_map_correlation_DT(place_fields,cluster_id,Task_info,Behaviour,options); % Roughly 6-7 mins for shuffle and plotting

            if options.probe_id == options.probe_V1 
                save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_V1.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
            elseif options.probe_id == options.probe_MEC 
                 save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_HVA.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
            end
            end
        end

%         if length(session_info(n).probe) > 1
%             % plot populational map and PV correlation
%             %%%%%% HPC
%             options.probe_combined = 1;
% 
%             options.region = 'HPC';
%             metric_param = create_cluster_selection_params('sorting_option',sorting_option);
%             metric_param.unstable_ids = @(x) x==0;
%             metric_param.region = @(x) contains(x,'HPC');
%             [~,cluster_id] = select_clusters(clusters_combined,metric_param);
% 
%             % Unique merged cell id
%             cluster_id = unique(clusters_combined.merged_cluster_id(cluster_id));
% 
%             [~,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = ...
%                 plot_place_cell_map_correlation(place_fields,cluster_id,Task_info,Behaviour,options); % Roughly 6-7 mins for shuffle and plotting
% 
%             if clusters(nprobe).probe_hemisphere == 1
%                 save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_HPC_combined.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
%             elseif clusters(nprobe).probe_hemisphere == 2
%                 save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_HPC_combined.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
%             end
% 
%             % plot populational map and PV correlation
%             %%%%%% V1
%             options.region = 'V1';
%             metric_param = create_cluster_selection_params('sorting_option',sorting_option);
%             metric_param.unstable_ids = @(x) x==0;
%             metric_param.region = @(x) contains(x,'V1');
% 
%             [~,cluster_id] = select_clusters(clusters_combined,metric_param);
% 
%             % Unique merged cell id
%             cluster_id = unique(clusters_combined.merged_cluster_id(cluster_id));
% 
%             [~,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = ...
%                 plot_place_cell_map_correlation(place_fields,cluster_id,Task_info,Behaviour,options); % Roughly 6-7 mins for shuffle and plotting
% 
%             if clusters(nprobe).probe_hemisphere == 1
%                 save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_V1_combined.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
%             elseif clusters(nprobe).probe_hemisphere == 2
%                 save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_V1_combined.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
%             end
%             options = rmfield(options,'probe_combined');
%         end

        if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','populational_map'))== 0
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','populational_map'))
        end

        save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','populational_map'),[])
        close all

        %         %GLM analysis
        %         spatial_modulation_GLM_analysis(V1_clusters_L,[],Behaviour,Task_info);
        %
        %         spatial_modulation_GLM_analysis(V1_clusters_L,place_fields_V1_L(1).cluster_id(place_fields_V1_L(1).all_good_cells_LIBERAL),Behaviour,Task_info);

end


%% Spatial bayesian decoding
% 
% clear all
% SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
% option = 'bilateral';
% experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Stimulus_type = 'RUN';
% % [1 2 3 4 9 10 12 14]
% 
% for nsession = [9 10 12 14]
%     session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%     stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%     load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));
% 
%     for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
%         options = session_info(n).probe(1);
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
% 
%         if exist(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
%             load(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%             clusters = merged_clusters;
%             sorting_option = 'spikeinterface';
%         elseif exist(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
%             load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%             clusters = clusters_ks3;
%             sorting_option = 'spikeinterface';
%         else
%             load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%             sorting_option = 'old';
%         end
% 
% 
%         load(fullfile(options.ANALYSIS_DATAPATH,'extracted_place_fields.mat'));
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
% 
%         if length(clusters) > 1
%             clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
%         else
%             clusters_combined = merged_clusters;
%         end
% 
%         estimated_position_lap_CV_HPC= [];
%         probability_ratio_RUN_lap_HPC= [];
%         probability_ratio_RUN_lap_V1= [];
%         estimated_position_lap_CV_V1 = [];
% 
%         for nprobe = 1:length(clusters)
%             options = session_info(n).probe(nprobe);
% 
%             % Bayesian decoding 10 fold cross validated
%             DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'RUN Bayesian decoding'));
%             if isempty(DIR)
%                 mkdir(fullfile(options.ANALYSIS_DATAPATH,'RUN Bayesian decoding'))
%             end
% 
%             %
%             metric_param = create_cluster_selection_params('sorting_option',sorting_option);
%             metric_param.unstable_ids = @(x) x==0;
%             if clusters(nprobe).probe_hemisphere == 1
%                 metric_param.region = @(x) contains(x,'HPC_L');
%                 options.region = 'HPC Left';
%             elseif clusters(nprobe).probe_hemisphere == 2
%                 metric_param.region = @(x) contains(x,'HPC_R');
%                 options.region = 'HPC Right';
%             end
% 
%             [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
%             [probability_ratio_RUN_lap_HPC{nprobe},estimated_position_lap_CV_HPC(nprobe).track] = bayesian_decoding_RUN_lap_cross_validation(selected_clusters,place_fields,Behaviour,Task_info,options);
% 
% 
%             if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'))== 0
%                 mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'))
%             end
% 
%             save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])
% 
%             metric_param = create_cluster_selection_params('sorting_option',sorting_option);
%             metric_param.unstable_ids = @(x) x==0;
%             if clusters(nprobe).probe_hemisphere == 1
%                 metric_param.region = @(x) contains(x,'V1_L');
%                 options.region = 'V1 Left';
%             elseif clusters(nprobe).probe_hemisphere == 2
%                 metric_param.region = @(x) contains(x,'V1_R');
%                 options.region = 'V1 Right';
%             end
% 
%             [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
%             [probability_ratio_RUN_lap_V1{nprobe},estimated_position_lap_CV_V1(nprobe).track] = bayesian_decoding_RUN_lap_cross_validation(selected_clusters,place_fields,Behaviour,Task_info,options);
% 
%             save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])
%         end
% 
%          if length(session_info(n).probe) > 1
%              options.region = 'HPC Combined';
%              metric_param = create_cluster_selection_params('sorting_option',sorting_option);
%              metric_param.unstable_ids = @(x) x==0;
%              metric_param.region = @(x) contains(x,'HPC');
% 
%              [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
%              [probability_ratio_RUN_lap_HPC_combined,estimated_position_lap_CV_HPC_combined.track] = bayesian_decoding_RUN_lap_cross_validation(selected_clusters,place_fields,Behaviour,Task_info,options);
% 
%              save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])
% 
%              options.region = 'V1 Combined';
%              metric_param = create_cluster_selection_params('sorting_option',sorting_option);
%              metric_param.unstable_ids = @(x) x==0;
%              metric_param.region = @(x) contains(x,'V1');
% 
%              [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
%              [probability_ratio_RUN_lap_V1_combined,estimated_position_lap_CV_V1_combined.track] = bayesian_decoding_RUN_lap_cross_validation(selected_clusters,place_fields,Behaviour,Task_info,options);
% 
%              save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])
% 
%          else
%              probability_ratio_RUN_lap_V1_combined = [];
%              estimated_position_lap_CV_V1_combined = [];
% 
%              probability_ratio_RUN_lap_HPC_combined = [];
%              estimated_position_lap_CV_HPC_combined = [];
%          end
% 
%          save(fullfile(options.ANALYSIS_DATAPATH,'estimated_position_lap_CV_HPC.mat'),'estimated_position_lap_CV_HPC_combined','estimated_position_lap_CV_HPC')
%          save(fullfile(options.ANALYSIS_DATAPATH,'probability_ratio_RUN_lap_HPC.mat'),'probability_ratio_RUN_lap_HPC_combined','probability_ratio_RUN_lap_HPC')
% 
%          save(fullfile(options.ANALYSIS_DATAPATH,'estimated_position_lap_CV_V1.mat'),'estimated_position_lap_CV_V1',"estimated_position_lap_CV_V1_combined")
%          save(fullfile(options.ANALYSIS_DATAPATH,'probability_ratio_RUN_lap_V1.mat'),'probability_ratio_RUN_lap_V1','probability_ratio_RUN_lap_V1_combined')
% 
%     end
% end


