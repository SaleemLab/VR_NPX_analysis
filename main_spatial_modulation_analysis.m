%% Main spatial modulation analysis code

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))
%% Spatial raster plot and spatial tuning curves & spatial modulation analysis

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]

for nsession = [1 2 3 4 9 10]
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


        probability_ratio_RUN_lap_HPC_combined= [];
        probability_ratio_RUN_lap_V1_combined= [];

        estimated_position_lap_CV_shuffled_HPC_combined.track = [];
        estimated_position_lap_CV_shuffled_V1_combined.track = [];
        estimated_position_lap_CV_HPC_combined.track = [];
        estimated_position_lap_CV_V1_combined.track = [];

        for nprobe = 1:length(clusters)
           
            %
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.unstable_ids = @(x) x==0;
            if clusters(nprobe).probe_hemisphere == 1
                metric_param.region = @(x) contains(x,'HPC_L');
                options.region = 'HPC Left';
            elseif clusters(nprobe).probe_hemisphere == 2
                metric_param.region = @(x) contains(x,'HPC_R');
                options.region = 'HPC Right';
            end

            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);


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


