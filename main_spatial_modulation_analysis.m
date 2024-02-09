%% Main place cell and V1 spatial modulation analysis code

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))


%% Spatial raster plot and spatial receptive field
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
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        clusters = clusters_ks3;

        metric_param = create_cluster_selection_params;
        clusters_R= [];
        clusters_L= [];
        place_field = struct();


        for nprobe = 1:length(clusters)
            
            if clusters(nprobe).probe_hemisphere == 1
                clusters_L = clusters(nprobe);
                V1_channels_L = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by probe');
                HPC_channels_L = determine_region_channels(best_channels{nprobe},options,'region','HPC','group','by probe');

            elseif clusters(nprobe).probe_hemisphere == 2
                clusters_R = clusters(nprobe);
                V1_channels_R = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by probe');
                HPC_channels_R = determine_region_channels(best_channels{nprobe},options,'region','HPC','group','by probe');
            end

        end

        % V1

        metric_param.peak_channel = @(x) ismember(x,V1_channels_L); % metric_param.depth_range = [] -- full range?
        [V1_clusters_L] = select_clusters(clusters_L,metric_param);

        metric_param.peak_channel = @(x) ismember(x,HPC_channels_L); % metric_param.depth_range = [] -- full range?
        [HPC_clusters_L] = select_clusters(clusters_L,metric_param);

        metric_param.peak_channel = @(x) ismember(x,V1_channels_R); % metric_param.depth_range = [] -- full range?
        [V1_clusters_R] = select_clusters(clusters_R,metric_param);

        metric_param.peak_channel = @(x) ismember(x,HPC_channels_R); % metric_param.depth_range = [] -- full range?
        [HPC_clusters_R] = select_clusters(clusters_R,metric_param);

        % HPC
%         plot_raster_both_track(selected_clusters_probe0,Task_info,Behaviour,[5 1])

        plot_raster_both_track(V1_clusters_L,Task_info,Behaviour,[5 1],[0 140],5)
        plot_raster_both_track(V1_clusters_R,Task_info,Behaviour,[5 1],[0 140],5)

        plot_raster_both_track(HPC_clusters_L,Task_info,Behaviour,[5 1],[0 140],5)
        plot_raster_both_track(HPC_clusters_R,Task_info,Behaviour,[5 1],[0 140],5)

        plot_raster_single_track(HPC_clusters_L,Task_info,Behaviour,[5 1],[0 140],5)

        plot_raster_both_track_extended(HPC_clusters_L,Task_info,Behaviour,[3 1],[0 140],2)


        %plot speed of each lap
        no_lap = size(Task_info.start_time_all,1);
        no_bin = 70;
        psth_speed = zeros(no_lap,no_bin);
        for iLap = 1:no_lap
            lap_index = Behaviour.tvec >= Task_info.start_time_all(iLap) & Behaviour.tvec <= Task_info.end_time_all(iLap);
            timevec_lap = Behaviour.tvec(lap_index);
            position_lap = Behaviour.position(lap_index);
            speed_lap = Behaviour.speed(lap_index);
            [N,edges,bin] = histcounts(position_lap,no_bin);
            for iBin = 1:no_bin
                psth_speed(iLap,iBin) = median(speed_lap(bin == iBin));
            end
        end
        figure;hold on;
        for iBlock = 1:5

            plot(1:2:139,mean(psth_speed((iBlock-1)*40+1:iBlock*40,:),'omitnan'))
        end



        spatial_modulation_GLM_analysis(clusters_L,Behaviour,Task_info);

    end
end



%% Spatial modulation GLM analysis
clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(options.ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        
        clusters = clusters_ks3;

        for nprobe = 1:length(session_info(n).probe)
            metric_param = [];
            metric_param.isi_violations_ratio = @(x) x<=0.1;
            metric_param.amplitude_cutoff = @(x) x<=0.1;
            metric_param.amplitude_median = @(x) x>50;
            metric_param.peak_depth = @(x) x>max(clusters(nprobe).peak_depth)-2000 & x<max(clusters(nprobe).peak_depth-1200); % metric_param.depth_range = [] -- full range?
            [selected_clusters good_cell_index] = select_clusters(clusters(nprobe),metric_param);

            spatial_modulation_GLM_analysis(selected_clusters(nprobe),Behaviour,Task_info);
        end
    end
end



% SUBJECTS = {'M23017'}
SUBJECT = 'M23028';
SESSION = '20230706';
Stimulus_type = 'Masa2tracks';
if contains(Stimulus_type,'Masa2tracks')
    session_files = dir(fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT,'analysis',SESSION,Stimulus_type,'session_info*.mat'));
    for n = 1:length(session_files) % May have PRE RUN and POST sessions rather than just one
        load(fullfile(session_files(n).folder, session_files(n).name))
        extract_and_preprocess_NPX(session_info,Stimulus_type)
    end
else
    load(fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT,'analysis',SESSION,Stimulus_type,'session_info.mat'))
    extract_and_preprocess_NPX(session_info,Stimulus_type)

end


