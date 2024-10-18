%% Set the data folders and processing parameters

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))

%% Check spatial response stability

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
            %             load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'))
            session_clusters_RUN1 = session_clusters;
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN2.mat'))
            session_clusters_RUN2 = session_clusters;
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
        end

        % reconstruct cluster structure and place field structure from
        % session clusters

        clusters_combined1= session_clusters_RUN1;
        clusters_combined1.spike_id=vertcat(session_clusters_RUN1.spike_id{:});
        clusters_combined1.spike_times=vertcat(session_clusters_RUN1.spike_times{:});
        [clusters_combined1.spike_times,index] =sort(clusters_combined1.spike_times);
        clusters_combined1.spike_id=clusters_combined1.spike_id(index);
        clusters_combined=clusters_combined1;

        clusters_combined2= session_clusters_RUN2;
        clusters_combined2.spike_id=vertcat(session_clusters_RUN2.spike_id{:});
        clusters_combined2.spike_times=vertcat(session_clusters_RUN2.spike_times{:});
        [clusters_combined2.spike_times,index] =sort(clusters_combined2.spike_times);
        clusters_combined2.spike_id=clusters_combined2.spike_id(index);


        % Cell with spatial tuning
        % ib = find((clusters_combined.peak_percentile(:,1)>0.95&clusters_combined.odd_even_stability(:,1)>0.95) ...
        %     | (clusters_combined.peak_percentile(:,2)>0.95&clusters_combined.odd_even_stability(:,2)>0.95));
        % ia = find(contains(clusters_combined.region,'HPC'))
        %         [C,ia,ic] = unique(clusters_combined.cluster_id);
        ia = find(clusters_combined1.odd_even_stability(:,1)>0.95 ...
            | clusters_combined1.odd_even_stability(:,2)>0.95);
        %         ia = find((clusters_combined1.odd_even_stability(:,1)>0.95 ...
        %             | clusters_combined1.odd_even_stability(:,2)>0.95)&contains(clusters_combined1.region,'HPC'));
        C = clusters_combined1.cluster_id(ia);


        ia = find((clusters_combined1.odd_even_stability(:,1)>0.95 ...
            | clusters_combined1.odd_even_stability(:,2)>0.95)&contains(clusters_combined1.region,'V1_L'));

        speed = clusters_combined.speed{1};
        speed(isnan(speed))=0;
        w = gausswin(9);
        w = w / sum(w);
        speed = filtfilt(w,1,speed')';
        lick_speed = interp1(clusters_combined.sglxTime_uncorrected{1},speed,clusters_combined.lick_time{1},'nearest');

        x_window = [0 140];
        x_bin_size =2;
        place_fields_all = calculate_spatial_cells(clusters_combined,clusters_combined.tvec{1},...
            clusters_combined.position{1},speed,clusters_combined.track_ID_all{1},clusters_combined.start_time_all{1},clusters_combined.end_time_all{1},x_window,x_bin_size);
        %
        %         ia = find((clusters_combined1.odd_even_stability(:,1)>0.95 ...
        %             | clusters_combined1.odd_even_stability(:,2)>0.95)&contains(clusters_combined1.region,'V1_L'));
        %         plot_place_cell_map_correlation(place_fields_all,ia,[],[],options)

        %         ia = find((clusters_combined1.odd_even_stability(:,1)>0.95 ...
        %             | clusters_combined1.odd_even_stability(:,2)>0.95));

        %         colorbar
        ia = find(clusters_combined1.odd_even_stability(:,1)>0.95 ...
            | clusters_combined1.odd_even_stability(:,2)>0.95);
%         ia=1:length(clusters_combined1.region);
        cluster_id= intersect(ia,find(contains(clusters_combined1.region,'HPC')));
        [~,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = ...
            plot_place_cell_map_correlation(place_fields_all,cluster_id,Task_info,Behaviour,options); % Roughly 6-7 mins for shuffle and plotting

        plot_spatial_map_stability(place_fields_all); 
        if  contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','spatial_map_stability',sprintf(erase(stimulus_name{n},'Masa2tracks_'))))
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','spatial_map_stability',sprintf(erase(stimulus_name{n},'Masa2tracks_'))),[])
        else
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','spatial_map_stability'),[])
        end

        close all

        %         %GLM analysis
        %         spatial_modulation_GLM_analysis(V1_clusters_L,[],Behaviour,Task_info);
        %
        %         spatial_modulation_GLM_analysis(V1_clusters_L,place_fields_V1_L(1).cluster_id(place_fields_V1_L(1).all_good_cells_LIBERAL),Behaviour,Task_info);
    end
end