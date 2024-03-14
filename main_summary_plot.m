%% Script for plotting key summary figures of the experiment 


addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))


%% Cluster density plots

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 6 7 8 9 10 12 14]

for nsession =[1 2 3 4 6 7 8 9 10 12 14]
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','extracted_PSD.mat'),'PSD');

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

        if exist(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_place_fields.mat'));
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
%             best_channels{nprobe}.xcoord
            xcoord_avaliable = best_channels{nprobe}.xcoord(~isnan(best_channels{nprobe}.surface_depth));

            power = [];
            for nchannel = 1:size(chan_config,1)
                power(nchannel,:) = PSD{nprobe}(nchannel).mean_power;
                xcoord(nchannel) = PSD{nprobe}(nchannel).xcoord;
                ycoord(nchannel) = PSD{nprobe}(nchannel).ycoord;
            end
%             xcoord_avaliable = unique(xcoord);

            % sort channel according to y coordinate
            [ycoord idx] = sort(ycoord,'ascend');
            xcoord = xcoord(idx);
            power = power(idx,:);
            chan_config = chan_config(idx,:);

            plot_cluster_density_profile(power(xcoord == xcoord_avaliable(1),:),chan_config,chan_config(xcoord == xcoord_avaliable(1),:),best_channels{nprobe},clusters(nprobe),options);


            spatial_cell_index = unique([find(place_fields(1).peak_percentile>0.95 & place_fields(1).odd_even_stability>0.95)...
                find(place_fields(2).peak_percentile>0.95 & place_fields(2).odd_even_stability>0.95)]);
             
            if nprobe == 1
                good_cell_index = intersect(spatial_cell_index,find(ismember(place_fields(1).cluster_id,clusters(nprobe).merged_cluster_id)));
                if isempty(good_cell_index)
                    continue
                end
                options.spatial_cell_id = place_fields(1).cluster_id(good_cell_index);

            elseif nprobe == 2
                good_cell_index = intersect(spatial_cell_index,find(ismember(place_fields(1).cluster_id-10000,clusters(nprobe).merged_cluster_id)));
                if isempty(good_cell_index)
                    continue
                end
                options.spatial_cell_id = place_fields(1).cluster_id(good_cell_index)-10000;
            end

            plot_cluster_density_profile(power(xcoord == xcoord_avaliable(1),:),chan_config,chan_config(xcoord == xcoord_avaliable(1),:),best_channels{nprobe},clusters(nprobe),options);
            options = rmfield(options,'spatial_cell_id');
        end

        if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','probe_depth_profile'))== 0
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','probe_depth_profile'))
        end

        save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','probe_depth_profile'),[])
      

        %         %GLM analysis
        %         spatial_modulation_GLM_analysis(V1_clusters_L,[],Behaviour,Task_info);
        %
        %         spatial_modulation_GLM_analysis(V1_clusters_L,place_fields_V1_L(1).cluster_id(place_fields_V1_L(1).all_good_cells_LIBERAL),Behaviour,Task_info);
    end
end

%% Across session spatial cell population summary

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 6 7 8 9 10 12 14]

for nsession =[1 2 3 4 6 7 8 9 10 12 14]
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','extracted_PSD.mat'),'PSD');

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

        if exist(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_place_fields.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'theta_modulation.mat'));
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
            
        end

        if exist(fullfile(options.ANALYSIS_DATAPATH,'..','..','figures','population_map'))== 0
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','..','figures','population_map'))
        end

        save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','..','figures','population_map'))

        %         %GLM analysis
        %         spatial_modulation_GLM_analysis(V1_clusters_L,[],Behaviour,Task_info);
        %
        %         spatial_modulation_GLM_analysis(V1_clusters_L,place_fields_V1_L(1).cluster_id(place_fields_V1_L(1).all_good_cells_LIBERAL),Behaviour,Task_info);
    end
end



%% Across sessions spatial map, theta modulation and ripple population summary

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 6 7 8 9 10 12 14]

place_fields_all = [];

for track_id = 1:2

    place_fields_all(track_id).session_id = [];
    place_fields_all(track_id).session = [];
    place_fields_all(track_id).cluster_id = [];
    place_fields_all(track_id).cell_type = [];
    place_fields_all(track_id).peak_channel_waveform = [];
    place_fields_all(track_id).dwell_map = [];
    place_fields_all(track_id).relative_depth = [];


    place_fields_all(track_id).x_bin_edges = [];
    place_fields_all(track_id).x_bin_centres =[];
    place_fields_all(track_id).average_map = [];
    place_fields_all(track_id).raw = [];

    place_fields_all(track_id).within_track_corr = [];
    place_fields_all(track_id).odd_even_stability = [];

    place_fields_all(track_id).t1_t2_remapping = [];
    place_fields_all(track_id).across_tracks_correlation = [];

    place_fields_all(pre_ripple_activation).ripple_PSTH = [];
    place_fields_all(track_id).ripple_modulation_percentile = [];
    place_fields_all(track_id).pre_ripple_activation = [];
    place_fields_all(track_id).pre_ripple_inhibition = [];
    place_fields_all(track_id).post_ripple_activation = [];
    place_fields_all(track_id).post_ripple_inhibition = [];

    place_fields_all(track_id).theta_modulation_percentile = [];
    place_fields_all(track_id).theta_phase_map = [];
    place_fields_all(track_id).position_phase_xcorr_map = [];
    place_fields_all(track_id).position_phase_map = [];
    place_fields_all(track_id).phase_coherence = [];
end

place_fields_all_L = place_fields_all;
place_fields_all_R = place_fields_all;

for nsession =[1 2 3 4 6 7 8 9 10 12 14]
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','extracted_PSD.mat'),'PSD');

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

        if exist(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_place_fields.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'theta_modulation.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'ripple_modulation.mat'));
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
        
        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end


        spatial_cell_index = unique([find(place_fields(1).peak_percentile>0.95 & place_fields(1).odd_even_stability>0.95)...
            find(place_fields(2).peak_percentile>0.95 & place_fields(2).odd_even_stability>0.95)]);

        if isfield(clusters_combined,'merged_cluster_id')
            good_cell_index = intersect(spatial_cell_index,find(ismember(place_fields(1).cluster_id,clusters_combined.merged_cluster_id)));
        else
            % Currently support just merged cluster id
        end

        ratemap_matrix = [];
        good_cell_index
        place_fields_all(1).session_id = [place_fields_all(1).session_id nsession*ones(1,length(good_cell_index))];
        place_fields_all.session = ;
        for track_id = 1:2
            ratemap_matrix = [place_fields(track_id).raw{good_cell_index}];
            ratemap_matrix = reshape(ratemap_matrix,size(place_fields(track_id).raw{1},1),[],length(good_cell_index));%laps X position bins X cells

            place_fields_all.average_map = squeeze(mean(ratemap_matrix,1));
            place_fields_all.raw = ratemap_matrix;

            all_ratemaps = [place_fields(1).raw{good_cell_index}];
            place_fields(1).within_track_corr(good_cell_index)
        end


        cluster_id = place_fields(1).cluster_id(good_cell_index);

        for nprobe = 1:length(session_info(n).probe) % ripple or theta
            options = session_info(n).probe(nprobe);
            options.importMode = 'KS';
            place_fields(nprobe)

            if options.probe_hemisphere == 1


            elseif options.probe_hemisphere == 2
                ripple_modulation_R
                theta_modulation_R

            end

        end


        %         %GLM analysis
        %         spatial_modulation_GLM_analysis(V1_clusters_L,[],Behaviour,Task_info);
        %
        %         spatial_modulation_GLM_analysis(V1_clusters_L,place_fields_V1_L(1).cluster_id(place_fields_V1_L(1).all_good_cells_LIBERAL),Behaviour,Task_info);
    end
end


if exist(fullfile(options.ANALYSIS_DATAPATH,'..','..','figures','population_map'))== 0
    mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','..','figures','population_map'))
end

save(fullfile(options.ANALYSIS_DATAPATH,'..','..','figures','population_map'))
save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','..','figures','population_map'))

%%


%%
clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';

for nsession =[1 2 3 4 6 7 8 9 10 12 14]
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


    end
end

