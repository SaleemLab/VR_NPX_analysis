%% Main spatial modulation analysis code

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))
%%  spatial tuning curves & spatial modulation analysis

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


        [place_fields] = spatial_modulation_calculation(place_fields,Task_info);
        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_place_fields.mat'),'place_fields');
    end
end
%% Raster plot lap start

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]

for nsession = [7 8 9 10 12 14]
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'power');
        load(fullfile(options.ANALYSIS_DATAPATH,'..','extracted_PSD.mat'));

        if contains(stimulus_name{n},'Masa2tracks')
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
        else
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP')
        end

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
%         load(fullfile(options.ANALYSIS_DATAPATH,'ripple_modulation.mat'));
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end
        
        % Plotting raster plot

        [C,ia,ic] = unique(clusters_combined.merged_cluster_id);

        plot_perievent_spiketimes_vs_spatial_response(clusters_combined.spike_times,clusters_combined.merged_spike_id,Task_info,Behaviour,[5 1],[-1 2],0.02,...
            'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',Task_info.start_time_all,...
            'event_label','Lap start','event_id',Task_info.track_ID_all,'place_fields',place_fields);

        if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','lap start PSTH'))== 0
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','lap start PSTH'))
        end

        save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','lap start PSTH'),[])



        plot_perievent_spiketimes_vs_spatial_response(clusters_combined.spike_times,clusters_combined.merged_spike_id,Task_info,Behaviour,[5 1],[-1 2],0.02,...
            'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',Task_info.end_time_all',...
            'event_label','Lap end','event_id',Task_info.track_ID_all,'place_fields',place_fields);

        if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','lap end PSTH'))== 0
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','lap end PSTH'))
        end

        save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','lap end PSTH'),[])



        plot_raster_both_track(clusters_combined.spike_times,clusters_combined.merged_spike_id,Task_info,Behaviour,[5 1],[0 140],2,...
            'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'place_fields',place_fields);

        
        if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','Spatial PSTH'))== 0
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','Spatial PSTH'))
        end

        save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','Spatial PSTH'),[])

        plot_raster_both_track_simple(clusters_combined.spike_times,clusters_combined.merged_spike_id,Task_info,Behaviour,[5 1],[0 140],2,...
            'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'place_fields',place_fields);

        
        if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','Spatial PSTH simple'))== 0
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','Spatial PSTH simple'))
        end

        save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','Spatial PSTH simple'),[])
    end

end
