%% Main ripple analysis codes

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))

%% Ripple modulation

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]

for nsession = [1 2 3 4 6 7 8 9 10 12 14]
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
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end

        ripple_modulation_L = [];
        ripple_modulation_R = [];
        ripple_modulation_combined = [];

        for nprobe = 1:length(clusters)

            options = session_info(n).probe(nprobe);
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
            [C,ia,ic] = unique(clusters_combined.merged_cluster_id);
            
            event_id = [ones(1,length(ripples(nprobe).T1_onset)) 2*ones(1,length(ripples(nprobe).T2_onset))];
            [event_times,index] = sort([ripples(nprobe).T1_onset ripples(nprobe).T2_onset]);

            if options.probe_hemisphere == 1
                [ripple_modulation_L]= ripple_modulation_analysis(clusters_combined.spike_times,clusters_combined.merged_spike_id,Task_info,Behaviour,[-2 2],0.02,...
                'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',event_times',...
                'event_label',{'Track 1','Track 2'},'event_id',event_id(index)','place_fields',place_fields);

            elseif options.probe_hemisphere == 2
                [ripple_modulation_R]= ripple_modulation_analysis(clusters_combined.spike_times,clusters_combined.merged_spike_id,Task_info,Behaviour,[-2 2],0.02,...
                'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',event_times',...
                'event_label',{'Track 1','Track 2'},'event_id',event_id(index)','place_fields',place_fields);

            end

        end

        if length(session_info(n).probe) <= 1
            if session_info(n).probe.probe_hemisphere == 1
                ripple_modulation_combined = ripple_modulation_L;
            elseif session_info(n).probe.probe_hemisphere == 2
                ripple_modulation_combined = ripple_modulation_R;
            end
        else
% 
%             event_times = sort([ripples(1).onset ripples(2).onset]);
            [C,ia,ic] = unique(clusters_combined.merged_cluster_id);

            event_id = [ones(1,length(ripples(1).T1_onset)) ones(1,length(ripples(2).T1_onset)) 2*ones(1,length(ripples(1).T2_onset)) 2*ones(1,length(ripples(2).T2_onset))];
            [event_times,index] = sort([ripples(1).T1_onset ripples(2).T1_onset ripples(1).T2_onset ripples(2).T2_onset]);

            [ripple_modulation_combined]= ripple_modulation_analysis(clusters_combined.spike_times,clusters_combined.merged_spike_id,Task_info,Behaviour,[-2 2],0.02,...
                'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',event_times',...
                'event_label',{'Track 1','Track 2'},'event_id',event_id(index)','place_fields',place_fields);
        end
        save(fullfile(options.ANALYSIS_DATAPATH,'ripple_modulation.mat'),"ripple_modulation_L","ripple_modulation_R","ripple_modulation_combined")

    end
end

%% Ripple spike time raster plot

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]

load(fullfile('D:\corticohippocampal_replay\summary','place_fields_all.mat'))
% 
% selected_cells_L = unique([find(contains(place_fields_all_combined(1).region,'V1') & ...
%     place_fields_all_combined(1).ripple_modulation_percentile>0.95) find(contains(place_fields_all_combined(1).region,'V1') & ...
%     place_fields_all_combined(2).ripple_modulation_percentile>0.95)]);
% 
% place_fields_all_combined(1).session_id(selected_cells_L)
% % selected_cells_L = unique([find(contains(place_fields_all_combined(1).region,'V1_L') & ...
% %     place_fields_all_combined(1).ripple_modulation_percentile>0.95) find(contains(place_fields_all_combined(1).region,'V1_L') & ...
% %     place_fields_all_combined(2).ripple_modulation_percentile>0.95)]);
for nsession = [1 2 3 4 6 7 8 9 10 12 14]
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
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end
        
        % Plotting raster plot
        metric_param = create_cluster_selection_params('sorting_option',sorting_option);
        metric_param.unstable_ids = @(x) x==0;


        for nprobe = 1:length(merged_clusters)
            options = session_info(n).probe(nprobe);
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

            %             metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            %             metric_param.unstable_ids = @(x) x==0;
            %             metric_param.region = @(x) contains(x,'V1');
            %             [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            %             clusters_combined.merged_cluster_id
            event_id = [ones(1,length(ripples(nprobe).T1_onset)) 2*ones(1,length(ripples(nprobe).T2_onset))];
            [event_times,index] = sort([ripples(nprobe).T1_onset ripples(nprobe).T2_onset]);

            if  options.probe_hemisphere == 1
                %                 ripple_cells = unique([find(place_fields_all_L(1).ripple_modulation_percentile>0.95 & place_fields_all_L(1).session_id == nsession),...
                %                     find(place_fields_all_L(2).ripple_modulation_percentile>0.95&place_fields_all_L(2).session_id == nsession)]);
                %                 ripple_cells = place_fields_all_L(1).cluster_id(ripple_cells);
                %                 [C,ia,ic] = unique(clusters_combined.merged_cluster_id);
                %                 [Lia,Locb] = ismember(C,ripple_cells');
                %                 ia = ia(Lia);

                metric_param = create_cluster_selection_params('sorting_option',sorting_option);
                metric_param.unstable_ids = @(x) x==0;
                %             metric_param.region = @(x) contains(x,'V1');
                %                 metric_param.merged_cluster_id = @(x) ismember(x,clusters_combined.merged_cluster_id(ia));
                metric_param.merged_cluster_id = @(x) ismember(x,clusters_combined.merged_cluster_id);
                [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
                [C,ia,ic] = unique(selected_clusters.merged_cluster_id);

                %                 [C,ia,ic] = unique(merged_clusters(nprobe).merged_cluster_id);
                plot_perievent_spiketimes_vs_spatial_response(selected_clusters.spike_times,selected_clusters.merged_spike_id,Task_info,Behaviour,[5 1],[-2 2],0.02,...
                    'unit_depth',selected_clusters.peak_depth(ia),'unit_region',selected_clusters.region(ia),'unit_id',C,'event_times',event_times',...
                    'event_id',event_id(index),'event_label','L ripple','place_fields',place_fields,'speed_filtering',1);

                if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','ripple PSTH','Left ripples'))== 0
                    mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','ripple PSTH','Left ripples'))
                end
                save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','ripple PSTH','Left ripples'),[])

            elseif  options.probe_hemisphere == 2
%                 ripple_cells = unique([find(place_fields_all_R(1).ripple_modulation_percentile>0.95 & place_fields_all_R(1).session_id == nsession),...
%                     find(place_fields_all_R(2).ripple_modulation_percentile>0.95&place_fields_all_R(2).session_id == nsession)]);
%                 ripple_cells = place_fields_all_R(1).cluster_id(ripple_cells);
%                 [C,ia,ic] = unique(clusters_combined.merged_cluster_id);
%                 [Lia,Locb] = ismember(C,ripple_cells');
%                 ia = ia(Lia);

                metric_param = create_cluster_selection_params('sorting_option',sorting_option);
                metric_param.unstable_ids = @(x) x==0;
                %                             metric_param.region = @(x) contains(x,'V1');
%                 metric_param.merged_cluster_id = @(x) ismember(x,clusters_combined.merged_cluster_id(ia));
                metric_param.merged_cluster_id = @(x) ismember(x,clusters_combined.merged_cluster_id);
                [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
                [C,ia,ic] = unique(selected_clusters.merged_cluster_id);


                %                 [C,ia,ic] = unique(merged_clusters(nprobe).merged_cluster_id);
                plot_perievent_spiketimes_vs_spatial_response(selected_clusters.spike_times,selected_clusters.merged_spike_id,Task_info,Behaviour,[5 1],[-2 2],0.02,...
                    'unit_depth',selected_clusters.peak_depth(ia),'unit_region',selected_clusters.region(ia),'unit_id',C,'event_times',event_times',...
                    'event_id',event_id(index),'event_label','R ripple','place_fields',place_fields,'speed_filtering',1);

                if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','ripple PSTH','Right ripples'))== 0
                    mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','ripple PSTH','Right ripples'))
                end

                save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','ripple PSTH','Right ripples'),[])
            end

        end

        if length(session_info(n).probe) <= 1

        else
%             ripple_cells = unique([find(place_fields_all_combined(1).ripple_modulation_percentile>0.95 & place_fields_all_combined(1).session_id == nsession),...
%                 find(place_fields_all_combined(2).ripple_modulation_percentile>0.95&place_fields_all_combined(2).session_id == nsession)]);
%             ripple_cells = place_fields_all_combined(1).cluster_id(ripple_cells);
%             [C,ia,ic] = unique(clusters_combined.merged_cluster_id);
%             [Lia,Locb] = ismember(C,ripple_cells');
%             ia = ia(Lia);

            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.unstable_ids = @(x) x==0;
            %             metric_param.region = @(x) contains(x,'V1');
            metric_param.merged_cluster_id = @(x) ismember(x,clusters_combined.merged_cluster_id);
%             metric_param.merged_cluster_id = @(x) ismember(x,clusters_combined.merged_cluster_id(ia));
            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            [C,ia,ic] = unique(selected_clusters.merged_cluster_id);

            event_id = [ones(1,length(ripples(1).T1_onset)) ones(1,length(ripples(2).T1_onset)) 2*ones(1,length(ripples(1).T2_onset)) 2*ones(1,length(ripples(2).T2_onset))];
            [event_times,index] = sort([ripples(1).T1_onset ripples(2).T1_onset ripples(1).T2_onset ripples(2).T2_onset]);

            plot_perievent_spiketimes_vs_spatial_response(selected_clusters.spike_times,selected_clusters.merged_spike_id,Task_info,Behaviour,[5 1],[-2 2],0.02,...
                'unit_depth',selected_clusters.peak_depth(ia),'unit_region',selected_clusters.region(ia),'unit_id',C,'event_times',event_times',...
                'event_id',event_id(index),'event_label','combined ripple','place_fields',place_fields,'speed_filtering',1);

            if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','ripple PSTH','combined'))== 0
                mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','ripple PSTH','combined'))
            end

            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','ripple PSTH','combined'),[])

        end

%         for nprobe = 1:length(merged_clusters)
%             options = session_info(n).probe(nprobe);
%             probe_no = session_info(n).probe(nprobe).probe_id + 1;
%             options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
%             [C,ia,ic] = unique(clusters_combined.merged_cluster_id);
% 
%             event_id = [ones(1,length(ripples(nprobe).T1_onset)) 2*ones(1,length(ripples(nprobe).T2_onset))];
%             [event_times,index] = sort([ripples(nprobe).T1_onset ripples(nprobe).T2_onset]);
%             
%             if  options.probe_hemisphere == 1
%                 %                 [C,ia,ic] = unique(merged_clusters(nprobe).merged_cluster_id);
%                 plot_perievent_spiketimes(clusters_combined.spike_times,clusters_combined.merged_spike_id,Task_info,Behaviour,[5 1],[-2 2],0.02,...
%                     'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',event_times',...
%                     'event_id',event_id(index),'event_label','L ripple','place_fields',place_fields);
% 
%             elseif  options.probe_hemisphere == 2
%                 %                 [C,ia,ic] = unique(merged_clusters(nprobe).merged_cluster_id);
%                 plot_perievent_spiketimes(clusters_combined.spike_times,clusters_combined.merged_spike_id,Task_info,Behaviour,[5 1],[-2 2],0.02,...
%                     'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',event_times',...
%                     'event_id',event_id(index),'event_label','R ripple','place_fields',place_fields);
%             end
% 
%         end

    end
end

%% Populational level V1-HPC interaction during context-selective ripples CCA

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]
load(fullfile('D:\corticohippocampal_replay\summary','place_fields_all.mat'),'place_fields_all_L','place_fields_all_R','place_fields_all_combined')

psthBinSize = 0.001; %create 1ms bin
% Define Gaussian window for smoothing
gaussianWindow = gausswin(0.2*1/psthBinSize);
% Normalize to have an area of 1 (i.e., to be a probability distribution)
gaussianWindow = gaussianWindow / sum(gaussianWindow);


for nsession = [1 2 3 4 6 7 8 9 10 12 14]
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
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%         tic
%         load(fullfile(options.ANALYSIS_DATAPATH,'ripple_modulation.mat'))
%         toc

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end

        ripple_modulation_L = [];
        ripple_modulation_R = [];
        ripple_modulation_combined = [];

        metric_param = create_cluster_selection_params('sorting_option',sorting_option);
        metric_param.unstable_ids = @(x) x==0;
        %             metric_param.region = @(x) contains(x,'V1');
        metric_param.merged_cluster_id = @(x) ismember(x,clusters_combined.merged_cluster_id);
        %             metric_param.merged_cluster_id = @(x) ismember(x,clusters_combined.merged_cluster_id(ia));
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        [C,ia,ic] = unique(selected_clusters.merged_cluster_id);

        event_id = [ones(1,length(ripples(1).T1_onset)) ones(1,length(ripples(2).T1_onset)) 2*ones(1,length(ripples(1).T2_onset)) 2*ones(1,length(ripples(2).T2_onset))];
        
        if length(merged_clusters)>1
            probe_hemisphere = [session_info(n).probe.probe_hemisphere];
            temp = [ones(1,length(ripples(1).T1_onset)) 2*ones(1,length(ripples(2).T1_onset)) ones(1,length(ripples(1).T2_onset)) 2*ones(1,length(ripples(2).T2_onset))];
        else
            probe_hemisphere = [session_info(n).probe.probe_hemisphere];
            temp = ones(1,length(event_id));
        end

        % whole session timevec variables
        timevec = Behaviour.tvec';
        timevec_edge = (timevec(1)-(psthBinSize)/2....
            :psthBinSize:...
            timevec(end)+(psthBinSize)/2)';

        probe_id = temp;
        for nprobe = 1:length(merged_clusters)
            probe_id(temp==nprobe) = probe_hemisphere(nprobe);
        end

%         [event_times,index] = sort([ripples(1).T1_onset ripples(2).T1_onset ripples(1).T2_onset ripples(2).T2_onset]);
        [event_times,index] = sort([Task_info.start_time_all]);
        event_id = Task_info.track_ID_all;

        event_id = event_id(index);
        probe_id = probe_id(index);

        all_spike_counts = [];
        %         event_id = [ones(1,size(ripple_modulation_combined(1).spike_count{ncell},1))  2*ones(1,size(ripple_modulation_combined(2).spike_count{ncell},1))];
%                 ripple_cells = unique([find(ripple_modulation_combined(1).ripple_modulation_percentile>=0.95) find(ripple_modulation_combined(2).ripple_modulation_percentile>=0.95)]);
        ripple_cells = 1:length(ripple_modulation_combined(1).ripple_modulation_percentile);

        for ncell = 1:length(ripple_cells)

            spike_times = clusters_combined.spike_times(clusters_combined.merged_spike_id == ripple_modulation_combined(1).cluster_id(ripple_cells(ncell)));

            %             % Convolve spike count time series with Gaussian window
            %             y = histcounts(spike_times, timevec_edge)';
            %             y = conv(y, gaussianWindow, 'same')/psthBinSize; % all spikes

            spike_speed = interp1(Behaviour.tvec,Behaviour.speed,spike_times,'nearest');
            spike_times = spike_times(spike_speed<5);

            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times,...
                event_times, [-2 2], psthBinSize);
            all_spike_counts(ncell,:,:) = binnedArray(:,bins<0.5&bins>-0.5)'; % ncell X time x events

%             if size(binnedArray,1)==1
%                 psth = conv(binnedArray/psthBinSize,gaussianWindow,'same');
%             else
%                 for nevent = 1:size(binnedArray,1)
%                     psth(nevent,:) = conv(binnedArray(nevent,:)/psthBinSize,gaussianWindow,'same');
%                 end
                %                 psth = (psth-mean(y))/std(y);
                %                 psth_track1 = mean(psth_track1,'omitnan');
%             end

%             all_spike_counts(ncell,:,:) = psth(:,bins<0.5&bins>-0.5)'; % ncell X time x events
            %                 all_spike_counts{track_id}(ncell,:,:) = [ripple_modulation_combined(1).spike_count{ncell}(:,bins<0.5&bins>-0.5);...
            %                     ripple_modulation_combined(2).spike_count{ncell}(:,bins<0.5&bins>-0.5);]; % ncell X time x events
        end

        % V1 L and R
        ripple_spike_count = [];
        cell_regions = ripple_modulation_combined(1).region(ripple_cells);
        for nprobe = 1:length(probe_hemisphere)

            if session_info(n).probe(nprobe).probe_hemisphere == 1
                ripple_spike_count{1} = all_spike_counts(contains(cell_regions,'V1_L'),:,:);
            elseif session_info(n).probe(nprobe).probe_hemisphere == 2
                ripple_spike_count{2} = all_spike_counts(contains(cell_regions,'V1_R'),:,:);
            end

        end

        % HPC
        ripple_spike_count{3} = all_spike_counts(contains(cell_regions,'HPC'),:,:);


        %%% Example computation of a cross-correlation map
        % Figs. 3 and 4

        % The units of the arguments are with respect to the binning window used
        % to bin spikes.
        pairs = [1 2;3 1;3 2];
        pair_name = {'V1 L - V1 R','V1 L - HPC','V1 R - HPC'};
        argOut_all = [];
        for npair = 1:size(pairs,1)
            argIn=[];
            argIn.BinWidth = 2;     % 20ms
            argIn.MaxDelay = 50;    % 500ms
            argIn.TimeStep = 4;     % 40ms
            argIn.WindowLength = 8; % 80ms

            %argIn.NumWorkers = Inf; % Requires Parallel Processing Toolbox

            disp(argIn)
                        argOut = ComputeCorrMap({ripple_spike_count{pairs(npair,1)} ripple_spike_count{pairs(npair,2)}}, event_id, argIn);
%             argOut = ComputeCorrMap({ripple_spike_count{2}(2:2:size(ripple_spike_count{2},1),:,:) ripple_spike_count{2}(1:2:size(ripple_spike_count{2},1),:,:)}, event_id, argIn);

            %%%
            CANONICAL_PAIR_IDX = 1;
            mapDim = size(argOut.CorrMap, 2);
            delays = (-argIn.MaxDelay:argIn.MaxDelay)*10; % Convert to ms
            t = (-0.5*100:argIn.TimeStep:0.5*100)*10; % Convert to ms

            figure(npair);
            subplot(2,2,1)

            imagesc( delays, t, argOut.CorrMap(:,:,CANONICAL_PAIR_IDX)' )

            ax = gca;
            ax.YDir = 'Normal';

            xlabel('Delay')
            ylabel('Time')
            colorbar
            clim([0.1 0.5])
            xlim([-200 200])
            ylim([-500 500])
%             colormap(flip(gray))
            title(pair_name{npair})

            subplot(2,2,2)

            imagesc( delays, t, argOut.CorrMap(:,:,2)' )

            ax = gca;
            ax.YDir = 'Normal';

            xlabel('Delay')
            ylabel('Time')
            colorbar
            clim([0.1 0.5])
            xlim([-200 200])
            ylim([-500 500])
%             colormap(flip(gray))
            title('2nd axis')

            subplot(2,2,3)

            imagesc( delays, t, argOut.FrMap(:,:)' )

            ax = gca;
            ax.YDir = 'Normal';

            xlabel('Delay')
            ylabel('Time')
            colorbar
%             clim([0.2 0.4])
            xlim([-200 200])
            ylim([-500 500])
            %             colormap(flip(gray))
            title('geomertic mean FR')

            subplot(2,2,4)

            plot( delays,mean(argOut.CorrMap(:,t>=-40 & t<=40,CANONICAL_PAIR_IDX),2),'k')
            hold on
            plot( delays,mean(argOut.CorrMap(:,t<=150 & t>=50,CANONICAL_PAIR_IDX),2) ,'r')
            plot( delays,mean(argOut.CorrMap(:,t>=-150 & t<=-50,CANONICAL_PAIR_IDX),2),'b' )
            argOut_all{npair} = argOut;
        end
        %%% Example computation of the interaction structure analysis
        % Fig. 6

        clear argIn

        argIn.TimePeriods = [...
            (  0:20:40)' ( 20:20:60)'; ...
            (128:20:168)' (148:20:188)'] + 5;

        % Can take up to 15min due to the 10-fold cross-validation
        argOut = CovStabilityAcrossTimeAnalysis(spikes, expCond, argIn);

        %%%
        RANK_TO_PLOT = 2;

        normFactor = diag(argOut.CvR(:,:,RANK_TO_PLOT));
        numTimePeriods = size(argIn.TimePeriods, 1);

        figure(2);

        imagesc(argOut.CvR(:,:,RANK_TO_PLOT)./repmat(normFactor', numTimePeriods, 1))

        ax = gca;
        ax.YDir = 'Normal';

        axis square

        xlabel('Time Used For Correlation')
        ylabel('Time Used For Fitting')


    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LFP analysis during ripples %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Peri ripple LFP correlation

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]

for nsession = 10
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
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end
        
        speed_interp = interp1(Behaviour.tvec,Behaviour.speed,LFP(1).tvec','linear');


        % Distribution of ripples
        hemisphere_text = {'Left ripples','Right ripples'};
        probe_color = {'b','r'};
        position_edges = 0:5:139.5;
        bin_centres = position_edges(1:end-1) + diff(position_edges)/2;

        T1_immob_occupancy = histcounts(Behaviour.position(Behaviour.speed<5 & Behaviour.track_ID==1),position_edges);
        T2_immob_occupancy = histcounts(Behaviour.position(Behaviour.speed<5 & Behaviour.track_ID==2),position_edges);
        T1_running_occupancy = histcounts(Behaviour.position(Behaviour.speed>5 & Behaviour.track_ID==1),position_edges);
        T2_running_occupancy = histcounts(Behaviour.position(Behaviour.speed>5 & Behaviour.track_ID==2),position_edges);


        ripples_position = [];
        ripples_track_id = [];
        for nprobe = 1:length(clusters)
            options = session_info(n).probe(nprobe);
            ripples_position{options.probe_hemisphere} = Behaviour.position(ismember(Behaviour.tvec',Restrict(Behaviour.tvec,[ripples(nprobe).peaktimes; ripples(nprobe).peaktimes+mean(diff(Behaviour.tvec))]')));
            ripples_track_id{options.probe_hemisphere} = Behaviour.track_ID(ismember(Behaviour.tvec',Restrict(Behaviour.tvec,[ripples(nprobe).peaktimes; ripples(nprobe).peaktimes+mean(diff(Behaviour.tvec))]')));
        end
        
        fig = figure;
        fig.Position = [160 70 1550 850];
        fig.Name = sprintf('Distribution of ripples on track')
        subplot(3,4,1)
        bar(bin_centres,T1_running_occupancy/max(T1_running_occupancy),'r','FaceAlpha',0.5)
        hold on;
        bar(bin_centres,T2_running_occupancy/max(T2_running_occupancy),'b','FaceAlpha',0.5)

        legend('Track Left','Track Right','Color','none')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        title('Occupancy map when speed > 5')

        subplot(3,4,2)
        bar(bin_centres,T1_immob_occupancy/max(T1_immob_occupancy),'r','FaceAlpha',0.5)
        hold on;
        bar(bin_centres,T2_immob_occupancy/max(T2_immob_occupancy),'b','FaceAlpha',0.5)
        legend('Track Left','Track Right','Color','none')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        title('Occupancy map when speed < 5')


        subplot(3,4,3)
        for nprobe = 1:length(clusters)
            options = session_info(n).probe(nprobe);
            h(options.probe_hemisphere)= histogram(ripples_position{options.probe_hemisphere}(ripples_track_id{options.probe_hemisphere}==1),35,'FaceColor',probe_color{options.probe_hemisphere},...
                'Normalization','probability');hold on
        end
        if length(clusters) ==2
            legend(h,hemisphere_text,'Color','none')
        else
            legend(h,hemisphere_text(options.probe_hemisphere),'Color','none')
        end
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        title('Track Left ripple distribution')

        subplot(3,4,4)
        for nprobe = 1:length(clusters)
            options = session_info(n).probe(nprobe);
            h(options.probe_hemisphere)= histogram(ripples_position{options.probe_hemisphere}(ripples_track_id{options.probe_hemisphere}==2),35,'FaceColor',probe_color{options.probe_hemisphere},...
                'Normalization','probability');hold on
        end
        if length(clusters) ==2
            legend(h,hemisphere_text,'Color','none')
        else
            legend(h,hemisphere_text(options.probe_hemisphere),'Color','none')
        end
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        title('Track Right ripple distribution')

        cortex_LFP = [];
        CA1_LFP = [];
        
        % Grabbing L and R probe CA1 and cortex LFP during T1 and T2. (Probe_hemisphere as 1 and 2)
        for nprobe = 1:length(clusters)
            options = session_info(n).probe(nprobe);

            %             LFP_ripple_tidx = ismember(LFP(nprobe).tvec',Restrict(LFP(nprobe).tvec,[ripples(nprobe).onset; ripples(nprobe).offset]'));
            LFP_T1_tidx = ismember(LFP(nprobe).tvec',Restrict(LFP(nprobe).tvec',[Task_info.start_time_all(Task_info.track_ID_all==1) Task_info.end_time_all(Task_info.track_ID_all==1)']));
            LFP_T2_tidx = ismember(LFP(nprobe).tvec',Restrict(LFP(nprobe).tvec',[Task_info.start_time_all(Task_info.track_ID_all==2) Task_info.end_time_all(Task_info.track_ID_all==2)']));
            LFP_T1_tidx(speed_interp<5) = 0;
            LFP_T2_tidx(speed_interp<5) = 0;

            if isfield(LFP(nprobe),'L4')
                cortex_LFP{1}{options.probe_hemisphere} = LFP(nprobe).L4(LFP_T1_tidx);
                cortex_LFP{2}{options.probe_hemisphere} = LFP(nprobe).L4(LFP_T2_tidx);
            elseif isfield(LFP(nprobe),'L5')
                cortex_LFP{1}{options.probe_hemisphere} = LFP(nprobe).L5(LFP_T1_tidx);
                cortex_LFP{2}{options.probe_hemisphere} = LFP(nprobe).L5(LFP_T2_tidx);
            elseif isfield(LFP(nprobe),'MEC')

            end

            CA1_LFP{1}{options.probe_hemisphere} = LFP(nprobe).CA1(LFP_T1_tidx);
            CA1_LFP{2}{options.probe_hemisphere} = LFP(nprobe).CA1(LFP_T2_tidx);
        end
        
        if length(clusters)==2 & ~isempty(CA1_LFP{1}) & ~isempty(CA1_LFP{2}) % If L and R probes

            %%%%% Left and Right CA1 interaction on T1
            lfp.timestamps = LFP(1).tvec;
            lfp.data = [CA1_LFP{1}{1}' CA1_LFP{1}{2}'];
            [ comod ] = LFP_comodulogram(lfp,[],'winsize',2,'space','log')
            corrcolor= [makeColorMap([1 1 1],[0 0 0.8],[0 0 0]);...
                makeColorMap([0 0 0],[0.8 0 0],[1 1 1])];

            fig = figure;
            fig.Position = [160 50 1300 950];
            fig.Name = sprintf('Left CA1 vs Right CA1 LFP correlation')
            % colormap(corrcolor)
            imagesc(log2(comod.freqs),log2(comod.freqs),comod.corrs)
            colorbar
            ColorbarWithAxis([-0.4 0.4],'Power-Power Correlation (rho)')
            LogScale('xy',2)
            xlabel('f (Hz)');ylabel('f (Hz)')
            colormap('copper')
            %     colorbar
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            LFP_comod.freqs = comod.freqs;
            LFP_corr.track.CA1_corr_L{1} = comod.corrs;
            %     NiceSave(['Comodulogram',figparms.plotname],figparms.figfolder,figparms.baseName)

            %%%%% Left and Right CA1 interaction on T2
            lfp.timestamps = LFP(1).tvec;
            lfp.data = [CA1_LFP{1}{1}' CA1_LFP{1}{2}'];
            [ comod ] = LFP_comodulogram(lfp,[],'winsize',2,'space','log')
            corrcolor= [makeColorMap([1 1 1],[0 0 0.8],[0 0 0]);...
                makeColorMap([0 0 0],[0.8 0 0],[1 1 1])];

            fig = figure;
            fig.Position = [160 50 1300 950];
            fig.Name = sprintf('Left CA1 vs Right CA1 LFP correlation')
            % colormap(corrcolor)
            imagesc(log2(comod.freqs),log2(comod.freqs),comod.corrs)
            colorbar
            ColorbarWithAxis([-0.4 0.4],'Power-Power Correlation (rho)')
            LogScale('xy',2)
            xlabel('f (Hz)');ylabel('f (Hz)')
            colormap('copper')
            %     colorbar
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            LFP_corr.freqs = comod.freqs;
            LFP_corr.CA1_L_R{1} = comod.corrs;

            %%%%% Left CA1 and Right V1 interaction on T1
            lfp.timestamps = LFP(1).tvec;
            lfp.data = [CA1_LFP{1}{1}' cortex_LFP{1}{2}'];
            [ comod ] = LFP_comodulogram(lfp,[],'winsize',2,'space','log')

            % Left and Right CA1 interaction on T1
            lfp.timestamps = LFP(1).tvec;
            lfp.data = [CA1_LFP{1}{1}' cortex_LFP{1}{1}'];
            [ comod ] = LFP_comodulogram(lfp,[],'winsize',2,'space','log')

            % Left and Right CA1 interaction on T2
             lfp.timestamps = LFP(1).tvec;
            lfp.data = [cortex_LFP{1}{2}' CA1_LFP{1}{2}'];
            [ comod ] = LFP_comodulogram(lfp,[],'winsize',2,'space','log')

            % Left and Right CA1 interaction
        end


        if length(clusters)==2 & ~isempty(CA1_LFP{1}) & ~isempty(CA1_LFP{2}) % If L and R probes

        end


    end
end

%% Peri ripple LFP and CSD

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]

for nsession = [3 4 6 7 8 9 10 12 14]
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
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end

        clear lfpAvg csd
        lfpAvg_R = [];
        csd_R = [];
        lfpAvg_L = [];
        csd_L = [];

        for nprobe = 1:length(clusters)

            options = session_info(n).probe(nprobe);
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

            options.importMode = 'KS';
            [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,1);% Since it is LF

            shank_id_avaliable = ceil(best_channels{nprobe}.xcoord./250); % based on xcoord in best_channels
%             if unique(shank_id_avaliable)>1 
                normalisation = 0;
%             else
%                 normalisation = 1;
%             end

            x_col = find(shank_id_avaliable==1);
            x_col =x_col(1);

            clear lfpAvg csd
            % Shank 1 peri ripple
            for iprobe = 1:length(clusters)
                probe_hemisphere = session_info(n).probe(iprobe).probe_hemisphere;

                if session_info(n).probe(iprobe).probe_hemisphere == 1
                    [lfpAvg(iprobe).track(1),csd(iprobe).track(1)] = perievent_LFP_profile('L Ripple T1',ripples(iprobe).T1_onset,[-0.2 0.5],PSD,clusters,options,...
                        'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col,'CSD_V1_CA1_normalisation',normalisation);

                    [lfpAvg(iprobe).track(2),csd(iprobe).track(2)] = perievent_LFP_profile('L Ripple T2',ripples(iprobe).T2_onset,[-0.2 0.5],PSD,clusters,options,...
                        'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col,'CSD_V1_CA1_normalisation',normalisation);
            
                elseif session_info(n).probe(iprobe).probe_hemisphere == 2
                    [lfpAvg(iprobe).track(1),csd(iprobe).track(1)] = perievent_LFP_profile('R Ripple T1',ripples(iprobe).T1_onset,[-0.2 0.5],PSD,clusters,options,...
                        'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col,'CSD_V1_CA1_normalisation',normalisation);

                    [lfpAvg(iprobe).track(2),csd(iprobe).track(2)] = perievent_LFP_profile('R Ripple T2',ripples(iprobe).T2_onset,[-0.2 0.5],PSD,clusters,options,...
                        'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col,'CSD_V1_CA1_normalisation',normalisation);
                end

            end

            if session_info(n).probe(nprobe).probe_hemisphere == 1
                lfpAvg_L = lfpAvg;
                csd_L = csd;

            else
                lfpAvg_R = lfpAvg;
                csd_R = csd;
            end

            if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','ripples'))== 0
                mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','ripples'))
            end
        
            if unique(shank_id_avaliable)==1
                save(fullfile(options.ANALYSIS_DATAPATH,'ripple_LFP_response.mat'),'lfpAvg_L','lfpAvg_R','csd_L','csd_R');
            end

            save_all_large_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','ripples'),[])
            close all


            if unique(shank_id_avaliable)>1 % if more than 1 shank
                clear lfpAvg csd
                shank_id_avaliable = ceil(best_channels{nprobe}.xcoord./250); % based on xcoord in best_channels
                x_col = find(shank_id_avaliable==3);
                x_col =x_col(1);

                % Shank 3 peri ripple
                for iprobe = 1:length(clusters)
                    probe_hemisphere = session_info(n).probe(iprobe).probe_hemisphere;

                    if session_info(n).probe(iprobe).probe_hemisphere == 1
                        [lfpAvg(iprobe).track(1),csd(iprobe).track(1)] = perievent_LFP_profile('L Ripple T1',ripples(iprobe).T1_onset,[-0.2 0.5],PSD,clusters,options,...
                            'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col);

                        [lfpAvg(iprobe).track(2),csd(iprobe).track(2)] = perievent_LFP_profile('L Ripple T2',ripples(iprobe).T2_onset,[-0.2 0.5],PSD,clusters,options,...
                            'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col);

                    elseif session_info(n).probe(iprobe).probe_hemisphere == 2
                        [lfpAvg(iprobe).track(1),csd(iprobe).track(1)] = perievent_LFP_profile('R Ripple T1',ripples(iprobe).T1_onset,[-0.2 0.5],PSD,clusters,options,...
                            'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col);

                        [lfpAvg(iprobe).track(2),csd(iprobe).track(2)] = perievent_LFP_profile('R Ripple T2',ripples(iprobe).T2_onset,[-0.2 0.5],PSD,clusters,options,...
                            'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col);
                    end

                end

                if session_info(n).probe(nprobe).probe_hemisphere == 1
                    lfpAvg_L(3:4) = lfpAvg;
                    csd_L(3:4) = csd;

                else
                    lfpAvg_R(3:4) = lfpAvg;
                    csd_R(3:4) = csd;
                end

                save(fullfile(options.ANALYSIS_DATAPATH,'ripple_LFP_response.mat'),'lfpAvg_L','lfpAvg_R','csd_L','csd_R');

                save_all_large_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','ripples'),[])
                close all
            end
            
        end


    end
end