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
theta_modulation_R = [];
theta_modulation_L = [];
ripple_modulation_R = [];
ripple_modulation_L = [];
ripple_modulation_combined = [];
V1_event_modulation_L = [];
V1_event_modulation_R = [];
V1_event_modulation_combined = [];

place_fields_all = [];

for track_id = 1:2
    place_fields_all(track_id).session_id = [];
    place_fields_all(track_id).session = [];
    place_fields_all(track_id).subject = [];
    place_fields_all(track_id).cluster_id = [];
    place_fields_all(track_id).cell_type = [];
    place_fields_all(track_id).region = [];
    %     place_fields_all(track_id).peak_channel_waveform = [];
    %     place_fields_all(track_id).dwell_map = [];
    place_fields_all(track_id).relative_depth = [];


    place_fields_all(track_id).x_bin_edges = [];
    place_fields_all(track_id).x_bin_centres =[];
    place_fields_all(track_id).average_map = [];
    place_fields_all(track_id).raw = [];
    place_fields_all(track_id).odd_lap_map = [];
    place_fields_all(track_id).even_lap_map = [];

    place_fields_all(track_id).odd_even_corr = [];
    place_fields_all(track_id).odd_even_stability = [];

    place_fields_all(track_id).t1_t2_corr = [];
    place_fields_all(track_id).t1_t2_remapping = [];

    place_fields_all(track_id).ripple_spike_count = [];
    place_fields_all(track_id).ripple_PSTH = [];
    place_fields_all(track_id).ripple_PSTH_zscored = [];
    place_fields_all(track_id).ripple_modulation_percentile = [];
    place_fields_all(track_id).pre_ripple_activation = [];
    place_fields_all(track_id).pre_ripple_inhibition = [];
    place_fields_all(track_id).post_ripple_activation = [];
    place_fields_all(track_id).post_ripple_inhibition = [];

    place_fields_all(track_id).V1_event_spike_count = [];
    place_fields_all(track_id).V1_event_PSTH = [];
    place_fields_all(track_id).V1_event_PSTH_zscored = [];
    place_fields_all(track_id).V1_event_modulation_percentile = [];
    place_fields_all(track_id).pre_V1_event_activation = [];
    place_fields_all(track_id).pre_V1_event_inhibition = [];
    place_fields_all(track_id).post_V1_event_activation = [];
    place_fields_all(track_id).post_V1_event_inhibition = [];

    place_fields_all(track_id).theta_modulation_percentile = [];
    place_fields_all(track_id).theta_phase_map = [];
    place_fields_all(track_id).position_phase_xcorr_map = [];
    place_fields_all(track_id).position_phase_map = [];

    place_fields_all(track_id).spatial_cell_correlation = [];
    place_fields_all(track_id).ripple_cell_correlation = [];

    place_fields_all(track_id).SMI = [];
    place_fields_all(track_id).lap_SMI = [];
    place_fields_all(track_id).within_block_SMI = [];
end

place_fields_all_L = place_fields_all;
place_fields_all_R = place_fields_all;
place_fields_all_combined = place_fields_all;

for nsession =[1 2 3 4 6 7 8 9 10 12 14]
    tic
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
            load(fullfile(options.ANALYSIS_DATAPATH,'V1_event_modulation.mat'));
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
        
        for track_id = 1:2
            place_fields_all_combined(track_id).session_id = [place_fields_all_combined(track_id).session_id nsession*ones(1,length(good_cell_index))];
            place_fields_all_combined(track_id).subject{nsession} = [options.SUBJECT];
            place_fields_all_combined(track_id).session{nsession} = [options.SESSION];

            ratemap_matrix = [place_fields(track_id).raw{good_cell_index}];
            ratemap_matrix = reshape(ratemap_matrix,size(place_fields(track_id).raw{1},1),[],length(good_cell_index));%laps X position bins X cells

            all_fields = fieldnames(place_fields_all_combined);
            all_fields(strcmp(all_fields,'raw'))=[];

            pf_fields = fieldnames(place_fields);
            ripple_fields = fieldnames(ripple_modulation_combined);
%             theta_fields = fieldnames(theta_modulation_L);
            V1_event_fields = fieldnames(V1_event_modulation_combined);

            for iField = 1:length(all_fields)

                if sum(strcmp(ripple_fields,all_fields{iField})) == 1
                    place_fields_all_combined(track_id).(all_fields{iField}) = [place_fields_all_combined(track_id).(all_fields{iField}) ...
                        ripple_modulation_combined(track_id).(ripple_fields{strcmp(ripple_fields,all_fields{iField})})];
                    continue
                end
% 
%                 if sum(strcmp(theta_fields,all_fields{iField})) == 1
%                     place_fields_all_combined(track_id).(all_fields{iField}) = [place_fields_all_combined(track_id).(all_fields{iField}) ...
%                         theta_modulation_L(track_id).(theta_fields{strcmp(theta_fields,all_fields{iField})})];
%                     continue
%                 end
                if sum(strcmp(V1_event_fields,all_fields{iField})) == 1
                    place_fields_all_combined(track_id).(all_fields{iField}) = [place_fields_all_combined(track_id).(all_fields{iField}) ...
                        V1_event_modulation_combined(track_id).(V1_event_fields{strcmp(V1_event_fields,all_fields{iField})})];
                    continue
                end

                if sum(strcmp(pf_fields,all_fields{iField})) == 1
                    if contains(all_fields{iField},'x_bin')==1
                        place_fields_all_combined(track_id).(all_fields{iField}) = [place_fields_all_combined(track_id).(all_fields{iField}) ...
                            place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})})];
                    elseif contains(all_fields{iField},'lap_SMI')==1
                        temp = place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})})(good_cell_index,:)';
                        for t = 1:length(good_cell_index)
                            place_fields_all_combined(track_id).(all_fields{iField}) = [place_fields_all_combined(track_id).(all_fields{iField}) ...
                                {temp(:,t)}];
                        end
                    else
                        if size(place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})}),1)~=1
                            temp = place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})})(good_cell_index,:)';
                        else
                            temp = place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})})(:,good_cell_index);
                        end

                        place_fields_all_combined(track_id).(all_fields{iField}) = [place_fields_all_combined(track_id).(all_fields{iField}) ...
                            temp];
                    end
                    continue
                end

                if strcmp(all_fields{iField},'ripple_PSTH')
                    PSTH_matrix = [ripple_modulation_combined(track_id).PSTH{:}];
                    PSTH_matrix = reshape(PSTH_matrix,[],length(good_cell_index));
                    place_fields_all_combined(track_id).ripple_PSTH = [place_fields_all_combined(track_id).ripple_PSTH PSTH_matrix];
                end

                if strcmp(all_fields{iField},'ripple_PSTH_zscored')
                    PSTH_matrix = [ripple_modulation_combined(track_id).PSTH_zscored{:}];
                    PSTH_matrix = reshape(PSTH_matrix,[],length(good_cell_index));
                    place_fields_all_combined(track_id).ripple_PSTH_zscored = [place_fields_all_combined(track_id).ripple_PSTH_zscored PSTH_matrix];
                end

                if strcmp(all_fields{iField},'ripple_spike_count')
                    place_fields_all_combined(track_id).ripple_spike_count = [place_fields_all_combined(track_id).ripple_spike_count ripple_modulation_combined(track_id).spike_count];
                end

                if strcmp(all_fields{iField},'V1_event_PSTH')
                    PSTH_matrix = [ripple_modulation_combined(track_id).PSTH{:}];
                    PSTH_matrix = reshape(PSTH_matrix,[],length(good_cell_index));
                    place_fields_all_combined(track_id).V1_event_PSTH = [place_fields_all_combined(track_id).V1_event_PSTH PSTH_matrix];
                end

                if strcmp(all_fields{iField},'V1_PSTH_zscored')
                    PSTH_matrix = [ripple_modulation_combined(track_id).PSTH_zscored{:}];
                    PSTH_matrix = reshape(PSTH_matrix,[],length(good_cell_index));
                    place_fields_all_combined(track_id).V1_event_PSTH_zscored = [place_fields_all_combined(track_id).V1_event_PSTH_zscored PSTH_matrix];
                end

                if strcmp(all_fields{iField},'V1_event_spike_count')
                    place_fields_all_combined(track_id).V1_event_spike_count = [place_fields_all_combined(track_id).V1_event_spike_count ripple_modulation_combined(track_id).spike_count];
                end
            end

            % Relative distance from the surface for clusters from
            % different shanks (mainly for cortex)
            relative_depth = [];

            for iprobe=1:length(best_channels)
                options = session_info(n).probe(iprobe);
                options.importMode = 'KS';
                [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(options,[]);

                if iprobe == 1
                    good_cell_index_this_probe = good_cell_index(place_fields(track_id).cluster_id(good_cell_index)<10000);
                elseif iprobe == 2
                    good_cell_index_this_probe = good_cell_index(place_fields(track_id).cluster_id(good_cell_index)>10000);
                end

                [Lia,Locb] = ismember(place_fields(track_id).peak_channel(good_cell_index_this_probe),chan_config.Channel);

                shank_id_avaliable = ceil(best_channels{iprobe}.xcoord./250); % based on xcoord in best_channels
                shanks_avaliable = unique(shank_id_avaliable);

                for nshank = shanks_avaliable
                    surface_depth(nshank) = median(best_channels{iprobe}.surface_depth(shank_id_avaliable == shanks_avaliable(nshank)),'omitnan');
                end


                [Li1,Loc2] = ismember(place_fields(track_id).cluster_id(good_cell_index_this_probe),ripple_modulation_combined(track_id).cluster_id);

                relative_depth = [relative_depth surface_depth(chan_config.Shank(Locb))'-ripple_modulation_combined(track_id).peak_depth(Loc2)];
            end

            place_fields_all_combined(track_id).relative_depth = [place_fields_all_combined(track_id).relative_depth  ...
                relative_depth];

            ratemap_matrix = [place_fields(track_id).raw{good_cell_index}];
            ratemap_matrix = reshape(ratemap_matrix,size(place_fields(track_id).raw{1},1),[],length(good_cell_index));%laps X position bins X cells
            place_fields_all_combined(track_id).raw = [place_fields_all_combined(track_id).raw;place_fields(track_id).raw(good_cell_index)];
            place_fields_all_combined(track_id).average_map = [place_fields_all_combined(track_id).average_map squeeze(mean(ratemap_matrix,1))];

        end

        % Per probe
        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);

            if options.probe_hemisphere == 1
                for track_id = 1:2
                    place_fields_all_L(track_id).session_id = [place_fields_all_L(track_id).session_id nsession*ones(1,length(good_cell_index))];
                    place_fields_all_L(track_id).subject{nsession} = [options.SUBJECT];
                    place_fields_all_L(track_id).session{nsession} = [options.SESSION];

                    ratemap_matrix = [place_fields(track_id).raw{good_cell_index}];
                    ratemap_matrix = reshape(ratemap_matrix,size(place_fields(track_id).raw{1},1),[],length(good_cell_index));%laps X position bins X cells

                    all_fields = fieldnames(place_fields_all_L);
                    all_fields(strcmp(all_fields,'raw'))=[];

                    pf_fields = fieldnames(place_fields);
                    ripple_fields = fieldnames(ripple_modulation_L);
                    theta_fields = fieldnames(theta_modulation_L);
                    V1_event_fields = fieldnames(V1_event_modulation_L);

                    for iField = 1:length(all_fields)

                        if sum(strcmp(ripple_fields,all_fields{iField})) == 1
                            place_fields_all_L(track_id).(all_fields{iField}) = [place_fields_all_L(track_id).(all_fields{iField}) ...
                                ripple_modulation_L(track_id).(ripple_fields{strcmp(ripple_fields,all_fields{iField})})];
                            continue
                        end

                        if sum(strcmp(theta_fields,all_fields{iField})) == 1
                            place_fields_all_L(track_id).(all_fields{iField}) = [place_fields_all_L(track_id).(all_fields{iField}) ...
                                theta_modulation_L(track_id).(theta_fields{strcmp(theta_fields,all_fields{iField})})];
                            continue
                        end

                        if sum(strcmp(V1_event_fields,all_fields{iField})) == 1
                            place_fields_all_L(track_id).(all_fields{iField}) = [place_fields_all_L(track_id).(all_fields{iField}) ...
                                V1_event_modulation_L(track_id).(V1_event_fields{strcmp(V1_event_fields,all_fields{iField})})];
                            continue
                        end

                        if sum(strcmp(pf_fields,all_fields{iField})) == 1
                            if contains(all_fields{iField},'x_bin')==1
                                place_fields_all_L(track_id).(all_fields{iField}) = [place_fields_all_L(track_id).(all_fields{iField}) ...
                                    place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})})];
                            elseif contains(all_fields{iField},'lap_SMI')==1
                                temp = place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})})(good_cell_index,:)';
                                for t = 1:length(good_cell_index)
                                    place_fields_all_L(track_id).(all_fields{iField}) = [place_fields_all_L(track_id).(all_fields{iField}) ...
                                        {temp(:,t)}];
                                end
                            else
                                if size(place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})}),1)~=1
                                    temp = place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})})(good_cell_index,:)';
                                else
                                    temp = place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})})(:,good_cell_index);
                                end

                                place_fields_all_L(track_id).(all_fields{iField}) = [place_fields_all_L(track_id).(all_fields{iField}) ...
                                    temp];
                            end
                            continue
                        end

                        if strcmp(all_fields{iField},'ripple_PSTH')
                            PSTH_matrix = [ripple_modulation_L(track_id).PSTH{:}];
                            PSTH_matrix = reshape(PSTH_matrix,[],length(good_cell_index));
                            place_fields_all_L(track_id).ripple_PSTH = [place_fields_all_L(track_id).ripple_PSTH PSTH_matrix];
                        end

                        if strcmp(all_fields{iField},'ripple_PSTH_zscored')
                            PSTH_matrix = [ripple_modulation_L(track_id).PSTH_zscored{:}];
                            PSTH_matrix = reshape(PSTH_matrix,[],length(good_cell_index));
                            place_fields_all_L(track_id).ripple_PSTH_zscored = [place_fields_all_L(track_id).ripple_PSTH_zscored PSTH_matrix];
                        end

                        if strcmp(all_fields{iField},'ripple_spike_count')
                            place_fields_all_L(track_id).ripple_spike_count = [place_fields_all_L(track_id).ripple_spike_count ripple_modulation_L(track_id).spike_count];
                        end

                        if strcmp(all_fields{iField},'V1_event_PSTH')
                            PSTH_matrix = [V1_event_modulation_L(track_id).PSTH{:}];
                            PSTH_matrix = reshape(PSTH_matrix,[],length(good_cell_index));
                            place_fields_all_L(track_id).V1_event_PSTH = [place_fields_all_L(track_id).V1_event_PSTH PSTH_matrix];
                        end

                        if strcmp(all_fields{iField},'V1_PSTH_zscored')
                            PSTH_matrix = [V1_event_modulation_L(track_id).PSTH_zscored{:}];
                            PSTH_matrix = reshape(PSTH_matrix,[],length(good_cell_index));
                            place_fields_all_L(track_id).V1_event_PSTH_zscored = [place_fields_all_L(track_id).V1_event_PSTH_zscored PSTH_matrix];
                        end

                        if strcmp(all_fields{iField},'V1_event_spike_count')
                            place_fields_all_L(track_id).V1_event_spike_count = [place_fields_all_L(track_id).V1_event_spike_count V1_event_modulation_L(track_id).spike_count];
                        end
                    end

                    % Relative distance from the surface for clusters from
                    % different shanks (mainly for cortex)
                    relative_depth = [];
                    
                    for iprobe=1:length(best_channels)
                        options = session_info(n).probe(iprobe);
                        options.importMode = 'KS';
                        [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(options,[]);

                        if iprobe == 1
                            good_cell_index_this_probe = good_cell_index(place_fields(track_id).cluster_id(good_cell_index)<10000);
                        elseif iprobe == 2
                            good_cell_index_this_probe = good_cell_index(place_fields(track_id).cluster_id(good_cell_index)>10000);
                        end
                        
                        [Lia,Locb] = ismember(place_fields(track_id).peak_channel(good_cell_index_this_probe),chan_config.Channel);

                        shank_id_avaliable = ceil(best_channels{iprobe}.xcoord./250); % based on xcoord in best_channels
                        shanks_avaliable = unique(shank_id_avaliable);

                        for nshank = shanks_avaliable
                            surface_depth(nshank) = median(best_channels{iprobe}.surface_depth(shank_id_avaliable == shanks_avaliable(nshank)),'omitnan');
                        end


                        [Li1,Loc2] = ismember(place_fields(track_id).cluster_id(good_cell_index_this_probe),ripple_modulation_L(track_id).cluster_id);

                        relative_depth = [relative_depth surface_depth(chan_config.Shank(Locb))'-ripple_modulation_L(track_id).peak_depth(Loc2)];
                    end

                    place_fields_all_L(track_id).relative_depth = [place_fields_all_L(track_id).relative_depth  ...
                        relative_depth];

                    ratemap_matrix = [place_fields(track_id).raw{good_cell_index}];
                    ratemap_matrix = reshape(ratemap_matrix,size(place_fields(track_id).raw{1},1),[],length(good_cell_index));%laps X position bins X cells
                    place_fields_all_L(track_id).raw = [place_fields_all_L(track_id).raw;place_fields(track_id).raw(good_cell_index)];
                    place_fields_all_L(track_id).average_map = [place_fields_all_L(track_id).average_map squeeze(mean(ratemap_matrix,1))];

                end

            elseif options.probe_hemisphere == 2
                for track_id = 1:2
                    place_fields_all_R(track_id).session_id = [place_fields_all_R(track_id).session_id nsession*ones(1,length(good_cell_index))];
                    place_fields_all_R(track_id).subject{nsession} = [options.SUBJECT];
                    place_fields_all_R(track_id).session{nsession} = [options.SESSION];

                    all_fields = fieldnames(place_fields_all_R);
                    all_fields(strcmp(all_fields,'raw'))=[];

                    pf_fields = fieldnames(place_fields);
                    ripple_fields = fieldnames(ripple_modulation_R);
                    theta_fields = fieldnames(theta_modulation_R);
                    V1_event_fields = fieldnames(V1_event_modulation_R);

                    for iField = 1:length(all_fields)

                        if sum(strcmp(ripple_fields,all_fields{iField})) == 1
                            place_fields_all_R(track_id).(all_fields{iField}) = [place_fields_all_R(track_id).(all_fields{iField}) ...
                                ripple_modulation_R(track_id).(ripple_fields{strcmp(ripple_fields,all_fields{iField})})];
                            continue
                        end

                        if sum(strcmp(theta_fields,all_fields{iField})) == 1
                            place_fields_all_R(track_id).(all_fields{iField}) = [place_fields_all_R(track_id).(all_fields{iField}) ...
                                theta_modulation_R(track_id).(theta_fields{strcmp(theta_fields,all_fields{iField})})];
                            continue
                        end

                        if sum(strcmp(V1_event_fields,all_fields{iField})) == 1
                            place_fields_all_R(track_id).(all_fields{iField}) = [place_fields_all_R(track_id).(all_fields{iField}) ...
                                V1_event_modulation_R(track_id).(V1_event_fields{strcmp(V1_event_fields,all_fields{iField})})];
                            continue
                        end


                        if sum(strcmp(pf_fields,all_fields{iField})) == 1
                            if contains(all_fields{iField},'x_bin')==1
                                place_fields_all_R(track_id).(all_fields{iField}) = place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})});
                            elseif contains(all_fields{iField},'lap_SMI')==1
                                temp = place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})})(good_cell_index,:)';
                                for t = 1:length(good_cell_index)
                                    place_fields_all_R(track_id).(all_fields{iField}) = [place_fields_all_R(track_id).(all_fields{iField}) ...
                                        {temp(:,t)}];
                                end
                            else
                                if size(place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})}),1)~=1
                                    temp = place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})})(good_cell_index,:)';
                                else
                                    temp = place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})})(:,good_cell_index);
                                end

                                place_fields_all_R(track_id).(all_fields{iField}) = [place_fields_all_R(track_id).(all_fields{iField}) ...
                                    temp];
                            end

                            continue
                        end

                        if strcmp(all_fields{iField},'ripple_PSTH')
                            PSTH_matrix = [ripple_modulation_R(track_id).PSTH{:}];
                            PSTH_matrix = reshape(PSTH_matrix,[],length(good_cell_index));
                            place_fields_all_R(track_id).ripple_PSTH = [place_fields_all_R(track_id).ripple_PSTH PSTH_matrix];
                        end

                        if strcmp(all_fields{iField},'ripple_PSTH_zscored')
                            PSTH_matrix = [ripple_modulation_R(track_id).PSTH_zscored{:}];
                            PSTH_matrix = reshape(PSTH_matrix,[],length(good_cell_index));
                            place_fields_all_R(track_id).ripple_PSTH_zscored = [place_fields_all_R(track_id).ripple_PSTH_zscored PSTH_matrix];
                        end
                        if strcmp(all_fields{iField},'ripple_spike_count')
                            place_fields_all_R(track_id).ripple_spike_count = [place_fields_all_R(track_id).ripple_spike_count ripple_modulation_R(track_id).spike_count];
                        end

                        if strcmp(all_fields{iField},'V1_event_PSTH')
                            PSTH_matrix = [V1_event_modulation_R(track_id).PSTH{:}];
                            PSTH_matrix = reshape(PSTH_matrix,[],length(good_cell_index));
                            place_fields_all_R(track_id).V1_event_PSTH = [place_fields_all_R(track_id).V1_event_PSTH PSTH_matrix];
                        end

                        if strcmp(all_fields{iField},'V1_event_PSTH_zscored')
                            PSTH_matrix = [V1_event_modulation_R(track_id).PSTH_zscored{:}];
                            PSTH_matrix = reshape(PSTH_matrix,[],length(good_cell_index));
                            place_fields_all_R(track_id).V1_event_PSTH_zscored = [place_fields_all_R(track_id).V1_event_PSTH_zscored PSTH_matrix];
                        end

                        if strcmp(all_fields{iField},'V1_event_spike_count')
                            place_fields_all_R(track_id).V1_event_spike_count = [place_fields_all_R(track_id).V1_event_spike_count V1_event_modulation_R(track_id).spike_count];
                        end
                    end

                    % Relative distance from the surface for clusters from
                    % different shanks (mainly for cortex)
                    relative_depth = [];

                    for iprobe=1:length(best_channels)
                        options = session_info(n).probe(iprobe);
                        options.importMode = 'KS';
                        [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(options,[]);

                        if iprobe == 1
                            good_cell_index_this_probe = good_cell_index(place_fields(track_id).cluster_id(good_cell_index)<10000);
                        elseif iprobe == 2
                            good_cell_index_this_probe = good_cell_index(place_fields(track_id).cluster_id(good_cell_index)>10000);
                        end

                        [Lia,Locb] = ismember(place_fields(track_id).peak_channel(good_cell_index_this_probe),chan_config.Channel);

                        shank_id_avaliable = ceil(best_channels{iprobe}.xcoord./250); % based on xcoord in best_channels
                        shanks_avaliable = unique(shank_id_avaliable);

                        for nshank = shanks_avaliable
                            surface_depth(nshank) = median(best_channels{iprobe}.surface_depth(shank_id_avaliable == shanks_avaliable(nshank)),'omitnan');
                        end


                        [Li1,Loc2] = ismember(place_fields(track_id).cluster_id(good_cell_index_this_probe),ripple_modulation_R(track_id).cluster_id);

                        relative_depth = [relative_depth surface_depth(chan_config.Shank(Locb))'-ripple_modulation_R(track_id).peak_depth(Loc2)];
                    end

                    place_fields_all_R(track_id).relative_depth = [place_fields_all_R(track_id).relative_depth  ...
                        relative_depth];


                    ratemap_matrix = [place_fields(track_id).raw{good_cell_index}];
                    ratemap_matrix = reshape(ratemap_matrix,size(place_fields(track_id).raw{1},1),[],length(good_cell_index));%laps X position bins X cells
                    place_fields_all_R(track_id).raw = [place_fields_all_R(track_id).raw;place_fields(track_id).raw(good_cell_index)];
                    place_fields_all_R(track_id).average_map = [place_fields_all_R(track_id).average_map squeeze(mean(ratemap_matrix,1))];

                    % spatial modulation
                end

            end

        end

        % Plot individual cells
    end
    toc
end

save(fullfile('D:\corticohippocampal_replay\summary','place_fields_all.mat'),'place_fields_all_L','place_fields_all_R','place_fields_all_combined')
% 
% save(fullfile(options.ANALYSIS_DATAPATH,'..','..','figures','population_map'))
% save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','..','figures','population_map'))

%% V1 poulation bursting number
clear all 
load(fullfile('D:\corticohippocampal_replay\summary','place_fields_all.mat'),'place_fields_all_L','place_fields_all_R','place_fields_all_combined')
counter = 1;
V1_events_no_L = [];
V1_events_no_R = [];
for nsession = unique(place_fields_all_R(1).session_id)

    for track_id = 1:2
        V1_events_no_L(counter,track_id) = nan;
        V1_events_no_R(counter,track_id) = nan;

        if sum(place_fields_all_L(1).session_id==nsession)>0
            this_session = find(place_fields_all_L(1).session_id==nsession);
            this_session = this_session(1);
            V1_events_no_L(counter,track_id) = size(place_fields_all_L(track_id).V1_event_spike_count{this_session},1);
        end

        if sum(place_fields_all_R(1).session_id==nsession)>0
            this_session = find(place_fields_all_R(1).session_id==nsession);
            this_session = this_session(1);
            V1_events_no_R(counter,track_id) = size(place_fields_all_R(track_id).V1_event_spike_count{this_session},1);
        end
    end
    counter = counter + 1;
end

subplot(2,2,1)
scatter(ones(1,length(V1_events_no_L(:,1))),V1_events_no_L(:,1),'blue')
hold on
scatter(2*ones(1,length(V1_events_no_L(:,2))),V1_events_no_L(:,2),'red')

for i = 1:length(V1_events_no_L(:,1))
    plot([ones(1,length(V1_events_no_L(i,1))); 2*ones(1,length(V1_events_no_L(i,2)))],[V1_events_no_L(i,1)';V1_events_no_L(i,2)'],'k')
end

scatter(3*ones(1,length(V1_events_no_R(:,1))),V1_events_no_R(:,1),'blue')

scatter(4*ones(1,length(V1_events_no_R(:,2))),V1_events_no_R(:,2),'red')
for i = 1:length(V1_events_no_R(:,1))
    plot([3*ones(1,length(V1_events_no_R(i,1))); 4*ones(1,length(V1_events_no_R(i,2)))],[V1_events_no_R(i,1)';V1_events_no_R(i,2)'],'k')
end
xticks(1:4)
xticklabels({'Left V1 event Track Left','Left V1 event Track Right','Right V1 event Track Left','Right V1 event Track Right'})
ylabel('No of V1 events')

subplot(2,2,3)
scatter(ones(1,length(V1_events_no_L(:,1))),V1_events_no_L(:,1),'blue')
hold on
scatter(3*ones(1,length(V1_events_no_L(:,2))),V1_events_no_L(:,2),'red')

for i = 1:length(V1_events_no_L(:,1))
    plot([ones(1,length(V1_events_no_L(i,1))); 2*ones(1,length(V1_events_no_L(i,2)))],[V1_events_no_L(i,1)';V1_events_no_R(i,1)'],'k')
end

scatter(2*ones(1,length(V1_events_no_R(:,1))),V1_events_no_R(:,1),'blue')

scatter(4*ones(1,length(V1_events_no_R(:,2))),V1_events_no_R(:,2),'red')
for i = 1:length(V1_events_no_R(:,1))
    plot([3*ones(1,length(V1_events_no_R(i,1))); 4*ones(1,length(V1_events_no_R(i,2)))],[V1_events_no_L(i,2)';V1_events_no_R(i,2)'],'k')
end
xticks(1:4)
xticklabels({'Left V1 event Track Left','Right V1 event Track Left','Left V1 event Track Right','Right V1 event Track Right'})
ylabel('No of V1 events')


%% ripple number
counter = 1;
ripple_no_L = [];
ripple_no_R = [];
for nsession = unique(place_fields_all_R(1).session_id)

    for track_id = 1:2
        ripple_no_L(counter,track_id) = nan;
        ripple_no_R(counter,track_id) = nan;

        if sum(place_fields_all_L(1).session_id==nsession)>0
            this_session = find(place_fields_all_L(1).session_id==nsession);
            this_session = this_session(1);
            ripple_no_L(counter,track_id) = size(place_fields_all_L(track_id).ripple_spike_count{this_session},1);
        end

        if sum(place_fields_all_R(1).session_id==nsession)>0
            this_session = find(place_fields_all_R(1).session_id==nsession);
            this_session = this_session(1);
            ripple_no_R(counter,track_id) = size(place_fields_all_R(track_id).ripple_spike_count{this_session},1);
        end
    end
    counter = counter + 1;
end

subplot(2,2,1)
scatter(ones(1,length(ripple_no_L(:,1))),ripple_no_L(:,1),'blue')
hold on
scatter(2*ones(1,length(ripple_no_L(:,2))),ripple_no_L(:,2),'red')

for i = 1:length(ripple_no_L(:,1))
    plot([ones(1,length(ripple_no_L(i,1))); 2*ones(1,length(ripple_no_L(i,2)))],[ripple_no_L(i,1)';ripple_no_L(i,2)'],'k')
end

scatter(3*ones(1,length(ripple_no_R(:,1))),ripple_no_R(:,1),'blue')

scatter(4*ones(1,length(ripple_no_R(:,2))),ripple_no_R(:,2),'red')
for i = 1:length(ripple_no_R(:,1))
    plot([3*ones(1,length(ripple_no_R(i,1))); 4*ones(1,length(ripple_no_R(i,2)))],[ripple_no_R(i,1)';ripple_no_R(i,2)'],'k')
end
xticks(1:4)
xticklabels({'Left ripple Track Left','Left ripple Track Right','Right ripple Track Left','Right ripple Track Right'})
ylabel('No of ripples')

subplot(2,2,3)
scatter(ones(1,length(ripple_no_L(:,1))),ripple_no_L(:,1),'blue')
hold on
scatter(3*ones(1,length(ripple_no_L(:,2))),ripple_no_L(:,2),'red')

for i = 1:length(ripple_no_L(:,1))
    plot([ones(1,length(ripple_no_L(i,1))); 2*ones(1,length(ripple_no_L(i,2)))],[ripple_no_L(i,1)';ripple_no_R(i,1)'],'k')
end

scatter(2*ones(1,length(ripple_no_R(:,1))),ripple_no_R(:,1),'blue')

scatter(4*ones(1,length(ripple_no_R(:,2))),ripple_no_R(:,2),'red')
for i = 1:length(ripple_no_R(:,1))
    plot([3*ones(1,length(ripple_no_R(i,1))); 4*ones(1,length(ripple_no_R(i,2)))],[ripple_no_L(i,2)';ripple_no_R(i,2)'],'k')
end
xticks(1:4)
xticklabels({'Left ripple Track Left','Right ripple Track Left','Left ripple Track Right','Right ripple Track Right'})
ylabel('No of ripples')

%% ripple combined
clear all
load(fullfile('D:\corticohippocampal_replay\summary','place_fields_all.mat'))
% V1 L cells
selected_cells_R = [];
selected_cells_L = [];

selected_cells = unique([find(contains(place_fields_all_combined(1).region,'V1') & ...
    place_fields_all_combined(1).ripple_modulation_percentile>0.99) find(contains(place_fields_all_combined(1).region,'V1') & ...
    place_fields_all_combined(2).ripple_modulation_percentile>0.99)]);

x_edges = -2:0.02:2;
x_bins = x_edges(2:end)-diff(x_edges)/2;
[coeff,score,latent,tsquared,explained]=pca(normalize(place_fields_all_combined(1).ripple_PSTH(x_bins>-1&x_bins<1,selected_cells),'range'));

plot(cumsum(explained))
xlabel('number of components')
ylabel('cumulative explained variance')
yline(90)
min_components = find(cumsum(explained)>=90);
min_components = min_components(1);

PCA_data = coeff(:,1:min_components);
max_clusters = min_components;
co_clustering = cell(max_clusters, 1);
total_prob = zeros(max_clusters, 1);
sumd_all = zeros(max_clusters, 100);

for k = 2:max_clusters
    co_clustering{k} = zeros(size(PCA_data, 1));
    for iter = 1:100
        rng(iter+k*1000)
        [idx,~,sumd,~] = kmeans(PCA_data,k);
        sumd_all(k,iter) = mean(sumd);
        for i = 1:size(PCA_data, 1)
            for j = 1:size(PCA_data, 1)
                if idx(i) == idx(j)
                    co_clustering{k}(i, j) = co_clustering{k}(i, j) + 1;
                end
            end
        end
    end
    co_clustering{k} = co_clustering{k} / 100; % Convert to probability
    mean_prob(k) = mean(co_clustering{k}(:)); % Sum all probabilities

    [~,index] = sort(idx);

    subplot(2,5,k)


    %     imagesc(co_clustering{k}(index,index))
%     imagesc(co_clustering{k}(index,index))
end

[~, optimal_clusters] = max(mean_prob); % Find the number of clusters with the highest sum of probabilities


selected_cells_L = unique([find(contains(place_fields_all_combined(1).region,'V1_L'))]);


x_edges = -2:0.02:2;
x_bins = x_edges(2:end)-diff(x_edges)/2;

sum_sqr_diff = sum(place_fields_all_combined(2).ripple_PSTH_zscored(x_bins>-0.5&x_bins<0.5,selected_cells_L)-place_fields_all_combined(1).ripple_PSTH_zscored(x_bins>-0.5&x_bins<0.5,selected_cells_L)).^2;
[~,index] = sort(sum_sqr_diff);


SMI = place_fields_all_combined(2).SMI(selected_cells_L);
[~,index] = sort(SMI);
selected_cells_L = selected_cells_L(index);
selected_cells_R = selected_cells_L;


plot_spatial_theta_ripple_population('All V1 L cells',selected_cells_L,selected_cells_R,place_fields_all_combined,place_fields_all_combined)


%%%%%%%%% V1 R
selected_cells_R = [];
selected_cells_L = [];

% selected_cells_L = unique([find(contains(place_fields_all_L(1).region,'V1_L') & ...
%     place_fields_all_L(1).ripple_modulation_percentile>0.95) find(contains(place_fields_all_L(1).region,'V1_L') & ...
%     place_fields_all_L(2).ripple_modulation_percentile>0.95)]);
selected_cells_L = unique([find(contains(place_fields_all_combined(1).region,'V1_R'))]);


x_edges = -2:0.02:2;
x_bins = x_edges(2:end)-diff(x_edges)/2;

sum_sqr_diff = sum(place_fields_all_combined(2).ripple_PSTH_zscored(x_bins>-0.5&x_bins<0.5,selected_cells_L)-place_fields_all_combined(1).ripple_PSTH_zscored(x_bins>-0.5&x_bins<0.5,selected_cells_L)).^2;
[~,index] = sort(sum_sqr_diff);
selected_cells_L = selected_cells_L(index);
selected_cells_R = selected_cells_L;


plot_spatial_theta_ripple_population('All V1 R cells',selected_cells_L,selected_cells_R,place_fields_all_combined,place_fields_all_combined)


%% ripple by magnitude 
clear all
load(fullfile('D:\corticohippocampal_replay\summary','place_fields_all.mat'))
% V1 L cells
selected_cells_R = [];
selected_cells_L = [];

% selected_cells_L = unique([find(contains(place_fields_all_L(1).region,'V1_L') & ...
%     place_fields_all_L(1).ripple_modulation_percentile>0.95) find(contains(place_fields_all_L(1).region,'V1_L') & ...
%     place_fields_all_L(2).ripple_modulation_percentile>0.95)]);

selected_cells_L = unique([find(contains(place_fields_all_L(1).region,'V1_L'))]);


x_edges = -2:0.02:2;
x_bins = x_edges(2:end)-diff(x_edges)/2;

sum_sqr_diff = sum(place_fields_all_L(2).ripple_PSTH_zscored(x_bins>-1&x_bins<1,selected_cells_L)-place_fields_all_L(1).ripple_PSTH_zscored(x_bins>-1&x_bins<1,selected_cells_L)).^2;
% [~,index] = sort(sum_sqr_diff);
% selected_cells_L = selected_cells_L(index);

SMI = place_fields_all_L(2).SMI(selected_cells_L);
[~,index] = sort(SMI);
selected_cells_L = selected_cells_L(index);



cluster_id = place_fields_all_L(1).cluster_id(selected_cells_L);
session_id = place_fields_all_L(1).session_id(selected_cells_L);
% 
counter = 1;
for ncell = 1:length(cluster_id)
    if sum(place_fields_all_R(1).cluster_id==cluster_id(ncell)&place_fields_all_R(1).session_id==session_id(ncell))>0
        selected_cells_R(counter) = find(place_fields_all_R(1).cluster_id==cluster_id(ncell)&place_fields_all_R(1).session_id==session_id(ncell))
        counter = counter + 1;
    end
end


plot_spatial_theta_ripple_population('All V1 L cells',selected_cells_L,selected_cells_R,place_fields_all_L,place_fields_all_R)


%%%%%%%%% V1 R
selected_cells_R = [];
selected_cells_L = [];

selected_cells_R = unique([find(contains(place_fields_all_R(1).region,'V1_R') & ...
    place_fields_all_R(1).ripple_modulation_percentile>0.90) find(contains(place_fields_all_R(1).region,'V1_R') & ...
    place_fields_all_R(2).ripple_modulation_percentile>0.90)]);

selected_cells_R = unique([find(contains(place_fields_all_R(1).region,'V1_L')& place_fields_all_R(1).session_id == 8)]);


x_edges = -2:0.02:2;
x_bins = x_edges(2:end)-diff(x_edges)/2;

sum_sqr_diff = sum(place_fields_all_R(2).ripple_PSTH_zscored(x_bins>-1&x_bins<1,selected_cells_R)-place_fields_all_R(1).ripple_PSTH_zscored(x_bins>-1&x_bins<1,selected_cells_R)).^2;
[~,index] = sort(sum_sqr_diff);
selected_cells_R = selected_cells_R(index);
% SMI = place_fields_all_R(1).SMI(selected_cells_R);
% [~,index] = sort(SMI);
% selected_cells_R = selected_cells_R(index);

cluster_id = place_fields_all_R(1).cluster_id(selected_cells_R);
session_id = place_fields_all_R(1).session_id(selected_cells_R);
% 
counter = 1;
for ncell = 1:length(cluster_id)
    if sum(place_fields_all_L(1).cluster_id==cluster_id(ncell)&place_fields_all_L(1).session_id==session_id(ncell))>0
        selected_cells_L(counter) = find(place_fields_all_L(1).cluster_id==cluster_id(ncell)&place_fields_all_L(1).session_id==session_id(ncell))
        counter = counter + 1;
    end
end



plot_spatial_theta_ripple_population('All V1 R cells',selected_cells_L,selected_cells_R,place_fields_all_L,place_fields_all_R)
plot_spatial_ripple_cells('All V1 R cells',selected_cells_L,selected_cells_R,place_fields_all_L,place_fields_all_R)


%% population map

clear all
load(fullfile('P:\corticohippocampal_replay\summary','place_fields_all.mat'))


selected_cells = unique([find(contains(place_fields_all_combined(1).region,'V1_R'))]);
options.region = 'V1 R';
normalised_raw_matrix_V1_R = plot_population_map_all_sessions(place_fields_all_combined,selected_cells,options); % Roughly 6-7 mins for shuffle and plotting

selected_cells = unique([find(contains(place_fields_all_combined(1).region,'V1_L'))]);
options.region = 'V1 L';
normalised_raw_matrix_V1_L = plot_population_map_all_sessions(place_fields_all_combined,selected_cells,options); % Roughly 6-7 mins for shuffle and plotting

selected_cells = unique([find(contains(place_fields_all_combined(1).region,'HPC'))]);
options.region = 'HPC';
normalised_raw_matrix_HPC = plot_population_map_all_sessions(place_fields_all_combined,selected_cells,options); % Roughly 6-7 mins for shuffle and plotting

if exist(fullfile('Z:\ibn-vision\DATA\PROJECTS\Masa2tracks','figures','population map all sessions'))== 0
    mkdir(fullfile('Z:\ibn-vision\DATA\PROJECTS\Masa2tracks','figures','population map all sessions'))
end

save_all_figures(fullfile('Z:\ibn-vision\DATA\PROJECTS\Masa2tracks','figures','population map all sessions'),[])

%%

theta_phase_map = reshape([place_fields_all_R(1).theta_phase_map{1:2}],[],size(place_fields_all_R(1).theta_phase_map{1},2))

plot_spatial_theta_ripple_individual_cell(selected_cluster_id,place_fields_all_L,place_fields_all_R)
plot_spatial_theta_ripple_population(selected_cluster_id,place_fields_all_L,place_fields_all_R)

% V1
plot_spatial_theta_ripple_population(selected_cluster_id,place_fields_all_L,place_fields_all_R)

%% Plot cell summary


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

