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
V1_event_modulation_L = [];
V1_event_modulation_R = [];

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
    place_fields_all(track_id).ripple_modulation_percentile = [];
    place_fields_all(track_id).pre_ripple_activation = [];
    place_fields_all(track_id).pre_ripple_inhibition = [];
    place_fields_all(track_id).post_ripple_activation = [];
    place_fields_all(track_id).post_ripple_inhibition = [];

    place_fields_all(track_id).V1_event_spike_count = [];
    place_fields_all(track_id).V1_event_PSTH = [];
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
                            elseif contains(all_fields{iField},'within_block_SMI')==1
                                temp = place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})})(good_cell_index,:)';
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


                        if strcmp(all_fields{iField},'ripple_spike_count')
                            place_fields_all_L(track_id).ripple_spike_count = [place_fields_all_L(track_id).ripple_spike_count ripple_modulation_L(track_id).spike_count];
                        end

                        if strcmp(all_fields{iField},'V1_event_PSTH')
                            PSTH_matrix = [V1_event_modulation_L(track_id).PSTH{:}];
                            PSTH_matrix = reshape(PSTH_matrix,[],length(good_cell_index));
                            place_fields_all_L(track_id).V1_event_PSTH = [place_fields_all_L(track_id).V1_event_PSTH PSTH_matrix];
                        end


                        if strcmp(all_fields{iField},'V1_event_spike_count')
                            place_fields_all_L(track_id).V1_event_spike_count = [place_fields_all_L(track_id).V1_event_spike_count V1_event_modulation_L(track_id).spike_count];
                        end
                    end

                    % Relative distance from the surface for clusters from
                    % different shanks (mainly for cortex)
                    options.importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP
                    [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(options,[]);

                    [Lia,Locb] = ismember(place_fields(track_id).peak_channel(good_cell_index),chan_config.Channel);
                    
                    shank_id_avaliable = ceil(best_channels{nprobe}.xcoord./250); % based on xcoord in best_channels
                    shanks_avaliable = unique(shank_id_avaliable);

                    for nshank = shanks_avaliable
                        surface_depth(nshank) = median(best_channels{nprobe}.surface_depth(shank_id_avaliable == shanks_avaliable(nshank)),'omitnan');
                    end

                    place_fields_all_L(track_id).relative_depth = [place_fields_all_L(track_id).relative_depth  ...
                        surface_depth(chan_config.Shank(Locb))'-ripple_modulation_L(track_id).peak_depth];

                    ratemap_matrix = [place_fields(track_id).raw{good_cell_index}];
                    ratemap_matrix = reshape(ratemap_matrix,size(place_fields(track_id).raw{1},1),[],length(good_cell_index));%laps X position bins X cells
                    place_fields_all_L(track_id).raw = [place_fields_all_L(track_id).raw;place_fields(track_id).raw(good_cell_index)];
                    place_fields_all_L(track_id).average_map = [place_fields_all_L(track_id).average_map squeeze(mean(ratemap_matrix,1))];

                    % spatial modulation
                    place_fields_all(track_id).SMI = [];
                    place_fields_all(track_id).lap_SMI = [];
                    place_fields_all(track_id).within_block_SMI = [];
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
                            place_fields_all_L(track_id).(all_fields{iField}) = [place_fields_all_L(track_id).(all_fields{iField}) ...
                                V1_event_modulation_R(track_id).(V1_event_fields{strcmp(V1_event_fields,all_fields{iField})})];
                            continue
                        end


                        if sum(strcmp(pf_fields,all_fields{iField})) == 1
                            if contains(all_fields{iField},'x_bin')==1
                                place_fields_all_R(track_id).(all_fields{iField}) = place_fields(track_id).(pf_fields{strcmp(pf_fields,all_fields{iField})});
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

                        if strcmp(all_fields{iField},'ripple_spike_count')
                            place_fields_all_R(track_id).ripple_spike_count = [place_fields_all_R(track_id).ripple_spike_count ripple_modulation_R(track_id).spike_count];
                        end

                        if strcmp(all_fields{iField},'V1_event_PSTH')
                            PSTH_matrix = [V1_event_modulation_R(track_id).PSTH{:}];
                            PSTH_matrix = reshape(PSTH_matrix,[],length(good_cell_index));
                            place_fields_all_R(track_id).V1_event_PSTH = [place_fields_all_R(track_id).V1_event_PSTH PSTH_matrix];
                        end


                        if strcmp(all_fields{iField},'V1_event_spike_count')
                            place_fields_all_R(track_id).V1_event_spike_count = [place_fields_all_R(track_id).V1_event_spike_count V1_event_modulation_R(track_id).spike_count];
                        end
                    end
                    
                    % Relative distance from the surface for clusters from
                    % different shanks (mainly for cortex)
                    options.importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP
                    [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(options,[]);

                    [Lia,Locb] = ismember(place_fields(track_id).peak_channel(good_cell_index),chan_config.Channel);
                    
                    shank_id_avaliable = ceil(best_channels{nprobe}.xcoord./250); % based on xcoord in best_channels
                    shanks_avaliable = unique(shank_id_avaliable);

                    for nshank = shanks_avaliable
                        surface_depth(nshank) = median(best_channels{nprobe}.surface_depth(shank_id_avaliable == shanks_avaliable(nshank)),'omitnan');
                    end

                    place_fields_all_R(track_id).relative_depth = [place_fields_all_R(track_id).relative_depth  ...
                        surface_depth(chan_config.Shank(Locb))'-ripple_modulation_R(track_id).peak_depth];


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

save(fullfile('D:\corticohippocampal_replay\summary','place_fields_all.mat'),'place_fields_all_L','place_fields_all_R')


% 
% save(fullfile(options.ANALYSIS_DATAPATH,'..','..','figures','population_map'))
% save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','..','figures','population_map'))

%% plot population map or individual cell with all information

% HPC
selected_cells_R = [];
selected_cells_L = [];

selected_cells_L = unique([find(contains(place_fields_all_L(1).region,'V1_L') & ...
    place_fields_all_L(1).ripple_modulation_percentile>0.95) find(contains(place_fields_all_L(1).region,'V1_L') & ...
    place_fields_all_L(2).ripple_modulation_percentile>0.95)]);

average_map = normalize([place_fields_all_L(1).average_map(:,selected_cells_L);place_fields_all_L(2).average_map(:,selected_cells_L)],'range')';
average_map = reshape(average_map,size(average_map,1),size(place_fields_all_L(1).average_map,1),[]);

[~,index]=max(squeeze(average_map(:,:,1)'));
[~,index]=sort(index);
selected_cells_L = selected_cells_L(index);
% selected_cells_L = find(contains(place_fields_all_L(1).region,'V1_R'));

cluster_id = place_fields_all_L(1).cluster_id(selected_cells_L);
session_id = place_fields_all_L(1).session_id(selected_cells_L);

counter = 1;
for ncell = 1:length(cluster_id)
    if sum(place_fields_all_R(1).cluster_id==cluster_id(ncell)&place_fields_all_R(1).session_id==session_id(ncell))>0
        selected_cells_R(counter) = find(place_fields_all_R(1).cluster_id==cluster_id(ncell)&place_fields_all_R(1).session_id==session_id(ncell))
        counter = counter + 1;
    end
end





selected_cells_R = [];
selected_cells_L = [];

selected_cells_R = unique([find(contains(place_fields_all_R(1).region,'V1_L') & ...
    place_fields_all_R(1).ripple_modulation_percentile<0.95) find(contains(place_fields_all_R(1).region,'V1_L') & ...
    place_fields_all_R(2).ripple_modulation_percentile<0.95)]);

average_map = normalize([place_fields_all_R(1).average_map(:,selected_cells_R);place_fields_all_R(2).average_map(:,selected_cells_R)],'range')';
average_map = reshape(average_map,size(average_map,1),size(place_fields_all_R(1).average_map,1),[]);

[~,index]=max(squeeze(average_map(:,:,1)'));
[~,index]=sort(index);
selected_cells_L = selected_cells_L(index);

cluster_id = place_fields_all_R(1).cluster_id(selected_cells_R);
session_id = place_fields_all_R(1).session_id(selected_cells_R);

counter = 1;
for ncell = 1:length(cluster_id)
    if sum(place_fields_all_L(1).cluster_id==cluster_id(ncell)&place_fields_all_L(1).session_id==session_id(ncell))>0
        selected_cells_L(counter) = find(place_fields_all_L(1).cluster_id==cluster_id(ncell)&place_fields_all_L(1).session_id==session_id(ncell))
        counter = counter + 1;
    end
end


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

