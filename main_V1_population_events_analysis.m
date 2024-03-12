%% Main V1 reactivation analysis codes

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))

%% V1 populational events modulation

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

        for nprobe = 1:length(clusters)

            options = session_info(n).probe(nprobe);
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
            [C,ia,ic] = unique(clusters_combined.merged_cluster_id);
            
            event_id = [ones(1,length(ripples(nprobe).T1_onset)) 2*ones(1,length(ripples(nprobe).T2_onset))];
            [event_times,index] = sort([ripples(nprobe).T1_onset ripples(nprobe).T2_onset]);

            if options.probe_hemisphere == 1
                [ripple_modulation_L]= ripple_modulation_analysis(clusters_combined.spike_times,clusters_combined.merged_spike_id,Task_info,Behaviour,[-0.3 0.3],0.02,...
                'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',event_times',...
                'event_label',{'Track 1','Track 2'},'event_id',event_id(index)','place_fields',place_fields);

            elseif options.probe_hemisphere == 2
                [ripple_modulation_R]= ripple_modulation_analysis(clusters_combined.spike_times,clusters_combined.merged_spike_id,Task_info,Behaviour,[-0.3 0.3],0.02,...
                'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',event_times',...
                'event_label',{'Track 1','Track 2'},'event_id',event_id(index)','place_fields',place_fields);

            end

        end

        save(fullfile(options.ANALYSIS_DATAPATH,'ripple_modulation.mat'),"ripple_modulation_L","ripple_modulation_R")

    end
end

%% Visualise peri ripple spike times

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

        for nprobe = 1:length(clusters)

            options = session_info(n).probe(nprobe);
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
            [C,ia,ic] = unique(clusters_combined.merged_cluster_id);
            
            event_id = [ones(1,length(ripples(nprobe).T1_onset)) 2*ones(1,length(ripples(nprobe).T2_onset))];
            [event_times,index] = sort([ripples(nprobe).T1_onset ripples(nprobe).T2_onset]);

            if options.probe_hemisphere == 1
                [ripple_modulation_L]= ripple_modulation_analysis(clusters_combined.spike_times,clusters_combined.merged_spike_id,Task_info,Behaviour,[-0.3 0.3],0.02,...
                'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',event_times',...
                'event_label',{'Track 1','Track 2'},'event_id',event_id(index)','place_fields',place_fields);

            elseif options.probe_hemisphere == 2
                [ripple_modulation_R]= ripple_modulation_analysis(clusters_combined.spike_times,clusters_combined.merged_spike_id,Task_info,Behaviour,[-0.3 0.3],0.02,...
                'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',event_times',...
                'event_label',{'Track 1','Track 2'},'event_id',event_id(index)','place_fields',place_fields);

            end

        end

        save(fullfile(options.ANALYSIS_DATAPATH,'ripple_modulation.mat'),"ripple_modulation_L","ripple_modulation_R")

    end
end

%% Peri V1 populational events

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
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_candidate_events_V1%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

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
                    [lfpAvg(iprobe).track(1),csd(iprobe).track(1)] = perievent_LFP_profile('L V1 events T1',V1_reactivations(iprobe).T1_onset,[-0.2 0.5],PSD,clusters,options,...
                        'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col,'CSD_V1_CA1_normalisation',normalisation);

                    [lfpAvg(iprobe).track(2),csd(iprobe).track(2)] = perievent_LFP_profile('L V1 events T2',V1_reactivations(iprobe).T2_onset,[-0.2 0.5],PSD,clusters,options,...
                        'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col,'CSD_V1_CA1_normalisation',normalisation);
            
                elseif session_info(n).probe(iprobe).probe_hemisphere == 2
                    [lfpAvg(iprobe).track(1),csd(iprobe).track(1)] = perievent_LFP_profile('R V1 events T1',V1_reactivations(iprobe).T1_onset,[-0.2 0.5],PSD,clusters,options,...
                        'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col,'CSD_V1_CA1_normalisation',normalisation);

                    [lfpAvg(iprobe).track(2),csd(iprobe).track(2)] = perievent_LFP_profile('R V1 events T2',V1_reactivations(iprobe).T2_onset,[-0.2 0.5],PSD,clusters,options,...
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

            if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','V1_events'))== 0
                mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','V1_events'))
            end
        
            if unique(shank_id_avaliable)==1
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('V1_events_LFP_response%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'lfpAvg_L','lfpAvg_R','csd_L','csd_R');

            end

            save_all_large_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','V1_events'),[])
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
                        [lfpAvg(iprobe).track(1),csd(iprobe).track(1)] = perievent_LFP_profile('L V1 events T1',V1_reactivations(iprobe).T1_onset,[-0.2 0.5],PSD,clusters,options,...
                            'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col,'CSD_V1_CA1_normalisation',normalisation);
                        
                        [lfpAvg(iprobe).track(2),csd(iprobe).track(2)] = perievent_LFP_profile('L V1 events T2',V1_reactivations(iprobe).T2_onset,[-0.2 0.5],PSD,clusters,options,...
                            'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col,'CSD_V1_CA1_normalisation',normalisation);
                        
                    elseif session_info(n).probe(iprobe).probe_hemisphere == 2
                        [lfpAvg(iprobe).track(1),csd(iprobe).track(1)] = perievent_LFP_profile('R V1 events T1',V1_reactivations(iprobe).T1_onset,[-0.2 0.5],PSD,clusters,options,...
                            'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col,'CSD_V1_CA1_normalisation',normalisation);
                        
                        [lfpAvg(iprobe).track(2),csd(iprobe).track(2)] = perievent_LFP_profile('R V1 events T2',V1_reactivations(iprobe).T2_onset,[-0.2 0.5],PSD,clusters,options,...
                            'filter_freq',[0.5 3; 9 17; 125 300],'filter_type',{'SO','spindles','ripples'},'x_col',x_col,'CSD_V1_CA1_normalisation',normalisation);
                    end
                    
                end

                if session_info(n).probe(nprobe).probe_hemisphere == 1
                    lfpAvg_L(3:4) = lfpAvg;
                    csd_L(3:4) = csd;

                else
                    lfpAvg_R(3:4) = lfpAvg;
                    csd_R(3:4) = csd;
                end
                
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('V1_events_LFP_response%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'lfpAvg_L','lfpAvg_R','csd_L','csd_R');
                
                
                save_all_large_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','V1_events'),[])
                close all
            end
            
        end


    end
end
