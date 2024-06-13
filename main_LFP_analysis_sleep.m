%% Main pipeline for extracting and saving and analysing LFP including ripple event and theta cycle detection
% This is a higher-order multi-purpose pipeline that calls dependent functions for
% different analysis piepline

% For mapping of visual receptive field, Please refer to
% Sparse_noise_RF_mapping_masa.mat (subject to change)


%% Set path
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

%% Extract and save LFP from noise channel, cortical L5 channel
clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
SUBJECTS={'M24017'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% experiment_info = experiment_info(4);
% Stimulus_type = 'RUN';
Stimulus_type = 'SleepChronic';

for nsession = [10]
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    if isempty(stimulus_name)
        continue
    end

    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
         load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'power');
%         load(fullfile(options.ANALYSIS_DATAPATH,'..','extracted_PSD.mat'),'power');
        raw_LFP = [];
        LFP = [];
        
        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);
            selected_channels = [];
            channel_regions = [];
            shank_id = [];

            all_fields = fieldnames(best_channels{nprobe});
            all_fields = {all_fields{contains(all_fields,'depth')}};
            all_fields = erase(all_fields,'_depth');
            for nregion = 1:length(all_fields)
                region_name = all_fields{nregion};
                [~,channels_temp] = determine_region_channels(best_channels{nprobe},options,'region',region_name,'group','by shank');
                %                 noise_channel{options.probe_hemisphere} = channels_temp;
                selected_channels = [selected_channels channels_temp];
                channel_regions = [channel_regions nregion*ones(1,length(channels_temp))];
                shank_id = [shank_id unique(ceil(best_channels{nprobe}.xcoord/250))];
            end

            channel_regions(isnan(selected_channels)) = []; % remove nan channel (Missing best channels for some shanks e.g. only 3 shanks with CA1)
            shank_id(isnan(selected_channels)) = [];
            selected_channels(isnan(selected_channels)) = [];
                

%             plot(power{1}(:,5))
%             hold on;
%             xline(selected_channels(4))
%             scatter(selected_channels(4),power{1}(selected_channels(4),5))

%                 column = 1;
            [raw_LFP{nprobe},tvec,SR,chan_config,~] = load_LFP_NPX(options,[],'selected_channels',selected_channels);
            
            selected_chan_config = chan_config(selected_channels,:);
            [PSD{nprobe},power{nprobe}] = calculate_channel_PSD(raw_LFP{nprobe},SR,selected_chan_config,options,'plot_option',0);
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'),'PSD','power')% save PSD for the sleep session

            % Save downsampled LFP from key channels
            options.importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP
            [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(options,[]);
            if isfield(options,'probe_hemisphere')
                LFP(nprobe).probe_hemisphere = options.probe_hemisphere;
                LFP(nprobe).probe_id = options.probe_id;
            else
                LFP(nprobe).probe_id = options.probe_id;
            end

            LFP(nprobe).tvec = tvec;
            for nregion = 1:length(all_fields)
                if sum(channel_regions == nregion)>0
                    LFP(nprobe).(all_fields{nregion}) = raw_LFP{nprobe}(channel_regions == nregion,:);
                    LFP(nprobe).(sprintf('%s_shank_id',all_fields{nregion})) = shank_id(channel_regions == nregion); % only avaliable shanks
                    LFP(nprobe).(sprintf('%s_channel',all_fields{nregion})) = selected_channels(channel_regions == nregion);
                    LFP(nprobe).(sprintf('%s_depth',all_fields{nregion})) = chan_config.Ks_ycoord(selected_channels(channel_regions == nregion));
                    LFP(nprobe).(sprintf('%s_power',all_fields{nregion})) = power{nprobe}(channel_regions == nregion,:);
                else
                    LFP(nprobe).(all_fields{nregion}) = [];
                    LFP(nprobe).(sprintf('%s_shank_id',all_fields{nregion})) = []; % only avaliable shanks
                    LFP(nprobe).(sprintf('%s_channel',all_fields{nregion})) = [];
                    LFP(nprobe).(sprintf('%s_depth',all_fields{nregion})) = [];
                    LFP(nprobe).(sprintf('%s_power',all_fields{nregion})) = [];
                end
                %                 LFP(nprobe).(sprintf('%s_depth',all_fields{nregion})) = chan_config.Ks_ycoord(selected_channels(channel_regions == nregion));
            end
        end

        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP_sleep.mat'),'LFP','-v7.3')
    end
end


%% Detect ripples and reactivation event

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';

% 1:length(experiment_info)
% [1 2 3 4 6 7 8 9 10 12 14]

for nsession =[5]
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        if ~contains(stimulus_name{n},'sleep')
            continue
        end
        
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'power');
        load(fullfile(options.ANALYSIS_DATAPATH,'..','extracted_PSD.mat'),'power');

        if contains(stimulus_name{n},'Masa2tracks')
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
        else
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP')
        end

        clear replay reactivations ripples slow_waves raw_LFP CA1_clusters V1_clusters V1_replay V1_reactivations replay_combined replay_combined

        if contains(stimulus_name{n},'Masa')
            %             load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            %             merged_clusters = clusters_ks3;
            if exist(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))) == 2
                load(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_place_fields.mat'));
            else exist(fullfile(options.ANALYSIS_DATAPATH,sprintf('across_session_merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))) == 2
                load(fullfile(options.ANALYSIS_DATAPATH,sprintf('across_session_merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_across_session_place_fields.mat'));
            end
        else
            load(fullfile(options.ANALYSIS_DATAPATH,'merged_clusters.mat'))
        end

        if length(merged_clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end


        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
            %                 Behavioural state detection

            abs(Behaviour.mobility-movmean(Behaviour.mobility,60*60))
            speed_interp = interp1(Behaviour.tvec,Behaviour.mobility,LFP(probe_no).tvec','linear');
%             mobility = movmean(Behaviour.mobility,120);
            mobility = abs([0 diff(movmean(Behaviour.mobility,30))])>1000;
            speedTreshold = 1;

            if isfield(LFP(probe_no),'L5')
                if ~isempty(LFP(probe_no).L5)
                    cortex_LFP = LFP(probe_no).L5;
                else
                     cortex_LFP = LFP(probe_no).L4;
                end

            elseif isfield(LFP(probe_no),'L4')
                if ~isempty(LFP(probe_no).L4)
                    cortex_LFP = LFP(probe_no).L4;
                else
                     cortex_LFP = []; 
                     disp('cortex LFP is missing')
                end

            elseif isfield(LFP(probe_no),'MEC')

            end

            if isfield(LFP(probe_no),'CA1')
                CA1_LFP = LFP(probe_no).CA1;
                %                 CA1_LFP = raw_LFP{nprobe}(229,:);
                [freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
                    [LFP(probe_no).tvec' cortex_LFP'],[LFP(probe_no).tvec' CA1_LFP'],...
                    [LFP(probe_no).tvec' speed_interp],speedTreshold);
            end
            %             behavioural_state.freezing = freezing;
            behavioural_state(probe_no).quietWake = quietWake;
            behavioural_state(probe_no).SWS = SWS;
            behavioural_state(probe_no).REM = REM;
            behavioural_state(probe_no).movement = movement;

            %%%%%%%%%%%%%%%%%%
            % Ripple and candidate reactivation events detection
            %%%%%%%%%%%%%%%%%%

            zscore_min = 0;
            zscore_max = 3;

            HPC_channels = determine_region_channels(best_channels{probe_no},options,'region','CA1','group','by probe');
            %             merged_clusters.region
            sorting_option = 'spikeinterface';
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.peak_channel = @(x) ismember(x,HPC_channels);
            %             metric_param.cell_type = @(x) x==1;
            %             merged_clusters(nprobe).region
            %             metric_param.merged_cluster_id = @(x) ismember(x,place_fields(nprobe).cluster_id(place_fields(nprobe).all_good_cells_LIBERAL));

            CA1_clusters(probe_no) = select_clusters(merged_clusters(probe_no),metric_param);

            if ~isempty(CA1_clusters(probe_no).cluster_id)
                [replay(probe_no),reactivations(probe_no)] = detect_candidate_events_masa(LFP(probe_no).tvec,CA1_LFP,...
                    [CA1_clusters(probe_no).spike_id CA1_clusters(probe_no).spike_times],Behaviour,zscore_min,zscore_max,options);
                event_no = length(reactivations(probe_no).onset);
            else
                event_no = 0;
            end


            if event_no < 50
                HPC_channels = determine_region_channels(best_channels{probe_no},options,'region','HPC','group','by probe');
                %             merged_clusters.region
                sorting_option = 'spikeinterface';
                metric_param = create_cluster_selection_params('sorting_option',sorting_option);
                metric_param.peak_channel = @(x) ismember(x,HPC_channels);
                metric_param.cell_type = @(x) x==1;
                CA1_clusters(probe_no) = select_clusters(merged_clusters(probe_no),metric_param);
                [replay(probe_no),reactivations(probe_no)] = detect_candidate_events_masa(LFP(probe_no).tvec,CA1_LFP,...
                    [CA1_clusters(probe_no).spike_id CA1_clusters(probe_no).spike_times],Behaviour,zscore_min,zscore_max,options);
                
            end

            % Detect V1 populational bursting events (Candidate events)
            V1_channels = determine_region_channels(best_channels{probe_no},options,'region','V1','group','by probe');
            %             merged_clusters.region
            sorting_option = 'spikeinterface';
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.peak_channel = @(x) ismember(x,V1_channels);
            metric_param.cell_type = @(x) x==1;

            V1_clusters(probe_no) = select_clusters(merged_clusters(probe_no),metric_param);

%             % Slow wave detections
%             if isfield(LFP(nprobe),'L5')==1
%                 slow_waves(nprobe) = DetectSlowWaves_masa('time',LFP(nprobe).tvec,'lfp',LFP(nprobe).L5,'spikes',V1_clusters(nprobe));
%             else
%                 slow_waves(nprobe) = DetectSlowWaves_masa('time',LFP(nprobe).tvec,'lfp',LFP(nprobe).L4,'spikes',V1_clusters(nprobe));
%             end

            % Select spatially tuned cells (reliable visual response to the landmark)
            spatial_cell_index = unique([find(place_fields(1).peak_percentile>0.95 & place_fields(1).odd_even_stability>0.95)...
                find(place_fields(2).peak_percentile>0.95 & place_fields(2).odd_even_stability>0.95)]);

            if nprobe ==1
                this_probe_cluster_id = place_fields(1).cluster_id;
            else
                this_probe_cluster_id = place_fields(1).cluster_id-10000;
            end

            % Select only if there are more than 5 spatial cells
            if sum(ismember(V1_clusters(probe_no).merged_cluster_id,this_probe_cluster_id(spatial_cell_index))) > 5
                metric_param.merged_cluster_id = @(x) ismember(x,this_probe_cluster_id(spatial_cell_index));
            else
                disp('less than 5 spatially tuned V1 cells in this session. Using all V1 cells for candidate event detection')
            end
               
%             metric_param.merged_cluster_id = @(x) ismember(x,place_fields(nprobe).cluster_id(place_fields(nprobe).all_good_cells_LIBERAL));

            V1_clusters(probe_no) = select_clusters(V1_clusters(probe_no),metric_param);

            [V1_replay(probe_no),V1_reactivations(probe_no)] = detect_candidate_events_masa(LFP(probe_no).tvec,CA1_LFP,...
                [V1_clusters(probe_no).spike_id V1_clusters(probe_no).spike_times],Behaviour,zscore_min,zscore_max,options);

            %             [reactivations.probe(nprobe).awake_offset,reactivations.probe(nprobe).awake_index] = RestrictInts(reactivations.probe(nprobe).offset,behavioural_state.quietWake);
            %             reactivations.probe(nprobe).awake_onset = reactivations.probe(nprobe).onset(reactivations.probe(nprobe).awake_index);
            %             reactivations.probe(nprobe).awake_peaktimes = reactivations.probe(nprobe).peaktimes(reactivations.awake_index);

            % Detect CA1 ripple events
            [ripples(probe_no)] = FindRipples_masa(LFP(probe_no).CA1',LFP(probe_no).tvec','behaviour',Behaviour,'minDuration',20,'durations',[30 200],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
                'noise',LFP(probe_no).surface','passband',[125 300],'thresholds',[3 5],'show','off');

        end


        %%%%%%%%%%%%%%%%%%
        % Candidate reactivation events detection (probe combined)
        %%%%%%%%%%%%%%%%%%
        reactivations_combined= [];
        replay_combined = [];
        clear CA1_clusters

        if length(session_info(n).probe)>1
            zscore_min = 0;
            zscore_max = 3;

            for nprobe = 1: length(session_info(n).probe)
                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

                HPC_channels = determine_region_channels(best_channels{probe_no},options,'region','HPC','group','by probe');
                %             merged_clusters.region
                sorting_option = 'spikeinterface';
                metric_param = create_cluster_selection_params('sorting_option',sorting_option);
                metric_param.peak_channel = @(x) ismember(x,HPC_channels);
                metric_param.cell_type = @(x) x==1;
                CA1_clusters(probe_no) = select_clusters(merged_clusters(probe_no),metric_param);
            end

            spike_times = [CA1_clusters(1).spike_times; CA1_clusters(2).spike_times];
            spike_id=  [CA1_clusters(1).spike_id; CA1_clusters(2).spike_id];
            [spike_times,index]= sort(spike_times);
            spike_id = spike_id(index);

            [replay_combined,reactivations_combined] = detect_candidate_events_masa(LFP(nprobe).tvec,CA1_LFP,...
                [spike_id spike_times],Behaviour,zscore_min,zscore_max,options);
        end

        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

            for event = 1:length(ripples(probe_no).onset)
                ripples(probe_no).speed(event) = mean(Behaviour.speed(find(Behaviour.sglxTime >= ripples(probe_no).onset(event) & Behaviour.sglxTime <= ripples(probe_no).offset(event))));
            end

            if ~contains(stimulus_name{n},'RUN') % If reactivation events during lap running
                continue
            end

            lap_times(1).start = Task_info.start_time_all(Task_info.track_ID_all==1);
            lap_times(1).end = Task_info.end_time_all(Task_info.track_ID_all==1)';

            lap_times(2).start = Task_info.start_time_all(Task_info.track_ID_all==2);
            lap_times(2).end = Task_info.end_time_all(Task_info.track_ID_all==2)';

            if contains(stimulus_name{n},'RUN') % If reactivation events during lap running
                [reactivations(probe_no).T1_offset,reactivations(probe_no).T1_index] = RestrictInts(reactivations(probe_no).offset',[lap_times(1).start lap_times(1).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                reactivations(probe_no).T1_onset = reactivations(probe_no).onset(reactivations(probe_no).T1_index);
                reactivations(probe_no).T1_midpoint = reactivations(probe_no).midpoint(reactivations(probe_no).T1_index);

                [reactivations(probe_no).T2_offset,reactivations(probe_no).T2_index] = RestrictInts(reactivations(probe_no).offset',[lap_times(2).start lap_times(2).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                reactivations(probe_no).T2_onset = reactivations(probe_no).onset(reactivations(probe_no).T2_index);
                reactivations(probe_no).T2_midpoint = reactivations(probe_no).midpoint(reactivations(probe_no).T2_index);

                [V1_reactivations(probe_no).T1_offset,V1_reactivations(probe_no).T1_index] = RestrictInts(V1_reactivations(probe_no).offset',[lap_times(1).start lap_times(1).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                V1_reactivations(probe_no).T1_onset = V1_reactivations(probe_no).onset(V1_reactivations(probe_no).T1_index);
                V1_reactivations(probe_no).T1_midpoint = V1_reactivations(probe_no).midpoint(V1_reactivations(probe_no).T1_index);

                [V1_reactivations(probe_no).T2_offset,V1_reactivations(probe_no).T2_index] = RestrictInts(V1_reactivations(probe_no).offset',[lap_times(2).start lap_times(2).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                V1_reactivations(probe_no).T2_onset = V1_reactivations(probe_no).onset(V1_reactivations(probe_no).T2_index);
                V1_reactivations(probe_no).T2_midpoint = V1_reactivations(probe_no).midpoint(V1_reactivations(probe_no).T2_index);
            end

            if contains(stimulus_name{n},'RUN') % If reactivation events during lap running
                [ripples(probe_no).T1_offset,ripples(probe_no).T1_index] = RestrictInts(ripples(probe_no).offset',[lap_times(1).start lap_times(1).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                ripples(probe_no).T1_onset = ripples(probe_no).onset(ripples(probe_no).T1_index);
                ripples(probe_no).T1_peaktimes = ripples(probe_no).peaktimes(ripples(probe_no).T1_index);

                [ripples(probe_no).T2_offset,ripples(probe_no).T2_index] = RestrictInts(ripples(probe_no).offset',[lap_times(2).start lap_times(2).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                ripples(probe_no).T2_onset = ripples(probe_no).onset(ripples(probe_no).T2_index);
                ripples(probe_no).T2_peaktimes = ripples(probe_no).peaktimes(ripples(probe_no).T2_index);
            end
        end
        
        if contains(stimulus_name{n},'RUN')
            if length(session_info(n).probe)>1
                [reactivations_combined.T1_offset,reactivations_combined.T1_index] = RestrictInts(reactivations_combined.offset',[lap_times(1).start lap_times(1).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                reactivations_combined.T1_onset = reactivations_combined.onset(reactivations_combined.T1_index);
                reactivations_combined.T1_midpoint = reactivations_combined.midpoint(reactivations_combined.T1_index);

                [reactivations_combined.T2_offset,reactivations_combined.T2_index] = RestrictInts(reactivations_combined.offset',[lap_times(2).start lap_times(2).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                reactivations_combined.T2_onset = reactivations_combined.onset(reactivations_combined.T2_index);
                reactivations_combined.T2_midpoint = reactivations_combined.midpoint(reactivations_combined.T2_index);
            end
        end

        clear lap_times

        save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'replay','reactivations','replay_combined','reactivations_combined')
        save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_candidate_events_V1%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'V1_replay','V1_reactivations')
        save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'ripples')
        save(fullfile(options.ANALYSIS_DATAPATH,sprintf('behavioural_state%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'behavioural_state')
       
        close all
    end
end
