% %% Main pipeline for extracting and saving and analysing LFP including ripple event and theta cycle detection
% % This is a higher-order multi-purpose pipeline that calls dependent functions for
% % different analysis piepline
% 
% % For mapping of visual receptive field, Please refer to
% % Sparse_noise_RF_mapping_masa.mat (subject to change)
% 
% 
% %% Set path
% addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
% addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
% addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))
% %% Extract and save LFP from noise channel, cortical L5 channel
% clear all
% SUBJECTS = {'M23032','M23034','M23037','M23038'};
% option = 'V1-MEC';
% experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Stimulus_type = 'Track';
% 
% for nsession =1:15
%     session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%     stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%     if isempty(stimulus_name)
%         continue
%     end
% 
%     load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));
% 
%     for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
%         options = session_info(n).probe(1);
%         load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
%         load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
%         %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'power');
%         load(fullfile(options.ANALYSIS_DATAPATH,'..','extracted_PSD.mat'),'power');
%         raw_LFP = [];
%         LFP = [];
%         
%         for nprobe = 1:length(session_info(n).probe)
%             options = session_info(n).probe(nprobe);
%             selected_channels = [];
%             channel_regions = [];
%             shank_id = [];
% 
%             all_fields = fieldnames(best_channels{nprobe});
%             all_fields = {all_fields{contains(all_fields,'depth')}};
%             all_fields = erase(all_fields,'_depth');
%             for nregion = 1:length(all_fields)
%                 region_name = all_fields{nregion};
%                 [~,channels_temp] = determine_region_channels(best_channels{nprobe},options,'region',region_name,'group','by shank');
%                 %                 noise_channel{options.probe_hemisphere} = channels_temp;
%                 selected_channels = [selected_channels channels_temp];
%                 channel_regions = [channel_regions nregion*ones(1,length(channels_temp))];
%                 shank_id = [shank_id unique(ceil(best_channels{nprobe}.xcoord/250))];
%             end
%             channel_regions(isnan(selected_channels)) = []; % remove nan channel (Missing best channels for some shanks e.g. only 3 shanks with CA1)
%             shank_id(isnan(selected_channels)) = [];
%             selected_channels(isnan(selected_channels)) = [];
%                 
% 
% %             plot(power{1}(:,5))
% %             hold on;
% %             xline(selected_channels(4))
% %             scatter(selected_channels(4),power{1}(selected_channels(4),5))
% 
%             %     column = 1;
%             if nprobe ~= 1
%                 % Information about AP band probe 1 sample size to align probe
%                 % 2 LFP traces.
%                 session_info(n).probe(1).importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP
%                 [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(session_info(n).probe(1),[]);
%                 [raw_LFP{nprobe} tvec SR chan_config ~] = load_LFP_NPX(options,[],'probe_no',nprobe,'probe_1_total_sample',imecMeta.nFileSamp,'selected_channels',selected_channels);
%             else
%                 %         session_info.probe(1).importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP
%                 %         [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(session_info.probe(1),[]);
%                 [raw_LFP{nprobe} tvec SR chan_config ~] = load_LFP_NPX(options,[],'selected_channels',selected_channels);
%                 %                 [raw_LFP{2} tvec SR chan_config ~] = load_LFP_NPX(options,[],'selected_channels',selected_channels,'duration_sec',10);
%             end
% 
%             % Save downsampled LFP from key channels
%             options.importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP
%             [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(options,[]);
%             if isfield(options,'probe_hemisphere')
%                 LFP(nprobe).probe_hemisphere = options.probe_hemisphere;
%                 LFP(nprobe).probe_id = options.probe_id;
%             else
%                 LFP(nprobe).probe_id = options.probe_id;
%                 LFP(nprobe).probe_V1 = options.probe_V1;
%                 LFP(nprobe).probe_MEC = options.probe_MEC;
%             end
% 
%             LFP(nprobe).tvec = tvec;
%             for nregion = 1:length(all_fields)
%                 LFP(nprobe).(all_fields{nregion}) = raw_LFP{nprobe}(channel_regions == nregion,:);
%                 LFP(nprobe).(sprintf('%s_shank_id',all_fields{nregion})) = shank_id(channel_regions == nregion); % only avaliable shanks
%                 LFP(nprobe).(sprintf('%s_channel',all_fields{nregion})) = selected_channels(channel_regions == nregion);
%                 LFP(nprobe).(sprintf('%s_depth',all_fields{nregion})) = chan_config.Ks_ycoord(selected_channels(channel_regions == nregion));
%                 LFP(nprobe).(sprintf('%s_power',all_fields{nregion})) = power{nprobe}(selected_channels(channel_regions == nregion),:);
%                 %                 LFP(nprobe).(sprintf('%s_depth',all_fields{nregion})) = chan_config.Ks_ycoord(selected_channels(channel_regions == nregion));
%             end
%         end
% 
%         save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP')
%     end
% end
% 

%% Detect ripples and reactivation event

clear all
SUBJECTS = {'M23038'};
option = 'V1-MEC';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'Track';

% 1:length(experiment_info)
% [1 2 3 4 6 7 8 9 10 12 14]
[1 2 3 4 6 7 8]

for nsession =1:2
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%     if isempty(stimulus_name)
%         continue
%     end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)

        
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'power');
        load(fullfile(options.ANALYSIS_DATAPATH,'..','extracted_PSD.mat'),'power');
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_place_fields.mat'));
        
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'))


        clear replay reactivations ripples slow_waves raw_LFP CA1_clusters V1_clusters V1_replay V1_reactivations


            load(fullfile(options.ANALYSIS_DATAPATH,'merged_clusters.mat'))
        

        if length(merged_clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end


        for nprobe = 2
            options = session_info(n).probe(nprobe);
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
            %                 Behavioural state detection
            speed_interp = interp1(Behaviour.tvec,Behaviour.speed,LFP(nprobe).tvec','linear');
            speedTreshold = 1;

            if isfield(LFP(nprobe),'L5')
                cortex_LFP = LFP(nprobe).L5;
            elseif isfield(LFP(nprobe),'L4')
                cortex_LFP = LFP(nprobe).L4;
            elseif isfield(LFP(nprobe),'MEC')
                cortex_LFP = LFP(nprobe).MEC_theta;
            end

            if isfield(LFP(nprobe),'CA1')
                CA1_LFP = LFP(nprobe).CA1;
                %                 CA1_LFP = raw_LFP{nprobe}(229,:);
                [freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
                    [LFP(nprobe).tvec' cortex_LFP'],[LFP(nprobe).tvec' CA1_LFP'],...
                    [LFP(nprobe).tvec' speed_interp],speedTreshold);
            end
            %             behavioural_state.freezing = freezing;
            behavioural_state(nprobe-1).quietWake = quietWake;
            behavioural_state(nprobe-1).SWS = SWS;
            behavioural_state(nprobe-1).REM = REM;
            behavioural_state(nprobe-1).movement = movement;

            %%%%%%%%%%%%%%%%%%
            % Ripple and candidate reactivation events detection
            %%%%%%%%%%%%%%%%%%

            zscore_min = 0;
            zscore_max = 3;

            HPC_channels = determine_region_channels(best_channels{nprobe},options,'region','CA1','group','by probe');
            %             merged_clusters.region
            sorting_option = 'spikeinterface';
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.peak_channel = @(x) ismember(x,HPC_channels);
            %             metric_param.cell_type = @(x) x==1;
            %             merged_clusters(nprobe).region
            %             metric_param.merged_cluster_id = @(x) ismember(x,place_fields(nprobe).cluster_id(place_fields(nprobe).all_good_cells_LIBERAL));

            CA1_clusters(nprobe-1) = select_clusters(merged_clusters(nprobe),metric_param);

            if ~isempty(CA1_clusters(nprobe-1).cluster_id)
                [replay(nprobe-1),reactivations(nprobe-1)] = detect_candidate_events_masa(LFP(nprobe).tvec,CA1_LFP,...
                    [CA1_clusters(nprobe-1).spike_id CA1_clusters(nprobe-1).spike_times],Behaviour,zscore_min,zscore_max,options);
                event_no = length(reactivations(nprobe-1).onset);
            else
                event_no = 0;
            end


            if event_no < 50
                HPC_channels = determine_region_channels(best_channels{nprobe},options,'region','HPC','group','by probe');
                %             merged_clusters.region
                sorting_option = 'spikeinterface';
                metric_param = create_cluster_selection_params('sorting_option',sorting_option);
                metric_param.peak_channel = @(x) ismember(x,HPC_channels);
                metric_param.cell_type = @(x) x==1;
                CA1_clusters(nprobe-1) = select_clusters(merged_clusters(nprobe),metric_param);
                [replay(nprobe-1),reactivations(nprobe-1)] = detect_candidate_events_masa(LFP(nprobe-1).tvec,CA1_LFP,...
                    [CA1_clusters(nprobe-1).spike_id CA1_clusters(nprobe-1).spike_times],Behaviour,zscore_min,zscore_max,options);
                
            end

            % Detect V1 populational bursting events (Candidate events)
            V1_channels = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by probe');
            %             merged_clusters.region
            sorting_option = 'spikeinterface';
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.peak_channel = @(x) ismember(x,V1_channels);
            metric_param.cell_type = @(x) x==1;

            V1_clusters(nprobe-1) = select_clusters(merged_clusters(nprobe),metric_param);

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
            if sum(ismember(V1_clusters(nprobe-1).merged_cluster_id,this_probe_cluster_id(spatial_cell_index))) > 5
                metric_param.merged_cluster_id = @(x) ismember(x,this_probe_cluster_id(spatial_cell_index));
            else
                disp('less than 5 spatially tuned V1 cells in this session. Using all V1 cells for candidate event detection')
            end
               
%             metric_param.merged_cluster_id = @(x) ismember(x,place_fields(nprobe).cluster_id(place_fields(nprobe).all_good_cells_LIBERAL));

            V1_clusters(nprobe-1) = select_clusters(V1_clusters(nprobe-1),metric_param);

            [V1_replay(nprobe-1),V1_reactivations(nprobe-1)] = detect_candidate_events_masa(LFP(nprobe).tvec,CA1_LFP,...
                [V1_clusters(nprobe-1).spike_id V1_clusters(nprobe-1).spike_times],Behaviour,zscore_min,zscore_max,options);

            %
            %             [reactivations.probe(nprobe).awake_offset,reactivations.probe(nprobe).awake_index] = RestrictInts(reactivations.probe(nprobe).offset,behavioural_state.quietWake);
            %             reactivations.probe(nprobe).awake_onset = reactivations.probe(nprobe).onset(reactivations.probe(nprobe).awake_index);
            %             reactivations.probe(nprobe).awake_peaktimes = reactivations.probe(nprobe).peaktimes(reactivations.awake_index);

            % Detect CA1 ripple events
            [ripples(nprobe-1)] = FindRipples_masa(LFP(nprobe).CA1',LFP(nprobe).tvec','behaviour',Behaviour,'minDuration',20,'durations',[30 200],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
                'noise',LFP(nprobe).surface','passband',[125 300],'thresholds',[3 5],'show','off');

           % [ripples(nprobe-1)] = DetectSlowWaves_masa(LFP(nprobe).CA1',LFP(nprobe).tvec','behaviour',Behaviour,'minDuration',20,'durations',[30 200],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
            %    'noise',LFP(nprobe).surface','passband',[125 300],'thresholds',[3 5],'show','off');


        end

        for nprobe = 2
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

            for event = 1:length(ripples(probe_no-1).onset)
                ripples(probe_no-1).speed(event) = mean(Behaviour.speed(find(Behaviour.sglxTime >= ripples(probe_no-1).onset(event) & Behaviour.sglxTime <= ripples(probe_no-1).offset(event))));
            end



            lap_times(1).start = Task_info.start_time_all(Task_info.track_ID_all==1);
            lap_times(1).end = Task_info.end_time_all(Task_info.track_ID_all==1)';

            lap_times(2).start = Task_info.start_time_all(Task_info.track_ID_all==2);
            lap_times(2).end = Task_info.end_time_all(Task_info.track_ID_all==2)';


                [reactivations(probe_no-1).T1_offset,reactivations(probe_no-1).T1_index] = RestrictInts(reactivations(probe_no-1).offset',[lap_times(1).start lap_times(1).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                reactivations(probe_no-1).T1_onset = reactivations(probe_no-1).onset(reactivations(probe_no-1).T1_index);
                reactivations(probe_no-1).T1_midpoint = reactivations(probe_no-1).midpoint(reactivations(probe_no-1).T1_index);

                [reactivations(probe_no-1).T2_offset,reactivations(probe_no-1).T2_index] = RestrictInts(reactivations(probe_no-1).offset',[lap_times(2).start lap_times(2).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                reactivations(probe_no-1).T2_onset = reactivations(probe_no-1).onset(reactivations(probe_no-1).T2_index);
                reactivations(probe_no-1).T2_midpoint = reactivations(probe_no-1).midpoint(reactivations(probe_no-1).T2_index);

                [V1_reactivations(probe_no-1).T1_offset,V1_reactivations(probe_no-1).T1_index] = RestrictInts(V1_reactivations(probe_no-1).offset',[lap_times(1).start lap_times(1).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                V1_reactivations(probe_no-1).T1_onset = V1_reactivations(probe_no-1).onset(V1_reactivations(probe_no-1).T1_index);
                V1_reactivations(probe_no-1).T1_midpoint = V1_reactivations(probe_no-1).midpoint(V1_reactivations(probe_no-1).T1_index);

                [V1_reactivations(probe_no-1).T2_offset,V1_reactivations(probe_no-1).T2_index] = RestrictInts(V1_reactivations(probe_no-1).offset',[lap_times(2).start lap_times(2).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                V1_reactivations(probe_no-1).T2_onset = V1_reactivations(probe_no-1).onset(V1_reactivations(probe_no-1).T2_index);
                V1_reactivations(probe_no-1).T2_midpoint = V1_reactivations(probe_no-1).midpoint(V1_reactivations(probe_no-1).T2_index);


                [ripples(probe_no-1).T1_offset,ripples(probe_no-1).T1_index] = RestrictInts(ripples(probe_no-1).offset',[lap_times(1).start lap_times(1).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                ripples(probe_no-1).T1_onset = ripples(probe_no-1).onset(ripples(probe_no-1).T1_index);
                ripples(probe_no-1).T1_peaktimes = ripples(probe_no-1).peaktimes(ripples(probe_no-1).T1_index);

                [ripples(probe_no-1).T2_offset,ripples(probe_no-1).T2_index] = RestrictInts(ripples(probe_no-1).offset',[lap_times(2).start lap_times(2).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                ripples(probe_no-1).T2_onset = ripples(probe_no-1).onset(ripples(probe_no-1).T2_index);
                ripples(probe_no-1).T2_peaktimes = ripples(probe_no-1).peaktimes(ripples(probe_no-1).T2_index);

        end
        clear lap_times

        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'),'replay','reactivations')
        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'),'V1_replay','V1_reactivations')
        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'),'ripples')
        save(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state.mat'),'behavioural_state')
       
        close all
    end
end



%%
%%


