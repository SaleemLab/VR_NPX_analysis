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
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'Masa2tracks';

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        raw_LFP = [];
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

            %     column = 1;
            if nprobe ~= 1
                % Information about AP band probe 1 sample size to align probe
                % 2 LFP traces.
                session_info(n).probe(1).importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP
                [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(session_info(n).probe(1),[]);
                [raw_LFP{nprobe} tvec SR chan_config ~] = load_LFP_NPX(options,[],'probe_no',nprobe,'probe_1_total_sample',imecMeta.nFileSamp,'selected_channels',selected_channels);
            else
                %         session_info.probe(1).importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP
                %         [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(session_info.probe(1),[]);
                [raw_LFP{nprobe} tvec SR chan_config ~] = load_LFP_NPX(options,[],'probe_no',nprobe,'probe_1_total_sample',imecMeta.nFileSamp,'selected_channels',selected_channels);
            end

            % Save downsampled LFP from key channels
            options.importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP
            [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(options,[]);
            LFP(nprobe).probe_hemisphere = options.probe_hemisphere;
            LFP(nprobe).tvec = tvec;
            for nregion = 1:length(all_fields)
                LFP(nprobe).(all_fields{nregion}) = raw_LFP{nprobe}(channel_regions == nregion,:);
                LFP(nprobe).(sprintf('%s_shank_id',all_fields{nregion})) = shank_id(channel_regions == nregion); % only avaliable shanks
                LFP(nprobe).(sprintf('%s_channel',all_fields{nregion})) = selected_channels(channel_regions == nregion);
                LFP(nprobe).(sprintf('%s_depth',all_fields{nregion})) = chan_config.Ks_ycoord(selected_channels(channel_regions == nregion));
                LFP(nprobe).(sprintf('%s_power',all_fields{nregion})) = power{nprobe}(selected_channels(channel_regions == nregion),:);
                %                 LFP(nprobe).(sprintf('%s_depth',all_fields{nregion})) = chan_config.Ks_ycoord(selected_channels(channel_regions == nregion));
            end

        end
        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP')


        clear replay reactivations ripples raw_LFP

        if contains(stimulus_name{n},'Masa')
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        else
            load(fullfile(options.ANALYSIS_DATAPATH,'merged_clusters.mat'))
        end

        if length(merged_clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end

        %%%%%%%%%%%%%%%%%%
        % Ripple and candidate events detection
        %%%%%%%%%%%%%%%%%%
        for nprobe = 1:length(session_info(n).probe)

            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
            %                 Behavioural state detection
            speed_interp = interp1(Behaviour.tvec,Behaviour.speed,LFP(nprobe).tvec','linear');
            speedTreshold = 1;
            if isfield(LFP(nprobe),'L4')
                cortex_LFP = LFP(nprobe).L4;
            elseif isfield(LFP(nprobe),'L5')
                cortex_LFP = LFP(nprobe).L5;
            elseif isfield(LFP(nprobe),'MEC')
            end

            if isfield(LFP(nprobe),'CA1')
                CA1_LFP = LFP(nprobe).CA1;
                cortex_channel = best_channels{nprobe}.L5_channel;
                CA1_channel = best_channels{nprobe}.CA1_channel;
                [freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
                    [LFP(nprobe).tvec' cortex_LFP'],[LFP(nprobe).tvec' CA1_LFP'],...
                    [LFP(nprobe).tvec' speed_interp],speedTreshold);
            end
            %             behavioural_state.freezing = freezing;
            behavioural_state.quietWake = quietWake;
            behavioural_state.SWS = SWS;
            behavioural_state.REM = REM;
            behavioural_state.movement = movement;


            % Detect CA1 populational bursting events (Candidate events)
            zscore_min = 0;
            zscore_max = 3;
            
            HPC_channels = determine_region_channels(best_channels{nprobe},options,'region','CA1','group','by probe');
%             merged_clusters.region
            sorting_option = 'spikeinterface';
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.peak_channel = @(x) ismember(x,HPC_channels);
            CA1_clusters = select_clusters(merged_clusters(nprobe),metric_param);            
            [replay(probe_no),reactivations(probe_no)] = detect_candidate_events_masa(LFP(nprobe).tvec,CA1_LFP,...
                [CA1_clusters.merged_spike_id CA1_clusters.spike_times],Behaviour,zscore_min,zscore_max,options);

            if length(reactivations.probe(probe_no).onset) < 50
                HPC_channels = determine_region_channels(best_channels{nprobe},options,'region','HPC','group','by probe');
                %             merged_clusters.region
                sorting_option = 'spikeinterface';
                metric_param = create_cluster_selection_params('sorting_option',sorting_option);
                metric_param.peak_channel = @(x) ismember(x,HPC_channels);
                CA1_clusters(nprobe) = select_clusters(merged_clusters(nprobe),metric_param);
                [replay(probe_no),reactivations(probe_no)] = detect_candidate_events_masa(LFP(nprobe).tvec,CA1_LFP,...
                    [CA1_clusters.merged_spike_id CA1_clusters.spike_times],Behaviour,zscore_min,zscore_max,options);
            end


            V1_channels = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by probe')
            %             merged_clusters.region
            sorting_option = 'spikeinterface';
            metric_param = create_cluster_selection_params('sorting_option',sorting_option);
            metric_param.peak_channel = @(x) ismember(x,V1_channels);
            V1_clusters = select_clusters(merged_clusters(nprobe),metric_param);
            [replay(probe_no),reactivations(probe_no)] = detect_candidate_events_masa(LFP(nprobe).tvec,CA1_LFP,...
                [V1_clusters.merged_spike_id V1_clusters.spike_times],Behaviour,zscore_min,zscore_max,options);

            %
            %             [reactivations.probe(nprobe).awake_offset,reactivations.probe(nprobe).awake_index] = RestrictInts(reactivations.probe(nprobe).offset,behavioural_state.quietWake);
            %             reactivations.probe(nprobe).awake_onset = reactivations.probe(nprobe).onset(reactivations.probe(nprobe).awake_index);
            %             reactivations.probe(nprobe).awake_peaktimes = reactivations.probe(nprobe).peaktimes(reactivations.awake_index);

            % Detect CA1 ripple events
            [ripples(nprobe)] = FindRipples_masa(CA1_LFP',LFP(nprobe).tvec','behaviour',Behaviour,'minDuration',20,'durations',[30 200],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
                'noise',LFP(nprobe).surface','passband',[125 300],'thresholds',[3 5],'show','on')
            
   
            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(CA1_clusters(nprobe).spike_times, ripples.probe(nprobe).peaktimes, [-0.2 0.2], 0.001);
            plot(rasterX,rasterY); hold on
            plot(bins,zscore(psth));
            yyaxis right

            close all

        end

        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

            for event = 1:length(ripples.probe(probe_no).onset)
                ripples.probe(probe_no).speed(event) = mean(Behaviour.speed(find(Behaviour.sglxTime >= ripples.probe(probe_no).onset(event) & Behaviour.sglxTime <= ripples.probe(probe_no).offset(event))));
            end

            if contains(Stimulus_type,'RUN') % If reactivation events during lap running
                [reactivations.probe(probe_no).T1_offset,reactivations.probe(probe_no).T1_index] = RestrictInts(reactivations.probe(probe_no).offset',[lap_times(1).start' lap_times(1).end']); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                reactivations.probe(probe_no).T1_onset = reactivations.probe(probe_no).onset(reactivations.probe(probe_no).T1_index);
                reactivations.probe(probe_no).T1_midpoint = reactivations.probe(probe_no).midpoint(reactivations.probe(probe_no).T1_index);

                [reactivations.probe(probe_no).T2_offset,reactivations.probe(probe_no).T2_index] = RestrictInts(reactivations.probe(probe_no).offset',[lap_times(2).start' lap_times(2).end']); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                reactivations.probe(probe_no).T2_onset = reactivations.probe(probe_no).onset(reactivations.probe(probe_no).T2_index);
                reactivations.probe(probe_no).T2_midpoint = reactivations.probe(probe_no).midpoint(reactivations.probe(probe_no).T2_index);
            end

            if contains(Stimulus_type,'RUN') % If reactivation events during lap running
                [ripples.probe(probe_no).T1_offset,ripples.probe(probe_no).T1_index] = RestrictInts(ripples.probe(probe_no).offset',[lap_times(1).start' lap_times(1).end']); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                ripples.probe(probe_no).T1_onset = ripples.probe(probe_no).onset(ripples.probe(probe_no).T1_index);
                ripples.probe(probe_no).T1_peaktimes = ripples.probe(probe_no).peaktimes(ripples.probe(probe_no).T1_index);

                [ripples.probe(probe_no).T2_offset,ripples.probe(probe_no).T2_index] = RestrictInts(ripples.probe(probe_no).offset',[lap_times(2).start' lap_times(2).end']); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                ripples.probe(probe_no).T2_onset = ripples.probe(probe_no).onset(ripples.probe(probe_no).T2_index);
                ripples.probe(probe_no).T2_peaktimes = ripples.probe(probe_no).peaktimes(ripples.probe(probe_no).T2_index);
            end
        end

        save(sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks')),'replay','reactivations')

        delete(sprintf('extracted_ripples_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        save(sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks')),'ripples')


    end
end

%% Ripple and Candidate event detection
clear all
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))


SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN'; % extract LFP during RUN
Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';

Stimulus_types_all = {'RUN','POST'};
% Stimulus_type = 'RUN'; % extract LFP during RUN
% Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
c = 1;

for epoch = 1:length(Stimulus_types_all)
    Stimulus_type= Stimulus_types_all{epoch};

    for nsession =1:length(experiment_info)
        tic
        session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        if isempty(session_info)
            continue
        end
        gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

        for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
            cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
            load best_channels
            load extracted_PSD
            load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('bonsai_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load('extracted_laps')
            column = 1;

            raw_LFP = [];
            for nprobe = 1:length(session_info(n).probe) % For each session, how many probes
                column = 1;
                session_info(n).probe(nprobe).task_type = stimulus_name{n};
                options = session_info(n).probe(nprobe);
                options.importMode = 'LF';
                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

                if nprobe ~= 1
                    session_info.probe(1).importMode = 'KS';

                    [~, imecMeta, ~, ~] = extract_NPX_channel_config(session_info.probe(1),column);
                    [raw_LFP{probe_no} tvec SR chan_config sorted_config] = load_LFP_NPX1(options,column,'probe_no',nprobe,'probe_1_total_sample',imecMeta.nFileSamp);
                else
                    [raw_LFP{probe_no} tvec SR chan_config sorted_config] = load_LFP_NPX1(options,column);
                end

                if ~isempty(tvec)
                    LFP_tvec = tvec;
                else
                    LFP_tvec = [];
                end
            end


            for nprobe = 1:length(session_info(n).probe) % For each session, how many probes
                options = session_info(n).probe(nprobe);
                options.importMode = 'KS';
                options.gFileNum = gFileNum(n);
                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
                options.ROOTPATH = ROOTPATH;
                % Load all spike data sorted according to the channel position

                if ~isempty(best_channels{probe_no}.CA1_channel)
                    [CA1_clusters.probe(probe_no) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{probe_no}.CA1_channel-10 best_channels{probe_no}.CA1_channel+10],'group','by region','cell_exporer','on');
                    %                 CA1_place_fields.probe(nprobe) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,CA1_clusters.probe(nprobe),[]);
                end
                % 'Whole' HPC from 100 micron above CA1 cell layer to 1000 micron
                % below
                [HPC_clusters.probe(probe_no) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{probe_no}.CA1_channel-100 best_channels{probe_no}.CA1_channel+10],'group','by region','cell_exporer','on');
                %             HPC_place_fields.probe(nprobe) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,HPC_clusters.probe(nprobe),[]);
                % all V1 spike data
                if best_channels{probe_no}.first_in_brain_channel-100 > best_channels{probe_no}.CA1_channel+10
                    [V1_clusters.probe(probe_no) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{probe_no}.first_in_brain_channel-100 best_channels{probe_no}.first_in_brain_channel],'group','by region','cell_exporer','on');
                else
                    [V1_clusters.probe(probe_no) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{probe_no}.CA1_channel+11 best_channels{probe_no}.first_in_brain_channel],'group','by region','cell_exporer','on');
                end
            end
            save(sprintf('extracted_CA1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'CA1_clusters');
            save(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'HPC_clusters');
            save(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'V1_clusters');


            clear replay reactivations ripples
            %%%%%%%%%%%%%%%%%%
            % Ripple and candidate events detection
            %%%%%%%%%%%%%%%%%%
            for nprobe = 1:length(session_info(n).probe)

                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
                %                 Behavioural state detection
                speed_interp = interp1(Behaviour.tvec,Behaviour.speed,tvec','linear');
                speedTreshold = 1;
                if isfield(LFP(nprobe),'L4')
                    cortex_LFP = LFP(nprobe).L4;
                elseif isfield(LFP(nprobe),'L5')
                    cortex_LFP = LFP(nprobe).L5;
                elseif isfield(LFP(nprobe),'MEC')
                end
                if isfield(LFP(nprobe),'CA1')
                    CA1_LFP = LFP(nprobe).CA1;
                    cortex_channel = best_channels{nprobe}.L5_channel;
                    CA1_channel = best_channels{nprobe}.CA1_channel;
                    [freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
                        [tvec' cortex_LFP'],[tvec' CA1_LFP'],...
                        [tvec' speed_interp],speedTreshold);
                end
                %             behavioural_state.freezing = freezing;
                behavioural_state.quietWake = quietWake;
                behavioural_state.SWS = SWS;
                behavioural_state.REM = REM;
                behavioural_state.movement = movement;


                % Detect CA1 populational bursting events (Candidate events)
                zscore_min = 0;
                zscore_max = 3;

                channel_to_use = find(sorted_config.Channel == best_channels{probe_no}.CA1_channel);
                [replay.probe(probe_no),reactivations.probe(probe_no)] = detect_candidate_events_masa(tvec,raw_LFP{probe_no}(channel_to_use,:),...
                    CA1_clusters.probe(probe_no).MUA_zscore,[CA1_clusters.probe(probe_no).spike_id CA1_clusters.probe(probe_no).spike_times],Behaviour,zscore_min,zscore_max,options);

                if length(reactivations.probe(probe_no).onset) < 50
                    [replay.probe(probe_no),reactivations.probe(probe_no)] = detect_candidate_events_masa(tvec,raw_LFP{probe_no}(channel_to_use,:),...
                        HPC_clusters.probe(probe_no).MUA_zscore,[HPC_clusters.probe(probe_no).spike_id HPC_clusters.probe(probe_no).spike_times],Behaviour,zscore_min,zscore_max,options);
                end

                %
                %             [reactivations.probe(nprobe).awake_offset,reactivations.probe(nprobe).awake_index] = RestrictInts(reactivations.probe(nprobe).offset,behavioural_state.quietWake);
                %             reactivations.probe(nprobe).awake_onset = reactivations.probe(nprobe).onset(reactivations.probe(nprobe).awake_index);
                %             reactivations.probe(nprobe).awake_peaktimes = reactivations.probe(nprobe).peaktimes(reactivations.awake_index);

                % Detect CA1 ripple events
                channel_to_use = find(sorted_config.Channel == best_channels{probe_no}.CA1_channel);
                [ripples.probe(probe_no)] = FindRipples_masa(raw_LFP{probe_no}(channel_to_use,:)',tvec','behaviour',Behaviour,'minDuration',20,'durations',[30 200],'frequency',SR,...
                    'noise',raw_LFP{probe_no}(2,:)','passband',[125 300],'thresholds',[3 5],'show','off')

                close all

            end

            for nprobe = 1:length(session_info(n).probe)
                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

                for event = 1:length(ripples.probe(probe_no).onset)
                    ripples.probe(probe_no).speed(event) = mean(Behaviour.speed(find(Behaviour.sglxTime >= ripples.probe(probe_no).onset(event) & Behaviour.sglxTime <= ripples.probe(probe_no).offset(event))));
                end

                if contains(Stimulus_type,'RUN') % If reactivation events during lap running
                    [reactivations.probe(probe_no).T1_offset,reactivations.probe(probe_no).T1_index] = RestrictInts(reactivations.probe(probe_no).offset',[lap_times(1).start' lap_times(1).end']); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                    reactivations.probe(probe_no).T1_onset = reactivations.probe(probe_no).onset(reactivations.probe(probe_no).T1_index);
                    reactivations.probe(probe_no).T1_midpoint = reactivations.probe(probe_no).midpoint(reactivations.probe(probe_no).T1_index);

                    [reactivations.probe(probe_no).T2_offset,reactivations.probe(probe_no).T2_index] = RestrictInts(reactivations.probe(probe_no).offset',[lap_times(2).start' lap_times(2).end']); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                    reactivations.probe(probe_no).T2_onset = reactivations.probe(probe_no).onset(reactivations.probe(probe_no).T2_index);
                    reactivations.probe(probe_no).T2_midpoint = reactivations.probe(probe_no).midpoint(reactivations.probe(probe_no).T2_index);
                end

                if contains(Stimulus_type,'RUN') % If reactivation events during lap running
                    [ripples.probe(probe_no).T1_offset,ripples.probe(probe_no).T1_index] = RestrictInts(ripples.probe(probe_no).offset',[lap_times(1).start' lap_times(1).end']); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                    ripples.probe(probe_no).T1_onset = ripples.probe(probe_no).onset(ripples.probe(probe_no).T1_index);
                    ripples.probe(probe_no).T1_peaktimes = ripples.probe(probe_no).peaktimes(ripples.probe(probe_no).T1_index);

                    [ripples.probe(probe_no).T2_offset,ripples.probe(probe_no).T2_index] = RestrictInts(ripples.probe(probe_no).offset',[lap_times(2).start' lap_times(2).end']); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                    ripples.probe(probe_no).T2_onset = ripples.probe(probe_no).onset(ripples.probe(probe_no).T2_index);
                    ripples.probe(probe_no).T2_peaktimes = ripples.probe(probe_no).peaktimes(ripples.probe(probe_no).T2_index);
                end
            end

            save(sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks')),'replay','reactivations')

            delete(sprintf('extracted_ripples_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            save(sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks')),'ripples')
        end
    end
end


%% Just HPC decoding
clear all
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))


SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN'; % extract LFP during RUN
Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';

Stimulus_types_all = {'RUN','POST'};
% Stimulus_type = 'RUN'; % extract LFP during RUN
% Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
c = 1;

for epoch = 1:length(Stimulus_types_all)
    Stimulus_type= Stimulus_types_all{epoch};

    for nsession =1:length(experiment_info)
        tic
        session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        if isempty(session_info)
            continue
        end
        gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

        for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
            cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
            load best_channels
            load extracted_PSD
            load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('bonsai_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load('extracted_laps')
            column = 1;

            %             load(sprintf('extracted_CA1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'CA1_clusters');
            %             load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            %             save(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'V1_clusters');
            load(sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))

            % Reactivation decoding
            if length(session_info(n).probe) > 1
                %         load(sprintf('extracted_HPC_place_fields_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')))

                load('extracted_HPC_place_fields_combined')
                if exist(sprintf('extracted_HPC_clusters_combined%s.mat',erase(stimulus_name{n},'Masa2tracks'))) == 0
                    HPC_clusters_combined = combine_clusters_from_multiple_probes(HPC_clusters);
                    save(sprintf('extracted_HPC_clusters_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')),'HPC_clusters_combined')
                else
                    load(sprintf('extracted_HPC_clusters_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')))
                end
                clusters = HPC_clusters_combined;
                place_fields_BAYESIAN = HPC_place_fields_combined;
            else
                load('extracted_HPC_place_fields')
                load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
                probe_no = session_info(n).probe(1).probe_id + 1;

                clusters = HPC_clusters.probe(probe_no);
                place_fields_BAYESIAN = HPC_place_fields.probe(probe_no);
            end

            % HPC
            decoded_ripple_events_shuffled = [];
            decoded_ripple_events = [];
            timebin = 0.05;
            tic
            % Reactivation log odds decoding
            for mprobe = 1:length(session_info(n).probe)
                % probe for ripples
                ripple_probe_no = session_info(n).probe(mprobe).probe_id + 1;
                %                 ripple_probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere;

                replay = ripples.probe(ripple_probe_no);
                %                 replay.offset = replay.onset(event) + 0.1;

                %                     [place_fields_BAYESIAN,decoded_events,probability_ratio_original] = log_odds_reactivation_analysis(clusters,place_fields_BAYESIAN,replay,position,[],stimulus_name{n},timebin);

                for event = 1:length(replay.onset)
                    replay1.offset = replay.onset(event) + 0.6;
                    replay1.onset = replay.onset(event) - 0.6;
                    [place_fields_BAYESIAN,decoded_events,~,decoded_events_shuffled] = log_odds_reactivation_analysis(clusters,place_fields_BAYESIAN,replay1,position,[],stimulus_name{n},timebin);

                    decoded_ripple_events(ripple_probe_no).track(1).replay_events(event).replay = decoded_events(1).replay_events.replay;
                    decoded_ripple_events(ripple_probe_no).track(2).replay_events(event).replay = decoded_events(2).replay_events.replay;

                    decoded_ripple_events(ripple_probe_no).track(1).replay_events(event).summed_probability = sum(decoded_events(1).replay_events.replay);
                    decoded_ripple_events(ripple_probe_no).track(2).replay_events(event).summed_probability = sum(decoded_events(2).replay_events.replay);

                    decoded_ripple_events(ripple_probe_no).track(1).replay_events(event).spikes = decoded_events(1).replay_events.spikes;
                    decoded_ripple_events(ripple_probe_no).track(2).replay_events(event).spikes = decoded_events(2).replay_events.spikes;

                    decoded_ripple_events(ripple_probe_no).track(1).replay_events(event).timebins_edges = decoded_events(1).replay_events.timebins_edges;
                    decoded_ripple_events(ripple_probe_no).track(2).replay_events(event).timebins_edges = decoded_events(2).replay_events.timebins_edges;

                    for nshuffle = 1:1000
                        decoded_ripple_events_shuffled(ripple_probe_no).track(1).replay_events(event).summed_probability(nshuffle,:) = sum(decoded_events_shuffled{nshuffle}(1).replay);
                        decoded_ripple_events_shuffled(ripple_probe_no).track(2).replay_events(event).summed_probability(nshuffle,:) = sum(decoded_events_shuffled{nshuffle}(2).replay);
                    end

                end


                for event = 1:length(replay.onset)
                    replay1.offset = replay.onset(event) + 0.6;
                    replay1.onset = replay.onset(event) - 0.6;
                    [place_fields_BAYESIAN,decoded_events,~,decoded_events_shuffled] = log_odds_reactivation_analysis(clusters,place_fields_BAYESIAN,replay1,position,[],stimulus_name{n},timebin);

                    decoded_ripple_events_global_remapped(ripple_probe_no).track(1).replay_events(event).replay = decoded_events(1).replay_events.replay;
                    decoded_ripple_events_global_remapped(ripple_probe_no).track(2).replay_events(event).replay = decoded_events(2).replay_events.replay;

                    decoded_ripple_events_global_remapped(ripple_probe_no).track(1).replay_events(event).summed_probability = sum(decoded_events(1).replay_events.replay);
                    decoded_ripple_events_global_remapped(ripple_probe_no).track(2).replay_events(event).summed_probability = sum(decoded_events(2).replay_events.replay);

                    decoded_ripple_events_global_remapped(ripple_probe_no).track(1).replay_events(event).spikes = decoded_events(1).replay_events.spikes;
                    decoded_ripple_events_global_remapped(ripple_probe_no).track(2).replay_events(event).spikes = decoded_events(2).replay_events.spikes;

                    decoded_ripple_events_global_remapped(ripple_probe_no).track(1).replay_events(event).timebins_edges = decoded_events(1).replay_events.timebins_edges;
                    decoded_ripple_events_global_remapped(ripple_probe_no).track(2).replay_events(event).timebins_edges = decoded_events(2).replay_events.timebins_edges;

                    for nshuffle = 1:1000
                        decoded_ripple_events_shuffled_global_remapped(ripple_probe_no).track(1).replay_events(event).summed_probability(nshuffle,:) = sum(decoded_events_shuffled{nshuffle}(1).replay);
                        decoded_ripple_events_shuffled_global_remapped(ripple_probe_no).track(2).replay_events(event).summed_probability(nshuffle,:) = sum(decoded_events_shuffled{nshuffle}(2).replay);
                    end

                end
            end
            toc


            save(sprintf('decoded_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks')),'decoded_ripple_events');
            save(sprintf('decoded_ripple_events_global_remapped%s.mat',erase(stimulus_name{n},'Masa2tracks')),'decoded_ripple_events_global_remapped');
            save(sprintf('decoded_ripple_events_shuffled%s.mat',erase(stimulus_name{n},'Masa2tracks')),'decoded_ripple_events_shuffled');
            save(sprintf('decoded_ripple_events_shuffled_global_remapped%s.mat',erase(stimulus_name{n},'Masa2tracks')),'decoded_ripple_events_shuffled_global_remapped');
            clear decoded_ripple_events_global_remapped decoded_ripple_events_shuffled_global_remapped decoded_ripple_events decoded_ripple_events_shuffled
        end

    end
end

%% V1 reactivation
clear all
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))


SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN'; % extract LFP during RUN
Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';

Stimulus_types_all = {'RUN','POST'};
% Stimulus_type = 'RUN'; % extract LFP during RUN
% Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
c = 1;

for epoch = 1:length(Stimulus_types_all)
    Stimulus_type= Stimulus_types_all{epoch};

    for nsession =1:length(experiment_info)
        tic
        session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        if isempty(session_info)
            continue
        end
        gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

        for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
            cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
            load best_channels
            load extracted_PSD
            load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('bonsai_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load('extracted_laps')
            column = 1;

            load(sprintf('extracted_CA1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            load(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            %             load(sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            load(sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load('extracted_V1_place_fields.mat')


            % V1
            decoded_ripple_events_V1_shuffled = [];
            decoded_ripple_events_V1 = [];
            timebin = 0.05;
            tic
            % Reactivation log odds decoding
            for mprobe = 1:length(session_info(n).probe)
                % probe for ripples
                ripple_probe_no = session_info(n).probe(mprobe).probe_id + 1;
                %                 ripple_probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere;

                replay = ripples.probe(ripple_probe_no);
                %                 replay.offset = replay.onset(event) + 0.1;

                for nprobe = 1:length(session_info(n).probe)
                    options = session_info(n).probe(nprobe);
                    options.importMode = 'KS';
                    options.gFileNum = gFileNum(n);
                    probe_no = session_info(n).probe(nprobe).probe_id + 1;
                    options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
                    options.ROOTPATH = ROOTPATH;
                    clusters = V1_clusters.probe(probe_no);
                    place_fields_BAYESIAN = V1_place_fields.probe(probe_no);

                    %                     [place_fields_BAYESIAN,decoded_events,probability_ratio_original] = log_odds_reactivation_analysis(clusters,place_fields_BAYESIAN,replay,position,[],stimulus_name{n},timebin);

                    for event = 1:length(replay.onset)
                        replay1.offset = replay.onset(event) + 0.6;
                        replay1.onset = replay.onset(event) - 0.6;
                        [place_fields_BAYESIAN,decoded_events,~,decoded_events_shuffled] = log_odds_reactivation_analysis(clusters,place_fields_BAYESIAN,replay1,position,[],stimulus_name{n},timebin);

                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).replay = decoded_events(1).replay_events.replay;
                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).replay = decoded_events(2).replay_events.replay;

                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability = sum(decoded_events(1).replay_events.replay);
                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability = sum(decoded_events(2).replay_events.replay);

                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).spikes = decoded_events(1).replay_events.spikes;
                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).spikes = decoded_events(2).replay_events.spikes;

                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).timebins_edges = decoded_events(1).replay_events.timebins_edges;
                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).timebins_edges = decoded_events(2).replay_events.timebins_edges;

                        for nshuffle = 1:1000
                            decoded_ripple_events_V1_shuffled(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(nshuffle,:) = sum(decoded_events_shuffled{nshuffle}(1).replay);
                            decoded_ripple_events_V1_shuffled(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability(nshuffle,:) = sum(decoded_events_shuffled{nshuffle}(2).replay);
                        end

                    end

                end
            end
            toc


            save(sprintf('decoded_ripple_events_V1%s.mat',erase(stimulus_name{n},'Masa2tracks')),'decoded_ripple_events_V1');
            %             save(sprintf('probability_ratio_global_remapped_V1%s.mat',erase(stimulus_name{n},'Masa2tracks')),'probability_ratio_global_remapped_V1');
            save(sprintf('decoded_ripple_events_V1_shuffled%s.mat',erase(stimulus_name{n},'Masa2tracks')),'decoded_ripple_events_V1_shuffled');
            clear decoded_ripple_events_V1_shuffled decoded_ripple_events_V1
        end
    end
end




%% combined V1 reactivation
clear all
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))


SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN'; % extract LFP during RUN
Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';

Stimulus_types_all = {'RUN','POST'};
% Stimulus_type = 'RUN'; % extract LFP during RUN
% Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
c = 1;

for epoch = 1:length(Stimulus_types_all)
    Stimulus_type= Stimulus_types_all{epoch};

    for nsession =5:length(experiment_info)
        tic
        session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        if isempty(session_info)
            continue
        end
        gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

        for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
            cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
            load best_channels
            load extracted_PSD
            load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('bonsai_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load('extracted_laps')
            column = 1;

            %             load(sprintf('extracted_CA1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            %             load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            if exist(sprintf('extracted_V1_clusters_combined%s.mat',erase(stimulus_name{n},'Masa2tracks'))) ~= 0

                load(sprintf('extracted_V1_clusters_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')));
                load('extracted_V1_place_fields_combined.mat')
            else
                %             load(sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks')));
                continue
            end
            load(sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))

            % V1
            decoded_ripple_events_V1_combined_shuffled = [];
            decoded_ripple_events_V1_combined = [];
            timebin = 0.05;
            tic

            if length(session_info(n).probe)>1
                % Reactivation log odds decoding
                for mprobe = 1:length(session_info(n).probe)
                    % probe for ripples
                    ripple_probe_no = session_info(n).probe(mprobe).probe_id + 1;
                    %                 ripple_probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere;

                    replay = ripples.probe(ripple_probe_no);
                    %                 replay.offset = replay.onset(event) + 0.1;

                    options = session_info(n).probe(1);
                    options.importMode = 'KS';
                    options.gFileNum = gFileNum(n);
                    options.ROOTPATH = ROOTPATH;
                    clusters = V1_clusters_combined;
                    place_fields_BAYESIAN = V1_place_fields_combined;

                    %                     [place_fields_BAYESIAN,decoded_events,probability_ratio_original] = log_odds_reactivation_analysis(clusters,place_fields_BAYESIAN,replay,position,[],stimulus_name{n},timebin);

                    for event = 1:length(replay.onset)
                        replay1.offset = replay.onset(event) + 0.6;
                        replay1.onset = replay.onset(event) - 0.6;
                        [place_fields_BAYESIAN,decoded_events,~,decoded_events_shuffled] = log_odds_reactivation_analysis(clusters,place_fields_BAYESIAN,replay1,position,[],stimulus_name{n},timebin);

                        decoded_ripple_events_V1_combined(ripple_probe_no).track(1).replay_events(event).replay = decoded_events(1).replay_events.replay;
                        decodeddecoded_ripple_events_V1_combined_ripple_events(ripple_probe_no).track(2).replay_events(event).replay = decoded_events(2).replay_events.replay;

                        decoded_ripple_events_V1_combined(ripple_probe_no).track(1).replay_events(event).summed_probability = sum(decoded_events(1).replay_events.replay);
                        decoded_ripple_events_V1_combined(ripple_probe_no).track(2).replay_events(event).summed_probability = sum(decoded_events(2).replay_events.replay);

                        decoded_ripple_events_V1_combined(ripple_probe_no).track(1).replay_events(event).spikes = decoded_events(1).replay_events.spikes;
                        decoded_ripple_events_V1_combined(ripple_probe_no).track(2).replay_events(event).spikes = decoded_events(2).replay_events.spikes;

                        decoded_ripple_events_V1_combined(ripple_probe_no).track(1).replay_events(event).timebins_edges = decoded_events(1).replay_events.timebins_edges;
                        decoded_ripple_events_V1_combined(ripple_probe_no).track(2).replay_events(event).timebins_edges = decoded_events(2).replay_events.timebins_edges;

                        for nshuffle = 1:1000
                            decoded_ripple_events_V1_combined_shuffled(ripple_probe_no).track(1).replay_events(event).summed_probability(nshuffle,:) = sum(decoded_events_shuffled{nshuffle}(1).replay);
                            decoded_ripple_events_V1_combined_shuffled(ripple_probe_no).track(2).replay_events(event).summed_probability(nshuffle,:) = sum(decoded_events_shuffled{nshuffle}(2).replay);
                        end
                    end
                end
                toc


                save(sprintf('decoded_ripple_events_V1_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')),'decoded_ripple_events_V1_combined');
                %             save(sprintf('probability_ratio_global_remapped_V1%s.mat',erase(stimulus_name{n},'Masa2tracks')),'probability_ratio_global_remapped_V1');
                save(sprintf('decoded_ripple_events_V1_combined_shuffled%s.mat',erase(stimulus_name{n},'Masa2tracks')),'decoded_ripple_events_V1_combined_shuffled');
                clear decoded_ripple_events_V1_combined_shuffled decoded_ripple_events_V1_combined
            else
                disp('Less than 2 probes')
            end


        end
    end
end

%% Superficial reactivation
clear all
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))


SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN'; % extract LFP during RUN
Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';

Stimulus_types_all = {'RUN','POST'};
% Stimulus_type = 'RUN'; % extract LFP during RUN
% Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
c = 1;
x_bins_width = 10;

for epoch = 1:length(Stimulus_types_all)
    Stimulus_type= Stimulus_types_all{epoch};

    for nsession =1:length(experiment_info)
        tic
        session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        if isempty(session_info)
            continue
        end
        gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

        for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
            cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
            load best_channels
            load extracted_PSD
            load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('bonsai_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load('extracted_laps')
            column = 1;

            load(sprintf('extracted_CA1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            load(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            %             load(sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            load(sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load('extracted_V1_place_fields.mat')


            % sV1
            decoded_ripple_events_V1_shuffled = [];
            decoded_ripple_events_V1 = [];
            timebin = 0.05;
            tic
            for nprobe = 1:length(session_info(n).probe)
                options = session_info(n).probe(nprobe);
                options.importMode = 'KS';
                options.gFileNum = gFileNum(n);
                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
                options.ROOTPATH = ROOTPATH;

                if isempty(best_channels{probe_no}.L4_channel)
                    if isempty(best_channels{probe_no}.L5_channel)
                        channel_boundary = [best_channels{probe_no}.first_in_brain_channel-45 best_channels{probe_no}.first_in_brain_channel]; % 450 micron from the surface
                    else
                        channel_boundary = [best_channels{probe_no}.L5_channel+20 best_channels{probe_no}.first_in_brain_channel]; % 200 micron above L5
                    end
                else
                    if  best_channels{probe_no}.first_in_brain_channel- best_channels{probe_no}.L4_channel < 50
                        channel_boundary = [best_channels{probe_no}.L4_channel best_channels{probe_no}.first_in_brain_channel];
                    else
                        channel_boundary = [best_channels{probe_no}.first_in_brain_channel-45 best_channels{probe_no}.first_in_brain_channel]; % 450 micron from the surface
                    end
                end

                LFP_tvec = [];
                [sV1_clusters.probe(probe_no) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[channel_boundary(1) channel_boundary(2)],'group','by region','cell_exporer','on');

                if contains(stimulus_name{n},'RUN')
                    sV1_place_fields.probe(probe_no) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,sV1_clusters.probe(probe_no),[]);
                    save extracted_sV1_place_fields sV1_place_fields
                else
                    load('extracted_sV1_place_fields')
                end
                %

                clusters = sV1_clusters.probe(probe_no);
                place_fields_BAYESIAN = sV1_place_fields.probe(probe_no);


                % Reactivation log odds decoding
                for mprobe = 1:length(session_info(n).probe)
                    % probe for ripples
                    ripple_probe_no = session_info(n).probe(mprobe).probe_id + 1;
                    %                 ripple_probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere;

                    replay = ripples.probe(ripple_probe_no);
                    %                 replay.offset = replay.onset(event) + 0.1;


                    if isempty(place_fields_BAYESIAN.good_place_cells_LIBERAL)
                        decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track = [];
                        decoded_ripple_events_sV1_shuffled(ripple_probe_no).probe(probe_no).track  = [];

                        continue
                    end
                    %                     [place_fields_BAYESIAN,decoded_events,probability_ratio_original] = log_odds_reactivation_analysis(clusters,place_fields_BAYESIAN,replay,position,[],stimulus_name{n},timebin);

                    for event = 1:length(replay.onset)
                        replay1.offset = replay.onset(event) + 0.6;
                        replay1.onset = replay.onset(event) - 0.6;
                        [place_fields_BAYESIAN,decoded_events,~,decoded_events_shuffled] = log_odds_reactivation_analysis(clusters,place_fields_BAYESIAN,replay1,position,[],stimulus_name{n},timebin);

                        decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).replay = decoded_events(1).replay_events.replay;
                        decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).replay = decoded_events(2).replay_events.replay;

                        decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability = sum(decoded_events(1).replay_events.replay);
                        decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability = sum(decoded_events(2).replay_events.replay);

                        decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).spikes = decoded_events(1).replay_events.spikes;
                        decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).spikes = decoded_events(2).replay_events.spikes;

                        decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).timebins_edges = decoded_events(1).replay_events.timebins_edges;
                        decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).timebins_edges = decoded_events(2).replay_events.timebins_edges;

                        for nshuffle = 1:1000
                            decoded_ripple_events_sV1_shuffled(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(nshuffle,:) = sum(decoded_events_shuffled{nshuffle}(1).replay);
                            decoded_ripple_events_sV1_shuffled(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability(nshuffle,:) = sum(decoded_events_shuffled{nshuffle}(2).replay);
                        end

                    end

                end
            end
            toc


            save(sprintf('decoded_ripple_events_sV1%s.mat',erase(stimulus_name{n},'Masa2tracks')),'decoded_ripple_events_sV1');
            save(sprintf('decoded_ripple_events_sV1_shuffled%s.mat',erase(stimulus_name{n},'Masa2tracks')),'decoded_ripple_events_sV1_shuffled');
            save(sprintf('extracted_sV1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'sV1_clusters')

            clear decoded_ripple_events_sV1_shuffled decoded_ripple_events_sV1 sV1_clusters sV1_place_fields
        end
    end
end

%% deep V1 reactivation
clear all
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))


SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN'; % extract LFP during RUN
Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';

Stimulus_types_all = {'RUN','POST'};
% Stimulus_type = 'RUN'; % extract LFP during RUN
% Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
c = 1;
x_bins_width = 10;


for epoch = 1:length(Stimulus_types_all)
    Stimulus_type= Stimulus_types_all{epoch};

    for nsession =1:length(experiment_info)
        tic
        session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        if isempty(session_info)
            continue
        end
        gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

        for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
            cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
            load best_channels
            load extracted_PSD
            load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('bonsai_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load('extracted_laps')
            column = 1;

            load(sprintf('extracted_CA1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            load(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            %             load(sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            load(sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load('extracted_V1_place_fields.mat')


            % dV1
            decoded_ripple_events_dV1_shuffled = [];
            decoded_ripple_events_dV1 = [];
            timebin = 0.05;
            tic
            for nprobe = 1:length(session_info(n).probe)
                options = session_info(n).probe(nprobe);
                options.importMode = 'KS';
                options.gFileNum = gFileNum(n);
                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
                options.ROOTPATH = ROOTPATH;


                channel_boundary = [];
                if ~isempty(best_channels{probe_no}.L5_channel)
                    if ~isempty(best_channels{probe_no}.L4_channel)
                        if best_channels{probe_no}.L5_channel + 10 < best_channels{probe_no}.L4_channel-6
                            channel_boundary = [best_channels{probe_no}.CA1_channel+15 best_channels{probe_no}.L5_channel+10];
                        else
                            channel_boundary = [best_channels{probe_no}.CA1_channel+15 best_channels{probe_no}.L4_channel-6];
                        end
                    else
                        channel_boundary = [best_channels{probe_no}.CA1_channel+15 best_channels{probe_no}.L5_channel+10];
                    end
                else
                    if  ~isempty(best_channels{probe_no}.L4_channel)
                        channel_boundary = [best_channels{probe_no}.L4_channel-6 best_channels{probe_no}.CA1_channel+15];
                    else
                        channel_boundary = [best_channels{probe_no}.CA1_channel+15 best_channels{probe_no}.CA1_channel+55]; % 400 from 15 micron from the HPC
                    end
                end

                LFP_tvec = [];
                [dV1_clusters.probe(probe_no) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[channel_boundary(1) channel_boundary(2)],'group','by region','cell_exporer','on');

                if contains(stimulus_name{n},'RUN')
                    dV1_place_fields.probe(probe_no) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,dV1_clusters.probe(probe_no),[]);
                    save extracted_dV1_place_fields dV1_place_fields
                else
                    load('extracted_dV1_place_fields')
                end

                clusters = dV1_clusters.probe(probe_no);
                place_fields_BAYESIAN = dV1_place_fields.probe(probe_no);

                % Reactivation log odds decoding
                for mprobe = 1:length(session_info(n).probe)
                    % probe for ripples
                    ripple_probe_no = session_info(n).probe(mprobe).probe_id + 1;
                    %                 ripple_probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere;

                    replay = ripples.probe(ripple_probe_no);
                    %                 replay.offset = replay.onset(event) + 0.1;



                    if isempty(place_fields_BAYESIAN.good_place_cells_LIBERAL)
                        decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track = [];
                        decoded_ripple_events_dV1_shuffled(ripple_probe_no).probe(probe_no).track = [];

                        continue
                    end

                    %                     [place_fields_BAYESIAN,decoded_events,probability_ratio_original] = log_odds_reactivation_analysis(clusters,place_fields_BAYESIAN,replay,position,[],stimulus_name{n},timebin);

                    for event = 1:length(replay.onset)
                        replay1.offset = replay.onset(event) + 0.6;
                        replay1.onset = replay.onset(event) - 0.6;
                        [place_fields_BAYESIAN,decoded_events,~,decoded_events_shuffled] = log_odds_reactivation_analysis(clusters,place_fields_BAYESIAN,replay1,position,[],stimulus_name{n},timebin);

                        decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).replay = decoded_events(1).replay_events.replay;
                        decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).replay = decoded_events(2).replay_events.replay;

                        decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability = sum(decoded_events(1).replay_events.replay);
                        decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability = sum(decoded_events(2).replay_events.replay);

                        decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).spikes = decoded_events(1).replay_events.spikes;
                        decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).spikes = decoded_events(2).replay_events.spikes;

                        decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).timebins_edges = decoded_events(1).replay_events.timebins_edges;
                        decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).timebins_edges = decoded_events(2).replay_events.timebins_edges;

                        for nshuffle = 1:1000
                            decoded_ripple_events_dV1_shuffled(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(nshuffle,:) = sum(decoded_events_shuffled{nshuffle}(1).replay);
                            decoded_ripple_events_dV1_shuffled(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability(nshuffle,:) = sum(decoded_events_shuffled{nshuffle}(2).replay);
                        end

                    end

                end
            end
            toc

            save(sprintf('decoded_ripple_events_dV1%s.mat',erase(stimulus_name{n},'Masa2tracks')),'decoded_ripple_events_dV1');
            save(sprintf('decoded_ripple_events_dV1_shuffled%s.mat',erase(stimulus_name{n},'Masa2tracks')),'decoded_ripple_events_dV1_shuffled');

            save(sprintf('extracted_dV1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'dV1_clusters')

            clear decoded_ripple_events_dV1_shuffled decoded_ripple_events_dV1 dV1_clusters dV1_place_fields
        end
    end
end


%% log odds caluclation (with sV1 and dV1) for ripple events
clear all
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))

addpath(genpath('P:\corticohippocampal_replay\code\buzcode\externalPackages'))
addpath(genpath('P:\corticohippocampal_replay\code\spikes'))

SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_types_all = {'RUN','POST'};
% Stimulus_type = 'RUN'; % extract LFP during RUN
% Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
c = 1;
%
% time_windows = [1:10;11:20;21:30;...
%     31:40;41:50;[51:59 nan]]; % for V1 log odds
time_windows = [1:4;5:8;9:12;...
    13:16;17:20;[21:23 nan]]; % for V1 log odds


for epoch = 1:length(Stimulus_types_all)
    Stimulus_type= Stimulus_types_all{epoch};
    for nsession =1:length(experiment_info)
        tic
        session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

        if isempty(session_info)
            continue
        end

        for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
            cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
            load best_channels
            load extracted_PSD
            column = 1;
            load('extracted_laps.mat') % remember to run the behavioural analysis first
            load(sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))

            load(sprintf('decoded_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('decoded_ripple_events_shuffled%s.mat',erase(stimulus_name{n},'Masa2tracks')))

            load(sprintf('decoded_ripple_events_global_remapped%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('decoded_ripple_events_shuffled_global_remapped%s.mat',erase(stimulus_name{n},'Masa2tracks')))

            load(sprintf('decoded_ripple_events_V1%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('decoded_ripple_events_V1_shuffled%s.mat',erase(stimulus_name{n},'Masa2tracks')))

            load(sprintf('decoded_ripple_events_sV1%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('decoded_ripple_events_sV1_shuffled%s.mat',erase(stimulus_name{n},'Masa2tracks')))

            load(sprintf('decoded_ripple_events_dV1%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('decoded_ripple_events_dV1_shuffled%s.mat',erase(stimulus_name{n},'Masa2tracks')))

            if length(session_info(n).probe) > 1
                load(sprintf('decoded_ripple_events_V1_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')))
                load(sprintf('decoded_ripple_events_V1_combined_shuffled%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            end

            load('extracted_V1_place_fields.mat')
            load('extracted_CA1_place_fields.mat')
            load('extracted_HPC_place_fields.mat')
            load('extracted_sV1_place_fields.mat')
            load('extracted_dV1_place_fields.mat')

            load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('extracted_CA1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('extracted_dV1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('extracted_sV1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))

            % global remapped for subtract mean log odds bias
            z_log_odds = [];
            z_log_odds_global_remapped = [];
            V1_z_log_odds = [];
            sV1_z_log_odds = [];
            dV1_z_log_odds = [];

            HPC_bayesian_bias = [];
            HPC_bayesian_bias_global_remapped = [];
            V1_bayesian_bias = [];
            sV1_bayesian_bias = [];
            dV1_bayesian_bias = [];

            for mprobe = 1:length(session_info(n).probe)
                ripple_probe_no = session_info(n).probe(mprobe).probe_id + 1;
                ripple_probe_hemisphere = session_info(n).probe(mprobe).probe_hemisphere;

                for event = 1:length(decoded_ripple_events(ripple_probe_no).track(1).replay_events)

                    for time_bin = 1:23;% sum of 100ms after ripple onset

                        data = log(sum(decoded_ripple_events(ripple_probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                            /sum(decoded_ripple_events(ripple_probe_no).track(2).replay_events(event).summed_probability(:,time_bin)));

                        shuffled_data = log(sum(decoded_ripple_events_shuffled(ripple_probe_no).track(1).replay_events(event).summed_probability(:,time_bin),2)...
                            ./sum(decoded_ripple_events_shuffled(ripple_probe_no).track(2).replay_events(event).summed_probability(:,time_bin),2));


                        HPC_bayesian_bias{ripple_probe_hemisphere}(event,time_bin) = nansum(decoded_ripple_events(ripple_probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                            /(nansum(decoded_ripple_events(ripple_probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                            + nansum(decoded_ripple_events(ripple_probe_no).track(2).replay_events(event).summed_probability(:,time_bin)));

                        z_log_odds{ripple_probe_hemisphere}(event,time_bin) = (data-mean(shuffled_data))/std(shuffled_data);


                        % global remapped

                        data = log(sum(decoded_ripple_events_global_remapped(ripple_probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                            /sum(decoded_ripple_events_global_remapped(ripple_probe_no).track(2).replay_events(event).summed_probability(:,time_bin)));

                        shuffled_data = log(sum(decoded_ripple_events_shuffled(ripple_probe_no).track(1).replay_events(event).summed_probability(:,time_bin),2)...
                            ./sum(decoded_ripple_events_shuffled(ripple_probe_no).track(2).replay_events(event).summed_probability(:,time_bin),2));


                        HPC_bayesian_bias_global_remapped{ripple_probe_hemisphere}(event,time_bin) = nansum(decoded_ripple_events_global_remapped(ripple_probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                            /(nansum(decoded_ripple_events_global_remapped(ripple_probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                            + nansum(decoded_ripple_events_global_remapped(ripple_probe_no).track(2).replay_events(event).summed_probability(:,time_bin)));

                        z_log_odds_global_remapped{ripple_probe_hemisphere}(event,time_bin) = (data-mean(shuffled_data))/std(shuffled_data);

                        if length(session_info(n).probe) > 1

                            % V1 combined log odds
                            data = log(sum(decoded_ripple_events_V1_combined(ripple_probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                                /sum(decoded_ripple_events_V1_combined(ripple_probe_no).track(2).replay_events(event).summed_probability(:,time_bin)));

                            shuffled_data = log(sum(decoded_ripple_events_shuffled(ripple_probe_no).track(1).replay_events(event).summed_probability(:,time_bin),2)...
                                ./sum(decoded_ripple_events_shuffled(ripple_probe_no).track(2).replay_events(event).summed_probability(:,time_bin),2));


                            V1_combined_bayesian_bias{ripple_probe_hemisphere}(event,time_bin) = nansum(decoded_ripple_events_V1_combined(ripple_probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                                /(nansum(decoded_ripple_events_V1_combined(ripple_probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                                + nansum(decoded_ripple_events_V1_combined(ripple_probe_no).track(2).replay_events(event).summed_probability(:,time_bin)));

                            V1_combined_z_log_odds{ripple_probe_hemisphere}(event,time_bin) = (data-mean(shuffled_data))/std(shuffled_data);
                        end
                    end

                    % V1 log odds for each probe
                    for nprobe = 1:length(session_info(n).probe)
                        probe_no = session_info(n).probe(nprobe).probe_id + 1;
                        probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere ;
                        for nbin = 1:23
                            time_bin = nbin;

                            % V1
                            data = log(sum(decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                                /sum(decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability(:,time_bin)));

                            shuffled_data = log(sum(decoded_ripple_events_V1_shuffled(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin),2)...
                                ./sum(decoded_ripple_events_V1_shuffled(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability(:,time_bin),2));


                            V1_bayesian_bias{ripple_probe_hemisphere}{probe_hemisphere}(event,nbin) = nansum(decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                                /(nansum(decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                                + nansum(decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability(:,time_bin)));

                            V1_z_log_odds{ripple_probe_hemisphere}{probe_hemisphere}(event,nbin) = (data-mean(shuffled_data))/std(shuffled_data);


                            % sV1
                            if ~isempty(decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track)
                                data = log(sum(decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                                    /sum(decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability(:,time_bin)));

                                shuffled_data = log(sum(decoded_ripple_events_sV1_shuffled(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin),2)...
                                    ./sum(decoded_ripple_events_sV1_shuffled(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability(:,time_bin),2));


                                sV1_bayesian_bias{ripple_probe_hemisphere}{probe_hemisphere}(event,nbin) = nansum(decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                                    /(nansum(decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                                    + nansum(decoded_ripple_events_sV1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability(:,time_bin)));

                                sV1_z_log_odds{ripple_probe_hemisphere}{probe_hemisphere}(event,nbin) = (data-mean(shuffled_data))/std(shuffled_data);
                            else
                                sV1_bayesian_bias{ripple_probe_hemisphere}{probe_hemisphere}(event,nbin) = nan;
                                sV1_z_log_odds{ripple_probe_hemisphere}{probe_hemisphere}(event,nbin) = nan;
                            end

                            % dV1
                            if ~isempty(decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track)
                                data = log(sum(decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                                    /sum(decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability(:,time_bin)));

                                shuffled_data = log(sum(decoded_ripple_events_dV1_shuffled(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin),2)...
                                    ./sum(decoded_ripple_events_dV1_shuffled(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability(:,time_bin),2));


                                dV1_bayesian_bias{ripple_probe_hemisphere}{probe_hemisphere}(event,nbin) = nansum(decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                                    /(nansum(decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                                    + nansum(decoded_ripple_events_dV1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability(:,time_bin)));

                                dV1_z_log_odds{ripple_probe_hemisphere}{probe_hemisphere}(event,nbin) = (data-mean(shuffled_data))/std(shuffled_data);
                            else
                                dV1_bayesian_bias{ripple_probe_hemisphere}{probe_hemisphere}(event,nbin) = nan;
                                dV1_z_log_odds{ripple_probe_hemisphere}{probe_hemisphere}(event,nbin) = nan;
                            end

                        end
                    end

                end
            end

            cell_index_L = [];
            cell_index_R = [];
            all_spike_data = [];
            for nprobe = 1:length(session_info(n).probe)
                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere ;

                % if left hemisphere find good cells with higher FR on T2
                if probe_hemisphere == 1
                    %                     [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                    %                         setdiff(V1_place_fields.probe(probe_no).track(2).good_cells_LIBERAL,V1_place_fields.probe(probe_no).track(1).good_cells_LIBERAL));
                    %                     [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                    %                         V1_place_fields.probe(probe_no).track(2).good_cells_LIBERAL);

                    cell_index_L = find((V1_place_fields.probe(probe_no).track(1).mean_rate_track - V1_place_fields.probe(probe_no).track(2).mean_rate_track)./...
                        (V1_place_fields.probe(probe_no).track(1).mean_rate_track + V1_place_fields.probe(probe_no).track(2).mean_rate_track) < 0);
                    [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                        cell_index_L);
                    all_spike_data{probe_hemisphere}{1} = [spike_id spike_times];

                    cell_index_L = find((sV1_place_fields.probe(probe_no).track(1).mean_rate_track - sV1_place_fields.probe(probe_no).track(2).mean_rate_track)./...
                        (sV1_place_fields.probe(probe_no).track(1).mean_rate_track + sV1_place_fields.probe(probe_no).track(2).mean_rate_track) < 0);
                    [spike_times,spike_id]= select_spikes_subset(sV1_clusters.probe(probe_no),...
                        cell_index_L);
                    all_spike_data{probe_hemisphere}{2} = [spike_id spike_times];

                    cell_index_L = find((dV1_place_fields.probe(probe_no).track(1).mean_rate_track - dV1_place_fields.probe(probe_no).track(2).mean_rate_track)./...
                        (dV1_place_fields.probe(probe_no).track(1).mean_rate_track + dV1_place_fields.probe(probe_no).track(2).mean_rate_track) < 0);
                    [spike_times,spike_id]= select_spikes_subset(dV1_clusters.probe(probe_no),...
                        cell_index_L);
                    all_spike_data{probe_hemisphere}{3} = [spike_id spike_times];


                    %                     [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                    %                         find((V1_place_fields.probe(probe_no).track(2).raw_peak - V1_place_fields.probe(probe_no).track(1).raw_peak)>0));


                    % %                                         spike_id = V1_clusters.probe(probe_no).spike_id;
                    % %                                         spike_times = V1_clusters.probe(probe_no).spike_times;

                    % if Right hemisphere find good cells with higher FR on T1
                elseif probe_hemisphere == 2
                    %
                    %                     [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                    %                         setdiff(V1_place_fields.probe(probe_no).track(1).good_cells_LIBERAL,V1_place_fields.probe(probe_no).track(2).good_cells_LIBERAL));
                    %                     [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                    %                         V1_place_fields.probe(probe_no).track(1).good_cells_LIBERAL);

                    cell_index_R = find((V1_place_fields.probe(probe_no).track(1).mean_rate_track - V1_place_fields.probe(probe_no).track(2).mean_rate_track)./...
                        (V1_place_fields.probe(probe_no).track(1).mean_rate_track + V1_place_fields.probe(probe_no).track(1).mean_rate_track) > 0);
                    [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                        cell_index_R);
                    all_spike_data{probe_hemisphere}{1} = [spike_id spike_times];

                    cell_index_R = find((sV1_place_fields.probe(probe_no).track(1).mean_rate_track - sV1_place_fields.probe(probe_no).track(2).mean_rate_track)./...
                        (sV1_place_fields.probe(probe_no).track(1).mean_rate_track + sV1_place_fields.probe(probe_no).track(1).mean_rate_track) > 0);
                    [spike_times,spike_id]= select_spikes_subset(sV1_clusters.probe(probe_no),...
                        cell_index_R);
                    all_spike_data{probe_hemisphere}{2} = [spike_id spike_times];

                    cell_index_R = find((dV1_place_fields.probe(probe_no).track(1).mean_rate_track - dV1_place_fields.probe(probe_no).track(2).mean_rate_track)./...
                        (dV1_place_fields.probe(probe_no).track(1).mean_rate_track + dV1_place_fields.probe(probe_no).track(1).mean_rate_track) > 0);
                    [spike_times,spike_id]= select_spikes_subset(dV1_clusters.probe(probe_no),...
                        cell_index_R);
                    all_spike_data{probe_hemisphere}{3} = [spike_id spike_times];
                    %        [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                    %                               find((V1_place_fields.probe(probe_no).track(1).raw_peak - V1_place_fields.probe(probe_no).track(2).raw_peak)>0));

                    %                                         spike_id = V1_clusters.probe(probe_no).spike_id;
                    %                                         spike_times = V1_clusters.probe(probe_no).spike_times;
                end

                [spike_times,spike_id]= select_spikes_subset(CA1_clusters.probe(probe_no),...
                    CA1_place_fields.probe(probe_no).good_place_cells_LIBERAL);
                all_spike_data{probe_hemisphere}{4} = [spike_id spike_times];

                [spike_times,spike_id]= select_spikes_subset(HPC_clusters.probe(probe_no),...
                    HPC_place_fields.probe(probe_no).good_place_cells_LIBERAL);
                all_spike_data{probe_hemisphere}{5} = [spike_id spike_times];
                group_name = {'V1','sV1','dV1','CA1','HPC'};
            end

            if length(session_info(n).probe) == 1 % if one probe, but it can be probe 2
                probe_no = session_info(n).probe(1).probe_id + 1;
                %                 replay = reactivations.probe(probe_no); % using probe 1 bursting for now
                probe_hemisphere = session_info(n).probe(1).probe_hemisphere;
                %                 events_probe_hemisphere = probe_hemisphere;
            else
                probe_hemisphere(1) = session_info(n).probe(1).probe_hemisphere;
                probe_hemisphere(2) = session_info(n).probe(2).probe_hemisphere;
            end

            for mprobe = 1:length(session_info(n).probe) %
                ripple_probe_no = session_info(n).probe(mprobe).probe_id + 1;
                ripple_probe_hemisphere = session_info(n).probe(mprobe).probe_hemisphere;
                replay = ripples.probe(ripple_probe_no);

                for event = 1:length(replay.onset) % For ripple events from one of the two probes

                    log_odds.index(c) = c;

                    %                     lap_times
                    log_odds.event_probe_hemisphere(c) = ripple_probe_hemisphere;
                    log_odds.experiment(c) = nsession;
                    log_odds.duration(c) = replay.offset(event) - replay.onset(event);
                    log_odds.onset(c) = replay.onset(event);
                    log_odds.offset(c) = replay.offset(event);
                    log_odds.speed(c) = replay.speed(event);

                    if isfield(replay,'T1_index') |isfield(replay,'T2_index')
                        if replay.T1_index(event) == 1
                            log_odds.behavioural_state(c) = 1;

                        elseif replay.T2_index(event) == 1
                            log_odds.behavioural_state(c) = 2;
                        else
                            log_odds.behavioural_state(c) = -1;
                        end
                    elseif contains(stimulus_name,'POST')
                        log_odds.behavioural_state(c) = 0;
                    elseif contains(stimulus_name,'PRE')
                        log_odds.behavioural_state(c) = -1;
                    end

                    if log_odds.behavioural_state(c) > 0
                        track_id = log_odds.behavioural_state(c);

                        log_odds.lap_id(c) = find(log_odds.onset(c)>=lap_times(track_id).start & log_odds.onset(c) <= lap_times(track_id).end);
                        log_odds.trial_type(c) = lap_times(track_id).lap(log_odds.lap_id(c)).trial_type;

                        if ~isempty(lap_times(track_id).lap(log_odds.lap_id(c)).lick_ratio)
                            log_odds.anticipatory_lick_ratio(c) = lap_times(track_id).lap(log_odds.lap_id(c)).lick_ratio;% how biased is the lick before reward zone towards the right direction
                        else
                            log_odds.anticipatory_lick_ratio(c) = 0;
                        end

                        if ~isempty(lap_times(track_id).lap(log_odds.lap_id(c)).lick_ratio_all)
                            log_odds.lick_ratio_all(c) = lap_times(track_id).lap(log_odds.lap_id(c)).lick_ratio_all;% overall lick bias including reward consumption
                        else
                            log_odds.lick_ratio_all(c) = 0; % 0 here means no lick at all
                        end

                    else % if PRE or POST
                        log_odds.lick_ratio_all(c) = nan; % no licking
                        log_odds.anticipatory_lick_ratio(c) = nan; % no licking
                        log_odds.trial_type(c) = nan;
                        log_odds.lap_id(c) = nan;
                    end


                    log_odds.HPC_bayesian_bias(c,:) = HPC_bayesian_bias{ripple_probe_hemisphere}(event,:);
                    log_odds.zscore(c,:) = z_log_odds{ripple_probe_hemisphere}(event,:);
                    log_odds.zscore_global_remapped(c,:) = z_log_odds_global_remapped{ripple_probe_hemisphere}(event,:);
                    log_odds.zscore_global_remapped_substracted(c,:) = z_log_odds_global_remapped{ripple_probe_hemisphere}(event,:)-mean(z_log_odds_global_remapped{ripple_probe_hemisphere});
                    log_odds.zscore_substracted(c,:) = z_log_odds{ripple_probe_hemisphere}(event)-mean(z_log_odds{ripple_probe_hemisphere});
                    log_odds.ripple_peak(c) = replay.peak_zscore (event);

                    if length(session_info(n).probe)>1
                        log_odds.V1_bayesian_bias_combined(c,:) = V1_combined_bayesian_bias{ripple_probe_hemisphere}(event,:);
                        log_odds.V1_log_odds_combined(c,:) = V1_combined_z_log_odds{ripple_probe_hemisphere}(event,:);
                    else
                        log_odds.V1_bayesian_bias_combined(c,:) = nan(1,23);
                        log_odds.V1_log_odds_combined(c,:) =  nan(1,23);
                    end

                    log_odds.V1_spike_count_L(c,:) = nan(1,100);
                    log_odds.sV1_spike_count_L(c,:) = nan(1,100);
                    log_odds.dV1_spike_count_L(c,:) = nan(1,100);

                    log_odds.CA1_spike_count_L(c,:) = nan(1,100);
                    log_odds.HPC_spike_count_L(c,:) = nan(1,100);

                    log_odds.V1_spike_count_R(c,:) = nan(1,100);
                    log_odds.sV1_spike_count_R(c,:) = nan(1,100);
                    log_odds.dV1_spike_count_R(c,:) = nan(1,100);

                    log_odds.CA1_spike_count_R(c,:) = nan(1,100);
                    log_odds.HPC_spike_count_R(c,:) = nan(1,100);
                    log_odds.HPC_spike_id_L{c} = [];
                    log_odds.HPC_spike_id_R{c} = [];

                    log_odds.V1_participation_L(c,:) = nan(1,10);
                    log_odds.V1_participation_R(c,:) = nan(1,10);
                    log_odds.V1_log_odds_L(c,:) = nan(1,23);
                    log_odds.V1_log_odds_R(c,:) = nan(1,23);
                    log_odds.V1_bayesian_bias_L(c,:) = nan(1,23);
                    log_odds.V1_bayesian_bias_R(c,:) = nan(1,23);
                    log_odds.V1_cell_spike_count_L{c} = [];
                    log_odds.V1_cell_spike_count_R{c} = [];
                    log_odds.V1_spike_id_L{c} = [];
                    log_odds.V1_spike_id_R{c} = [];


                    log_odds.sV1_participation_L(c,:) = nan(1,10);
                    log_odds.sV1_participation_R(c,:) = nan(1,10);
                    log_odds.sV1_log_odds_L(c,:) = nan(1,23);
                    log_odds.sV1_log_odds_R(c,:) = nan(1,23);
                    log_odds.sV1_bayesian_bias_L(c,:) = nan(1,23);
                    log_odds.sV1_bayesian_bias_R(c,:) = nan(1,23);
                    log_odds.sV1_cell_spike_count_L{c} = [];
                    log_odds.sV1_cell_spike_count_R{c} = [];
                    log_odds.sV1_spike_id_L{c} = [];
                    log_odds.sV1_spike_id_R{c} = [];

                    log_odds.dV1_participation_L(c,:) = nan(1,10);
                    log_odds.dV1_participation_R(c,:) = nan(1,10);
                    log_odds.dV1_log_odds_L(c,:) = nan(1,23);
                    log_odds.dV1_log_odds_R(c,:) = nan(1,23);
                    log_odds.dV1_bayesian_bias_L(c,:) = nan(1,23);
                    log_odds.dV1_bayesian_bias_R(c,:) = nan(1,23);
                    log_odds.dV1_cell_spike_count_L{c} = [];
                    log_odds.dV1_cell_spike_count_R{c} = [];
                    log_odds.dV1_spike_id_L{c} = [];
                    log_odds.dV1_spike_id_R{c} = [];



                    for nprobe = 1:length(session_info(n).probe)
                        %                     PSTH = plot_perievent_spiketime_histogram(all_spike_data{probe_hemisphere(nprobe)},replay.onset(event),'group','by cell zscore',...
                        %                         'group_name',group_name,'twin',[-1 1],'bin_size',0.1,'plot_option',0,'smooth_option',0);
                        probe_no = session_info(n).probe(nprobe).probe_id + 1;
                        probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere ;

                        PSTH = plot_perievent_spiketime_histogram(all_spike_data{probe_hemisphere},replay.onset(event),'group','by cell zscore',...
                            'group_name',group_name,'twin',[-1 1],'bin_size',0.02,'plot_option',0,'smooth_option',1);

                        spike_count_this_event = [];
                        for nregion = 1:3
                            spike_count_this_event{nregion} = squeeze(PSTH(nregion).spike_count);

                            if size(spike_count_this_event{nregion},2) == 1
                                spike_count_this_event{nregion} = spike_count_this_event{nregion}';
                            end
                        end

                        % Empty variables
                        cell_spike_count_this_event = [];
                        participation_this_event = [];

                        %                     for k = 1:length(PSTH)
                        %                         % Remove non-spiking cells for average spiking
                        %                         PSTH(k).cell_psth(find(sum( PSTH(k).cell_psth,2) == 0),:) = [];
                        %                     end

                        if probe_hemisphere == 1
                            %                     mean(zscore(squeeze(PSTH(1).spike_count),0,2));
                            log_odds.V1_total_spike_count_L(c,:) = sum(spike_count_this_event{1},1);

                            log_odds.V1_spike_count_L(c,:) = nanmean(normalize(PSTH(1).cell_psth,2,'zscore'));
                            log_odds.sV1_spike_count_L(c,:) = nanmean(normalize(PSTH(2).cell_psth,2,'zscore'));
                            log_odds.dV1_spike_count_L(c,:) = nanmean(normalize(PSTH(3).cell_psth,2,'zscore'));


                            log_odds.CA1_spike_count_L(c,:) = nanmean(normalize(PSTH(4).cell_psth,2,'zscore'));
                            log_odds.HPC_spike_count_L(c,:) = nanmean(normalize(PSTH(5).cell_psth,2,'zscore'));

                            spikes_this_event = find(all_spike_data{probe_hemisphere}{5}(:,2) >= log_odds.onset(c)-0.6 ...
                                & all_spike_data{probe_hemisphere}{5}(:,2) <= log_odds.onset(c)+0.6);
                            log_odds.HPC_spike_id_L{c} = all_spike_data{probe_hemisphere}{5}(spikes_this_event,:);


                            log_odds.V1_log_odds_L(c,:) = V1_z_log_odds{ripple_probe_hemisphere}{probe_hemisphere}(event,:);
                            log_odds.V1_bayesian_bias_L(c,:) = V1_bayesian_bias{ripple_probe_hemisphere}{probe_hemisphere}(event,:);

                            log_odds.sV1_log_odds_L(c,:) = sV1_z_log_odds{ripple_probe_hemisphere}{probe_hemisphere}(event,:);
                            log_odds.sV1_bayesian_bias_L(c,:) = sV1_bayesian_bias{ripple_probe_hemisphere}{probe_hemisphere}(event,:);

                            log_odds.dV1_log_odds_L(c,:) = dV1_z_log_odds{ripple_probe_hemisphere}{probe_hemisphere}(event,:);
                            log_odds.dV1_bayesian_bias_L(c,:) = dV1_bayesian_bias{ripple_probe_hemisphere}{probe_hemisphere}(event,:);

                            spikes_this_event = find(all_spike_data{probe_hemisphere}{1}(:,2) >= log_odds.onset(c)-0.6 ...
                                & all_spike_data{probe_hemisphere}{1}(:,2) <= log_odds.onset(c)+0.6);
                            log_odds.V1_spike_id_L{c} = all_spike_data{probe_hemisphere}{1}(spikes_this_event,:);

                            if ~isempty(all_spike_data{probe_hemisphere}{2})
                                spikes_this_event = find(all_spike_data{probe_hemisphere}{2}(:,2) >= log_odds.onset(c)-0.6 ...
                                    & all_spike_data{probe_hemisphere}{2}(:,2) <= log_odds.onset(c)+0.6);
                                log_odds.sV1_spike_id_L{c} = all_spike_data{probe_hemisphere}{2}(spikes_this_event,:);
                            else
                                log_odds.sV1_spike_id_L{c} = [];
                            end

                            if ~isempty(all_spike_data{probe_hemisphere}{3})
                                spikes_this_event = find(all_spike_data{probe_hemisphere}{3}(:,2) >= log_odds.onset(c)-0.6 ...
                                    & all_spike_data{probe_hemisphere}{3}(:,2) <= log_odds.onset(c)+0.6);
                                log_odds.dV1_spike_id_L{c} = all_spike_data{probe_hemisphere}{3}(spikes_this_event,:);
                            else
                                log_odds.dV1_spike_id_L{c} = [];
                            end


                            for nregion = 1:3
                                time_window_c = 1;
                                if isempty(spike_count_this_event{nregion})
                                    continue
                                end
                                for nwin = 1:20
                                    cell_spike_count_this_event{nregion}(:,nwin) = sum(spike_count_this_event{nregion}(:,time_window_c:time_window_c+4),2);
                                    time_window_c = time_window_c + 5;
                                end
                            end

                            log_odds.V1_cell_spike_count_L{c} = cell_spike_count_this_event{1};
                            if ~isempty(spike_count_this_event{2})
                                log_odds.sV1_cell_spike_count_L{c} = cell_spike_count_this_event{2};
                            end

                            if ~isempty(spike_count_this_event{3})
                                log_odds.dV1_cell_spike_count_L{c} = cell_spike_count_this_event{3};
                            end
                            %                          PSTH = plot_perievent_spiketime_histogram(all_spike_data{probe_hemisphere(nprobe)},replay.onset(event),'group','by cell zscore',...
                            %                              'group_name',group_name,'twin',[-1 1],'bin_size',0.1,'plot_option',0,'smooth_option',0);


                            for nregion = 1:3
                                time_window_c = 1;
                                participation_this_event{nregion} = nan;
                                if isempty(spike_count_this_event{nregion})
                                    continue
                                end

                                if size(spike_count_this_event{nregion},1) > 2
                                    for nwin = 1:10
                                        participation_this_event{nregion}(nwin) = sum(sum(spike_count_this_event{nregion}(:,time_window_c:time_window_c+9),2)>0)/size(spike_count_this_event{nregion},1);
                                        time_window_c = time_window_c + 10;
                                    end
                                end
                            end
                            log_odds.V1_participation_L(c,:) = participation_this_event{1};
                            log_odds.sV1_participation_L(c,:) = participation_this_event{2};
                            log_odds.dV1_participation_L(c,:) = participation_this_event{3};

                        elseif probe_hemisphere == 2

                            log_odds.V1_spike_count_R(c,:) = nanmean(normalize(PSTH(1).cell_psth,2,'zscore'));
                            log_odds.V1_total_spike_count_R(c,:) = sum(spike_count_this_event{1},1);
                            log_odds.sV1_spike_count_R(c,:) = nanmean(normalize(PSTH(2).cell_psth,2,'zscore'));
                            log_odds.dV1_spike_count_R(c,:) = nanmean(normalize(PSTH(3).cell_psth,2,'zscore'));

                            log_odds.CA1_spike_count_R(c,:) = nanmean(normalize(PSTH(4).cell_psth,2,'zscore'));
                            log_odds.HPC_spike_count_R(c,:) = nanmean(normalize(PSTH(5).cell_psth,2,'zscore'));

                            spikes_this_event = find(all_spike_data{probe_hemisphere}{5}(:,2) >= log_odds.onset(c)-0.6 ...
                                & all_spike_data{probe_hemisphere}{5}(:,2) <= log_odds.onset(c)+0.6);
                            log_odds.HPC_spike_id_R{c} = all_spike_data{probe_hemisphere}{5}(spikes_this_event,:);


                            log_odds.V1_log_odds_R(c,:) = V1_z_log_odds{ripple_probe_hemisphere}{probe_hemisphere}(event,:);
                            log_odds.V1_bayesian_bias_R(c,:) = V1_bayesian_bias{ripple_probe_hemisphere}{probe_hemisphere}(event,:);

                            log_odds.sV1_log_odds_R(c,:) = sV1_z_log_odds{ripple_probe_hemisphere}{probe_hemisphere}(event,:);
                            log_odds.sV1_bayesian_bias_R(c,:) = sV1_bayesian_bias{ripple_probe_hemisphere}{probe_hemisphere}(event,:);

                            log_odds.dV1_log_odds_R(c,:) = dV1_z_log_odds{ripple_probe_hemisphere}{probe_hemisphere}(event,:);
                            log_odds.dV1_bayesian_bias_R(c,:) = dV1_bayesian_bias{ripple_probe_hemisphere}{probe_hemisphere}(event,:);

                            spikes_this_event = find(all_spike_data{probe_hemisphere}{1}(:,2) >= log_odds.onset(c)-0.6 ...
                                & all_spike_data{probe_hemisphere}{1}(:,2) <= log_odds.onset(c)+0.6);
                            log_odds.V1_spike_id_R{c} = all_spike_data{probe_hemisphere}{1}(spikes_this_event,:);

                            if ~isempty(all_spike_data{probe_hemisphere}{3})
                                spikes_this_event = find(all_spike_data{probe_hemisphere}{2}(:,2) >= log_odds.onset(c)-0.6 ...
                                    & all_spike_data{probe_hemisphere}{2}(:,2) <= log_odds.onset(c)+0.6);
                                log_odds.sV1_spike_id_R{c} = all_spike_data{probe_hemisphere}{2}(spikes_this_event,:);
                            end

                            if ~isempty(all_spike_data{probe_hemisphere}{3})
                                spikes_this_event = find(all_spike_data{probe_hemisphere}{3}(:,2) >= log_odds.onset(c)-0.6 ...
                                    & all_spike_data{probe_hemisphere}{3}(:,2) <= log_odds.onset(c)+0.6);
                                log_odds.dV1_spike_id_R{c} = all_spike_data{probe_hemisphere}{3}(spikes_this_event,:);
                            end

                            for nregion = 1:3
                                time_window_c = 1;
                                if isempty(spike_count_this_event{nregion})
                                    continue
                                end
                                for nwin = 1:20
                                    cell_spike_count_this_event{nregion}(:,nwin) = sum(spike_count_this_event{nregion}(:,time_window_c:time_window_c+4),2);
                                    time_window_c = time_window_c + 5;
                                end
                            end

                            log_odds.V1_cell_spike_count_R{c} = cell_spike_count_this_event{1};
                            if ~isempty(spike_count_this_event{2})
                                log_odds.sV1_cell_spike_count_R{c} = cell_spike_count_this_event{2};
                            end

                            if ~isempty(spike_count_this_event{3})
                                log_odds.dV1_cell_spike_count_R{c} = cell_spike_count_this_event{3};
                            end
                            %                          PSTH = plot_perievent_spiketime_histogram(all_spike_data{probe_hemisphere(nprobe)},replay.onset(event),'group','by cell zscore',...
                            %                              'group_name',group_name,'twin',[-1 1],'bin_size',0.1,'plot_option',0,'smooth_option',0);


                            for nregion = 1:3
                                time_window_c = 1;
                                participation_this_event{nregion} = nan;
                                if isempty(spike_count_this_event{nregion})
                                    continue
                                end
                                if size(spike_count_this_event{nregion},1) > 2
                                    for nwin = 1:10
                                        participation_this_event{nregion}(nwin) = sum(sum(spike_count_this_event{nregion}(:,time_window_c:time_window_c+9),2)>0)/size(spike_count_this_event{nregion},1);
                                        time_window_c = time_window_c + 10;
                                    end
                                end
                            end
                            log_odds.V1_participation_R(c,:) = participation_this_event{1};
                            log_odds.sV1_participation_R(c,:) = participation_this_event{2};
                            log_odds.dV1_participation_R(c,:) = participation_this_event{3};

                        end
                    end
                    c = c + 1;
                end
            end


        end
        toc
    end
end

cd('Z:\ibn-vision\USERS\Masa\V1_HPC_project')
% save('log_odds_best_probe_events',"log_odds")
% save('log_odds_re_zscored',"log_odds")
% save('log_odds_50ms',"log_odds")
save('log_odds_ripples_50ms_all',"log_odds")
% save('log_odds_50ms_ratio_threshold_0.3',"log_odds")

%% Visualisation of V1 and HPC spiking and bayesian decoding and log odds around each ripple event with high and low z-score
clear all
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))
% addpath(genpath('P:\corticohippocampal_replay\code\spikes'))
% addpath(genpath('P:\corticohippocampal_replay\code\spikes'))
cd('Z:\ibn-vision\USERS\Masa\V1_HPC_project')
SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
% Stimulus_types_all = {'RUN'};
Stimulus_type = 'RUN'; % extract LFP during RUN
% Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
load('log_odds_ripples_50ms_all')

for nsession =1:10
    tic
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    if isempty(session_info)
        continue
    end
    %     for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
    for n = 1
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
        load best_channels
        column = 1;

        load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')));
        load('extracted_V1_place_fields.mat')
        load('extracted_CA1_place_fields.mat')
        load('extracted_HPC_place_fields.mat')
        load(sprintf('decoded_ripple_events_V1%s.mat',erase(stimulus_name{n},'Masa2tracks')));
        load(sprintf('decoded_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks')));
        load('extracted_laps')
        %         load('extracted_V1_clusters_PRE_RUN.mat')
        load(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))

        if length(session_info(n).probe) > 1
            load(sprintf('extracted_HPC_clusters_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('extracted_V1_clusters_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load('extracted_V1_place_fields_combined.mat')
            load('extracted_HPC_place_fields_combined.mat')
            if contains(stimulus_name{n},'RUN')
                load('estimated_position_lap_CV_HPC_combined.mat')
                load('estimated_position_lap_CV_V1_combined.mat')
                estimated_position_lap_CV_HPC = estimated_position_lap_CV_HPC_combined;
                estimated_position_lap_CV_V1 = estimated_position_lap_CV_V1_combined;
            end
            HPC_cell_index = [];
            HPC_cell_hemisphere = [];

            % Get sorted cell id based on spatial activity peak
            place_fields = HPC_place_fields_combined;
            normalised_raw_matrix = [];
            for track_id = 1:length(HPC_place_fields_combined.track)
                max_FR = place_fields.track(track_id).raw_peak(place_fields.good_place_cells_LIBERAL)';
                raw_matrix = cat(1,place_fields.track(track_id).raw{place_fields.good_place_cells_LIBERAL});
                normalised_raw_matrix{track_id} = bsxfun(@rdivide, raw_matrix, max_FR);
                normalised_raw_matrix{track_id}(isnan(normalised_raw_matrix{track_id})) = 0;
                [~,peak_location] = max(normalised_raw_matrix{track_id},[],2);
                [~,HPC_cell_index(track_id,:)] = sort(peak_location);
                HPC_cell_index(track_id,:) = place_fields.good_place_cells_LIBERAL(HPC_cell_index(track_id,:));

                for ncell = 1:length(HPC_cell_index)
                    if HPC_clusters_combined.id_conversion((HPC_cell_index(track_id,ncell) == HPC_clusters_combined.id_conversion(:,1)),2) < 10000
                        HPC_cell_hemisphere(track_id,ncell) = 1;
                    else
                        HPC_cell_hemisphere(track_id,ncell) = 2;
                    end
                end
            end

            place_fields = V1_place_fields_combined;
            V1_cell_index = [];
            V1_cell_hemisphere = [];

            normalised_raw_matrix = [];
            for track_id = 1:length(HPC_place_fields_combined.track)
                max_FR = place_fields.track(track_id).raw_peak(place_fields.good_place_cells_LIBERAL)';
                raw_matrix = cat(1,place_fields.track(track_id).raw{place_fields.good_place_cells_LIBERAL});
                normalised_raw_matrix{track_id} = bsxfun(@rdivide, raw_matrix, max_FR);
                normalised_raw_matrix{track_id}(isnan(normalised_raw_matrix{track_id})) = 0;
                [~,peak_location] = max(normalised_raw_matrix{track_id},[],2);
                [~,V1_cell_index(track_id,:)] = sort(peak_location);
                V1_cell_index(track_id,:) = place_fields.good_place_cells_LIBERAL(V1_cell_index(track_id,:));

                for ncell = 1:length(V1_cell_index)
                    if V1_clusters_combined.id_conversion((V1_cell_index(track_id,ncell) == V1_clusters_combined.id_conversion(:,1)),2) < 10000
                        V1_cell_hemisphere(track_id,ncell) = 1;
                    else
                        V1_cell_hemisphere(track_id,ncell) = 2;
                    end
                end
            end

            V1_clusters = V1_clusters_combined;
            HPC_clusters = HPC_clusters_combined;

        else
            if contains(stimulus_name{n},'RUN')
                load('estimated_position_lap_CV_HPC.mat')
                load('estimated_position_lap_CV_V1.mat')
            end

            HPC_cell_index = HPC_place_fields.probe(1).good_place_cells_LIBERAL;
            V1_cell_index = V1_place_fields.probe(1).good_place_cells_LIBERAL;

            % Get sorted cell id based on spatial activity peak
            place_fields = HPC_place_fields.probe(1);
            normalised_raw_matrix = [];
            for track_id = 1:length(place_fields.track)
                max_FR = place_fields.track(track_id).raw_peak(place_fields.good_place_cells_LIBERAL)';
                raw_matrix = cat(1,place_fields.track(track_id).raw{place_fields.good_place_cells_LIBERAL});
                normalised_raw_matrix{track_id} = bsxfun(@rdivide, raw_matrix, max_FR);
                normalised_raw_matrix{track_id}(isnan(normalised_raw_matrix{track_id})) = 0;
                [~,peak_location] = max(normalised_raw_matrix{track_id},[],2);
                [~,HPC_cell_index(track_id,:)] = sort(peak_location);
                HPC_cell_index(track_id,:) = place_fields.good_place_cells_LIBERAL(HPC_cell_index(track_id,:));
            end


            place_fields = V1_place_fields.probe(1);
            normalised_raw_matrix = [];
            for track_id = 1:length(HPC_place_fields_combined.track)
                max_FR = place_fields.track(track_id).raw_peak(place_fields.good_place_cells_LIBERAL)';
                raw_matrix = cat(1,place_fields.track(track_id).raw{place_fields.good_place_cells_LIBERAL});
                normalised_raw_matrix{track_id} = bsxfun(@rdivide, raw_matrix, max_FR);
                normalised_raw_matrix{track_id}(isnan(normalised_raw_matrix{track_id})) = 0;
                [~,peak_location] = max(normalised_raw_matrix{track_id},[],2);
                [~,V1_cell_index(track_id,:)] = sort(peak_location);
                V1_cell_index(track_id,:) = place_fields.good_place_cells_LIBERAL(V1_cell_index(track_id,:));
            end

            V1_clusters = V1_clusters.probe(1);
            HPC_clusters = HPC_clusters.probe(1);
        end

        %         fig.Name = sprintf('relationship between V1 SUA spike count (Right hemisphere) and log odds during %s Session %i',behavioural_epoches_text{epoch});
        if exist('ripples') == 0
            mkdir('ripples')
        end

        nfigure = 0;
        for mprobe = 1:length(V1_place_fields.probe)
            ripple_probe_hemisphere = session_info(n).probe(mprobe).probe_hemisphere;

            % Find ripple events with 'coherent' reactivation when log odds
            % in V1 and HPC are coherently > 0.5  or < -0.5.

            if length(session_info(n).probe) > 1
                index = find(log_odds.experiment == nsession & log_odds.ripple_peak >= 5 & log_odds.event_probe_hemisphere == ripple_probe_hemisphere&...
                    ((log_odds.V1_log_odds_combined(:,12)'>0.5 & log_odds.zscore(:,12)'>0.5)|(log_odds.V1_log_odds_combined(:,12)'<-0.5 & log_odds.zscore(:,12)'<-0.5)));
            elseif ripple_probe_hemisphere == 1
                index = find(log_odds.experiment == nsession & log_odds.ripple_peak >= 5 & log_odds.event_probe_hemisphere == ripple_probe_hemisphere&...
                    ((log_odds.V1_log_odds_L(:,12)'>0.5 & log_odds.zscore(:,12)'>0.5)|(log_odds.V1_log_odds_L(:,12)'<-0.5 & log_odds.zscore(:,12)'<-0.5)));
            elseif ripple_probe_hemisphere == 2
                index = find(log_odds.experiment == nsession & log_odds.ripple_peak >= 5 & log_odds.event_probe_hemisphere == ripple_probe_hemisphere&...
                    ((log_odds.V1_log_odds_R(:,12)'>0.5 & log_odds.zscore(:,12)'>0.5)|(log_odds.V1_log_odds_R(:,12)'<-0.5 & log_odds.zscore(:,12)'<-0.5)));
            end

            if isempty(index)
                continue
            end

            for event = 1:length(index)

                fig = figure(nfigure)
                fig.Position = [300 150 810 800];
                fig.Name = sprintf('%s %s SWR response probe %s SWR (%i)',options.SUBJECT,options.SESSION,probe_hemisphere_text{probe_hemisphere},nfigure);

                onset = log_odds.onset(index(event));
                offset = onset + 1;

                for track_id = 1:length(HPC_place_fields.probe(1).track)
                    subplot(4,2,track_id)
                    for ncell = 1:length(HPC_cell_index)
                        cluster_id = HPC_clusters.id_conversion(HPC_cell_index(track_id,ncell) == HPC_clusters.id_conversion(:,1),2);
                        spike_times_this_cell = HPC_clusters.spike_times(HPC_clusters.spike_id == cluster_id);
                        spike_times_this_cell = spike_times_this_cell((spike_times_this_cell>onset-1)&(spike_times_this_cell<offset));
                        speed_during_spike = interp1(position.t,position.v_cm,spike_times_this_cell,'nearest');
                        spike_times_this_cell = spike_times_this_cell(speed_during_spike>5);



                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cell,...
                            [zeros(1,ncell-1) onset], [-1 offset-onset], 0.001);
                        hold on
                        if HPC_cell_hemisphere(track_id,ncell) == 1
                            plot(rasterX,rasterY,'b')
                        else
                            plot(rasterX,rasterY,'r')
                        end
                        xlim([0 offset-onset])
                        %                         [psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(spikeTimes, eventTimes, window, psthBinSize)
                    end
                    title(sprintf('HPC spiking Track %i',track_id))
                end


                for track_id = 1:length(V1_place_fields.probe(1).track)
                    subplot(4,2,track_id+2)
                    for ncell = 1:length(V1_cell_index)
                        cluster_id = V1_clusters.id_conversion(V1_cell_index(track_id,ncell) == V1_clusters.id_conversion(:,1),2);
                        spike_times_this_cell = V1_clusters.spike_times(V1_clusters.spike_id == cluster_id);
                        spike_times_this_cell = spike_times_this_cell((spike_times_this_cell>onset-2)&(spike_times_this_cell<offset));
                        speed_during_spike = interp1(position.t,position.v_cm,spike_times_this_cell,'nearest');
                        spike_times_this_cell = spike_times_this_cell(speed_during_spike>5);

                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cell,...
                            [zeros(1,ncell-1) onset], [-1 offset-onset], 0.001);
                        hold on
                        if V1_cell_hemisphere(track_id,ncell) == 1
                            plot(rasterX,rasterY,'b')
                        else
                            plot(rasterX,rasterY,'r')
                        end
                        xlim([0 offset-onset])
                        %                         [psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(spikeTimes, eventTimes, window, psthBinSize)
                    end
                    title(sprintf('V1 spiking Track %i',track_id))
                end

                subplot(4,2,5)
                decoded_ripple_events(ripple_probe_hemisphere).track(1).replay_events(event).replay

                estimated_position_lap_CV = estimated_position_lap_CV_HPC.track;
                if ~isempty(estimated_position_lap_CV(track).lap(lap_id))
                    imagesc([decoded_ripple_events(ripple_probe_hemisphere).track(1).replay_events(event).replay;...
                        decoded_ripple_events(ripple_probe_hemisphere).track(2).replay_events(event).replay])
                    hold on
                    plot(estimated_position_lap_CV(track).lap(lap_id).track(1).run_actual_position/10,'r')
                    plot(estimated_position_lap_CV(track).lap(lap_id).track(2).run_actual_position/10 + 15,'b')
                    yticks([30 50 70 90 110 140 170 190 210 230 250 280]/10)
                    yline(14.5,'LineWidth',2,'Color','k','DisplayName','Track 2')
                    yticklabels([30 50 70 90 110 140 30 50 70 90 110 140])
                    run_time_edges = estimated_position_lap_CV(track).lap(lap_id).track(1).run_time_edges;

                    xticks(linspace(1,length(run_time_edges),5))
                    xticklabels(linspace(run_time_edges(1),run_time_edges(end),5))
                end
                set(gca,"TickDir","out",'box', 'off','Color','none')
                %                 colorbar
                colormap(flip(gray))
                title('HPC')

                subplot(4,2,6)
                estimated_position_lap_CV = estimated_position_lap_CV_V1.track;
                if ~isempty(estimated_position_lap_CV(track).lap(lap_id))
                    imagesc([estimated_position_lap_CV(track).lap(lap_id).track(1).run; estimated_position_lap_CV(track).lap(lap_id).track(2).run])
                    hold on
                    plot(estimated_position_lap_CV(track).lap(lap_id).track(1).run_actual_position/10,'r')
                    plot(estimated_position_lap_CV(track).lap(lap_id).track(2).run_actual_position/10 + 15,'b')
                    yticks([30 50 70 90 110 140 170 190 210 230 250 280]/10)
                    yline(14.5,'LineWidth',2,'Color','k','DisplayName','Track 2')
                    yticklabels([30 50 70 90 110 140 30 50 70 90 110 140])
                    run_time_edges = estimated_position_lap_CV(track).lap(lap_id).track(1).run_time_edges;

                    xticks(linspace(1,length(run_time_edges),5))
                    xticklabels(linspace(run_time_edges(1),run_time_edges(end),5))
                end
                set(gca,"TickDir","out",'box', 'off','Color','none')
                %                 colorbar
                colormap(flip(gray))
                title('V1')

                subplot(4,2,7)
                time_this_event = position.t(position.t>onset-1 & position.t<offset);
                plot( time_this_event-time_this_event(1)-1,position.x(position.t>onset-1 & position.t<offset))
                xlim([-1 offset-onset])
                title(sprintf('Track %i lap %i',track,lap_id))



            end

            index1 = intersect(index,find(log_odds.zscore > 0));
            [~,idx] = sort(log_odds.zscore(index1));
            index1 = index1(idx);

            index2 = intersect(index,find(log_odds.zscore < 0));
            [~,idx] = sort(log_odds.zscore(index2));
            index2 = index2(idx);

            index0 = setdiff(index,[index1 index2]);
            [~,idx] = sort(log_odds.zscore(index0));
            index0 = index0(idx);

            sorted_index = [index0 index1 index2];
            for nprobe = 1:length(session_info(n).probe)
                options = session_info(n).probe(nprobe);
                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere ;

            end
        end



    end
end


%% Individual cell response to candidate events response

clear all
cd('Z:\ibn-vision\USERS\Masa\V1_HPC_project')
% save('log_odds_best_probe_events',"log_odds")
% save('log_odds_re_zscored',"log_odds")
% save('log_odds_50ms',"log_odds")
% load('log_odds_50ms_ratio_threshold_0')
load('log_odds_ripples_50ms_ratio_threshold_0')
load V1_all_cells_FR_bias
load HPC_all_cells_FR_bias
load V1_all_cells
load V1_spatial_cells
load HPC_all_cells
load HPC_spatial_cells

SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
% Stimulus_types_all = {'RUN'};
Stimulus_type = 'RUN'; % extract LFP during RUN
% Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
c = 1;
%
% time_windows = [1:10;11:20;21:30;...
%     31:40;41:50;[51:59 nan]]; % for V1 log odds
time_windows = [1:4;5:8;9:12;...
    13:16;17:20;[21:23 nan]]; % for V1 log odds

behavioural_epoches = [-1,0,1,2];
behavioural_epoches_text = {'PRE','POST','Track 1','Track 2'};
probe_hemisphere_text = {'left','right'};
cell_slope_L = [];
cell_slope_R = [];
index_session = [];

for nsession =1:10
    tic
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    if isempty(session_info)
        continue
    end
    %     for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
    for n = 1
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
        load best_channels
        column = 1;

        load('extracted_V1_place_fields.mat')
        load('extracted_CA1_place_fields.mat')
        load('extracted_HPC_place_fields.mat')
        %         load('extracted_V1_clusters_PRE_RUN.mat')
        load(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        %         fig.Name = sprintf('relationship between V1 SUA spike count (Right hemisphere) and log odds during %s Session %i',behavioural_epoches_text{epoch});


        index = find(log_odds.experiment == nsession & log_odds.ripple_peak >= 3);

        for epoch = 1:4
            index_session{nsession}{epoch} = [];
            index_session{nsession}{epoch} = intersect(index,find(log_odds.behavioural_state == behavioural_epoches(epoch)));
            cell_slope_L{nsession}{epoch} = [];
            cell_slope_R{nsession}{epoch} = [];
        end

        spike_data = [];


        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere ;
            for mprobe = 1:length(V1_place_fields.probe)
                ripple_probe_hemisphere = session_info(n).probe(mprobe).probe_hemisphere;

                spike_data{ripple_probe_hemisphere}{probe_hemisphere} =[];
                spike_data{ripple_probe_hemisphere}{probe_hemisphere} =[];
            end

            cell_index{probe_hemisphere} = [];
            % if left hemisphere find good cells with higher FR on T2
            if probe_hemisphere == 1
                cell_index{probe_hemisphere} = find((V1_place_fields.probe(probe_no).track(1).mean_rate_track - V1_place_fields.probe(probe_no).track(2).mean_rate_track)./...
                    (V1_place_fields.probe(probe_no).track(1).mean_rate_track + V1_place_fields.probe(probe_no).track(1).mean_rate_track) < 0);
                % if left hemisphere find good cells with higher FR on T1
            elseif probe_hemisphere == 2
                cell_index{probe_hemisphere} = find((V1_place_fields.probe(probe_no).track(1).mean_rate_track - V1_place_fields.probe(probe_no).track(2).mean_rate_track)./...
                    (V1_place_fields.probe(probe_no).track(1).mean_rate_track + V1_place_fields.probe(probe_no).track(1).mean_rate_track) > 0);
            end
        end



        nfigure = 0;
        for mprobe = 1:length(V1_place_fields.probe)
            ripple_probe_hemisphere = session_info(n).probe(mprobe).probe_hemisphere;

            % For ripple events, they all have peak ripple above 5
            index = find(log_odds.experiment == nsession & log_odds.ripple_peak >= 5 & log_odds.event_probe_hemisphere == ripple_probe_hemisphere);
            if isempty(index)
                continue
            end
            for event = 1:length(index)
                spike_data{ripple_probe_hemisphere}{1} = [spike_data{ripple_probe_hemisphere}{1}; log_odds.V1_spike_id_L{index(event)}];
                spike_data{ripple_probe_hemisphere}{2} = [spike_data{ripple_probe_hemisphere}{2}; log_odds.V1_spike_id_R{index(event)}];
            end

            index1 = intersect(index,find(log_odds.zscore > 0));
            [~,idx] = sort(log_odds.zscore(index1));
            index1 = index1(idx);

            index2 = intersect(index,find(log_odds.zscore < 0));
            [~,idx] = sort(log_odds.zscore(index2));
            index2 = index2(idx);

            index0 = setdiff(index,[index1 index2]);
            [~,idx] = sort(log_odds.zscore(index0));
            index0 = index0(idx);

            sorted_index = [index0 index1 index2];
            for nprobe = 1:length(session_info(n).probe)
                options = session_info(n).probe(nprobe);
                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere ;
                region = 'V1';
                pcount = 1;
                nfigure = nfigure+ 1;

                for ncell = 1:length(cell_index{probe_hemisphere})
                    this_cell = V1_all_cells{nsession}{probe_hemisphere}(cell_index{probe_hemisphere}(ncell));
                    spike_index = find(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(:,1)==this_cell);
                    if isempty(spike_index)
                        continue
                    end


                    if pcount >= 20
                        nfigure = nfigure + 1;
                        pcount = 1;
                    end

                    fig = figure(nfigure)
                    fig.Position = [300 150 810 800];
                    if ischar(nprobe)
                        fig.Name = sprintf('%s %s %s cell SWR response log odds %s (%i)',options.SUBJECT,options.SESSION,region,nprobe,nfigure);
                    else
                        fig.Name = sprintf('%s %s %s cell SWR response probes %s SWR %s (%i)',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere},probe_hemisphere_text{ripple_probe_hemisphere},nfigure);
                    end

                    %                     if ~isempty(index0)
                    %                         [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(spike_index,:),...
                    %                             log_odds.onset(index0), [-0.6 0.6], 0.001);
                    %
                    %                         subplot(6,4,pcount)
                    %                         plot(rasterX,rasterY,'k');hold on
                    %                         ylim([1 length(index)])
                    % %                         [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(spike_index,:),...
                    % %                             log_odds.onset(index0), [-0.6 0.6], 0.04);
                    %                         subplot(6,4,pcount+1)
                    %                         plot(mean(binnedArray)/0.04,'k');hold on
                    %                     end


                    if ~isempty(index1)
                        event_index = [0.001*ones(1,1+length(index0)) log_odds.onset(index1)];
                        subplot(5,4,pcount)
                        for event = 1:length(event_index)
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(spike_index,:),...
                                [zeros(1,event-1) event_index(event)], [-0.5 0.5], 0.001);

                            %                     [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{probe_hemisphere}(spike_index,:),...
                            %                         [log_odds.onset(index1)], [-0.6 0.6], 0.001);
                            %                     subplot(6,4,pcount)
                            plot(rasterX,rasterY,'r');hold on
                        end
                        %                     xline(0,'r')
                        xlim([-0.5 0.5])
                        ylim([1 length(index)])

                        subplot(5,4,pcount+1)
                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(spike_index,:),...
                            [log_odds.onset(index1)], [-0.5 0.5],0.05);
                        binnedArray(binnedArray>25) = nan;
                        plot(nanmean(binnedArray)/0.05,'r');hold on
                        ymax = max(binnedArray);

                    end

                    if ~isempty(index2)
                        %                     [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{probe_hemisphere}(spike_index,:),...
                        %                         [zeros(1,length([index0 index1])) log_odds.onset(index2)], [-0.6 0.6], 0.001);

                        event_index = [0.001*ones(1,1+length(index0)+length(index1)) log_odds.onset(index2)];
                        subplot(5,4,pcount)
                        for event = 1:length(event_index)
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(spike_index,:),...
                                [zeros(1,event-1) event_index(event)], [-0.5 0.5], 0.001);

                            %                     [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{probe_hemisphere}(spike_index,:),...
                            %                         [log_odds.onset(index1)], [-0.6 0.6], 0.001);
                            %                     subplot(6,4,pcount)
                            plot(rasterX,rasterY,'b');hold on
                        end
                        xlim([-0.5 0.5])
                        ylim([1 length(index)])
                        %                     xline(0,'r')
                        subplot(5,4,pcount+1)
                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(spike_index,:),...
                            [log_odds.onset(index2)], [-0.5 0.5], 0.05);
                        binnedArray(binnedArray>25) = nan;
                        plot(nanmean(binnedArray)/0.05,'b');hold on
                        ymax = [ymax max(binnedArray)];
                    end

                    subplot(5,4,pcount)
                    xlabel('Perievent Time')
                    ylabel('SWR Events')
                    xline(0,'r')

                    subplot(5,4,pcount+1)
                    xlabel('Perievent Time')
                    ylabel('Firing Rate')
                    xlim([1 20])
                    xticks([1 5 10.5 15 20])
                    xticklabels([-0.5 -0.25 0 0.25 0.5])
                    xline(10.5,'r')
                    %                     ylim([0 ymax])

                    pcount = pcount + 2;
                    this_cell_channel = V1_clusters.probe(probe_no).peak_channel((V1_clusters.probe(probe_no).id_conversion(:,2) == this_cell));
                    if ~isempty(best_channels{probe_no}.L4_channel)
                        distance = 10*(best_channels{probe_no}.L4_channel-this_cell_channel);
                    elseif ~isempty(best_channels{probe_no}.L5_channel) % rough estimation of
                        distance = 10*(best_channels{probe_no}.L5_channel+10-this_cell_channel);
                    else
                        distance = 10*(best_channels{probe_no}.first_in_brain_channel-45-this_cell_channel);
                    end

                    if sum(this_cell == V1_spatial_cells{nsession}{probe_hemisphere}{1}) == 1 &...
                            sum(this_cell == V1_spatial_cells{nsession}{probe_hemisphere}{2}) == 1

                        title(sprintf('Spatial cell %i (%i from L4)',this_cell,distance),'Color','m')
                    elseif sum(this_cell == V1_spatial_cells{nsession}{probe_hemisphere}{1}) == 1
                        title(sprintf('Spatial cell %i (%i from L4)',this_cell,distance),'Color','r')
                    elseif sum(this_cell == V1_spatial_cells{nsession}{probe_hemisphere}{2}) == 1
                        title(sprintf('Spatial cell %i (%i from L4)',this_cell,distance),'Color','b')
                    else
                        title(sprintf('cell %i (%i from L4)',this_cell,distance),'Color','k')
                    end

                    if ischar(nprobe)
                        sgtitle(sprintf('%s %s %s cell SWR response %s (%i)',options.SUBJECT,options.SESSION,region,nprobe,nfigure));
                    else
                        sgtitle(sprintf('%s %s %s cell SWR response probes %s SWR %s (%i)',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere},probe_hemisphere_text{ripple_probe_hemisphere},nfigure));
                    end
                end

                save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project\V1 cell response to reactivation SWR',[])
            end
        end



        % HPC spiking
        spike_data = [];
        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere ;
            for mprobe = 1:length(V1_place_fields.probe)
                ripple_probe_hemisphere = session_info(n).probe(mprobe).probe_hemisphere;

                spike_data{ripple_probe_hemisphere}{probe_hemisphere} =[];
                spike_data{ripple_probe_hemisphere}{probe_hemisphere} =[];
            end
        end



        nfigure = 0;
        for mprobe = 1:length(V1_place_fields.probe)
            ripple_probe_hemisphere = session_info(n).probe(mprobe).probe_hemisphere;

            % For ripple events, they all have peak ripple above 5
            index = find(log_odds.experiment == nsession & log_odds.ripple_peak >= 5 & log_odds.event_probe_hemisphere == ripple_probe_hemisphere);
            if isempty(index)
                continue
            end
            for event = 1:length(index)
                spike_data{ripple_probe_hemisphere}{1} = [spike_data{ripple_probe_hemisphere}{1}; log_odds.HPC_spike_id_L{index(event)}];
                spike_data{ripple_probe_hemisphere}{2} = [spike_data{ripple_probe_hemisphere}{2}; log_odds.HPC_spike_id_R{index(event)}];
            end

            index1 = intersect(index,find(log_odds.zscore > 0));
            [~,idx] = sort(log_odds.zscore(index1));
            index1 = index1(idx);

            index2 = intersect(index,find(log_odds.zscore < 0));
            [~,idx] = sort(log_odds.zscore(index2));
            index2 = index2(idx);

            index0 = setdiff(index,[index1 index2]);
            [~,idx] = sort(log_odds.zscore(index0));
            index0 = index0(idx);

            sorted_index = [index0 index1 index2];
            for nprobe = 1:length(session_info(n).probe)
                options = session_info(n).probe(nprobe);
                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere ;
                region = 'HPC';
                pcount = 1;
                nfigure = nfigure+ 1;

                for ncell = 1:length(HPC_all_cells{nsession}{probe_hemisphere})
                    this_cell = HPC_all_cells{nsession}{probe_hemisphere}(ncell);
                    spike_index = find(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(:,1)==this_cell);
                    if isempty(spike_index)
                        continue
                    end


                    if pcount >= 20
                        nfigure = nfigure + 1;
                        pcount = 1;
                    end

                    fig = figure(nfigure)
                    fig.Position = [300 150 810 800];
                    if ischar(nprobe)
                        fig.Name = sprintf('%s %s %s cell SWR response log odds %s (%i)',options.SUBJECT,options.SESSION,region,nprobe,nfigure);
                    else
                        fig.Name = sprintf('%s %s %s cell SWR response probes %s SWR %s (%i)',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere},probe_hemisphere_text{ripple_probe_hemisphere},nfigure);
                    end

                    %                     if ~isempty(index0)
                    %                         [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(spike_index,:),...
                    %                             log_odds.onset(index0), [-0.6 0.6], 0.001);
                    %
                    %                         subplot(6,4,pcount)
                    %                         plot(rasterX,rasterY,'k');hold on
                    %                         ylim([1 length(index)])
                    % %                         [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(spike_index,:),...
                    % %                             log_odds.onset(index0), [-0.6 0.6], 0.04);
                    %                         subplot(6,4,pcount+1)
                    %                         plot(mean(binnedArray)/0.04,'k');hold on
                    %                     end


                    if ~isempty(index1)
                        event_index = [0.001*ones(1,1+length(index0)) log_odds.onset(index1)];
                        subplot(5,4,pcount)
                        for event = 1:length(event_index)
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(spike_index,:),...
                                [zeros(1,event-1) event_index(event)], [-0.5 0.5], 0.001);

                            %                     [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{probe_hemisphere}(spike_index,:),...
                            %                         [log_odds.onset(index1)], [-0.6 0.6], 0.001);
                            %                     subplot(6,4,pcount)
                            plot(rasterX,rasterY,'r');hold on
                        end
                        %                     xline(0,'r')
                        xlim([-0.5 0.5])
                        ylim([1 length(index)])

                        subplot(5,4,pcount+1)
                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(spike_index,:),...
                            [log_odds.onset(index1)], [-0.5 0.5],0.05);
                        binnedArray(binnedArray>25) = nan;
                        plot(nanmean(binnedArray)/0.05,'r');hold on
                        ymax = max(binnedArray);

                    end

                    if ~isempty(index2)
                        %                     [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{probe_hemisphere}(spike_index,:),...
                        %                         [zeros(1,length([index0 index1])) log_odds.onset(index2)], [-0.6 0.6], 0.001);

                        event_index = [0.001*ones(1,1+length(index0)+length(index1)) log_odds.onset(index2)];
                        subplot(5,4,pcount)
                        for event = 1:length(event_index)
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(spike_index,:),...
                                [zeros(1,event-1) event_index(event)], [-0.5 0.5], 0.001);

                            %                     [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{probe_hemisphere}(spike_index,:),...
                            %                         [log_odds.onset(index1)], [-0.6 0.6], 0.001);
                            %                     subplot(6,4,pcount)
                            plot(rasterX,rasterY,'b');hold on
                        end
                        xlim([-0.5 0.5])
                        ylim([1 length(index)])
                        %                     xline(0,'r')
                        subplot(5,4,pcount+1)
                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(spike_index,:),...
                            [log_odds.onset(index2)], [-0.5 0.5], 0.05);
                        binnedArray(binnedArray>25) = nan;
                        plot(nanmean(binnedArray)/0.05,'b');hold on
                        ymax = [ymax max(binnedArray)];
                    end

                    subplot(5,4,pcount)
                    xlabel('Perievent Time')
                    ylabel('SWR Events')
                    xline(0,'r')

                    subplot(5,4,pcount+1)
                    xlabel('Perievent Time')
                    ylabel('Firing Rate')
                    xlim([1 20])
                    xticks([1 5 10.5 15 20])
                    xticklabels([-0.5 -0.25 0 0.25 0.5])
                    xline(10.5,'r')
                    %                     ylim([0 ymax])

                    pcount = pcount + 2;
                    this_cell_channel = HPC_clusters.probe(probe_no).peak_channel((HPC_clusters.probe(probe_no).id_conversion(:,2) == this_cell));
                    if ~isempty(best_channels{probe_no}.CA1_channel)
                        distance = 10*(best_channels{probe_no}.CA1_channel-this_cell_channel);
                    end

                    if sum(this_cell == HPC_spatial_cells{nsession}{probe_hemisphere}{1}) == 1 &...
                            sum(this_cell == HPC_spatial_cells{nsession}{probe_hemisphere}{2}) == 1

                        title(sprintf('Spatial cell %i (%i from CA1)',this_cell,distance),'Color','m')
                    elseif sum(this_cell == HPC_spatial_cells{nsession}{probe_hemisphere}{1}) == 1
                        title(sprintf('Spatial cell %i (%i from CA1)',this_cell,distance),'Color','r')
                    elseif sum(this_cell == HPC_spatial_cells{nsession}{probe_hemisphere}{2}) == 1
                        title(sprintf('Spatial cell %i (%i from CA1)',this_cell,distance),'Color','b')
                    else
                        title(sprintf('cell %i (%i from CA1)',this_cell,distance),'Color','k')
                    end

                    if ischar(nprobe)
                        sgtitle(sprintf('%s %s %s cell SWR response %s (%i)',options.SUBJECT,options.SESSION,region,nprobe,nfigure));
                    else
                        sgtitle(sprintf('%s %s %s cell SWR response probes %s SWR %s (%i)',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere},probe_hemisphere_text{ripple_probe_hemisphere},nfigure));
                    end
                end

                save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project\V1 cell response to reactivation SWR',[])
            end
        end


    end
end

PSTH = plot_perievent_spiketime_histogram(all_spike_data{probe_hemisphere(nprobe)},replay.onset(event),'group','by cell zscore',...
    'group_name',group_name,'twin',[-1 1],'bin_size',0.02,'plot_option',0,'smooth_option',1);

%% V1 spiking vs HPC log odds
% time_windows = [31:35;36:40;41:45;46:50;...
%     51:55;56:60;61:65;66:70];
% time_windows_text = {'-300 to -400ms','-200 to -300ms','-100 to -200ms','0 to -100ms','0 to 100ms','100 to 200ms','200 to 300ms','300ms to 400ms'};
behavioural_epoches = [-1,0,1,2];
behavioural_epoches_text = {'PRE','POST','Track 1','Track 2'};
time_windows = [31:40;41:50;...
    51:60;61:70];
time_windows_text = {'-200 to -400ms','0 to -200ms','0 to 200ms','200 to 400ms'};


for epoch = 1:length(behavioural_epoches)
    c = 1;
    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between normalised V1 spike count (Right hemisphere) and HPC log odds during %s',behavioural_epoches_text{epoch});

    for nwin =1:size(time_windows,1)
        time_bin= time_windows(nwin,:);

        subplot(2,3,nwin)
        hold on
        %                 index = find( log_odds.behavioural_state == behavioural_epoches(epoch));
        index = find( log_odds.behavioural_state >0);
        %         index = intersect(index,find(log_odds.experiment >= 5))
        index = intersect(index,find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 8 ))
        %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
        index = intersect(index,find(log_odds.ripple_peak >=3));
        %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));


        scatter(log_odds.zscore_substracted(index),trapz(log_odds.V1_spike_count_R(index,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.zscore_substracted(index)',trapz(log_odds.V1_spike_count_R(index,time_bin)'));
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.zscore_substracted') max(log_odds.zscore_substracted')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.zscore_substracted < -1));
        scatter(log_odds.zscore_substracted(index2),trapz(log_odds.V1_spike_count_R(index2,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.zscore_substracted > 1));
        scatter(log_odds.zscore_substracted(index1),trapz(log_odds.V1_spike_count_R(index1,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        title(time_windows_text{nwin})

        set(gca,"TickDir","out",'box', 'off','Color','none')

    end

    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between normalised V1 spike count (Left hemisphere) and HPC log odds during %s',behavioural_epoches_text{epoch});
    for nwin =1:size(time_windows,1)
        time_bin= time_windows(nwin,:);

        subplot(2,2,nwin)
        hold on

        scatter(log_odds.zscore_substracted(index),trapz(log_odds.V1_spike_count_L(index,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.zscore_substracted(index)',trapz(log_odds.V1_spike_count_L(index,time_bin)'));
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.zscore_substracted') max(log_odds.zscore_substracted')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.zscore_substracted < -1));
        scatter(log_odds.zscore_substracted(index2),trapz(log_odds.V1_spike_count_L(index2,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.zscore_substracted > 1));
        scatter(log_odds.zscore_substracted(index1),trapz(log_odds.V1_spike_count_L(index1,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        set(gca,"TickDir","out",'box', 'off','Color','none')
        title(time_windows_text{nwin})
    end

end

%% V1 spiking vs HPC bias
% time_windows = [31:35;36:40;41:45;46:50;...
%     51:55;56:60;61:65;66:70];
% time_windows_text = {'-300 to -400ms','-200 to -300ms','-100 to -200ms','0 to -100ms','0 to 100ms','100 to 200ms','200 to 300ms','300ms to 400ms'};
behavioural_epoches = [-1,0,1,2];
behavioural_epoches_text = {'PRE','POST','Track 1','Track 2'};
time_windows = [31:40;41:50;...
    51:60;61:70];
time_windows_text = {'-200 to -400ms','0 to -200ms','0 to 200ms','200 to 400ms'};


for epoch = 2:length(behavioural_epoches)
    c = 1;
    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between normalised V1 spike count (Right hemisphere) and HPC log odds during %s',behavioural_epoches_text{epoch});

    for nwin =1:size(time_windows,1)
        time_bin= time_windows(nwin,:);

        subplot(2,2,nwin)
        hold on
        index = find( log_odds.behavioural_state == behavioural_epoches(epoch));
        %         index = find( log_odds.behavioural_state > 0);
        %         index = intersect(index,find(log_odds.experiment >= 5))
        index = intersect(index,find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 8  ));
        %         index = intersect(index,find(log_odds.experiment == 8));
        %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
        index = intersect(index,find(log_odds.ripple_peak >=3));
        %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));


        scatter(log_odds.HPC_bayesian_bias(index),mean(log_odds.V1_spike_count_R(index,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.HPC_bayesian_bias(index)',mean(log_odds.V1_spike_count_R(index,time_bin),2)');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.HPC_bayesian_bias') max(log_odds.HPC_bayesian_bias')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.4));
        scatter(log_odds.HPC_bayesian_bias(index2),mean(log_odds.V1_spike_count_R(index2,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.6));
        scatter(log_odds.HPC_bayesian_bias(index1),mean(log_odds.V1_spike_count_R(index1,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        title(time_windows_text{nwin})

        set(gca,"TickDir","out",'box', 'off','Color','none')

    end

    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between normalised V1 spike count (Left hemisphere) and HPC log odds during %s',behavioural_epoches_text{epoch});
    for nwin =1:size(time_windows,1)
        time_bin= time_windows(nwin,:);

        subplot(2,2,nwin)
        hold on


        scatter(log_odds.HPC_bayesian_bias(index),mean(log_odds.V1_spike_count_L(index,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.HPC_bayesian_bias(index)',mean(log_odds.V1_spike_count_L(index,time_bin)'));
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.HPC_bayesian_bias') max(log_odds.HPC_bayesian_bias')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.4));
        scatter(log_odds.HPC_bayesian_bias(index2),mean(log_odds.V1_spike_count_L(index2,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.6));
        scatter(log_odds.HPC_bayesian_bias(index1),mean(log_odds.V1_spike_count_L(index1,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        set(gca,"TickDir","out",'box', 'off','Color','none')
        title(time_windows_text{nwin})
    end

end

%% V1 total spiking vs HPC bias
% time_windows = [31:35;36:40;41:45;46:50;...
%     51:55;56:60;61:65;66:70];
% time_windows_text = {'-300 to -400ms','-200 to -300ms','-100 to -200ms','0 to -100ms','0 to 100ms','100 to 200ms','200 to 300ms','300ms to 400ms'};
behavioural_epoches = [-1,0,1,2];
behavioural_epoches_text = {'PRE','POST','Track 1','Track 2'};
time_windows = [31:40;41:50;...
    51:60;61:70];
time_windows_text = {'-200 to -400ms','0 to -200ms','0 to 200ms','200 to 400ms'};


for epoch = 2:length(behavioural_epoches)
    c = 1;
    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between normalised V1 spike count (Right hemisphere) and HPC log odds during %s',behavioural_epoches_text{epoch});

    for nwin =1:size(time_windows,1)
        time_bin= time_windows(nwin,:);

        subplot(2,2,nwin)
        hold on
        index = find( log_odds.behavioural_state == behavioural_epoches(epoch));
        %         index = find( log_odds.behavioural_state > 0);
        %         index = intersect(index,find(log_odds.experiment >= 5))
        %         index = intersect(index,find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 7  ));
        index = intersect(index,find(log_odds.experiment >= 8));
        %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
        index = intersect(index,find(log_odds.ripple_peak >=3));
        %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));


        scatter(log_odds.HPC_bayesian_bias(index),sum(log_odds.V1_total_spike_count_R(index,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.HPC_bayesian_bias(index)',sum(log_odds.V1_total_spike_count_R(index,time_bin),2)');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.HPC_bayesian_bias') max(log_odds.HPC_bayesian_bias')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.4));
        scatter(log_odds.HPC_bayesian_bias(index2),sum(log_odds.V1_total_spike_count_R(index2,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.6));
        scatter(log_odds.HPC_bayesian_bias(index1),sum(log_odds.V1_total_spike_count_R(index1,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        title(time_windows_text{nwin})

        set(gca,"TickDir","out",'box', 'off','Color','none')

    end

    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between normalised V1 spike count (Left hemisphere) and HPC log odds during %s',behavioural_epoches_text{epoch});
    for nwin =1:size(time_windows,1)
        time_bin= time_windows(nwin,:);

        subplot(2,2,nwin)
        hold on


        scatter(log_odds.HPC_bayesian_bias(index),sum(log_odds.V1_total_spike_count_L(index,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.HPC_bayesian_bias(index)',sum(log_odds.V1_total_spike_count_L(index,time_bin)'));
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.HPC_bayesian_bias') max(log_odds.HPC_bayesian_bias')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.4));
        scatter(log_odds.HPC_bayesian_bias(index2),sum(log_odds.V1_total_spike_count_L(index2,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.6));
        scatter(log_odds.HPC_bayesian_bias(index1),sum(log_odds.V1_total_spike_count_L(index1,time_bin),2),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        set(gca,"TickDir","out",'box', 'off','Color','none')
        title(time_windows_text{nwin})
    end

end

%% V1 log odds vs HPC log odds
behavioural_epoches = [-1,0,1,2];
behavioural_epoches_text = {'PRE','POST','Track 1','Track 2'};
% time_windows = [31:40;41:50;...
%     51:60;61:70];
time_windows_text = {'-600 to -400 ms','-400 to -200ms','-200 to 0ms','0 to 200ms','200 to 400ms','400 to 600ms'};


for epoch = 2:length(behavioural_epoches)
    c = 1;
    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between V1 log odds (Right hemisphere) and HPC log odds during %s',behavioural_epoches_text{epoch});

    for nbin =1:6

        subplot(2,3,nbin)
        hold on
        index = find( log_odds.behavioural_state == behavioural_epoches(epoch));
        %         index = find( log_odds.behavioural_state>0);
        index = intersect(index,find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 7  ));
        %         index = intersect(index,find(log_odds.experiment  >= 8 ));
        %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
        index = intersect(index,find(log_odds.ripple_peak >=3));
        %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));


        scatter(log_odds.zscore(index),log_odds.V1_log_odds_R(index,nbin)','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.zscore(index)',log_odds.V1_log_odds_R(index,nbin)');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.zscore') max(log_odds.zscore')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.zscore < -1));
        scatter(log_odds.zscore(index2),log_odds.V1_log_odds_R(index2,nbin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.zscore > 1));
        scatter(log_odds.zscore(index1),log_odds.V1_log_odds_R(index1,nbin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        title(time_windows_text{nbin})

        set(gca,"TickDir","out",'box', 'off','Color','none')

    end

    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between V1 log odds (Left hemisphere) and HPC log odds during %s',behavioural_epoches_text{epoch});

    for nbin =1:6

        subplot(2,3,nbin)
        hold on


        scatter(log_odds.zscore(index),log_odds.V1_log_odds_L(index,nbin)','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.zscore(index)',log_odds.V1_log_odds_L(index,nbin)');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.zscore') max(log_odds.zscore')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.zscore < -1));
        scatter(log_odds.zscore(index2),log_odds.V1_log_odds_L(index2,nbin)','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.zscore > 1));
        scatter(log_odds.zscore(index1),log_odds.V1_log_odds_L(index1,nbin)','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        set(gca,"TickDir","out",'box', 'off','Color','none')
        title(time_windows_text{nbin})
    end

end

%% dV1 log odds vs HPC log odds
behavioural_epoches = [-1,0,1,2];
behavioural_epoches_text = {'PRE','POST','Track 1','Track 2'};
% time_windows = [31:40;41:50;...
%     51:60;61:70];
time_windows_text = {'-600 to -400 ms','-400 to -200ms','-200 to 0ms','0 to 200ms','200 to 400ms','400 to 600ms'};


for epoch = 2:length(behavioural_epoches)
    c = 1;
    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between V1 log odds (Right hemisphere) and HPC log odds during %s',behavioural_epoches_text{epoch});

    for nbin =1:6

        subplot(2,3,nbin)
        hold on
        %         index = find( log_odds.behavioural_state == behavioural_epoches(epoch));
        index = find( log_odds.behavioural_state == behavioural_epoches(epoch) & log_odds.event_probe_hemisphere == 2);
        %         index = find( log_odds.behavioural_state>0);
        index = intersect(index,find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 8  ));
        %                 index = intersect(index,find(log_odds.experiment  == 8 ));
        %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
        index = intersect(index,find(log_odds.ripple_peak >=3));
        %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));


        scatter(log_odds.zscore(index),log_odds.dV1_log_odds_R(index,nbin)','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.zscore(index)',log_odds.dV1_log_odds_R(index,nbin)');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.zscore') max(log_odds.zscore')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.zscore < -1));
        scatter(log_odds.zscore(index2),log_odds.dV1_log_odds_R(index2,nbin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.zscore > 1));
        scatter(log_odds.zscore(index1),log_odds.dV1_log_odds_R(index1,nbin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        title(time_windows_text{nbin})

        set(gca,"TickDir","out",'box', 'off','Color','none')

    end


    index = find( log_odds.behavioural_state == behavioural_epoches(epoch) & log_odds.event_probe_hemisphere == 2);
    %         index = find( log_odds.behavioural_state>0);
    index = intersect(index,find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 8  ));
    %             index = intersect(index,find(log_odds.experiment  ==8 ));
    %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
    index = intersect(index,find(log_odds.ripple_peak >=3));
    %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));

    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between V1 log odds (Left hemisphere) and HPC log odds during %s',behavioural_epoches_text{epoch});

    for nbin =1:6

        subplot(2,3,nbin)
        hold on


        scatter(log_odds.zscore(index),log_odds.dV1_log_odds_L(index,nbin)','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.zscore(index)',log_odds.dV1_log_odds_L(index,nbin)');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.zscore') max(log_odds.zscore')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.zscore < -1));
        scatter(log_odds.zscore(index2),log_odds.dV1_log_odds_L(index2,nbin)','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.zscore > 1));
        scatter(log_odds.zscore(index1),log_odds.dV1_log_odds_L(index1,nbin)','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        set(gca,"TickDir","out",'box', 'off','Color','none')
        title(time_windows_text{nbin})
    end

end

%% V1 vs HPC bayesian bias
behavioural_epoches = [-1,0,1,2];
behavioural_epoches_text = {'PRE','POST','Track 1','Track 2'};
% time_windows = [31:40;41:50;...
%     51:60;61:70];
time_windows_text = {'-600 to -400 ms','-400 to -200ms','-200 to 0ms','0 to 200ms','200 to 400ms','400 to 600ms'};
for epoch = 3:4
    c = 1;
    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between V1 bayesian bias (Right hemisphere) and HPC bayesian bias during %s',behavioural_epoches_text{epoch});

    for nbin =1:6

        subplot(2,3,nbin)
        hold on
        index = find( log_odds.behavioural_state == behavioural_epoches(epoch));
        %         index = find( log_odds.behavioural_state > 0);
        index = intersect(index,find(log_odds.experiment <= 0|log_odds.experiment == 5  |log_odds.experiment >= 8  ))
        %         index = intersect(index,find(log_odds.experiment => 8 ));


        %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
        index = intersect(index,find(log_odds.ripple_peak >=3));
        %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));

        index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.5));
        %         scatter3(log_odds.V1_bayesian_bias_R(index2,nbin)',log_odds.V1_bayesian_bias_L(index2,nbin)',log_odds.HPC_bayesian_bias(index2),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        scatter3(log_odds.V1_bayesian_bias_R(index2,nbin)',log_odds.HPC_bayesian_bias(index2)',log_odds.V1_bayesian_bias_L(index2,nbin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')

        hold on
        index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.5));
        scatter3(log_odds.V1_bayesian_bias_R(index1,nbin)',log_odds.HPC_bayesian_bias(index1)',log_odds.V1_bayesian_bias_L(index1,nbin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %         scatter3(log_odds.V1_bayesian_bias_R(index1,nbin)',log_odds.V1_bayesian_bias_L(index1,nbin)',log_odds.HPC_bayesian_bias(index1),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')
        xlabel('V1 R')
        zlabel('V1 L')
        ylabel('HPC')
    end
end


% bayesian bias V1 vs HPC
behavioural_epoches = [-1,0,1,2];
behavioural_epoches_text = {'PRE','POST','Track 1','Track 2'};
% time_windows = [31:40;41:50;...
%     51:60;61:70];
time_windows_text = {'-600 to -400 ms','-400 to -200ms','-200 to 0ms','0 to 200ms','200 to 400ms','400 to 600ms'};
for epoch = 2:4
    c = 1;
    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between V1 bayesian bias (Right hemisphere) and HPC bayesian bias during %s',behavioural_epoches_text{epoch});

    index = find( log_odds.behavioural_state == behavioural_epoches(epoch) & log_odds.event_probe_hemisphere >= 1);
    %         index = find( log_odds.behavioural_state>0);
    index = intersect(index,find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 8 ));
    %             index = intersect(index,find(log_odds.experiment == 9 ));
    %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
    %     index = intersect(index,find(log_odds.ripple_peak >=3));
    %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));

    for nbin =1:6

        subplot(2,3,nbin)
        hold on


        scatter(log_odds.HPC_bayesian_bias(index),log_odds.V1_bayesian_bias_R(index,nbin)','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.V1_bayesian_bias_R(index,nbin)');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.HPC_bayesian_bias(index)') max(log_odds.HPC_bayesian_bias(index)')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.5));
        scatter(log_odds.HPC_bayesian_bias(index2),log_odds.V1_bayesian_bias_R(index2,nbin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.5));
        scatter(log_odds.HPC_bayesian_bias(index1),log_odds.V1_bayesian_bias_R(index1,nbin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')


        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        title(time_windows_text{nbin})

        set(gca,"TickDir","out",'box', 'off','Color','none')

    end


    index = find( log_odds.behavioural_state == behavioural_epoches(epoch)& log_odds.event_probe_hemisphere >= 1);
    %             index = find( log_odds.behavioural_state==2);
    index = intersect(index,find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 8 ));
    %             index = intersect(index,find(log_odds.experiment == 9 ));
    %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
    %     index = intersect(index,find(log_odds.ripple_peak >=3));
    %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));

    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between V1 bayesian bias (Left hemisphere) and HPC bayesian bias during %s',behavioural_epoches_text{epoch});

    for nbin =1:6

        subplot(2,3,nbin)
        hold on

        scatter(log_odds.HPC_bayesian_bias(index),log_odds.V1_bayesian_bias_L(index,nbin)','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.V1_bayesian_bias_L(index,nbin)');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.HPC_bayesian_bias(index)') max(log_odds.HPC_bayesian_bias(index)')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.5));
        scatter(log_odds.HPC_bayesian_bias(index2),log_odds.V1_bayesian_bias_L(index2,nbin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.5));
        scatter(log_odds.HPC_bayesian_bias(index1),log_odds.V1_bayesian_bias_L(index1,nbin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        title(time_windows_text{nbin})

        set(gca,"TickDir","out",'box', 'off','Color','none')
    end
end


% dV1 vs HPC bayesian bias
behavioural_epoches = [-1,0,1,2];
behavioural_epoches_text = {'PRE','POST','Track 1','Track 2'};
% time_windows = [31:40;41:50;...
%     51:60;61:70];
time_windows_text = {'-600 to -400 ms','-400 to -200ms','-200 to 0ms','0 to 200ms','200 to 400ms','400 to 600ms'};
for epoch = 2:4
    c = 1;
    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between dV1 bayesian bias (Right hemisphere) and HPC bayesian bias during %s',behavioural_epoches_text{epoch});


    index = find( log_odds.behavioural_state == behavioural_epoches(epoch)& log_odds.event_probe_hemisphere >= 2);
    %         index = find( log_odds.behavioural_state>0);
    index = intersect(index,find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 8 ));
    %             index = intersect(index,find(log_odds.experiment >= 8 ));
    %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
    index = intersect(index,find(log_odds.ripple_peak >=3));
    %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));

    for nbin =1:6

        subplot(2,3,nbin)
        hold on


        scatter(log_odds.HPC_bayesian_bias(index),log_odds.dV1_bayesian_bias_R(index,nbin)','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.dV1_bayesian_bias_R(index,nbin)');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.HPC_bayesian_bias(index)') max(log_odds.HPC_bayesian_bias(index)')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.5));
        scatter(log_odds.HPC_bayesian_bias(index2),log_odds.dV1_bayesian_bias_R(index2,nbin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.5));
        scatter(log_odds.HPC_bayesian_bias(index1),log_odds.dV1_bayesian_bias_R(index1,nbin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')


        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        title(time_windows_text{nbin})

        set(gca,"TickDir","out",'box', 'off','Color','none')

    end


    index = find( log_odds.behavioural_state == behavioural_epoches(epoch)& log_odds.event_probe_hemisphere >= 2);
    %             index = find( log_odds.behavioural_state>0& log_odds.event_probe_hemisphere >= 2);
    index = intersect(index,find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 8 ));
    %         index = intersect(index,find(log_odds.experiment >= 8 ));
    %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
    index = intersect(index,find(log_odds.ripple_peak >=3));
    %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));


    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between dV1 bayesian bias (Left hemisphere) and HPC bayesian bias during %s',behavioural_epoches_text{epoch});

    for nbin =1:6

        subplot(2,3,nbin)
        hold on

        scatter(log_odds.HPC_bayesian_bias(index),log_odds.dV1_bayesian_bias_L(index,nbin)','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.dV1_bayesian_bias_L(index,nbin)');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.HPC_bayesian_bias(index)') max(log_odds.HPC_bayesian_bias(index)')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.5));
        scatter(log_odds.HPC_bayesian_bias(index2),log_odds.dV1_bayesian_bias_L(index2,nbin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.5));
        scatter(log_odds.HPC_bayesian_bias(index1),log_odds.dV1_bayesian_bias_L(index1,nbin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        title(time_windows_text{nbin})

        set(gca,"TickDir","out",'box', 'off','Color','none')
    end
end



% sV1 vs HPC bayesian bias
behavioural_epoches = [-1,0,1,2];
behavioural_epoches_text = {'PRE','POST','Track 1','Track 2'};
% time_windows = [31:40;41:50;...
%     51:60;61:70];
time_windows_text = {'-600 to -400 ms','-400 to -200ms','-200 to 0ms','0 to 200ms','200 to 400ms','400 to 600ms'};
for epoch = 2:4
    c = 1;
    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between sV1 bayesian bias (Right hemisphere) and HPC bayesian bias during %s',behavioural_epoches_text{epoch});


    index = find( log_odds.behavioural_state == behavioural_epoches(epoch)& log_odds.event_probe_hemisphere >= 2);
    %         index = find( log_odds.behavioural_state>0);
    index = intersect(index,find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 8 ));
    %             index = intersect(index,find(log_odds.experiment >= 8 ));
    %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
    index = intersect(index,find(log_odds.ripple_peak >=3));
    %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));

    for nbin =1:6

        subplot(2,3,nbin)
        hold on


        scatter(log_odds.HPC_bayesian_bias(index),log_odds.sV1_bayesian_bias_R(index,nbin)','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.sV1_bayesian_bias_R(index,nbin)');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.HPC_bayesian_bias(index)') max(log_odds.HPC_bayesian_bias(index)')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.5));
        scatter(log_odds.HPC_bayesian_bias(index2),log_odds.sV1_bayesian_bias_R(index2,nbin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.5));
        scatter(log_odds.HPC_bayesian_bias(index1),log_odds.sV1_bayesian_bias_R(index1,nbin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')


        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        title(time_windows_text{nbin})

        set(gca,"TickDir","out",'box', 'off','Color','none')

    end


    index = find( log_odds.behavioural_state == behavioural_epoches(epoch)& log_odds.event_probe_hemisphere >= 2);
    %         index = find( log_odds.behavioural_state>0);
    index = intersect(index,find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 8 ));
    %         index = intersect(index,find(log_odds.experiment >= 8 ));
    %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
    index = intersect(index,find(log_odds.ripple_peak >=3));
    %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));


    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between sV1 bayesian bias (Left hemisphere) and HPC bayesian bias during %s',behavioural_epoches_text{epoch});

    for nbin =1:6

        subplot(2,3,nbin)
        hold on

        scatter(log_odds.HPC_bayesian_bias(index),log_odds.sV1_bayesian_bias_L(index,nbin)','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.sV1_bayesian_bias_L(index,nbin)');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.HPC_bayesian_bias(index)') max(log_odds.HPC_bayesian_bias(index)')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.5));
        scatter(log_odds.HPC_bayesian_bias(index2),log_odds.sV1_bayesian_bias_L(index2,nbin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.5));
        scatter(log_odds.HPC_bayesian_bias(index1),log_odds.sV1_bayesian_bias_L(index1,nbin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        title(time_windows_text{nbin})

        set(gca,"TickDir","out",'box', 'off','Color','none')
    end
end
save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project\Log odds analysis',[])


%% On track awake reactivation correlation
% On track awake replay

index_resampled = [];
index11 = [];
index22 = [];
for nsession = [1 2 5 7 8 9 10]

    index = find(log_odds.experiment ==nsession );
    index = intersect(index,find(log_odds.ripple_peak >=3));

    s = RandStream('mrg32k3a','Seed',nsession); % Set random seed for resampling
    temp = intersect(index,find(log_odds.behavioural_state == 1));
    temp = datasample(s,temp,length(temp),'Replace',false);
    index_resampled = [index_resampled temp];
    index11 = [index11 temp];

    s = RandStream('mrg32k3a','Seed',nsession+100); % Set random seed for resampling
    temp = intersect(index,find(log_odds.behavioural_state == 2));
    temp = datasample(s,temp,length(temp),'Replace',false);
    index_resampled = [index_resampled temp];
    index22 = [index22 temp];
end
index_resampled = [index11 index22];

% index = find(log_odds.experiment ==8 );
index = find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 7 );
index = intersect(index,find(log_odds.ripple_peak >=3));
index = intersect(index,find(log_odds.behavioural_state > 0));
index1 = intersect(index,find( log_odds.behavioural_state == 1));
index2 = intersect(index,find( log_odds.behavioural_state == 2));
index = [index1 index2];

timebin = 6;

fig = figure
fig.Position = [750 240 950 640];
fig.Name = 'On track awake reactivation HPC vs V1 bayesian bias';
subplot(2,2,1)
scatter(log_odds.HPC_bayesian_bias(index1),log_odds.V1_bayesian_bias_R(index1,timebin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.2')
hold on
scatter(log_odds.HPC_bayesian_bias(index2),log_odds.V1_bayesian_bias_R(index2,timebin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.2')
xlabel('HPC bayesian bias')
ylabel('V1 Right bayesian bias')
mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.V1_bayesian_bias_R(index,timebin)');
[pval,~,~] = coefTest(mdl);
x =[min(log_odds.HPC_bayesian_bias(index)') max(log_odds.HPC_bayesian_bias(index)')];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);

if pval <= 0.05
    plot(x,y_est,'r:')
    %         xlim([-30 120])

    %     title(sprintf('Session %i',s),'Color','red')
    %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
    text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
else
    plot(x,y_est,'k:')
    %         xlim([-30 120])
    %     title(sprintf('Session %i',s))
    %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
    text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
end
legend('Track 1','Track 2','Color','none')
set(gca,"TickDir","out",'box', 'off','Color','none')
title('HPC vs V1 right bayesian bias')

if sum(isnan(log_odds.V1_bayesian_bias_L(index,1))) ~= length(index)
    subplot(2,2,2)
    scatter(log_odds.HPC_bayesian_bias(index1),log_odds.V1_bayesian_bias_L(index1,timebin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.2')
    hold on
    scatter(log_odds.HPC_bayesian_bias(index2),log_odds.V1_bayesian_bias_L(index2,timebin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.2')
    xlabel('HPC bayesian bias')
    ylabel('V1 Left bayesian bias')
    mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.V1_bayesian_bias_L(index,timebin)');
    [pval,~,~] = coefTest(mdl);
    x =[min(log_odds.HPC_bayesian_bias(index)') max(log_odds.HPC_bayesian_bias(index)')];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);

    if pval <= 0.05
        plot(x,y_est,'r:')
        %         xlim([-30 120])

        %     title(sprintf('Session %i',s),'Color','red')
        %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
        text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
    else
        plot(x,y_est,'k:')
        %         xlim([-30 120])
        %     title(sprintf('Session %i',s))
        %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
        text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
    end
    legend('Track 1','Track 2','Color','none')
end
set(gca,"TickDir","out",'box', 'off','Color','none')
title('HPC vs V1 left bayesian bias')

subplot(2,2,3)
scatter(log_odds.HPC_bayesian_bias(index1),log_odds.V1_bayesian_bias_R(index11,timebin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.2')
hold on
scatter(log_odds.HPC_bayesian_bias(index2),log_odds.V1_bayesian_bias_R(index22,timebin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.2')
xlabel('HPC bayesian bias')
ylabel('V1 Right bayesian bias')
mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.V1_bayesian_bias_R(index_resampled,timebin)');
[pval,~,~] = coefTest(mdl);
x =[min(log_odds.HPC_bayesian_bias(index)') max(log_odds.HPC_bayesian_bias(index)')];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);

if pval <= 0.05
    plot(x,y_est,'r:')
    %         xlim([-30 120])

    %     title(sprintf('Session %i',s),'Color','red')
    %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
    text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
else
    plot(x,y_est,'k:')
    %         xlim([-30 120])
    %     title(sprintf('Session %i',s))
    %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
    text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
end
legend('Track 1','Track 2','Color','none')
set(gca,"TickDir","out",'box', 'off','Color','none')
title('HPC vs V1 right bayesian bias (within track resampled)')

if sum(isnan(log_odds.V1_bayesian_bias_L(index,1))) ~= length(index)
    subplot(2,2,4)
    scatter(log_odds.HPC_bayesian_bias(index1),log_odds.V1_bayesian_bias_L(index11,timebin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.2')
    hold on
    scatter(log_odds.HPC_bayesian_bias(index2),log_odds.V1_bayesian_bias_L(index22,timebin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.2')
    xlabel('HPC bayesian bias')
    ylabel('V1 Left bayesian bias')
    mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.V1_bayesian_bias_L(index_resampled,timebin)');
    [pval,~,~] = coefTest(mdl);
    x =[min(log_odds.HPC_bayesian_bias(index)') max(log_odds.HPC_bayesian_bias(index)')];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);

    if pval <= 0.05
        plot(x,y_est,'r:')
        %         xlim([-30 120])

        %     title(sprintf('Session %i',s),'Color','red')
        %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
        text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
    else
        plot(x,y_est,'k:')
        %         xlim([-30 120])
        %     title(sprintf('Session %i',s))
        %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
        text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
    end
    legend('Track 1','Track 2','Color','none')
end
set(gca,"TickDir","out",'box', 'off','Color','none')
title('HPC vs V1 left bayesian bias (within track resampled)')
save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project\Log odds analysis',[])

%% log odds


index_resampled = [];
index11 = [];
index22 = [];
for nsession = [1 2 5 7 8 9 10]

    index = find(log_odds.experiment ==nsession );
    index = intersect(index,find(log_odds.ripple_peak >=3));

    s = RandStream('mrg32k3a','Seed',nsession); % Set random seed for resampling
    temp = intersect(index,find(log_odds.behavioural_state == 1));
    temp = datasample(s,temp,length(temp),'Replace',false);
    index_resampled = [index_resampled temp];
    index11 = [index11 temp];

    s = RandStream('mrg32k3a','Seed',nsession+100); % Set random seed for resampling
    temp = intersect(index,find(log_odds.behavioural_state == 2));
    temp = datasample(s,temp,length(temp),'Replace',false);
    index_resampled = [index_resampled temp];
    index22 = [index22 temp];
end
index_resampled = [index11 index22];

% index = find(log_odds.experiment ==8 );
index = find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 7 );
index = intersect(index,find(log_odds.ripple_peak >=3));
index = intersect(index,find(log_odds.behavioural_state > 0));
index1 = intersect(index,find( log_odds.behavioural_state == 1));
index2 = intersect(index,find( log_odds.behavioural_state == 2));
index = [index1 index2];

timebin = 1;


fig = figure
fig.Position = [750 240 950 640];
fig.Name = 'On track awake reactivation HPC vs V1 log odds';
subplot(2,2,1)
scatter(log_odds.zscore(index1),log_odds.V1_log_odds_R(index1,4),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.2')
hold on
scatter(log_odds.zscore(index2),log_odds.V1_log_odds_R(index2,4),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.2')
xlabel('HPC log odds')
ylabel('V1 Right log odds')
mdl = fitlm(log_odds.zscore(index)',log_odds.V1_log_odds_R(index,4)');
[pval,~,~] = coefTest(mdl);
x =[min(log_odds.zscore(index)') max(log_odds.V1_log_odds_R(index)')];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);

if pval <= 0.05
    plot(x,y_est,'r:')
    %         xlim([-30 120])

    %     title(sprintf('Session %i',s),'Color','red')
    %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
    text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
else
    plot(x,y_est,'k:')
    %         xlim([-30 120])
    %     title(sprintf('Session %i',s))
    %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
    text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
end
legend('Track 1','Track 2','Color','none')
set(gca,"TickDir","out",'box', 'off','Color','none')
title('HPC vs V1 right log odds')

subplot(2,2,2)
scatter(log_odds.zscore(index1),log_odds.V1_log_odds_L(index1,4),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.2')
hold on
scatter(log_odds.zscore(index2),log_odds.V1_log_odds_L(index2,4),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.2')
xlabel('HPC log odds')
ylabel('V1 Left log odds')
mdl = fitlm(log_odds.zscore(index)',log_odds.V1_log_odds_L(index,4)');
[pval,~,~] = coefTest(mdl);
x =[min(log_odds.zscore(index)') max(log_odds.zscore(index)')];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);

if pval <= 0.05
    plot(x,y_est,'r:')
    %         xlim([-30 120])

    %     title(sprintf('Session %i',s),'Color','red')
    %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
    text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
else
    plot(x,y_est,'k:')
    %         xlim([-30 120])
    %     title(sprintf('Session %i',s))
    %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
    text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
end
legend('Track 1','Track 2','Color','none')
set(gca,"TickDir","out",'box', 'off','Color','none')
title('HPC vs V1 left log odds')

subplot(2,2,3)
scatter(log_odds.zscore(index1),log_odds.V1_log_odds_R(index11,4),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.2')
hold on
scatter(log_odds.zscore(index2),log_odds.V1_log_odds_R(index22,4),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.2')
xlabel('HPC log odds')
ylabel('V1 Right log odds')
mdl = fitlm(log_odds.zscore(index)',log_odds.V1_log_odds_R(index_resampled,4)');
[pval,~,~] = coefTest(mdl);
x =[min(log_odds.zscore(index)') max(log_odds.zscore(index)')];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);

if pval <= 0.05
    plot(x,y_est,'r:')
    %         xlim([-30 120])

    %     title(sprintf('Session %i',s),'Color','red')
    %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
    text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
else
    plot(x,y_est,'k:')
    %         xlim([-30 120])
    %     title(sprintf('Session %i',s))
    %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
    text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
end
legend('Track 1','Track 2','Color','none')
set(gca,"TickDir","out",'box', 'off','Color','none')
title('HPC vs V1 right log odds (within track resampled)')

subplot(2,2,4)
scatter(log_odds.zscore(index1),log_odds.V1_log_odds_L(index11,4),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.2')
hold on
scatter(log_odds.zscore(index2),log_odds.V1_log_odds_L(index22,4),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.2')
xlabel('HPC log odds')
ylabel('V1 Left log odds')
mdl = fitlm(log_odds.zscore(index)',log_odds.V1_log_odds_L(index_resampled,4)');
[pval,~,~] = coefTest(mdl);
x =[min(log_odds.zscore(index)') max(log_odds.HPC_bayesian_bias(index_resampled)')];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);

if pval <= 0.05
    plot(x,y_est,'r:')
    %         xlim([-30 120])

    %     title(sprintf('Session %i',s),'Color','red')
    %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
    text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
else
    plot(x,y_est,'k:')
    %         xlim([-30 120])
    %     title(sprintf('Session %i',s))
    %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
    text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
end
legend('Track 1','Track 2','Color','none')
set(gca,"TickDir","out",'box', 'off','Color','none')
title('HPC vs V1 left log odds (within track resampled)')
save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project\Log odds analysis',[])




%%
scatter3(log_odds.HPC_bayesian_bias(index1),log_odds.V1_bayesian_bias_L(index1,4),log_odds.V1_bayesian_bias_R(index1,4),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.2')
hold on
scatter3(log_odds.HPC_bayesian_bias(index2),log_odds.V1_bayesian_bias_L(index2,4),log_odds.V1_bayesian_bias_R(index2,4),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.2')
xlabel('HPC bayesian bias')
ylabel('V1 Left bayesian bias')
zlabel('V1 Right bayesian bias')



index = find(log_odds.experiment <= 0|log_odds.experiment == 5  |log_odds.experiment >= 8  );

% index = find(log_odds.experiment >0 );
index = intersect(index,find(log_odds.ripple_peak <=3));

index1 = intersect(index,find(log_odds.behavioural_state ==1));
scatter(log_odds.ripple_peak(index1),log_odds.HPC_bayesian_bias(index1),'r');hold on

index2 = intersect(index,find(log_odds.behavioural_state ==2));
scatter(log_odds.ripple_peak(index2),log_odds.HPC_bayesian_bias(index2),'b');hold on

% scatter(log_odds.V1_participation_R(index2,4),log_odds.HPC_bayesian_bias(index2),'b')

index2 = intersect(index,find(log_odds.behavioural_state ==2));
scatter(log_odds.ripple_peak(index),log_odds.V1_participation_R(index,4),'b')
% scatter(log_odds.V1_participation_R(index2,4),log_odds.HPC_bayesian_bias(index2),'b')


scatter(log_odds.V1_participation_R(index1,6),log_odds.V1_bayesian_bias_R(index1,4))
scatter(log_odds.V1_participation_R(index2,6),log_odds.V1_bayesian_bias_R(index2,4))


scatter(log_odds.V1_participation_L(index,6),log_odds.V1_bayesian_bias_L(index,4))

figure
for nwin = 1:6
    subplot(2,5,nwin)
    scatter(log_odds.ripple_peak(index),log_odds.V1_bayesian_bias_R(index,nwin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.2')
    hold on
    mdl = fitlm(log_odds.ripple_peak(index)',log_odds.V1_bayesian_bias_R(index,nwin)');
    [pval,~,~] = coefTest(mdl);
    x =[min(log_odds.ripple_peak(index)') max(log_odds.ripple_peak(index)')];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);

    if pval <= 0.05
        plot(x,y_est,'r:')
        %         xlim([-30 120])

        %     title(sprintf('Session %i',s),'Color','red')
        %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
        text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
    else
        plot(x,y_est,'k:')
        %         xlim([-30 120])
        %     title(sprintf('Session %i',s))
        %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
        text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
    end
end


figure
for nwin = 1:6
    subplot(2,5,nwin)
    scatter(log_odds.V1_participation_L(index2,nwin),log_odds.V1_bayesian_bias_L(index2,nwin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.2')
    hold on
    mdl = fitlm(log_odds.V1_participation_L(index2,nwin)',log_odds.V1_bayesian_bias_L(index2,nwin)');
    [pval,~,~] = coefTest(mdl);
    x =[min(log_odds.V1_participation_L(index2,nwin)') max(log_odds.V1_participation_L(index2,nwin)')];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);

    if pval <= 0.05
        plot(x,y_est,'r:')
        %         xlim([-30 120])

        %     title(sprintf('Session %i',s),'Color','red')
        %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
        text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
    else
        plot(x,y_est,'k:')
        %         xlim([-30 120])
        %     title(sprintf('Session %i',s))
        %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
        text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
    end
end
%
scatter(log_odds.V1_participation_L(index,6),log_odds.V1_bayesian_bias_R(index,4))


plot(log_odds.V1_bayesian_bias_R(index,4),'b');hold on
plot(log_odds.V1_participation_R(index,6),'r');hold on


index1 = intersect(index,find(log_odds.behavioural_state == 1));
% index1 = datasample(index1,length(index1),'Replace',false);
index2 = intersect(index,find(log_odds.behavioural_state == 2));
% index2 = datasample(index2,length(index2),'Replace',false);

index = [index1 index2];
% index1 = intersect(index,find( log_odds.behavioural_state == 0));hold on;
% index2 = intersect(index,find( log_odds.behavioural_state == 2));
scatter(log_odds.HPC_bayesian_bias(index),log_odds.V1_bayesian_bias_R(index,4),'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.2')
hold on
scatter(log_odds.HPC_bayesian_bias(index1),log_odds.V1_bayesian_bias_R(index1,3),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.2')
scatter(log_odds.HPC_bayesian_bias(index2),log_odds.V1_bayesian_bias_R(index2,3),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.2')
% scatter3(log_odds.HPC_bayesian_bias(index2),log_odds.V1_bayesian_bias_L(index2,nbin),log_odds.V1_bayesian_bias_R(index2,nbin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.2')
xlabel('HPC bayesian bias')
% ylabel('V1 Left bayesian bias')
ylabel('V1 Right participation')
mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.V1_bayesian_bias_R(index,4)');
[pval,~,~] = coefTest(mdl);
x =[min(log_odds.HPC_bayesian_bias(index)') max(log_odds.HPC_bayesian_bias(index)')];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);

if pval <= 0.05
    plot(x,y_est,'r:')
    %         xlim([-30 120])

    %     title(sprintf('Session %i',s),'Color','red')
    %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
    text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
else
    plot(x,y_est,'k:')
    %         xlim([-30 120])
    %     title(sprintf('Session %i',s))
    %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
    text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
end
set(gca,"TickDir","out",'box', 'off','Color','none')



index = find(log_odds.experiment == 1|log_odds.experiment == 5  |log_odds.experiment >= 8  );

% index = find(log_odds.experiment >0 );
index = intersect(index,find(log_odds.ripple_peak >=3));

index1 = intersect(index,find(log_odds.behavioural_state == 1));
% index1 = datasample(index1,length(index1),'Replace',false);
index2 = intersect(index,find(log_odds.behavioural_state == 2));
% index2 = datasample(index2,length(index2),'Replace',false);
index = intersect(index,find(log_odds.behavioural_state ==0));

histogram(log_odds.HPC_bayesian_bias(index),30,'FaceColor','k','FaceAlpha',0.2,'Normalization','probability');hold on;
histogram(log_odds.HPC_bayesian_bias(index1),30,'FaceColor','r','FaceAlpha',0.2,'Normalization','probability');hold on;
histogram(log_odds.HPC_bayesian_bias(index2),30,'FaceColor','b','FaceAlpha',0.2,'Normalization','probability');hold on;


% histogram(log_odds.zscore(index),30,'FaceColor','k','FaceAlpha',0.2,'Normalization','probability');hold on;
histogram(log_odds.zscore(index1(~isnan(log_odds.zscore(index1)))),20,'FaceColor','r','FaceAlpha',0.2,'Normalization','probability');hold on;
histogram(log_odds.zscore(index2(~isnan(log_odds.zscore(index2)))),20,'FaceColor','b','FaceAlpha',0.2,'Normalization','probability');hold on;

scatter(log_odds.zscore(index),log_odds.V1_bayesian_bias_R(index1,5))



scatter(log_odds.zscore(index1),log_odds.V1_log_odds_R(index1,3),'r','MarkerFaceAlpha',0.2)
hold on
scatter(log_odds.zscore(index2),log_odds.V1_log_odds_R(index2,3),'b','MarkerFaceAlpha',0.2)

histogram(log_odds.V1_bayesian_bias_L(index1(~isnan(log_odds.V1_bayesian_bias_L(index1,5))),5),20,'FaceColor','r','FaceAlpha',0.2,'Normalization','cdf');hold on;
histogram(log_odds.V1_bayesian_bias_L(index2(~isnan(log_odds.V1_bayesian_bias_L(index2,5))),5),20,'FaceColor','b','FaceAlpha',0.2,'Normalization','cdf');hold on;


histogram(log_odds.V1_bayesian_bias_R(index1(~isnan(log_odds.V1_bayesian_bias_R(index1,5))),5),20,'FaceColor','r','FaceAlpha',0.2,'Normalization','cdf');hold on;
histogram(log_odds.V1_bayesian_bias_R(index2(~isnan(log_odds.V1_bayesian_bias_R(index2,5))),5),20,'FaceColor','b','FaceAlpha',0.2,'Normalization','cdf');hold on;


histogram(log_odds.V1_log_odds_L(index1(~isnan(log_odds.V1_log_odds_L(index1,5))),5),20,'FaceColor','r','FaceAlpha',0.2,'Normalization','cdf');hold on;
histogram(log_odds.V1_log_odds_L(index2(~isnan(log_odds.V1_log_odds_L(index2,5))),5),20,'FaceColor','b','FaceAlpha',0.2,'Normalization','cdf');hold on;


histogram(log_odds.V1_log_odds_R(index1(~isnan(log_odds.V1_log_odds_R(index1,5))),5),20,'FaceColor','r','FaceAlpha',0.2,'Normalization','cdf');hold on;
histogram(log_odds.V1_log_odds_R(index2(~isnan(log_odds.V1_log_odds_R(index2,5))),5),20,'FaceColor','b','FaceAlpha',0.2,'Normalization','cdf');hold on;



histogram(log_odds.V1_participation_L(index1(~isnan(log_odds.V1_log_odds_L(index1,4))),4),20,'FaceColor','r','FaceAlpha',0.2,'Normalization','cdf');hold on;
histogram(log_odds.V1_participation_L(index2(~isnan(log_odds.V1_log_odds_L(index2,4))),4),20,'FaceColor','b','FaceAlpha',0.2,'Normalization','cdf');hold on;

histogram(log_odds.V1_participation_R(index1(~isnan(log_odds.V1_log_odds_R(index1,4))),4),20,'FaceColor','r','FaceAlpha',0.2,'Normalization','cdf');hold on;
histogram(log_odds.V1_participation_R(index2(~isnan(log_odds.V1_log_odds_R(index2,4))),4),20,'FaceColor','b','FaceAlpha',0.2,'Normalization','cdf');hold on;

histogram(log_odds.V1_participation_R(index1(~isnan(log_odds.V1_log_odds_R(index1,4))),4),20,'FaceColor','r','FaceAlpha',0.2,'Normalization','cdf');hold on;
histogram(log_odds.V1_participation_R(index2(~isnan(log_odds.V1_log_odds_R(index2,4))),4),20,'FaceColor','b','FaceAlpha',0.2,'Normalization','cdf');hold on;


%% V1 participation

% time_windows = [31:35;36:40;41:45;46:50;...
%     51:55;56:60;61:65;66:70];
% time_windows_text = {'-300 to -400ms','-200 to -300ms','-100 to -200ms','0 to -100ms','0 to 100ms','100 to 200ms','200 to 300ms','300ms to 400ms'};
behavioural_epoches = [-1,0,1,2];
behavioural_epoches_text = {'PRE','POST','Track 1','Track 2'};
time_windows = [3; 4; 5; 6; 7; 8];
time_windows_text = {'-600 to -400ms','-400 to -200ms','0 to -200ms','0 to 200ms','200 to 400ms','400 to 600ms'};


for epoch = 3:4
    c = 1;
    fig = figure(1)
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    %         fig.Name = sprintf('session %i V1 cell participation (Right hemisphere) vs HPC log odds during %s',nsession,behavioural_epoches_text{epoch})

    for nsession = [1 5 8 9 10]
        %         fig = figure
        %         fontsize(fig, 14, "points")
        %         fig.Position = [661 431 945 512];
        %         fig.Name = sprintf('session %i V1 cell participation (Right hemisphere) vs HPC log odds during %s',nsession,behavioural_epoches_text{epoch})

        for nwin =1:size(time_windows,1)
            time_bin= time_windows(nwin,:);

            subplot(2,5,nwin)
            hold on
            index = find( log_odds.behavioural_state == behavioural_epoches(epoch));
            %             index = find( log_odds.behavioural_state> 0);

            %         index = intersect(index,find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >=  7));
            index = intersect(index,find(log_odds.experiment == nsession));
            %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
            index = intersect(index,find(log_odds.ripple_peak >=3));
            %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));

            if length(index) < 10
                continue
            end

            scatter(log_odds.HPC_bayesian_bias(index),log_odds.V1_participation_R(index,time_bin),'MarkerFaceAlpha','0.1')
            %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
            mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.V1_participation_R(index,time_bin)');
            [pval,~,~] = coefTest(mdl);
            x =[min(log_odds.HPC_bayesian_bias') max(log_odds.HPC_bayesian_bias')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'r:')
                %         xlim([-30 120])

                %     title(sprintf('Session %i',s),'Color','red')
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
            else
                plot(x,y_est,'k:')
                %         xlim([-30 120])
                %     title(sprintf('Session %i',s))
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end

            %             index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.5));
            %             scatter(log_odds.HPC_bayesian_bias(index2),log_odds.V1_participation_R(index2,time_bin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
            %             hold on
            %             index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.5));
            %             scatter(log_odds.HPC_bayesian_bias(index1),log_odds.V1_participation_R(index1,time_bin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

            %     c = c + 1;
            %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
            title(time_windows_text{nwin})

            set(gca,"TickDir","out",'box', 'off','Color','none')

        end
    end

    fig = figure(2)
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('session %i V1 cell participation (Left hemisphere) vs HPC log odds during %s',nsession,behavioural_epoches_text{epoch})


    for nsession = [1 5 8 9 10]


        for nwin =1:size(time_windows,1)
            time_bin= time_windows(nwin,:);

            subplot(2,5,nwin)
            hold on

            if sum(isnan(log_odds.V1_participation_L(index,time_bin))) == length(index)
                continue
            end

            if length(index) < 10
                continue
            end



            scatter(log_odds.HPC_bayesian_bias(index),log_odds.V1_participation_L(index,time_bin),'MarkerFaceAlpha','0.1')
            %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
            mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.V1_participation_L(index,time_bin)');
            [pval,~,~] = coefTest(mdl);
            x =[min(log_odds.HPC_bayesian_bias') max(log_odds.HPC_bayesian_bias')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'r:')
                %         xlim([-30 120])

                %     title(sprintf('Session %i',s),'Color','red')
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
            else
                plot(x,y_est,'k:')
                %         xlim([-30 120])
                %     title(sprintf('Session %i',s))
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end
            %
            %             index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.5));
            %             scatter(log_odds.HPC_bayesian_bias(index2),log_odds.V1_participation_L(index2,time_bin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
            %             hold on
            %             index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.5));
            %             scatter(log_odds.HPC_bayesian_bias(index1),log_odds.V1_participation_L(index1,time_bin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

            %     c = c + 1;
            %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
            set(gca,"TickDir","out",'box', 'off','Color','none')
            title(time_windows_text{nwin})
        end

    end
end

save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project',[])


%% participation combined

% time_windows = [31:35;36:40;41:45;46:50;...
%     51:55;56:60;61:65;66:70];
% time_windows_text = {'-300 to -400ms','-200 to -300ms','-100 to -200ms','0 to -100ms','0 to 100ms','100 to 200ms','200 to 300ms','300ms to 400ms'};
behavioural_epoches = [-1,0,1,2];
behavioural_epoches_text = {'PRE','POST','Track 1','Track 2'};
% time_windows = [3; 4; 5; 6; 7; 8];
time_windows_text = {'-600 to -400ms','-400 to -200ms','0 to -200ms','0 to 200ms','200 to 400ms','400 to 600ms'};


time_windows = [1;2;3; 4; 5; 6; 7; 8;9;10];
time_windows_text = {'-1000 to -800ms','-800 to -600','-600 to -400ms','-400 to -200ms','0 to -200ms','0 to 200ms','200 to 400ms','400 to 600ms','600 to 800ms','800 to 1000ms'};


for epoch = 3:length(behavioural_epoches)
    c = 1;
    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between V1 cell participation (Right hemisphere) and HPC log odds during %s',behavioural_epoches_text{epoch})
    for nwin =1:size(time_windows,1)
        time_bin= time_windows(nwin,:);

        subplot(2,5,nwin)
        hold on
        index = find( log_odds.behavioural_state == behavioural_epoches(epoch));
        %         index = find( log_odds.behavioural_state> 0);
        index = intersect(index,find(log_odds.experiment ==1 |log_odds.experiment == 5  |log_odds.experiment >=8));
        %         index = intersect(index,find(log_odds.experiment == nsession));
        %     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
        index = intersect(index,find(log_odds.ripple_peak >=3));
        %     index = intersect(find( log_odds.behavioural_state == 0),find( log_odds.ripple_peak >=3));



        scatter(log_odds.HPC_bayesian_bias(index),(log_odds.V1_participation_R(index,time_bin)-mean(log_odds.V1_participation_R(index,1:9:10),2)),'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %         scatter(log_odds.HPC_bayesian_bias(index),log_odds.V1_participation_R(index,time_bin),'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %         mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.V1_participation_R(index,time_bin)');
        mdl = fitlm(log_odds.HPC_bayesian_bias(index)',(log_odds.V1_participation_R(index,time_bin)-mean(log_odds.V1_participation_R(index,:),2))');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.HPC_bayesian_bias(index)') max(log_odds.HPC_bayesian_biasindex')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.4));
        scatter(log_odds.HPC_bayesian_bias(index2),(log_odds.V1_participation_R(index2,time_bin)-mean(log_odds.V1_participation_R(index2,1:9:10),2)),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.6));
        scatter(log_odds.HPC_bayesian_bias(index1),(log_odds.V1_participation_R(index1,time_bin)-mean(log_odds.V1_participation_R(index1,1:9:10),2)),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        title(time_windows_text{nwin})

        set(gca,"TickDir","out",'box', 'off','Color','none')

    end


    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between V1 cell participation (Left hemisphere) and HPC log odds during %s',behavioural_epoches_text{epoch})

    for nwin =1:size(time_windows,1)
        time_bin= time_windows(nwin,:);
        if sum(isnan(log_odds.V1_participation_L(index,time_bin))) == length(log_odds.V1_participation_L(index,time_bin))
            continue
        end
        subplot(2,5,nwin)
        hold on
        scatter(log_odds.HPC_bayesian_bias(index),(log_odds.V1_participation_L(index,time_bin)-mean(log_odds.V1_participation_L(index,1:9:10),2)),'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.HPC_bayesian_bias(index)',(log_odds.V1_participation_L(index,time_bin)-mean(log_odds.V1_participation_L(index,1:9:10),2))');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.HPC_bayesian_bias(index)') max(log_odds.HPC_bayesian_bias(index)')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.4));
        scatter(log_odds.HPC_bayesian_bias(index2),(log_odds.V1_participation_L(index2,time_bin)-mean(log_odds.V1_participation_L(index2,1:9:10),2)),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.6));
        scatter(log_odds.HPC_bayesian_bias(index1),(log_odds.V1_participation_L(index1,time_bin)-mean(log_odds.V1_participation_L(index1,1:9:10),2)),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        set(gca,"TickDir","out",'box', 'off','Color','none')
        title(time_windows_text{nwin})
    end
end


save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project',[])



for epoch = 3:length(behavioural_epoches)
    c = 1;
    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between V1 cell participation (Right hemisphere) and HPC log odds during %s',behavioural_epoches_text{epoch})

    time_bin= time_windows(nwin,:);

    index_resampled = [];

    for nsession = [ 5 8 9 10]

        index = find(log_odds.experiment ==nsession );
        index = intersect(index,find(log_odds.ripple_peak >=3));

        s = RandStream('mrg32k3a','Seed',nsession); % Set random seed for resampling
        temp = intersect(index,find(log_odds.behavioural_state == behavioural_epoches(epoch)));
        temp = datasample(s,temp,length(temp),'Replace',false);
        index_resampled = [index_resampled temp];

    end

    % index = find(log_odds.experiment ==8 );
    index = find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 8 );
    index = intersect(index,find(log_odds.ripple_peak >=3));
    index = intersect(index,find(log_odds.behavioural_state== behavioural_epoches(epoch)));

    for nwin =1:size(time_windows,1)

        time_bin= time_windows(nwin,:);

        subplot(2,5,nwin)
        hold on

        %         scatter(log_odds.HPC_bayesian_bias(index),log_odds.V1_participation_R(index,time_bin)-mean(log_odds.V1_participation_R(index,:),2),'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        scatter(log_odds.HPC_bayesian_bias(index),log_odds.V1_participation_R(index_resampled,time_bin),'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.V1_participation_R(index_resampled,time_bin)');
        %         mdl = fitlm(log_odds.HPC_bayesian_bias(index)',(log_odds.V1_participation_R(index,time_bin)-mean(log_odds.V1_participation_R(index,:),2))');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.HPC_bayesian_bias') max(log_odds.HPC_bayesian_bias')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.4));
        index2_resampled = intersect(index_resampled,find(log_odds.HPC_bayesian_bias < 0.4));
        scatter(log_odds.HPC_bayesian_bias(index2),log_odds.V1_participation_R(index2_resampled,time_bin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.6));
        index1_resampled = intersect(index_resampled,find(log_odds.HPC_bayesian_bias > 0.6));
        scatter(log_odds.HPC_bayesian_bias(index1),log_odds.V1_participation_R(index1_resampled,time_bin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        title(time_windows_text{nwin})

        set(gca,"TickDir","out",'box', 'off','Color','none')

    end


    fig = figure
    fontsize(fig, 14, "points")
    fig.Position = [661 431 945 512];
    fig.Name = sprintf('relationship between V1 cell participation (Left hemisphere) and HPC log odds during %s',behavioural_epoches_text{epoch})

    for nwin =1:size(time_windows,1)
        time_bin= time_windows(nwin,:);
        if sum(isnan(log_odds.V1_participation_L(index,time_bin))) == length(log_odds.V1_participation_L(index,time_bin))
            continue
        end
        subplot(2,5,nwin)
        hold on
        scatter(log_odds.HPC_bayesian_bias(index),log_odds.V1_participation_L(index_resampled,time_bin),'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha','0.1')
        %     mdl = fitlm(log_odds.zscore_substracted(index)',log_odds.V1_spike_count_R(index,time_bin)');
        mdl = fitlm(log_odds.HPC_bayesian_bias(index)',log_odds.V1_participation_L(index_resampled,time_bin)');
        [pval,~,~] = coefTest(mdl);
        x =[min(log_odds.HPC_bayesian_bias') max(log_odds.HPC_bayesian_bias')];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);

        if pval <= 0.05
            plot(x,y_est,'r:')
            %         xlim([-30 120])

            %     title(sprintf('Session %i',s),'Color','red')
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
        else
            plot(x,y_est,'k:')
            %         xlim([-30 120])
            %     title(sprintf('Session %i',s))
            %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
            text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
        end

        %         index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.4));
        scatter(log_odds.HPC_bayesian_bias(index2),log_odds.V1_participation_L(index2_resampled,time_bin),'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
        hold on
        %         index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.6));
        scatter(log_odds.HPC_bayesian_bias(index1),log_odds.V1_participation_L(index1_resampled,time_bin),'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')

        %     c = c + 1;
        %     title(sprintf('%.2f ms relative to CA1 candidate events',(1/10*nbin-1)*1000))
        set(gca,"TickDir","out",'box', 'off','Color','none')
        title(time_windows_text{nwin})
    end
end




%%
index1 = find(log_odds.zscore_substracted > 1);
for n = 1:length(index1)
    plot(log_odds.V1_spike_count_R(index1(n),:),'r')
    hold on
end

for n = 1:length(index2)
    index2 = find(log_odds.zscore_substracted < -1);
    plot(log_odds.V1_spike_count_R(index2(n),:),'b')
    hold on
end


c = 1;
for nsession = 1:10
    index1 = find(log_odds.zscore > 1 & log_odds.experiment == nsession & log_odds.behavioural_state >0);
    index2 = find(log_odds.zscore < -1 & log_odds.experiment == nsession & log_odds.behavioural_state >0);
    subplot(4,6,c)
    plot(nanmean(log_odds.V1_spike_count_R(index1,:)),'r')
    hold on
    plot(nanmean(log_odds.V1_spike_count_R(index2,:)),'b')
    xlim([25 75])
    xticks(25:10:75)
    xticklabels(-500:200:500)
    xline(50,'k','LineWidth',2)
    ylabel('Averaged FR (z)')
    xlabel('Time (ms)')
    title(sprintf('Right session %i',nsession))

    c = c + 1;
    subplot(4,6,c)
    plot(nanmean(log_odds.V1_spike_count_L(index1,:)),'r')
    hold on
    plot(nanmean(log_odds.V1_spike_count_L(index2,:)),'b')
    title(sprintf('Left session %i',nsession))
    xlim([25 75])
    xticks(25:10:75)
    xticklabels(-500:200:500)
    xline(50,'k','LineWidth',2)
    ylabel('Averaged FR (z)')
    xlabel('Time (ms)')
    c = c + 1;
end
legend('log odds > 0.5','log odds < -0.5')
sgtitle('spiking')

figure
c = 1;
for nsession = 1:10
    index1 = find(log_odds.zscore > 0.5 & log_odds.experiment == nsession & log_odds.behavioural_state >= 0);
    index2 = find(log_odds.zscore < -0.5 & log_odds.experiment == nsession & log_odds.behavioural_state >= 0);
    subplot(4,6,c)
    plot(nanmean(log_odds.V1_participation_R(index1,:)-nanmean(log_odds.V1_participation_R(index1,:),2)),'r')
    hold on
    plot(nanmean(log_odds.V1_participation_R(index2,:)-nanmean(log_odds.V1_participation_R(index2,:),2)),'b')
    %     xlim([25 75])
    title(sprintf('Right session %i',nsession))

    c = c + 1;
    subplot(4,6,c)
    plot(nanmean(log_odds.V1_participation_L(index1,:)-nanmean(log_odds.V1_participation_L(index1,:),2)),'r')
    hold on
    plot(nanmean(log_odds.V1_participation_L(index2,:)-nanmean(log_odds.V1_participation_L(index2,:),2)),'b')
    title(sprintf('Left session %i',nsession))
    %     xlim([25 75])
    c = c + 1;
end
legend('log odds > 1','log odds < -1')
sgtitle('participation')


figure
c = 1;
for nsession = 1:10
    index1 = find(log_odds.HPC_bayesian_bias > 0.6 & log_odds.experiment == nsession & log_odds.behavioural_state >= 0);
    index2 = find(log_odds.HPC_bayesian_bias < 0.4 & log_odds.experiment == nsession & log_odds.behavioural_state >= 0);
    subplot(4,6,c)
    for n = 1:length(index1)
        plot(log_odds.V1_participation_R(index1(n),:)-mean(log_odds.V1_participation_R(index1(n),:),2),'r','LineWidth',0.2)
        hold on
    end

    for n = 1:length(index2)
        plot(log_odds.V1_participation_R(index2(n),:)-mean(log_odds.V1_participation_R(index2(n),:),2),'b','LineWidth',0.2)
        %     xlim([25 75])
    end

    title(sprintf('Right session %i',nsession))

    c = c + 1;
    subplot(4,6,c)
    for n = 1:length(index1)
        plot(log_odds.V1_participation_L(index1(n),:)-mean(log_odds.V1_participation_L(index1(n),:),2),'r','LineWidth',0.2)
        hold on
    end

    for n = 1:length(index2)
        plot(log_odds.V1_participation_L(index2(n),:)-mean(log_odds.V1_participation_L(index2(n),:),2),'b','LineWidth',0.2)
        %     xlim([25 75])
    end
    title(sprintf('Left session %i',nsession))
    %     xlim([25 75])
    c = c + 1;
end
legend('HPC_bayesian_bias > 0.5','HPC_bayesian_bias < 0.4')

legend('log odds > 1','log odds < -1')
sgtitle('participation')


for nsession = 1:10
    subplot(2,5,nsession)
    hold on
    scatter(track_label{1}  .* (rand(1,length(z_log_odds{1}))),z_log_odds{1},'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')
    hold on
    scatter(track_label{2}./2 .* (2+rand(1,length(z_log_odds{2}))),z_log_odds{2},'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
    xticks([0.5 2.5])
    xticklabels({'Track 1','Track 2'})
    xlabel('Track ID')
    ylabel('zscored log odds')
    set(gca,"TickDir","out",'box', 'off','Color','none')
end


index = find(log_odds.experiment <= 2|log_odds.experiment == 5  |log_odds.experiment >= 7);
%         index = intersect(index,find(log_odds.experiment == nsession));
%     index = find( log_odds.behavioural_state == 0 | log_odds.behavioural_state == 2);
index = intersect(index,find(log_odds.ripple_peak >=3));

index1 = intersect(index,find(log_odds.behavioural_state ==1));
index2 = intersect(index,find(log_odds.behavioural_state ==2));

% index1 = intersect(index,find(log_odds.HPC_bayesian_bias > 0.6));
% index2 = intersect(index,find(log_odds.HPC_bayesian_bias < 0.4));


scatter(log_odds.HPC_bayesian_bias(index),log_odds.V1_bayesian_bias_L(index,4))


scatter(log_odds.HPC_bayesian_bias(index),log_odds.V1_bayesian_bias_R(index,4))



figure
for nwin = 1:10
    subplot(2,5,nwin)
    histogram(log_odds.V1_participation_R(index1,nwin)-mean(log_odds.V1_participation_R(index1,1:9:10),2),-0.15:0.01:0.15,'FaceColor','r','FaceAlpha',0.4,"Normalization","cdf")
    hold on
    histogram(log_odds.V1_participation_R(index2,nwin)-mean(log_odds.V1_participation_R(index2,1:9:10),2),-0.15:0.01:0.15,'FaceColor','b','FaceAlpha',0.4,"Normalization","cdf")
end

figure
for nwin = 1:10
    subplot(2,5,nwin)
    histogram(log_odds.V1_participation_L(index1,nwin)-mean(log_odds.V1_participation_L(index1,1:9:10),2),-0.1:0.01:0.1,'FaceColor','r','FaceAlpha',0.4,"Normalization","cdf")
    hold on
    histogram(log_odds.V1_participation_L(index2,nwin)-mean(log_odds.V1_participation_L(index2,1:9:10),2),-0.1:0.01:0.1,'FaceColor','b','FaceAlpha',0.4,"Normalization","cdf")
end
%% Peri event event histogram

event1_name = [];
event2_name = [];
event1 = {ripples.probe(2).peaktimes',reactivations.probe(2).midpoint(reactivations.probe(2).ripple_peak >= 3)'};
event2 = {ripples.probe(1).peaktimes',reactivations.probe(1).midpoint(reactivations.probe(1).ripple_peak >= 3)'};
event1_name = {'Ripple probe 1','CA1 bursting events probe 1'};
event2_name = {'Ripple probe 2','CA1 bursting events probe 2'};

figure
for n = 1:length(event1)
    nexttile
    plot_perievent_event_histogram(event1{n},event2{n},'twin',[-1 1])

    title(sprintf('%s relative to %s',event1_name{n},event2_name{n}))
    set(gca,"TickDir","out",'box', 'off','Color','none')
end
sgtitle('Awake')


histogram(reactivations.probe(2).midpoint,50)
hold on
histogram(reactivations.probe(1).midpoint,50)

histogram(reactivations.probe(1).midpoint,50)
hold on
histogram(ripples.probe(1).peaktimes,50)
hold on
histogram(ripples.probe(2).peaktimes,50)


plot_perievent_event_histogram(ripples.probe(2).onset,ripples.probe(1).onset,'twin',[-1 1])
%%

% %% Detect ripple and cortical slow wave oscillation (SO) and cortical spindles
%
% % Detect CA1 populational bursting events (Candidate events)
% zscore_min = 0;
% zscore_max = 3;
%
% cd(options.EPHYS_DATAPATH)
% channel_to_use = find(sorted_config.Channel == best_channels.CA1_channel);
% [replay,reactivations] = detect_candidate_events_masa(tvec,raw_LFP(channel_to_use,:),...
%     CA1_clusters.MUA_zscore,[CA1_clusters.spike_id CA1_clusters.spike_times],peripherals,zscore_min,zscore_max,options)
% save extracted_candidate_events replay reactivations
%
% channel_to_use = find(sorted_config.Channel == best_channels.L5_channel);
% [~,L5_reactivations] = detect_candidate_events_masa(tvec,raw_LFP(channel_to_use,:),...
%     L5_clusters.MUA_zscore,[L5_clusters.spike_id L5_clusters.spike_times],peripherals,zscore_min,zscore_max,options)
% save extracted_L5_candidate_events L5_reactivations
%
% % Detect CA1 ripple events
% channel_to_use = find(sorted_config.Channel == best_channels.CA1_channel);
% [ripples] = FindRipples_masa(raw_LFP(channel_to_use,:)',tvec','minDuration',20,'durations',[30 200],'frequency',SR,...
%     'noise',raw_LFP(2,:)','passband',[125 300],'thresholds',[3 5])
%
% figure
% [ripples.SWS_offset,ripples.SWS_index] = RestrictInts(ripples.offset,behavioural_state.SWS);
% ripples.SWS_onset = ripples.onset(ripples.SWS_index);
% ripples.SWS_peaktimes = ripples.peaktimes(ripples.SWS_index);
%
% histogram(abs(ripples.SWS_onset-ripples.SWS_offset)',0:0.005:0.2,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r','Normalization','probability')
%
% [ripples.awake_offset,ripples.awake_index] = RestrictInts(ripples.offset,behavioural_state.quietWake);
% ripples.awake_onset = ripples.onset(ripples.awake_index);
% ripples.awake_peaktimes = ripples.peaktimes(ripples.awake_index);
% hold on
% histogram(abs(ripples.awake_onset-ripples.awake_offset)',0:0.005:0.2,'FaceColor','b','FaceAlpha',0.5,'EdgeColor','b','Normalization','probability')
% legend('NREM Ripples','awake Ripples')
% ylabel('Probability')
% xlabel('Duration (sec)')
%
% save extracted_ripples_events ripples
%
% % Detect Cortical ripple events
% channel_to_use = find(sorted_config.Channel == best_channels.first_in_brain_channel -12);% 240 micron
% [V1_ripples] = FindRipples_masa(raw_LFP(channel_to_use,:)',tvec','minDuration',20,'durations',[30 200],'frequency',SR,'noise',raw_LFP(1,:)','passband',[125 300])
% [V1_ripples.SWS_offset,V1_ripples.SWS_index] = RestrictInts(V1_ripples.offset,behavioural_state.SWS);
% V1_ripples.SWS_onset = V1_ripples.onset(V1_ripples.SWS_index);
% V1_ripples.SWS_peaktimes = V1_ripples.peaktimes(V1_ripples.SWS_index);
%
% [V1_ripples.awake_offset,V1_ripples.awake_index] = RestrictInts(V1_ripples.offset,behavioural_state.quietWake);
% V1_ripples.awake_onset = V1_ripples.onset(V1_ripples.awake_index);
% V1_ripples.awake_peaktimes = V1_ripples.peaktimes(V1_ripples.awake_index);
%
% V1_superficial_ripples = V1_ripples;
%
% channel_to_use = find(sorted_config.Channel == best_channels.L5_channel);% 240 micron
% [V1_ripples] = FindRipples_masa(raw_LFP(channel_to_use,:)',tvec','minDuration',20,'durations',[30 200],'frequency',SR,'noise',raw_LFP(1,:)','passband',[125 300])
% [V1_ripples.SWS_offset,V1_ripples.SWS_index] = RestrictInts(V1_ripples.offset,behavioural_state.SWS);
% V1_ripples.SWS_onset = V1_ripples.onset(V1_ripples.SWS_index);
% V1_ripples.SWS_peaktimes = V1_ripples.peaktimes(V1_ripples.SWS_index);
%
% [V1_ripples.awake_offset,V1_ripples.awake_index] = RestrictInts(V1_ripples.offset,behavioural_state.quietWake);
% V1_ripples.awake_onset = V1_ripples.onset(V1_ripples.awake_index);
% V1_ripples.awake_peaktimes = V1_ripples.peaktimes(V1_ripples.awake_index);
%
% V1_deep_ripples = V1_ripples;
% save extracted_V1_ripples_events V1_superficial_ripples V1_deep_ripples
%
% figure
% histogram(abs(V1_ripples.SWS_onset-V1_ripples.SWS_offset)',0:0.005:0.2,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r','Normalization','probability')
% hold on
% histogram(abs(V1_ripples.awake_onset-V1_ripples.awake_offset)',0:0.005:0.2,'FaceColor','b','FaceAlpha',0.5,'EdgeColor','b','Normalization','probability')
% legend('NREM cortical Ripples','awake cortical Ripples')
% ylabel('Probability')
% xlabel('Duration (sec)')
%
% % Detect cortical gamma events
% channel_to_use = find(sorted_config.Channel == best_channels.first_in_brain_channel -12);% 240 micron
% [gamma_events] = FindRipples_masa(raw_LFP(channel_to_use,:)',tvec','minDuration',20,'durations',[30 200],'frequency',SR,'noise',raw_LFP(1,:)','passband',[60 100])
% [gamma_events.SWS_offset,gamma_events.SWS_index] = RestrictInts(gamma_events.offset,behavioural_state.SWS);
% gamma_events.SWS_onset = gamma_events.onset(gamma_events.SWS_index);
% gamma_events.SWS_peaktimes = gamma_events.peaktimes(gamma_events.SWS_index);
%
% [gamma_events.awake_offset,gamma_events.awake_index] = RestrictInts(gamma_events.offset,behavioural_state.quietWake);
% gamma_events.awake_onset = gamma_events.onset(gamma_events.awake_index);
% gamma_events.awake_peaktimes = gamma_events.peaktimes(gamma_events.awake_index);
%
% V1_superficial_gamma_events = gamma_events;
%
% channel_to_use = find(sorted_config.Channel == best_channels.L5_channel);
% [gamma_events] = FindRipples_masa(raw_LFP(channel_to_use,:)',tvec','minDuration',20,'durations',[30 200],'frequency',SR,'noise',raw_LFP(1,:)','passband',[60 100])
% [gamma_events.SWS_offset,gamma_events.SWS_index] = RestrictInts(gamma_events.offset,behavioural_state.SWS);
% gamma_events.SWS_onset = gamma_events.onset(gamma_events.SWS_index);
% gamma_events.SWS_peaktimes = gamma_events.peaktimes(gamma_events.SWS_index);
%
% [gamma_events.awake_offset,gamma_events.awake_index] = RestrictInts(gamma_events.offset,behavioural_state.quietWake);
% gamma_events.awake_onset = gamma_events.onset(gamma_events.awake_index);
% gamma_events.awake_peaktimes = gamma_events.peaktimes(gamma_events.awake_index);
%
% V1_deep_gamma_events = gamma_events;
%
% save extracted_V1_gamma_events V1_superficial_gamma_events V1_deep_gamma_events
%
% % Detect Cortical spindle events
% [spindles] = FindSpindles_masa(raw_LFP(channel_to_use,:)',tvec','durations',[400 3000],'frequency',SR,'noise',raw_LFP(1,:)','passband',[9 17],'thresholds',[1.5 3])
% [spindles.SWS_offset,spindles.SWS_index] = RestrictInts(spindles.offset,SWS);
% spindles.SWS_onset = spindles.onset(spindles.SWS_index);
% spindles.SWS_peaktimes = spindles.peaktimes(spindles.SWS_index);
%
% [spindles.awake_offset,spindles.awake_index] = RestrictInts(spindles.offset,behavioural_state.quietWake);
% spindles.awake_onset = spindles.onset(spindles.awake_index);
% spindles.awake_peaktimes = spindles.peaktimes(spindles.awake_index);
% save extracted_spindles_events spindles
%
% % Detect Slow wave Up and Down states (Using all layer 5)
% options.importMode = 'KS';
% [spikes chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels.L5_channel-10 best_channels.L4_channel + 10] ,'group','Buz style');
%
% channel_to_use = find(sorted_config.Channel == best_channels.L5_channel);
% slow_waves= DetectSlowWaves_masa('time',tvec,'lfp',raw_LFP(channel_to_use,:)','NREMInts',behavioural_state.SWS,'spikes',spikes);
%
% % save best_channels best_spindle_channel best_spindle_channel best_ripple_channel
% save extracted_slow_waves slow_waves

%% Peri-event LFP amplitude and phase

% Stimulus
if strcmp(StimulusName,'replay_Masa2tracks') ;
    all_events = {MousePos.stimuli_onset(MousePos.stimuli_track == 1)',MousePos.stimuli_onset(MousePos.stimuli_track == 2)'};
    event_group = {'T1 stimuli','T2 stimuli'};


    lfpAvg = [];
    csd = [];
    lfpAvg.filter_type = {'SO','all'};
    lfpAvg.event_group = event_group;

    % Filtered at broad band (0.5Hz - 300Hz)
    [ csd.all, lfpAvg.all ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
        'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[0.2 0.8],'filter',[0.5 300]);

    % Filtered at slow wave oscilation band (0.5Hz - 4Hz)
    [ csd.SO, lfpAvg.SO ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
        'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[0.2 0.8],'filter',[0.5 4]);

    plot_perievent_CSD_LFP_amplitude_phase(lfpAvg,csd,power,chan_config,sorted_config,best_channels);
    save visual_scene_LFP lfpAvg csd
end

% Brain events
all_events = {ripples.SWS_peaktimes',slow_waves.ints.UP(:,1)',spindles.SWS_peaktimes,V1_superficial_ripples.SWS_peaktimes',V1_deep_ripples.SWS_peaktimes'};
event_group = {'Ripple','UP','Spindle','V1 Superficial Ripple','V1 Deep Ripple'};

all_events = {ripples.probe(2).onset(41)'};
event_group = {'Ripple'};

lfpAvg = [];
csd = [];
lfpAvg.event_group = event_group;
% lfpAvg.filter_type = {'SO','spindle','gamma','ripple','all'};
lfpAvg.filter_type = {'SO','spindle','gamma','ripple','all'};

% Filtered at slow wave oscilation band (0.5Hz - 4Hz)
[ csd.SO, lfpAvg.SO ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
    'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[0.5 0.5],'filter',[0.5 4]);

% [ csd1, lfp1 ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
%     'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[1 1],'filter',[]);

% Filtered at spindle range (9Hz - 17Hz)
[ csd.spindle, lfpAvg.spindle ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
    'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[0.5 0.5],'filter',[9 17]);

% Filtered at gamma range (30Hz - 100Hz)
[ csd.gamma, lfpAvg.gamma ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
    'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[0.5 0.5],'filter',[30 100]);


% Filtered at ripple range (125Hz - 300Hz)
[ csd.ripple, lfpAvg.ripple ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
    'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[0.5 0.5],'filter',[125 300]);

% Filtered at broad band (0.5Hz - 300Hz)
[ csd.all, lfpAvg.all ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
    'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[0.5 0.5],'filter',[0.5 300]);

save peri_event_LFP lfpAvg csd




all_events = {ripples.probe(2).onset(61)'};
event_group = {'Ripple'};
lfpAvg.filter_type = {'all'}
% Filtered at broad band (0.5Hz - 300Hz)
[ csd.all, lfpAvg.all ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP{2}',tvec',all_events,...
    'channels',1:1:size(raw_LFP{2},1),'samplingRate',SR,'twin',[0.5 0.5],'filter',[]);
plot_perievent_CSD_LFP_amplitude_phase(lfpAvg,csd,power{2},chan_config,sorted_config,best_channels{2});



%% Peri-event spike time histogram
clear all
SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_types_all = {'RUN','POST'};
% Stimulus_type = 'RUN'; % extract LFP during RUN
% Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
c = 1;

time_windows = [1:10;11:20;21:30;...
    31:40;41:50;[51:59 nan]]; % for V1 log odds

for epoch = 1:length(Stimulus_types_all)
    Stimulus_type= Stimulus_types_all{epoch};
    for nsession =1:length(experiment_info)
        tic
        session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

        if isempty(session_info)
            continue
        end

        for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
            cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
            load best_channels
            load extracted_PSD
            column = 1;
            load(sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('decoded_replay_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('decoded_replay_events_global_remapped%s.mat',erase(stimulus_name{n},'Masa2tracks')))

            load(sprintf('probability_ratio_original%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('probability_ratio_global_remapped%s.mat',erase(stimulus_name{n},'Masa2tracks')))

            load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('extracted_CA1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('decoded_replay_events_V1%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('decoded_replay_events_shuffled_V1%s.mat',erase(stimulus_name{n},'Masa2tracks')))

            load('extracted_V1_place_fields.mat')
            load('extracted_CA1_place_fields.mat')
            load('extracted_HPC_place_fields.mat')

            % global remapped for subtract mean log odds bias
            z_log_odds_global_remapped = [];
            V1_z_log_odds = [];
            HPC_bayesian_bias = [];
            V1_bayesian_bias = [];

            for event = 1:length(decoded_replay_events(1).replay_events)
                T1_T2_ratio_shuffled = [];
                % T1/T2 ratio
                for nshuffle = 1:1000
                    T1_T2_ratio_shuffled(nshuffle) = probability_ratio_global_remapped{2}{nshuffle}(1,event);
                end

                % Calculate and save zscored log odd
                data = log(probability_ratio_global_remapped{1}(1,event));
                shuffled_data = log(T1_T2_ratio_shuffled);
                z_log_odds_global_remapped(event) = (data-mean(shuffled_data))/std(shuffled_data);

                T1_T2_ratio_shuffled = [];
                % T1/T2 ratio
                for nshuffle = 1:1000
                    T1_T2_ratio_shuffled(nshuffle) = probability_ratio_original{2}{nshuffle}(1,event);
                end

                % Calculate and save zscored log odd
                data = log(probability_ratio_original{1}(1,event));
                shuffled_data = log(T1_T2_ratio_shuffled);
                z_log_odds(event) = (data-mean(shuffled_data))/std(shuffled_data);

                HPC_bayesian_bias(event) = nansum(sum(decoded_replay_events(1).replay_events(event).decoded_position))/...
                    (nansum(sum(decoded_replay_events(1).replay_events(event).decoded_position)) + nansum(sum(decoded_replay_events(2).replay_events(event).decoded_position)));

                % V1 log odds for each probe
                for nprobe = 1:length(session_info(n).probe)
                    probe_no = session_info(n).probe(nprobe).probe_id + 1;
                    probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere ;
                    for nbin = 1:6
                        time_bin = time_windows(nbin,~isnan(time_windows(nbin,:)));
                        data = log(sum(decoded_replay_events_V1.probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                            /sum(decoded_replay_events_V1.probe(probe_no).track(2).replay_events(event).summed_probability(:,time_bin)));

                        shuffled_data = log(sum(decoded_replay_events_shuffled_V1.probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin),2)...
                            ./sum(decoded_replay_events_shuffled_V1.probe(probe_no).track(2).replay_events(event).summed_probability(:,time_bin),2));


                        V1_bayesian_bias{probe_hemisphere}(event,nbin) = nansum(decoded_replay_events_V1.probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                            /(nansum(decoded_replay_events_V1.probe(probe_no).track(1).replay_events(event).summed_probability(:,time_bin))...
                            + nansum(decoded_replay_events_V1.probe(probe_no).track(2).replay_events(event).summed_probability(:,time_bin)));

                        V1_z_log_odds{probe_hemisphere}(event,nbin) = (data-mean(shuffled_data))/std(shuffled_data);
                    end
                end
            end

            all_spike_data = [];
            for nprobe = 1:length(session_info(n).probe)
                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere ;

                % if left hemisphere find good cells with higher FR on T2
                if probe_hemisphere == 1
                    cell_index = find((V1_place_fields.probe(probe_no).track(1).mean_rate_track - V1_place_fields.probe(probe_no).track(2).mean_rate_track)./...
                        (V1_place_fields.probe(probe_no).track(1).mean_rate_track + V1_place_fields.probe(probe_no).track(2).mean_rate_track) < 0);
                    [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                        cell_index);
                    %                     [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                    %                         V1_place_fields.probe(probe_no).track(2).good_cells_LIBERAL);

                    %                     [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                    %                         find(V1_place_fields.probe(probe_no).track(1).mean_rate_track < 0.1 & V1_place_fields.probe(probe_no).track(2).mean_rate_track > 0.1));
                    %                     [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                    %                         find((V1_place_fields.probe(probe_no).track(2).raw_peak - V1_place_fields.probe(probe_no).track(1).raw_peak)>0));


                    %                                         spike_id = V1_clusters.probe(probe_no).spike_id;
                    %                                         spike_times = V1_clusters.probe(probe_no).spike_times;

                    % if Right hemisphere find good cells with higher FR on T1
                elseif probe_hemisphere == 2
                    %
                    cell_index = find((V1_place_fields.probe(probe_no).track(1).mean_rate_track - V1_place_fields.probe(probe_no).track(2).mean_rate_track)./...
                        (V1_place_fields.probe(probe_no).track(1).mean_rate_track + V1_place_fields.probe(probe_no).track(2).mean_rate_track) > 0);
                    [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                        cell_index);
                    %                     [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                    %                         V1_place_fields.probe(probe_no).track(1).good_cells_LIBERAL);

                    %                     [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                    %                               find((V1_place_fields.probe(probe_no).track(1).mean_rate_track - V1_place_fields.probe(probe_no).track(2).mean_rate_track)>0));
                    %        [spike_times,spike_id]= select_spikes_subset(V1_clusters.probe(probe_no),...
                    %                               find((V1_place_fields.probe(probe_no).track(1).raw_peak - V1_place_fields.probe(probe_no).track(2).raw_peak)>0));


                    %
                    %                                         spike_id = V1_clusters.probe(probe_no).spike_id;
                    %                                         spike_times = V1_clusters.probe(probe_no).spike_times;
                end
                all_spike_data{probe_hemisphere}{1} = [spike_id spike_times];

                [spike_times,spike_id]= select_spikes_subset(CA1_clusters.probe(probe_no),...
                    CA1_place_fields.probe(probe_no).good_place_cells_LIBERAL);
                all_spike_data{probe_hemisphere}{2} = [spike_id spike_times];
                [spike_times,spike_id]= select_spikes_subset(HPC_clusters.probe(probe_no),...
                    HPC_place_fields.probe(probe_no).good_place_cells_LIBERAL);
                all_spike_data{probe_hemisphere}{3} = [spike_id spike_times];
                group_name = {'V1','CA1','HPC'};
            end

            if length(session_info(n).probe) == 1 % if one probe, but it can be probe 2
                probe_no = session_info(n).probe(1).probe_id + 1;
                replay = reactivations.probe(probe_no); % using probe 1 bursting for now
                probe_hemisphere = session_info(n).probe(1).probe_hemisphere;
            else
                probe_hemisphere(1) = session_info(n).probe(1).probe_hemisphere;
                probe_hemisphere(2) = session_info(n).probe(2).probe_hemisphere;

                %                 Use the probe with more events
                if sum(reactivations.probe(1).ripple_peak>=3) > sum(reactivations.probe(2).ripple_peak>=3)
                    replay = reactivations.probe(1);
                else
                    replay = reactivations.probe(2);
                end
            end


            %             z_log_odds = z_log_odds-z_log_odds_global_remapped;

            ripple_thresholded = replay.ripple_peak >= 2;

            if isempty(ripple_thresholded)
                continue
            end

            probe_hemisphere_text = {'Left','Right'}
            for nprobe = 1:length(session_info(n).probe)
                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere;
                options = session_info(n).probe(nprobe);


                track_1_events = replay.onset(HPC_bayesian_bias>prctile(HPC_bayesian_bias,70) & ripple_thresholded == 1);
                track_2_events = replay.onset(HPC_bayesian_bias>prctile(HPC_bayesian_bias,30) & ripple_thresholded == 1);


                PSTH = plot_perievent_spiketime_histogram(all_spike_data{probe_hemisphere},track_1_events,'group','by cell zscore','group_name',group_name,'event_name','T1 reactivation','twin',[-0.6 0.6])
                %                 plot_perievent_spiketime_histogram(all_spike_data{probe_hemisphere},track_1_events,'group','by region','group_name',group_name,'event_name','T1 reactivation','twin',[-1 1])
                sgtitle(sprintf('%s %s %s T1 biased reactivation probe %s',options.SUBJECT,options.SESSION,Stimulus_types_all{epoch},probe_hemisphere_text{probe_hemisphere}))
                set(gcf, 'Name', sprintf('%s %s %s T1 biased reactivation probe %s',options.SUBJECT,options.SESSION,Stimulus_types_all{epoch},probe_hemisphere_text{probe_hemisphere}))
                save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project\Log odds analysis',[])


                PSTH = plot_perievent_spiketime_histogram(all_spike_data{probe_hemisphere},track_2_events,'group','by cell zscore','group_name',group_name,'event_name','T2 reactivation','twin',[-0.6 0.6])
                %                 plot_perievent_spiketime_histogram(all_spike_data{probe_hemisphere},track_2_events,'group','by region','group_name',group_name,'event_name','T2 reactivation','twin',[-1 1])
                sgtitle(sprintf('%s %s %s T2 biased reactivation probe %s',options.SUBJECT,options.SESSION,Stimulus_types_all{epoch},probe_hemisphere_text{probe_hemisphere}))
                set(gcf, 'Name', sprintf('%s %s %s T2 biased reactivation probe %s',options.SUBJECT,options.SESSION,Stimulus_types_all{epoch},probe_hemisphere_text{probe_hemisphere}))
                save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project\Log odds analysis',[])

                for mprobe = 1:length(session_info(n).probe)
                    if session_info(n).probe(mprobe).probe_hemisphere == 1
                        if isempty(reactivations.probe(mprobe).onset(reactivations.probe(mprobe).ripple_peak>=2))
                            continue
                        end
                        PSTH = plot_perievent_spiketime_histogram(all_spike_data{probe_hemisphere},reactivations.probe(mprobe).onset(reactivations.probe(mprobe).ripple_peak>=2),'group','by cell zscore','group_name',group_name,'event_name','T2 reactivation','twin',[-0.6 0.6])
                        %                 plot_perievent_spiketime_histogram(all_spike_data{probe_hemisphere},track_2_events,'group','by region','group_name',group_name,'event_name','T2 reactivation','twin',[-1 1])
                        sgtitle(sprintf('%s %s %s candidate events probe %s',options.SUBJECT,options.SESSION,Stimulus_types_all{epoch},probe_hemisphere_text{probe_hemisphere}))
                        set(gcf, 'Name', sprintf('%s %s %s candidate events (left detected) probe %s',options.SUBJECT,options.SESSION,Stimulus_types_all{epoch},probe_hemisphere_text{probe_hemisphere}))
                        save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project\Log odds analysis',[])

                    else
                        if isempty(reactivations.probe(mprobe).onset(reactivations.probe(mprobe).ripple_peak>=2))
                            continue
                        end
                        PSTH = plot_perievent_spiketime_histogram(all_spike_data{probe_hemisphere},reactivations.probe(mprobe).onset(reactivations.probe(mprobe).ripple_peak>=2),'group','by cell zscore','group_name',group_name,'event_name','T2 reactivation','twin',[-0.6 0.6])
                        %                 plot_perievent_spiketime_histogram(all_spike_data{probe_hemisphere},track_2_events,'group','by region','group_name',group_name,'event_name','T2 reactivation','twin',[-1 1])
                        sgtitle(sprintf('%s %s %s candidate events probe %s',options.SUBJECT,options.SESSION,Stimulus_types_all{epoch},probe_hemisphere_text{probe_hemisphere}))
                        set(gcf, 'Name', sprintf('%s %s %s candidate events (right detected) probe %s',options.SUBJECT,options.SESSION,Stimulus_types_all{epoch},probe_hemisphere_text{probe_hemisphere}))
                        save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project\Log odds analysis',[])
                    end
                end
            end
        end
        toc
    end

    %     save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project\Log odds analysis',[])
end
