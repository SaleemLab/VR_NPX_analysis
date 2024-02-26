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
SUBJECTS = {'M23034','M23087'};
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