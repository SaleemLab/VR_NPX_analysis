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
% SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
SUBJECTS={'M24016','M24017','M24018'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
% experiment_info = experiment_info(4);
% Stimulus_type = 'RUN';
% Stimulus_type = 'SleepChronic';
% all_stimulus_type={'RUN'};
all_stimulus_type={'SleepChronic','RUN'};
for nstimuli = 1:length(all_stimulus_type)
    for nsession = 1:length(experiment_info)
        session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,all_stimulus_type{nstimuli}));
        stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,all_stimulus_type{nstimuli}));
        if isempty(stimulus_name)
            continue
        end

        load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

        for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
            options = session_info(n).probe(1);
            % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
            if contains(stimulus_name{n},'Masa2tracks')
                load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            else
                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
            end
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

            if contains(stimulus_name{n},'Masa2tracks')
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PSD')
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP','-v7.3')
            else
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'),'PSD','power')% save PSD for the sleep session
                % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP','-v7.3')
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP','-v7.3')
            end
        end
    end
end

%% Detect UP/DOWN and spindles ripples and reactivation event

clear all
% SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
SUBJECTS={'M24016','M24017','M24018'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
Stimulus_type = 'Sleep';
experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
% 1:length(experiment_info)
% [1 2 3 4 6 7 8 9 10 12 14]

all_stimulus_type={'SleepChronic','RUN'};
for nstimuli = 1:length(all_stimulus_type)
    for nsession = 1:length(experiment_info)
        session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,all_stimulus_type{nstimuli}));
        stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,all_stimulus_type{nstimuli}));

        SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT},option);
        % find right date number based on all experiment dates of the subject
        iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));
        if isempty(stimulus_name)
            continue
        end
        load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

        for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
            %         if ~contains(stimulus_name{n},'Sleep')
            %             continue
            %         end

            options = session_info(n).probe(1);

            DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));
            if isempty(DIR)
                continue
            end

            clear clusters
            if contains(stimulus_name{n},'Masa2tracks')
                % load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%                 save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP');

                load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP');
                load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
                load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
                DIR = dir(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3*')));
                if isempty(DIR)
                    continue
                end
                load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
                clusters=clusters_ks3;
            elseif contains(stimulus_name{n},'Sleep')
                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
                %             load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'));
                clusters=clusters_ks3;
            else
                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'));
                clusters=clusters_ks3;
            end
            clear CA1_clusters V1_clusters

            DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
            DIR1 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
            DIR2 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN2.mat'));

            DIR3 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters.mat'));

            session_clusters_RUN=[];
            session_clusters_RUN1=[];
            session_clusters_RUN2=[];

            if ~isempty(DIR)
                load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
                session_clusters_RUN=session_clusters;
            end

            if ~isempty(DIR1)
                load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
                session_clusters_RUN1=session_clusters;
                session_clusters_RUN=session_clusters_RUN1;
            end

            if ~isempty(DIR2)
                load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN2.mat'));
                session_clusters_RUN2=session_clusters;
            end

            session_clusters=[];
            if ~isempty(DIR3)
                load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters.mat'));
            end

            params = create_cluster_selection_params('sorting_option','masa');
            clear selected_clusters
            for nprobe = 1:length(clusters)
                selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters
            end

            for nprobe = 1:length(clusters)
                % Convert to unique spike/cluster id
                selected_clusters(nprobe).spike_id = selected_clusters(nprobe).spike_id + nprobe*10^4 + iDate* 10^6 + str2double(options.SUBJECT(2:end))*10^8;
                selected_clusters(nprobe).cluster_id = selected_clusters(nprobe).cluster_id + nprobe*10^4 + iDate* 10^6 + str2double(options.SUBJECT(2:end))*10^8;
                if session_info(n).probe(nprobe).probe_hemisphere==1
                    selected_clusters(nprobe).region = session_clusters_RUN.region(contains(session_clusters_RUN.region,'L'));
                    selected_clusters(nprobe).spatial_response = session_clusters_RUN.spatial_response(contains(session_clusters_RUN.region,'L'),:);
                    selected_clusters(nprobe).odd_even_stability = session_clusters_RUN.odd_even_stability(contains(session_clusters_RUN.region,'L'),:);
                    selected_clusters(nprobe).peak_percentile = session_clusters_RUN.peak_percentile(contains(session_clusters_RUN.region,'L'),:);
                elseif session_info(n).probe(nprobe).probe_hemisphere==2
                    selected_clusters(nprobe).region = session_clusters_RUN.region(contains(session_clusters_RUN.region,'R'));
                    selected_clusters(nprobe).spatial_response = session_clusters_RUN.spatial_response(contains(session_clusters_RUN.region,'R'),:);
                    selected_clusters(nprobe).odd_even_stability = session_clusters_RUN.odd_even_stability(contains(session_clusters_RUN.region,'R'),:);
                    selected_clusters(nprobe).peak_percentile = session_clusters_RUN.peak_percentile(contains(session_clusters_RUN.region,'R'),:);
                end
            end

            % spatial_cell_index = find((session_clusters_RUN.peak_percentile(:,1)>0.95&session_clusters_RUN.odd_even_stability(:,1)>0.95) ...
            %     | (session_clusters_RUN.peak_percentile(:,2)>0.95&session_clusters_RUN.odd_even_stability(:,2)>0.95));


            % As long as stable firing durng track running. It doesn't have to
            % be single-peaked place cell like cells
            spatial_cell_index = find(session_clusters_RUN.odd_even_stability(:,1)>0.95 ...
                | session_clusters_RUN.odd_even_stability(:,2)>0.95);

            clear replay reactivations ripples spindles slow_waves raw_LFP CA1_clusters V1_clusters
            clear V1_replay V1_reactivations replay_combined replay_combined behavioural_state
            for nprobe = 1:length(session_info(n).probe)
                clear cortex_LFP CA1_LFP
                options = session_info(n).probe(nprobe);
                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
                %                 Behavioural state detection


                %             mobility = movmean(Behaviour.mobility,120);

                if contains(stimulus_name{n},'Sleep')
                    mobility_thresholded = abs([0 diff(movmean(Behaviour.mobility,1/mean(diff(Behaviour.tvec))))])>2000;%1 second movemean windows

                    mob_index = find(mobility_thresholded==1);
                    for nindex = 1:length(mob_index)
                        if nindex<length(mob_index)
                            if Behaviour.tvec(mob_index(nindex+1))-Behaviour.tvec(mob_index(nindex))<1 % if immobility less than 1 second, treats that as movement period
                                mobility_thresholded(mob_index(nindex):mob_index(nindex+1))=1;
                            end
                        end
                    end
                    %             plot(Behaviour.tvec,Behaviour.mobility);hold on;
                    %             plot(Behaviour.tvec,mobility*500000);plot(Behaviour.tvec,abs([0 diff(movmean(Behaviour.mobility,1/mean(diff(Behaviour.tvec))))]))
                    %             plot(Behaviour.mobility_zscore)
                    mobility_thresholded = interp1(Behaviour.tvec,double(mobility_thresholded),LFP(probe_no).tvec,'linear');
                    Behaviour.mobility_zscore = abs([0 diff(movmean(Behaviour.mobility,1/mean(diff(Behaviour.tvec))))]);% Diff of pixel change
                    Behaviour.mobility_zscore(isnan(Behaviour.mobility_zscore))=mean(Behaviour.mobility_zscore,'omitnan');
                    Behaviour.mobility_zscore=zscore(Behaviour.mobility_zscore);

                    speed = mobility_thresholded;

                else
                    speed = Behaviour.speed;
                    speed(isnan(speed))=0;
                    w = gausswin(9);
                    w = w / sum(w);
                    speed = filtfilt(w,1,speed')';

                    speed = interp1(Behaviour.tvec,speed,LFP(probe_no).tvec,'linear');
                end

                if isfield(LFP(probe_no),'L5')
                    if ~isempty(LFP(probe_no).L5)
                        bad_channels=[];
                        all_shanks = 1:size(LFP(probe_no).L5_power,1);
                        for nshank = 1:size(LFP(probe_no).L5_power,1)
                            bad_channels(nshank,:) = LFP(probe_no).L5_power(nshank,:)>3*mean(LFP(probe_no).L5_power(all_shanks~=nshank,:));
                        end
                        bad_channels = sum(bad_channels,2)>4;% ignore channels if the mean power is 3 x the rest of channels (possibly will be nosiy)
                        [~,best_channel]=max(LFP(probe_no).L5_power(~bad_channels,7));
                        good_channels = find(~bad_channels);
                        cortex_LFP = LFP(probe_no).L5(good_channels(best_channel),:);
                    else
                        bad_channels=[];
                        all_shanks = 1:size(LFP(probe_no).L4_power,1);
                        for nshank = 1:size(LFP(probe_no).L4_power,1)
                            bad_channels(nshank,:) = LFP(probe_no).L4_power(nshank,:)>3*mean(LFP(probe_no).L4_power(all_shanks~=nshank,:));
                        end
                        bad_channels = sum(bad_channels,2)>4;
                        [~,best_channel]=max(LFP(probe_no).L4_power(~bad_channels,7));
                        good_channels = find(~bad_channels);
                        cortex_LFP = LFP(probe_no).L4(good_channels(best_channel),:);
                    end

                elseif isfield(LFP(probe_no),'L4')
                    if ~isempty(LFP(probe_no).L4)
                        bad_channels=[];
                        all_shanks = 1:size(LFP(probe_no).L4_power,1);
                        for nshank = 1:size(LFP(probe_no).L4_power,1)
                            bad_channels(nshank,:) = LFP(probe_no).L4_power(nshank,:)>3*mean(LFP(probe_no).L4_power(all_shanks~=nshank,:));
                        end
                        bad_channels = sum(bad_channels,2)>4;
                        [~,best_channel]=max(LFP(probe_no).L4_power(~bad_channels,7));
                        good_channels = find(~bad_channels);
                        cortex_LFP = LFP(probe_no).L4(good_channels(best_channel),:);
                    else
                        cortex_LFP = [];
                        disp('cortex LFP is missing')
                    end

                elseif isfield(LFP(probe_no),'MEC')

                end

                if isfield(LFP(probe_no),'CA1')
                    bad_channels=[];
                    all_shanks = 1:size(LFP(probe_no).CA1_power,1);
                    for nshank = 1:size(LFP(probe_no).CA1_power,1)
                        bad_channels(nshank,:) = LFP(probe_no).CA1_power(nshank,:)>3*mean(LFP(probe_no).CA1_power(all_shanks~=nshank,:));
                    end
                    bad_channels = sum(bad_channels,2)>4;

                    [~,best_channel]=max(LFP(probe_no).CA1_power(~bad_channels,6));
                    good_channels = find(~bad_channels);

                    CA1_LFP = LFP(probe_no).CA1(good_channels(best_channel),:);

                    speedTreshold = 1;
                    %                 CA1_LFP = raw_LFP{nprobe}(229,:);
                    [freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
                        [LFP(probe_no).tvec' cortex_LFP'],[LFP(probe_no).tvec' CA1_LFP'],...
                        [LFP(probe_no).tvec' speed'],speedTreshold);
                end
                %             behavioural_state.freezing = freezing;
                behavioural_state(probe_no).quietWake = quietWake;
                behavioural_state(probe_no).SWS = SWS;
                behavioural_state(probe_no).REM = REM;
                behavioural_state(probe_no).movement = movement;

                %%%%%%%%%%%%%%%%%%
                % UP/Down states and ripple and candidate reactivation events detection
                %%%%%%%%%%%%%%%%%%

                zscore_min = 0;
                zscore_max = 3;
                metric_param =[];
%                 metric_param.cluster_id = @(x) ismember(x,session_clusters_RUN.cluster_id(spatial_cell_index));

                if options.probe_hemisphere==1
                    metric_param.region = @(x) contains(x,'HPC_L');
                    CA1_clusters(probe_no) = select_clusters(selected_clusters(nprobe),metric_param);
                elseif options.probe_hemisphere==2
                    metric_param.region = @(x) contains(x,'HPC_R');
                    CA1_clusters(probe_no) = select_clusters(selected_clusters(nprobe),metric_param);
                end

                if ~isempty(CA1_clusters(probe_no).cluster_id)
                    if length(CA1_clusters(probe_no).cluster_id)>10
                        [replay(probe_no),reactivations(probe_no)] = detect_candidate_events_masa(LFP(probe_no).tvec,CA1_LFP,...
                            [CA1_clusters(probe_no).spike_id CA1_clusters(probe_no).spike_times],Behaviour,zscore_min,zscore_max,options);
                    end
                end

                % Detect V1 populational bursting events (Candidate events)
                metric_param =[];
                %             metric_param.cluster_id = @(x) ismember(x,session_clusters_RUN.cluster_id(spatial_cell_index));

                if options.probe_hemisphere==1
                    metric_param.region = @(x) contains(x,'V1_L');
                    V1_clusters(probe_no) = select_clusters(selected_clusters(nprobe),metric_param);
                elseif options.probe_hemisphere==2
                    metric_param.region = @(x) contains(x,'V1_R');
                    V1_clusters(probe_no) = select_clusters(selected_clusters(nprobe),metric_param);
                end

                %%%%%%% Slow wave detections

                if isfield(LFP(nprobe),'L5')==1 & ~isempty(behavioural_state(probe_no).SWS)
                    slow_waves(nprobe) = DetectSlowWaves_masa('time',LFP(nprobe).tvec,'lfp',cortex_LFP,'spikes',V1_clusters(nprobe),'NREMInts',behavioural_state(probe_no).SWS);
                elseif isfield(LFP(nprobe),'L4')==1 & ~isempty(behavioural_state(probe_no).SWS)
                    slow_waves(nprobe) = DetectSlowWaves_masa('time',LFP(nprobe).tvec,'lfp',cortex_LFP,'spikes',V1_clusters(nprobe),'NREMInts',behavioural_state(probe_no).SWS);
                else
                    if nprobe == length(session_info(n).probe)
                        if exist('slow_waves')==0
                            slow_waves(nprobe) = struct();
                        end
                    end
                end

                %%%%%% Populational burtsting events
                metric_param =[];
%                 metric_param.cluster_id = @(x) ismember(x,session_clusters_RUN.cluster_id(spatial_cell_index));

                if options.probe_hemisphere==1
                    metric_param.region = @(x) contains(x,'V1_L');
                    V1_clusters(probe_no) = select_clusters(selected_clusters(nprobe),metric_param);
                elseif options.probe_hemisphere==2
                    metric_param.region = @(x) contains(x,'V1_R');
                    V1_clusters(probe_no) = select_clusters(selected_clusters(nprobe),metric_param);
                end

                if ~isempty(V1_clusters(probe_no).cluster_id)
                    if length(V1_clusters(probe_no).cluster_id)>10
                        [V1_replay(probe_no),V1_reactivations(probe_no)] = detect_candidate_events_masa(LFP(probe_no).tvec,CA1_LFP,...
                            [V1_clusters(probe_no).spike_id V1_clusters(probe_no).spike_times],Behaviour,zscore_min,zscore_max,options);
                    end
                end

                % Detect V1 spindle events
                [spindles(probe_no)] = FindSpindles_masa(cortex_LFP,LFP(probe_no).tvec','behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
                    'noise',[],'passband',[9 17],'thresholds',[1 3],'show','off');

                % Detect CA1 ripple events
                [ripples(probe_no)] = FindRipples_masa(CA1_LFP,LFP(probe_no).tvec','behaviour',Behaviour,'minDuration',30,'durations',[30 200],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
                    'noise',[],'passband',[125 300],'thresholds',[2 5],'show','off');
            end

            for nprobe = 1:length(behavioural_state)
                if contains(stimulus_name{n},'Sleep')
                    if ~isempty(CA1_clusters(probe_no).cluster_id)
                        if length(CA1_clusters(probe_no).cluster_id)>10
                            [reactivations(nprobe).awake_offset,reactivations(nprobe).awake_index] = RestrictInts(reactivations(nprobe).offset',behavioural_state(nprobe).quietWake);
                            reactivations(nprobe).awake_onset = reactivations(nprobe).onset(reactivations(nprobe).awake_index)';
                        end
                    end

                    if ~isempty(V1_reactivations(nprobe).onset)
                        [V1_reactivations(nprobe).awake_offset,V1_reactivations(nprobe).awake_index] = RestrictInts(V1_reactivations(nprobe).offset',behavioural_state(nprobe).quietWake);
                        V1_reactivations(nprobe).awake_onset = V1_reactivations(nprobe).onset(V1_reactivations(nprobe).awake_index)';
                    end

                    [ripples(nprobe).awake_offset,ripples(nprobe).awake_index] = RestrictInts(ripples(nprobe).offset,behavioural_state(nprobe).quietWake);
                    ripples(nprobe).awake_onset = ripples(nprobe).onset(ripples(nprobe).awake_index);
                    ripples(nprobe).awake_peaktimes = ripples(nprobe).peaktimes(ripples(nprobe).awake_index);

                    [spindles(nprobe).awake_offset,spindles(nprobe).awake_index] = RestrictInts(spindles(nprobe).offset,behavioural_state(nprobe).quietWake);
                    spindles(nprobe).awake_onset = spindles(nprobe).onset(spindles(nprobe).awake_index);
                    spindles(nprobe).awake_peaktimes = spindles(nprobe).peaktimes(spindles(nprobe).awake_index);
                    %                 reactivations(nprobe).awake_peaktimes = reactivations(nprobe).peaktimes(reactivations.awake_index);

                    if ~isempty(behavioural_state(nprobe).SWS)
                        if length(CA1_clusters(probe_no).cluster_id)>10
                            [reactivations(nprobe).SWS_offset,reactivations(nprobe).SWS_index] = RestrictInts(reactivations(nprobe).offset',behavioural_state(nprobe).SWS);
                            reactivations(nprobe).SWS_onset = reactivations(nprobe).onset(reactivations(nprobe).SWS_index)';
                        end

                        [ripples(nprobe).SWS_offset,ripples(nprobe).SWS_index] = RestrictInts(ripples(nprobe).offset,behavioural_state(nprobe).SWS);
                        ripples(nprobe).SWS_onset = ripples(nprobe).onset(ripples(nprobe).SWS_index);
                        ripples(nprobe).SWS_peaktimes = ripples(nprobe).peaktimes(ripples(nprobe).SWS_index);

                        [spindles(nprobe).SWS_offset,spindles(nprobe).SWS_index] = RestrictInts(spindles(nprobe).offset,behavioural_state(nprobe).SWS);
                        spindles(nprobe).SWS_onset = spindles(nprobe).onset(spindles(nprobe).SWS_index);
                        spindles(nprobe).SWS_peaktimes = spindles(nprobe).peaktimes(spindles(nprobe).SWS_index);
                    end
                else

                end

            end


            %%%%%%%%%%%%%%%%%%
            % Candidate reactivation events detection (probe combined)
            %%%%%%%%%%%%%%%%%%
            reactivations_combined= [];
            replay_combined = [];
            clear CA1_clusters_combined

            if length(session_info(n).probe)>1
                zscore_min = 0;
                zscore_max = 3;

                metric_param =[];
                metric_param.cluster_id = @(x) ismember(x,session_clusters_RUN.cluster_id(spatial_cell_index));
                metric_param.region = @(x) contains(x,'HPC');
                clusters_combined = combine_clusters_from_multiple_probes(selected_clusters(1),selected_clusters(2))
                CA1_clusters_combined = select_clusters(clusters_combined,metric_param);

                [replay_combined,reactivations_combined] = detect_candidate_events_masa(LFP(nprobe).tvec,CA1_LFP,...
                    [CA1_clusters_combined.spike_id CA1_clusters_combined.spike_times],Behaviour,zscore_min,zscore_max,options);
            end

            for nprobe = 1:length(session_info(n).probe)
                probe_no = session_info(n).probe(nprobe).probe_id + 1;
                options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

                if isfield(Behaviour,'speed')
                    for event = 1:length(ripples(probe_no).onset)
                        ripples(probe_no).speed(event) = mean(Behaviour.speed(find(Behaviour.sglxTime >= ripples(probe_no).onset(event) & Behaviour.sglxTime <= ripples(probe_no).offset(event))));
                    end
                    if ~isempty(spindles(probe_no).onset)
                        for event = 1:length(spindles(probe_no).onset)
                            spindles(probe_no).speed(event) = mean(Behaviour.speed(find(Behaviour.sglxTime >= ripples(probe_no).onset(event) & Behaviour.sglxTime <= ripples(probe_no).offset(event))));
                        end
                    end
                else
                    for event = 1:length(ripples(probe_no).onset)
                        ripples(probe_no).mobility_zscore(event) = mean(Behaviour.mobility_zscore(find(Behaviour.sglxTime >= ripples(probe_no).onset(event) & Behaviour.sglxTime <= ripples(probe_no).offset(event))));
                    end

                    if ~isempty(spindles(probe_no).onset)
                        for event = 1:length(spindles(probe_no).onset)
                            spindles(probe_no).mobility_zscore(event) = mean(Behaviour.mobility_zscore(find(Behaviour.sglxTime >= ripples(probe_no).onset(event) & Behaviour.sglxTime <= ripples(probe_no).offset(event))));
                        end
                    end
                end

                if ~contains(stimulus_name{n},'RUN') % If reactivation events during lap running
                    continue
                end

                lap_times(1).start = session_clusters_RUN.start_time_all{1}(session_clusters_RUN.track_ID_all{1}==1);
                lap_times(1).end = session_clusters_RUN.end_time_all{1}(session_clusters_RUN.track_ID_all{1}==1)';

                lap_times(2).start = session_clusters_RUN.start_time_all{1}(session_clusters_RUN.track_ID_all{1}==2);
                lap_times(2).end = session_clusters_RUN.end_time_all{1}(session_clusters_RUN.track_ID_all{1}==2)';

                if contains(stimulus_name{n},'RUN') % If reactivation events during lap running
                    if length(CA1_clusters(probe_no).cluster_id)>10
                        if ~isempty(reactivations(probe_no).onset)
                            [reactivations(probe_no).T1_offset,reactivations(probe_no).T1_index] = RestrictInts(reactivations(probe_no).offset',[lap_times(1).start lap_times(1).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                            reactivations(probe_no).T1_onset = reactivations(probe_no).onset(reactivations(probe_no).T1_index);
                            reactivations(probe_no).T1_midpoint = reactivations(probe_no).midpoint(reactivations(probe_no).T1_index);

                            [reactivations(probe_no).T2_offset,reactivations(probe_no).T2_index] = RestrictInts(reactivations(probe_no).offset',[lap_times(2).start lap_times(2).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                            reactivations(probe_no).T2_onset = reactivations(probe_no).onset(reactivations(probe_no).T2_index);
                            reactivations(probe_no).T2_midpoint = reactivations(probe_no).midpoint(reactivations(probe_no).T2_index);
                        end
                    end

                    if length(V1_clusters(probe_no).cluster_id)>10
                        if ~isempty(V1_reactivations(probe_no).onset)
                            [V1_reactivations(probe_no).T1_offset,V1_reactivations(probe_no).T1_index] = RestrictInts(V1_reactivations(probe_no).offset',[lap_times(1).start lap_times(1).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                            V1_reactivations(probe_no).T1_onset = V1_reactivations(probe_no).onset(V1_reactivations(probe_no).T1_index);
                            V1_reactivations(probe_no).T1_midpoint = V1_reactivations(probe_no).midpoint(V1_reactivations(probe_no).T1_index);

                            [V1_reactivations(probe_no).T2_offset,V1_reactivations(probe_no).T2_index] = RestrictInts(V1_reactivations(probe_no).offset',[lap_times(2).start lap_times(2).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
                            V1_reactivations(probe_no).T2_onset = V1_reactivations(probe_no).onset(V1_reactivations(probe_no).T2_index);
                            V1_reactivations(probe_no).T2_midpoint = V1_reactivations(probe_no).midpoint(V1_reactivations(probe_no).T2_index);
                        end
                    end
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
            if contains(stimulus_name{n},'Masa2tracks')
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'replay','reactivations','replay_combined','reactivations_combined')
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_candidate_events_V1%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'V1_replay','V1_reactivations')
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_spindle_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'spindles')
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'ripples')
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_slow_wave_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'slow_waves')
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'),'slow_waves')
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('behavioural_state%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'behavioural_state')
            else
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'),'replay','reactivations','replay_combined','reactivations_combined')
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'),'V1_replay','V1_reactivations')
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'),'ripples')
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'),'spindles')
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'),'slow_waves')
                save(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state.mat'),'behavioural_state')
            end
            close all
        end
    end
end

%% UP DOWN state and ripple and spindle analysis

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

clear all
% SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
SUBJECTS={'M24016','M24017','M24018'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
Stimulus_type = 'Sleep';
experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
% 1:length(experiment_info)
% [1 2 3 4 6 7 8 9 10 12 14]

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT},option);
    % find right date number based on all experiment dates of the subject
    iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        
        DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));
        if isempty(DIR)
            continue
        end

        if contains(stimulus_name{n},'Masa2tracks')
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP');
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            clusters=clusters_ks3;
        elseif contains(stimulus_name{n},'Sleep')
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
%             load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
            
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state.mat'));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
%             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'));
            clusters=clusters_ks3;
        else
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
            
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'));
            clusters=clusters_ks3;
        end
        clear CA1_clusters V1_clusters 
        
        DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
        DIR1 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
        DIR2 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN2.mat'));
        
        DIR3 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters.mat'));
        
        session_clusters_RUN=[];
        session_clusters_RUN1=[];
        session_clusters_RUN2=[];
        
        if ~isempty(DIR)
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
            session_clusters_RUN=session_clusters;
        end

        if ~isempty(DIR1)
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
            session_clusters_RUN1=session_clusters;
            session_clusters_RUN=session_clusters_RUN1;
        end

        if ~isempty(DIR2)
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN2.mat'));
            session_clusters_RUN2=session_clusters;
        end
        
        session_clusters=[];
        if ~isempty(DIR3)
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters.mat'));
        end

        % From cell structure back to spike times and spike id
        session_clusters_RUN.spike_id=vertcat(session_clusters_RUN.spike_id{:});
        session_clusters_RUN.spike_times=vertcat(session_clusters_RUN.spike_times{:});
        [session_clusters_RUN.spike_times,index] =sort(session_clusters_RUN.spike_times);
        session_clusters_RUN.spike_id=session_clusters_RUN.spike_id(index);

        params = create_cluster_selection_params('sorting_option','masa');
        clear selected_clusters
        for nprobe = 1:length(clusters)
            selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters
        end

        for nprobe = 1:length(clusters)
            % Convert to unique spike/cluster id
            selected_clusters(nprobe).spike_id = selected_clusters(nprobe).spike_id + nprobe*10^4 + iDate* 10^6 + str2double(options.SUBJECT(2:end))*10^8;
            selected_clusters(nprobe).cluster_id = selected_clusters(nprobe).cluster_id + nprobe*10^4 + iDate* 10^6 + str2double(options.SUBJECT(2:end))*10^8;
            if session_info(n).probe(nprobe).probe_hemisphere==1
                selected_clusters(nprobe).region = session_clusters_RUN.region(contains(session_clusters_RUN.region,'L'));
                selected_clusters(nprobe).spatial_response = session_clusters_RUN.spatial_response(contains(session_clusters_RUN.region,'L'),:);
                selected_clusters(nprobe).odd_even_stability = session_clusters_RUN.odd_even_stability(contains(session_clusters_RUN.region,'L'),:);
                selected_clusters(nprobe).peak_percentile = session_clusters_RUN.peak_percentile(contains(session_clusters_RUN.region,'L'),:);
            elseif session_info(n).probe(nprobe).probe_hemisphere==2
                selected_clusters(nprobe).region = session_clusters_RUN.region(contains(session_clusters_RUN.region,'R'));
                selected_clusters(nprobe).spatial_response = session_clusters_RUN.spatial_response(contains(session_clusters_RUN.region,'R'),:);
                selected_clusters(nprobe).odd_even_stability = session_clusters_RUN.odd_even_stability(contains(session_clusters_RUN.region,'R'),:);
                selected_clusters(nprobe).peak_percentile = session_clusters_RUN.peak_percentile(contains(session_clusters_RUN.region,'R'),:);
            end
        end

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(selected_clusters(1),selected_clusters(2));
        else
            clusters_combined = selected_clusters;
        end
        clear selected_clusters
%         spatial_cell_index = find((session_clusters_RUN.peak_percentile(:,1)>0.95&session_clusters_RUN.odd_even_stability(:,1)>0.95) ...
%             | (session_clusters_RUN.peak_percentile(:,2)>0.95&session_clusters_RUN.odd_even_stability(:,2)>0.95));
        spatial_cell_index = find(session_clusters_RUN.odd_even_stability(:,1)>0.95 ...
            | session_clusters_RUN.odd_even_stability(:,2)>0.95);

        metric_param =[];
        metric_param.cluster_id = @(x) ismember(x,session_clusters_RUN.cluster_id(spatial_cell_index));
        [selected_clusters,cluster_id] = select_clusters(session_clusters_RUN,metric_param);
        
        x_window = [0 140];
        x_bin_width = 2;
        place_fields = calculate_spatial_cells(selected_clusters,selected_clusters.tvec{1},...
            selected_clusters.position{1},selected_clusters.speed{1},selected_clusters.track_ID_all{1},selected_clusters.start_time_all{1},selected_clusters.end_time_all{1},x_window,x_bin_width);

        
        for nprobe = 1:2
            if ~isempty(behavioural_state(nprobe).SWS)
                [V1_reactivations(nprobe).SWS_offset,V1_reactivations(nprobe).SWS_index] = RestrictInts(V1_reactivations(nprobe).offset',behavioural_state(nprobe).SWS);
                V1_reactivations(nprobe).SWS_onset = V1_reactivations(nprobe).onset(V1_reactivations(nprobe).SWS_index)';
            end
        end

        [reactivations_combined.SWS_offset,reactivations_combined.SWS_index] = RestrictInts(reactivations_combined.offset',behavioural_state(1).SWS);
        reactivations_combined.SWS_onset = reactivations_combined.onset(reactivations_combined.SWS_index)';


        figure
        nexttile
%         plot_perievent_event_histogram(ripples(2).SWS_peaktimes,ripples(1).SWS_peaktimes,'twin',[-1 1],'event_name','Right ripple')
        time_wondows = [-2 2];
        time_bin = 0.02;
        probabilities = calculate_event_probability(ripples(1).SWS_peaktimes,ripples(2).SWS_peaktimes, [time_wondows(1):time_bin:time_wondows(2)])
        plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,3))
        title('Probability of Left ripple during Right HPC')
        ylabel('probability')
        xlabel('Time(s)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
            
        nexttile
%         plot_perievent_event_histogram(slow_waves(2).ints.UP(:,1),slow_waves(1).ints.UP(:,1),'twin',[-1 1],'event_name','Right UP events')
        time_wondows = [-2 2];
        time_bin = 0.02;
        probabilities = calculate_event_probability(slow_waves(1).ints.UP(:,1),slow_waves(2).ints.UP(:,1), [time_wondows(1):time_bin:time_wondows(2)])
        plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,3))
        title('Probability of UP state from Left V1 during Right V1 Upstates')
        title('Up state from Left V1 relative to Up state from Right V1')
        ylabel('probability')
        xlabel('Time(s)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        nexttile
%         plot_perievent_event_histogram(spindles(2).SWS_peaktimes,spindles(1).SWS_peaktimes,'twin',[-1 1],'event_name','Right spindles events')
        time_wondows = [-3 3];
        time_bin = 0.05;
        probabilities = calculate_event_probability(spindles(1).SWS_peaktimes,spindles(2).SWS_peaktimes, [time_wondows(1):time_bin:time_wondows(2)])
        plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,3))
        title('Probability of left spindles during right spindles')
        title('Spindles from Left V1 relative to Spindles from Right V1')
        ylabel('probability')
        xlabel('Time(s)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        nexttile
%         plot_perievent_event_histogram(V1_reactivations(2).SWS_onset,V1_reactivations(1).SWS_onset,'twin',[-1 1],'event_name','Right V1 populational burtsting events')
        time_wondows = [-2 2];
        time_bin = 0.02;
        probabilities = calculate_event_probability(V1_reactivations(1).SWS_onset,V1_reactivations(2).SWS_onset, [time_wondows(1):time_bin:time_wondows(2)])
        plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,3))
        title('Probability of left V1 bursting during right V1 bursting')
        ylabel('probability')
        xlabel('Time(s)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        nexttile
%         plot_perievent_event_histogram(reactivations(2).SWS_onset,reactivations(1).SWS_onset,'twin',[-1 1],'event_name','Right HPC populational burtsting events')
        time_wondows = [-2 2];
        time_bin = 0.02;
        probabilities = calculate_event_probability(reactivations(2).SWS_onset,reactivations(1).SWS_onset, [time_wondows(1):time_bin:time_wondows(2)])
        plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,3))
        title('Probability of left HPC bursting during right HPC bursting')
        ylabel('probability')
        xlabel('Time(s)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        sgtitle('Left and Right events Synchronisation')

        
        figure
        nexttile
        time_wondows = [-1 1];
        time_bin = 0.02;
        probabilities = calculate_event_probability(ripples(1).SWS_peaktimes, slow_waves(1).ints.UP(:,1), [time_wondows(1):time_bin:time_wondows(2)])
        plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
        title('probability of Left ripples during UP states')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        nexttile
        time_wondows = [-1 1];
        time_bin = 0.02;
        probabilities = calculate_event_probability(ripples(1).SWS_peaktimes, slow_waves(1).ints.DOWN(:,1), [time_wondows(1):time_bin:time_wondows(2)])
        plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
        title('probability of Left ripples during DOWN states')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        nexttile
        time_wondows = [-1 1];
        time_bin = 0.02;
        probabilities = calculate_event_probability(ripples(2).SWS_peaktimes, slow_waves(1).ints.UP(:,1), [time_wondows(1):time_bin:time_wondows(2)])
        plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
        title('probability of Right ripples during UP states')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        nexttile
        time_wondows = [-1 1];
        time_bin = 0.02;
        probabilities = calculate_event_probability(ripples(2).SWS_peaktimes, slow_waves(1).ints.DOWN(:,1), [time_wondows(1):time_bin:time_wondows(2)])
        plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
        title('probability of Right ripples relative to DOWN states')
        sgtitle('Left V1')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        figure
        nexttile
        time_wondows = [-1 1];
        time_bin = 0.02;
        probabilities = calculate_event_probability(ripples(1).SWS_peaktimes, slow_waves(2).ints.UP(:,1), [time_wondows(1):time_bin:time_wondows(2)])
        plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
        title('probability of Left ripples relative to UP states')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        
        nexttile
        time_wondows = [-1 1];
        time_bin = 0.02;
        probabilities = calculate_event_probability(ripples(1).SWS_peaktimes, slow_waves(2).ints.DOWN(:,1), [time_wondows(1):time_bin:time_wondows(2)])
        plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
        title('probability of Left ripples relative to DOWN states')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        nexttile
        time_wondows = [-1 1];
        time_bin = 0.02;
        probabilities = calculate_event_probability(ripples(2).SWS_peaktimes, slow_waves(2).ints.UP(:,1), [time_wondows(1):time_bin:time_wondows(2)])
        plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
        title('probability of Right ripples relative to UP states')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        nexttile
        time_wondows = [-1 1];
        time_bin = 0.02;
        probabilities = calculate_event_probability(ripples(2).SWS_peaktimes, slow_waves(2).ints.DOWN(:,1), [time_wondows(1):time_bin:time_wondows(2)])
        plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
        title('probability of Right ripples relative to DOWN states')
        sgtitle('Right V1')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        sgtitle('Left V1 and HPC interaction')



        % Spike relative to events

        spatial_cell_index = find(clusters_combined.odd_even_stability(:,1)>0.95 ...
            | clusters_combined.odd_even_stability(:,2)>0.95);

        metric_param =[];
        metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
        metric_param.region = @(x) contains(x,'V1_L');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        all_spikes{1}=[selected_clusters.spike_id selected_clusters.spike_times];

        metric_param.region = @(x) contains(x,'V1_R');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        all_spikes{2}=[selected_clusters.spike_id selected_clusters.spike_times];

        metric_param.region = @(x) contains(x,'HPC');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        all_spikes{3}=[selected_clusters.spike_id selected_clusters.spike_times];
        group_name = {'V1_L','V1_R','HPC'};

        % UP
        plot_perievent_spiketime_histogram(all_spikes,slow_waves(1).ints.UP(:,1),'group','by cell zscore','group_name',group_name,'event_name','left UP onset','twin',[-1 1])
        plot_perievent_spiketime_histogram(all_spikes,slow_waves(1).ints.UP(:,1),'group','by region','group_name',group_name,'event_name','left UP onset','twin',[-1 1])

        plot_perievent_spiketime_histogram(all_spikes,slow_waves(2).ints.UP(:,1),'group','by cell zscore','group_name',group_name,'event_name','right UP onset','twin',[-1 1])
        plot_perievent_spiketime_histogram(all_spikes,slow_waves(2).ints.UP(:,1),'group','by region','group_name',group_name,'event_name','right UP onset','twin',[-1 1])

        % Ripples
        plot_perievent_spiketime_histogram(all_spikes,ripples(1).SWS_peaktimes,'group','by cell zscore','group_name',group_name,'event_name','Left Ripples','twin',[-1 1])
        plot_perievent_spiketime_histogram(all_spikes,ripples(1).SWS_peaktimes,'group','by region','group_name',group_name,'event_name','Left Ripples','twin',[-1 1])
        plot_perievent_spiketime_histogram(all_spikes,ripples(2).SWS_peaktimes,'group','by cell zscore','group_name',group_name,'event_name','Right Ripples','twin',[-1 1])
        plot_perievent_spiketime_histogram(all_spikes,ripples(2).SWS_peaktimes,'group','by region','group_name',group_name,'event_name','Right Ripples','twin',[-1 1])
        

        plot_perievent_spiketime_histogram(all_spikes,ripples(1).awake_peaktimes,'group','by cell zscore','group_name',group_name,'event_name','Left Ripples','twin',[-1 1])
        plot_perievent_spiketime_histogram(all_spikes,ripples(1).awake_peaktimes,'group','by region','group_name',group_name,'event_name','Left Ripples','twin',[-1 1])

        % Ripples
        metric_param =[];
        metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
        metric_param.region = @(x) contains(x,'V1');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);

        ia = find(cluster_id);
        C = clusters_combined.cluster_id(ia)
        event_id = [ones(1,length(ripples(1).SWS_peaktimes))];
        event_times = [ripples(1).SWS_peaktimes];
        plot_perievent_spiketimes(clusters_combined.spike_times,clusters_combined.spike_id,[],[],[5 1],[-1 1],0.02,...
            'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',event_times,...
            'event_id',event_id,'event_label','ripple','place_fields',place_fields,'plot_option','by time');

        % 
        
        V1_reactivations(nprobe).SWS_onset = V1_reactivations(nprobe).onset(V1_reactivations(nprobe).SWS_index)';

        zero_meaned_log_odds = zscore([decoded_ripple_events(1).track(1).replay_events(:).z_log_odds]);
        
        T2_events = ripples(2).peaktimes(find(zero_meaned_log_odds<-1));
        [T2_events,~] = RestrictInts(T2_events,behavioural_state(1).SWS);

        T1_events = ripples(2).peaktimes(find(zero_meaned_log_odds>1));
        [T1_events,~] = RestrictInts(T1_events,behavioural_state(1).SWS);
%         T1_events = ripples(1).peaktimes(find(zero_meaned_log_odds>0.5));
        plot_perievent_spiketime_histogram(all_spikes,T2_events,'group','by cell zscore','group_name',group_name,'event_name','Left Ripples','twin',[-1 1])
        plot_perievent_spiketime_histogram(all_spikes,T2_events,'group','by region','group_name',group_name,'event_name','Left Ripples','twin',[-1 1])


        metric_param =[];
        metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
        metric_param.region = @(x) contains(x,'HPC');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        ia = find(cluster_id);
%         ia = ia((41:60));
        C = clusters_combined.cluster_id(ia);

        event_id = [ones(1,length(T1_events)) 2*ones(1,length(T2_events))];
        event_times = [T1_events; T2_events];
        [event_times,index] = sort(event_times);
        event_id=event_id(index);

        plot_perievent_spiketimes(clusters_combined.spike_times,clusters_combined.spike_id,[],[],[5 1],[-1 1],0.02,...
            'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',event_times,...
            'event_id',event_id,'event_label','ripple T1','place_fields',place_fields,'plot_option','by time');

    end

end

zero_meaned_log_odds = [decoded_ripple_events(1).track(1).replay_events(:).z_log_odds] - mean([decoded_ripple_events(1).track(1).replay_events(:).z_log_odds]);
histogram(zero_meaned_log_odds,100)


% histogram([decoded_ripple_events(1).track(1).replay_events(:).z_log_odds],100)
find(zero_meaned_log_odds<-1)
find(zero_meaned_log_odds>1)

zero_meaned_log_odds = zscore([decoded_ripple_events(1).track(1).replay_events(:).z_log_odds]);
T2_events = find(zero_meaned_log_odds<-1);
figure
count=1;
for event = 1:5:length(T2_events)

    subplot(5,5,count)
    T1_data = decoded_ripple_events(1).track(1).replay_events(T2_events(event));
    T2_data = decoded_ripple_events(1).track(2).replay_events(T2_events(event));

    timebins = T1_data.timebins_edges(1:end-1);
    imagesc(timebins,...
        [1:2*size(T2_data.replay,1)],...
        [T1_data.replay; T2_data.replay])
    hold on
    colormap(flipud(gray))
    yticks([30 50 70 90 110 140 170 190 210 230 250 280]/5)
    yline(14.5,'LineWidth',2,'Color','k','DisplayName','Track 2')
    yticklabels([30 50 70 90 110 140 30 50 70 90 110 140])


    [b,index]=min(abs(T1_data.onset-timebins));
    xline(timebins(index)-mean(diff(timebins)/2))
    [b,index]=min(abs(T1_data.offset-timebins));
    xline(timebins(index)-mean(diff(timebins)/2))
    count = count + 1;
end


T1_events = find(zero_meaned_log_odds>1);
figure
count = 1;
for event = 1:length(T1_events)
    subplot(5,5,count)
    T1_data = decoded_ripple_events(1).track(1).replay_events(T1_events(event));
    T2_data = decoded_ripple_events(1).track(2).replay_events(T1_events(event));

    timebins = T1_data.timebins_edges(1:end-1);
    imagesc(timebins,...
        [1:2*size(size(T2_data.replay,1))],...
        [T1_data.replay; T2_data.replay])
    colormap(flipud(gray))
    yticks([30 50 70 90 110 140 170 190 210 230 250 280]/10)
    yline(14.5,'LineWidth',2,'Color','k','DisplayName','Track 2')
    yticklabels([30 50 70 90 110 140 30 50 70 90 110 140])

    hold on
    [b,index]=min(abs(T1_data.onset-timebins));
    xline(timebins(index)-mean(diff(timebins)/2))
    [b,index]=min(abs(T1_data.offset-timebins));
    xline(timebins(index)-mean(diff(timebins)/2))
    count = count + 1;
end




