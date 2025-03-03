%% Main pipeline for extracting and saving and analysing LFP including ripple event and theta cycle detection
% This is a higher-order multi-purpose pipeline that calls dependent functions for
% different analysis piepline

% For mapping of visual receptive field, Please refer to
% Sparse_noise_RF_mapping_masa.mat (subject to change)


%% Set path
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

%% Extract and save session clusters for sleep sessions
clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar 
% experiment_info=experiment_info([4 5 6 ]);
experiment_info=experiment_info([4 5 6 18 19 21 34 35 44 45 58 59 60 71]);
% experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
% experiment_info = experiment_info(4);
% Stimulus_type = 'RUN';
% Stimulus_type = 'SleepChronic';
% all_stimulus_type={'RUN'};
all_stimulus_type={'SleepChronic'};
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
            options = session_info(n).probe(1);
            % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
            DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));
            if isempty(DIR)
                continue
            end

            DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
            DIR1 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));

            if ~isempty(DIR)
                load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
                session_clusters_RUN=session_clusters;
                clear session_clusters
            end

            if ~isempty(DIR1)
                load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
                session_clusters_RUN=session_clusters;
                clear session_clusters
            end

            if contains(stimulus_name{n},'Masa2tracks')
                load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
                session_clusters= session_clusters_RUN;
                load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks4%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
                clusters = clusters_ks4;
            else               
                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
                clusters = clusters_ks4;
            end
            %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'power');
            %         load(fullfile(options.ANALYSIS_DATAPATH,'..','extracted_PSD.mat'),'power');


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
                Behaviour.mobility_thresholded = mobility_thresholded;
%                 mobility_thresholded = interp1(Behaviour.tvec,double(mobility_thresholded),LFP(probe_no).tvec,'linear');
                Behaviour.mobility_zscore = abs([0 diff(movmean(Behaviour.mobility,1/mean(diff(Behaviour.tvec))))]);% Diff of pixel change
                Behaviour.mobility_zscore(isnan(Behaviour.mobility_zscore))=mean(Behaviour.mobility_zscore,'omitnan');
                Behaviour.mobility_zscore=zscore(Behaviour.mobility_zscore);
%                 speed = mobility_thresholded;
            else
                speed = Behaviour.speed;
                speed(isnan(speed))=0;
                w = gausswin(9);
                w = w / sum(w);
                speed = filtfilt(w,1,speed')';
%                 speed = interp1(Behaviour.tvec,speed,LFP(probe_no).tvec,'linear');
            end

            metric_param =[];
            %                  metric_param.cluster_id = @(x) ismember(x,session_clusters_RUN.cluster_id(spatial_cell_index));
            params = create_cluster_selection_params('sorting_option','masa');
            clear selected_clusters V1_clusters session_clusters

            for nprobe = 1:length(clusters)
                options = session_info(n).probe(nprobe);
                probe_no =  session_info(n).probe(nprobe).probe_id+1;

                clusters(nprobe).region = strings(length(clusters(nprobe).cluster_id),1);
                if clusters(nprobe).probe_hemisphere == 1
                    clusters(nprobe).region(:) = 'n.a_L';
                    kClusters=kmeans(clusters(nprobe).peak_depth,2);
                    if mean(clusters(nprobe).peak_depth(kClusters==1))>mean(clusters(nprobe).peak_depth(kClusters==2))
                        % if mean ocation of cluster one is above cluster two, it is
                        % Cortex.
                        % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==1);
                        V1_cell_id = find(kClusters==1);
                        HPC_cell_id = find(kClusters==2);
                    else
                        % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==2);
                        V1_cell_id = find(kClusters==2);
                        HPC_cell_id = find(kClusters==1);
                    end
                    clusters(nprobe).region(V1_cell_id) = 'V1_L';
                    clusters(nprobe).region(HPC_cell_id) = 'HPC_L';

                elseif clusters(nprobe).probe_hemisphere == 2
                    clusters(nprobe).region(:) = 'n.a_R';

                    kClusters=kmeans(clusters(nprobe).peak_depth,2);
                    if mean(clusters(nprobe).peak_depth(kClusters==1))>mean(clusters(nprobe).peak_depth(kClusters==2))
                        % if mean ocation of cluster one is above cluster two, it is
                        % Cortex.
                        % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==1);
                        V1_cell_id = find(kClusters==1);
                        HPC_cell_id = find(kClusters==2);
                    else
                        % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==2);
                        V1_cell_id = find(kClusters==2);
                        HPC_cell_id = find(kClusters==1);
                    end
                    clusters(nprobe).region(V1_cell_id) = 'V1_R';
                    clusters(nprobe).region(HPC_cell_id) = 'HPC_R';
                end

                % Convert to unique spike/cluster id
                probe_clusters = select_clusters(clusters(nprobe),params); %only look at good clusters

                single_cluster_spike_id = cell(size(probe_clusters.cluster_id));
                single_cluster_spike_times = cell(size(probe_clusters.cluster_id));
                unique_spike_id = probe_clusters.spike_id + nprobe*10^4 + iDate* 10^6 + str2double(options.SUBJECT(2:end))*10^8;
                unique_cluster_id = probe_clusters.cluster_id + nprobe*10^4 + iDate* 10^6 + str2double(options.SUBJECT(2:end))*10^8;
                [unique_cluster_ids,first_index] = unique(unique_cluster_id);
                for ncell =1:length(unique_cluster_ids)
                    single_cluster_spike_id{ncell,1} = unique_spike_id(unique_spike_id == unique_cluster_ids(ncell));
                    single_cluster_spike_times{ncell,1} = probe_clusters.spike_times(unique_spike_id == unique_cluster_ids(ncell));
                end
                probe_clusters.cluster_id = unique_cluster_id(first_index);
                probe_clusters.spike_id = single_cluster_spike_id;
                probe_clusters.spike_times = single_cluster_spike_times;

                if session_info(n).probe(nprobe).probe_hemisphere==1
                    probe_clusters.region = session_clusters_RUN.region(contains(session_clusters_RUN.region,'L'));
                    probe_clusters.spatial_response = session_clusters_RUN.spatial_response(contains(session_clusters_RUN.region,'L'),:);
                    probe_clusters.odd_even_stability = session_clusters_RUN.odd_even_stability(contains(session_clusters_RUN.region,'L'),:);
                    probe_clusters.peak_percentile = session_clusters_RUN.peak_percentile(contains(session_clusters_RUN.region,'L'),:);
                elseif session_info(n).probe(nprobe).probe_hemisphere==2
                    probe_clusters.region = session_clusters_RUN.region(contains(session_clusters_RUN.region,'R'));
                    probe_clusters.spatial_response = session_clusters_RUN.spatial_response(contains(session_clusters_RUN.region,'R'),:);
                    probe_clusters.odd_even_stability = session_clusters_RUN.odd_even_stability(contains(session_clusters_RUN.region,'R'),:);
                    probe_clusters.peak_percentile = session_clusters_RUN.peak_percentile(contains(session_clusters_RUN.region,'R'),:);
                end

                session_clusters(nprobe)=probe_clusters;

                metric_param =[];
                if options.probe_hemisphere==1
                    metric_param.region = @(x) contains(x,'V1_L');
                    %                     metric_param.region = @(x) contains(x,'x');
                    V1_clusters(probe_no) = select_clusters(clusters(nprobe),metric_param);
                elseif options.probe_hemisphere==2
                    metric_param.region = @(x) contains(x,'V1_R');
                    V1_clusters(probe_no) = select_clusters(clusters(nprobe),metric_param);
                end
            end

            session_clusters= combine_clusters_from_multiple_probes(session_clusters(1),session_clusters(2));
            % add behaviour info for each session
            behaviour_fields = fieldnames(Behaviour);
            for iF = 1:length(behaviour_fields)
                session_clusters.(behaviour_fields{iF}) = {Behaviour.(behaviour_fields{iF})};
            end

            if contains(lower(stimulus_name{n}),'sleep') % 
                save(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
                clusters_ks4 = clusters;
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'),'clusters_ks4'); % save region info 
            else
                clusters_ks4 = clusters;
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks4%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'clusters_ks4');
            end
        end
    end
end


%% Extract LFP for sleep and RUN, analyse sleep and slow oscillation and then save LFP from selected channels
clear all
% SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar 
% experiment_info=experiment_info([4 5 6 ]);
experiment_info=experiment_info([4 5 6 18 19 21 34 35 44 45 58 59 60 71]);
% experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
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
            % options = session_info(n).probe(1);
            LFP_extraction_and_event_detection_pipeline(session_info(n),stimulus_name{n},best_channels)
        end
    end
end

%% LFP PSD slope and cortical wave direction
% pyversion('C:\Users\masahiro.takigawa\.conda\envs\fooof\python')
pyversion('C:\Users\masah\anaconda3\envs\fooof\python')

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
pyenv("ExecutionMode","OutOfProcess")

clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar 
% experiment_info=experiment_info([4 5 6 ]);
experiment_info=experiment_info([4 5 6 18 19 21 34 35 44 45 58 59 60 71]);
Stimulus_type = 'Sleep';

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

        DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
        DIR1 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));

        if ~isempty(DIR)
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
            session_clusters_RUN=session_clusters;
            clear session_clusters
        end

        if ~isempty(DIR1)
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
            session_clusters_RUN=session_clusters;
            clear session_clusters
        end
        tic
        if contains(stimulus_name{n},'Masa2tracks')
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP');
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks4%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            clusters=clusters_ks4;
        elseif contains(stimulus_name{n},'Sleep')
            % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
            %             load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
            % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state_merged.mat'));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
            % load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
            %             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
            clusters=clusters_ks4;
        else
            % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
            % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
            clusters=clusters_ks4;
        end
        toc
        
        params = create_cluster_selection_params('sorting_option','masa');

        %%% Test slow wave detection
        tvec = LFP(1).tvec;
        clear slow_waves_markov
        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id+1;
            options = session_info(n).probe(nprobe);

            %%%%% V1 MUA spikes
            metric_param =[]; % get all V1 spikes from left or right hemisphere
            if options.probe_hemisphere==1
                metric_param.region = @(x) contains(x,'HPC_L');
                HPC_clusters = select_clusters(clusters(nprobe),metric_param);

                V1_clusters=[];
                V1_clusters.cluster_id = session_clusters.cluster_id(session_clusters.region=="V1_L");
                V1_clusters.spike_times = vertcat(session_clusters.spike_times{session_clusters.region=="V1_L"});
                V1_clusters.spike_id = vertcat(session_clusters.spike_id{session_clusters.region=="V1_L"});

                if length(V1_clusters.cluster_id)<150
                    V1_clusters=[];
                    metric_param.region = @(x) contains(x,'V1_L');
                    %                 metric_param.amplitude_median = @(x) abs(x)>50; %IBL 50 but bombcell 20
                    V1_clusters = select_clusters(clusters(nprobe),metric_param);
                end
            elseif options.probe_hemisphere==2
                %                 metric_param.amplitude_median = @(x) abs(x)>50; %IBL 50 but bombcell 20
                metric_param.region = @(x) contains(x,'HPC_R');
                HPC_clusters = select_clusters(clusters(nprobe),metric_param);

                V1_clusters=[];
                V1_clusters.cluster_id = session_clusters.cluster_id(session_clusters.region=="V1_R");
                V1_clusters.spike_times = vertcat(session_clusters.spike_times{session_clusters.region=="V1_R"});
                V1_clusters.spike_id = vertcat(session_clusters.spike_id{session_clusters.region=="V1_R"});

                if length(V1_clusters.cluster_id)<150
                    V1_clusters=[];
                    metric_param.region = @(x) contains(x,'V1_R');
                    %                 metric_param.amplitude_median = @(x) abs(x)>50; %IBL 50 but bombcell 20
                    V1_clusters = select_clusters(clusters(nprobe),metric_param);
                end
            end
            
            
            best_channel = find(LFP(probe_no).best_V1_channel==slow_waves(nprobe).best_channel);

            if isempty(best_channel)
                [~,best_channel] = max(LFP(nprobe).best_V1_high_freq_power(:,7));
                slow_waves_markov(probe_no) = detect_UP_DOWN_markov(tvec,slow_waves(probe_no),[V1_clusters.spike_times V1_clusters.spike_id],behavioural_state_merged);
%                 temp = DetectSlowWaves_masa('time',tvec,'lfp',LFP(probe_no).best_V1_high_freq(best_channel,:),'spikes',V1_clusters(probe_no),'NREMInts',behavioural_state_merged.SWS,'sensitivity',0.5);
                temp = detect_UP_DOWN_markov(tvec,slow_waves(probe_no),[V1_clusters.spike_times V1_clusters.spike_id],behavioural_state_merged);
                temp.best_channel = LFP(nprobe).best_V1_high_freq_channel(best_channel);
                slow_waves_markov(probe_no) = temp;
            else
%                 temp = DetectSlowWaves_masa('time',tvec,'lfp',LFP(probe_no).best_V1(best_channel,:),'spikes',V1_clusters(probe_no),'NREMInts',behavioural_state_merged.SWS,'sensitivity',0.5);
                temp = detect_UP_DOWN_markov(tvec,slow_waves(probe_no),[V1_clusters.spike_times V1_clusters.spike_id],behavioural_state_merged);
                temp.best_channel = slow_waves(probe_no).best_channel;
                slow_waves_markov(probe_no) = temp;
            end 
        end

        % PSD slope quantification using fooof ()
        tvec = LFP(1).tvec;
        SR = round(1/mean(diff(tvec)));
        nfft_seconds= 2;
        nfft = 2^(nextpow2(SR*nfft_seconds));
        win  = hanning(nfft);

        clipDur = 5; % seconds
        timebin_edges = tvec(1):clipDur:tvec(end); % 5 seconds timebin edges for PSD slope
        nClipSamps = round(SR*clipDur);

        PSD_slope=[];

        disp('PSD slope for slow waves started')
        tic
        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id+1;

            slow_waves_markov(probe_no).power = [];
            slow_waves_markov(probe_no).frequency = [];
            slow_waves_markov(probe_no).PSD_slope = [];
            slow_waves_markov(probe_no).timebin_edges = [];
            slow_waves_markov(probe_no).DOWN_PSD_slope = [];
            slow_waves_markov(probe_no).UP_PSD_slope = [];

            if ~isempty(behavioural_state_merged.SWS)
                best_channel = find(LFP(probe_no).best_V1_channel==slow_waves_markov(probe_no).best_channel);
                V1_LFP = LFP(probe_no).best_V1(best_channel,:);

            else
                continue
            end

            % slow_waves(probe_no).power=[];
            % slow_waves(probe_no).powerdB=[];
            nClips = floor(length(V1_LFP)/nClipSamps);
            samples_to_pass = 0;

            for clip = 1:nClips
                tidx = [1+samples_to_pass:samples_to_pass+nClipSamps]; % in samples
                % timeWindow = [1+start_samp+samples_to_pass:start_samp+samples_to_pass+nClipSamps]; % in seconds
                % timebin(clip) = round((tvec(tidx(1))+ tvec(tidx(end)))/2);


                [pxx,fxx] = pwelch(V1_LFP(tidx),win,[],nfft,SR);

                slow_waves_markov(probe_no).power(clip,:) = pxx;
                % slow_waves(probe_no).powerdB(clip,:) = 10*log10(pxx);
                fxx = fxx';
                % FOOOF settings
                settings = struct();
                f_range = [1, 10];

                % Run FOOOF to fit the model and quantify the slope of power spectra
                fooof_results = fooof(fxx, pxx, f_range, settings);
                % fooof_plot(fooof_results)
                PSD_slope(clip) = fooof_results.aperiodic_params(2); % slope of aperiodic component

                samples_to_pass = samples_to_pass + nClipSamps;
            end

            slow_waves_markov(probe_no).frequency = fxx;
            slow_waves_markov(probe_no).PSD_slope = PSD_slope;
            slow_waves_markov(probe_no).timebin_edges = timebin_edges;
            
            delta_power=[];
            delta_power = zscore(mean(10*log10(slow_waves_markov(probe_no).power(:,slow_waves_markov(probe_no).frequency>0.5&slow_waves_markov(probe_no).frequency<4)),2));
%             delta_power = mean(10*log10(slow_waves(probe_no).power(:,slow_waves(probe_no).frequency>0.5&slow_waves(probe_no).frequency<4)),2);
            
%             timebin_centre = timebin_edges(1)+mean(diff(timebin_edges))/2:mean(diff(timebin_edges)):timebin_edges(end)-mean(diff(timebin_edges))/2;
%             NREM_delta_power = delta_power(find(ismember(timebin_centre,Restrict(timebin_centre,behavioural_state_merged.SWS))));
%             delta_power = (delta_power-mean(NREM_delta_power))/std(NREM_delta_power);

            event_midpoint = mean(slow_waves_markov(probe_no).UP_ints,2);
            % Find the bin indices for each event time
            [~, ~, bin_indices] = histcounts(event_midpoint, timebin_edges);
            slow_waves_markov(probe_no).UP_PSD_slope= PSD_slope(bin_indices);
            slow_waves_markov(probe_no).UP_delta_power= delta_power(bin_indices);


            event_midpoint = mean(slow_waves_markov(probe_no).DOWN_ints,2);
            % Find the bin indices for each event time
            [~, ~, bin_indices] = histcounts(event_midpoint, timebin_edges);
            slow_waves_markov(probe_no).DOWN_PSD_slope= PSD_slope(bin_indices);
            slow_waves_markov(probe_no).DOWN_delta_power= delta_power(bin_indices);
        end
        toc
        disp('PSD slope for slow waves finished')

        disp('cortical wave travleing direction analysis started')
        tic
        %%%%%%%%%%%% Cortical wave direction during DOWN state peak
        % -1 is posterior -> anterior, 0 is no delay or noisy delay and 1 is anterior -> posterior
        filterparms.deltafilter = [0.5 8];%heuristically defined.  room for improvement here.
        filterparms.gammafilter = [100 400];
        filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
        filterparms.gammanormwin = 20; %window for gamma normalization (s)

        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id+1;
            lfp.samplingRate = round(1/mean(diff(LFP(probe_no).tvec)));
            lfp.timestamps = LFP(probe_no).tvec;
            lfp.data=[];
            DOWN_peaks_shank = [];
            peaks_latency = [];
            DOWN_traveling = [];
            slow_waves_markov(probe_no).DOWN_peaks_shank = [];
            slow_waves_markov(probe_no).DOWN_peaks_latency = [];
            slow_waves_markov(probe_no).DOWN_traveling = [];

            if isfield(LFP(probe_no),'average_V1_xcoord') & ~isempty(behavioural_state_merged.SWS)% if exist best V1 channel for sleep

                lfp.data= LFP(probe_no).average_V1';
                deltaLFP = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);

                zscored_LFP = [];
                zscored_LFP = zscore(deltaLFP.data);
                % [ordered_xcoord,~]=sort(LFP(probe_no).best_V1_xcoord);

                for nevent = 1:size(slow_waves_markov(probe_no).DOWN_ints,1)
                    % for nevent = 600:640
                    midpoint = mean(slow_waves_markov(probe_no).DOWN_ints(nevent,:));
                    tidx = FindInInterval(tvec,[midpoint-0.1 midpoint+0.1]);
                    % tidx = FindInInterval(tvec,[slow_waves(probe_no).ints.UP(nevent,1)-0.3 slow_waves(probe_no).ints.UP(nevent,1)+0.3]);
                    tidx=tidx(1):tidx(end);
                    [~,idx]=min(abs(midpoint-tvec));

                    for nShank=1:length(LFP(probe_no).average_V1_shank_id)

                        [~,peak_id] = findpeaks(zscored_LFP(tidx(1):tidx(end),nShank));

                        [~,temp]=min(abs(midpoint-tvec(tidx(peak_id))));
                        if ~isempty(temp)
                            DOWN_peaks_shank(nShank,nevent) = tvec(tidx(peak_id(temp)));
                        else
                            DOWN_peaks_shank(nShank,nevent) = nan;
                        end
                    end

                    % (diff(ordered_xcoord)/1000000)'./diff(DOWN_peaks_shank(:,nevent))
                    %
                    % diff(DOWN_peaks_shank(:,425))

                    % Putatively
                    if sum(~isnan(DOWN_peaks_shank(:,nevent)))==4 % if delta peaks on four shanks
                        peaks_latency(nevent) = mean(diff(DOWN_peaks_shank(:,nevent)),'omitnan');
                    elseif sum(~isnan(DOWN_peaks_shank(:,nevent)))==3 % if delta peaks on three shanks
                        skipped_shank= diff(LFP(probe_no).average_V1_shank_id(~isnan(DOWN_peaks_shank(:,nevent))))>1;
                        if sum(skipped_shank)>0 % if delta peak skipped one shank

                            if skipped_shank(1)==1 % if shanks [1 2 4] then delay using 1 -> 2
                                peaks_latency(nevent) = DOWN_peaks_shank(1,nevent)-DOWN_peaks_shank(2,nevent);
                            elseif skipped_shank(2) == 1 % if shanks [1 3 4] then delay using 3 -> 4
                                peaks_latency(nevent) = DOWN_peaks_shank(2,nevent)-DOWN_peaks_shank(3,nevent);
                            end
                        else
                            peaks_latency(nevent) = mean(diff(DOWN_peaks_shank(:,nevent)),'omitnan');
                        end
                    elseif sum(~isnan(DOWN_peaks_shank(:,nevent)))==2 % if delta peaks on two shanks
                        peaks_latency(nevent) = diff(DOWN_peaks_shank(1:2,nevent));
                    else % only one peak. Can't calculate latency
                        peaks_latency(nevent) =nan;
                    end

                    % From https://www.jneurosci.org/content/34/26/8875#sec-2
                    % speed is ~40 milimeter per seconds
                    % 6.25 miliseconds to travel 250 micrometer (rough shank spacing)
                    % putatively set the minimum delay threshold to be 3
                    % miliseconds.
                    if session_clusters.probe_hemisphere(nprobe)==1

                        % if direction of mean latency is consistent with the latency more than half
                        % of the shank latency
                        % (i.e. if four shanks, at least 2 jumps between three shanks should be in the same direction as the mean latency)
                        % (if three shanks, then still 2 jumps needed)
                        if peaks_latency(nevent)>0.003 & sum(DOWN_peaks_shank(2:end,nevent)-DOWN_peaks_shank(1,nevent)>0)>=length(LFP(probe_no).average_V1_shank_id)/2
                            DOWN_traveling(nevent) = 1; % anterior to posterior
                        elseif peaks_latency(nevent)<-0.003 & sum(DOWN_peaks_shank(2:end,nevent)-DOWN_peaks_shank(1,nevent)<0)>=length(LFP(probe_no).average_V1_shank_id)/2
                            DOWN_traveling(nevent) = -1; % posterior to anterior
                        else
                            DOWN_traveling(nevent) = 0; % noisy or standing wave?
                        end
                        % DOWN_peaks_shank(:,nevent) - DOWN_peaks_shank(1,nevent)
                    elseif session_clusters.probe_hemisphere(nprobe)==2
                        if peaks_latency(nevent)>0.003 & sum(DOWN_peaks_shank(2:end,nevent)-DOWN_peaks_shank(1,nevent)>0)>=length(LFP(probe_no).average_V1_shank_id)/2
                            DOWN_traveling(nevent) = -1; % posterior to anterior
                        elseif peaks_latency(nevent)<-0.003 & sum(DOWN_peaks_shank(2:end,nevent)-DOWN_peaks_shank(1,nevent)<0)>=length(LFP(probe_no).average_V1_shank_id)/2
                            DOWN_traveling(nevent) = 1; % anterior to posterior
                        else
                            DOWN_traveling(nevent) = 0; % noisy or standing wave?
                        end
                    end

                    % nexttile
                    % hold on;xline(tvec(idx));
                    % plot(tvec(tidx),zscored_LFP(tidx,:));
                    % % hold on;xline(median(tvec(tidx(1):tidx(end))))
                    % xline(DOWN_peaks_shank(:,nevent)')
                end


            end

            slow_waves_markov(probe_no).DOWN_peaks_shank = DOWN_peaks_shank;
            slow_waves_markov(probe_no).DOWN_peaks_latency = peaks_latency;
            slow_waves_markov(probe_no).DOWN_traveling = DOWN_traveling;
        end
        disp('cortical wave travleing direction analysis finished')
        toc
        % 

        disp('Sharp wave amplitude started')
        tic
        %%%%%%%%%%%% Cortical wave direction during DOWN state peak
        % -1 is posterior -> anterior, 0 is no delay or noisy delay and 1 is anterior -> posterior
        filterparms.deltafilter = [0.5 9];%heuristically defined.  room for improvement here.
        filterparms.gammafilter = [100 400];
        filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
        filterparms.gammanormwin = 20; %window for gamma normalization (s)

        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id+1;
            lfp.samplingRate = round(1/mean(diff(LFP(probe_no).tvec)));
            lfp.timestamps = LFP(probe_no).tvec;
            lfp.data=[];
            % DOWN_peaks_shank = [];
            sharp_wave_peaks_shank=[];
            peaks_latency = [];
            wave_traveling = [];
            sharp_wave_zscore_shank=[];


            ripples(probe_no).sharp_wave_peaktimes = [];
            ripples(probe_no).sharp_wave_latency = [];
            ripples(probe_no).sharp_wave_zscore = sharp_wave_zscore_shank;
            ripples(probe_no).SW_traveling = [];

            if isfield(LFP(probe_no),'best_HPC') % if exist best HPC channels

                lfp.data= LFP(probe_no).best_HPC';
                deltaLFP = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);

                zscored_LFP = [];
                zscored_LFP = zscore(deltaLFP.amp);
                % [ordered_xcoord,~]=sort(LFP(probe_no).best_V1_xcoord);

                for nevent = 1:length(ripples(probe_no).onset)
                % for nevent = 600:640
                    tidx = FindInInterval(tvec,[ripples(probe_no).onset(nevent)-0.05 ripples(probe_no).offset(nevent)]);
                    % tidx = FindInInterval(tvec,[slow_waves(probe_no).ints.UP(nevent,1)-0.3 slow_waves(probe_no).ints.UP(nevent,1)+0.3]);
                    tidx=tidx(1):tidx(end);
                    [~,idx]=min(abs(ripples(probe_no).onset(nevent)-tvec));

                    for nShank=1:length(LFP(probe_no).best_HPC_shank_id)

                        [~,peak_id] = findpeaks(zscored_LFP(tidx(1):tidx(end),nShank));

                        if ~isempty(peak_id)
                            [~,temp]=min(abs(ripples(probe_no).onset(nevent)-tvec(tidx(peak_id))));
                            sharp_wave_peaks_shank(nShank,nevent) = tvec(tidx(peak_id(temp)));
                            sharp_wave_zscore_shank(nShank,nevent) = zscored_LFP(tidx(peak_id(temp)),nShank);
                        else
                            sharp_wave_peaks_shank(nShank,nevent) = nan;
                            sharp_wave_zscore_shank(nShank,nevent) = nan;
                        end
                    end

                    % (diff(ordered_xcoord)/1000000)'./diff(DOWN_peaks_shank(:,nevent))
                    %
                    % diff(DOWN_peaks_shank(:,425))

                    % Putatively
                    if sum(~isnan(sharp_wave_peaks_shank(:,nevent)))==4 % if delta peaks on four shanks
                        peaks_latency(nevent) = mean(diff(sharp_wave_peaks_shank(:,nevent)),'omitnan');
                    elseif sum(~isnan(sharp_wave_peaks_shank(:,nevent)))==3 % if delta peaks on three shanks
                        skipped_shank= diff(LFP(probe_no).average_V1_shank_id(~isnan(sharp_wave_peaks_shank(:,nevent))))>1;
                        if sum(skipped_shank)>0 % if delta peak skipped one shank

                            if skipped_shank(1)==1 % if shanks [1 2 4] then delay using 1 -> 2
                                peaks_latency(nevent) = sharp_wave_peaks_shank(1,nevent)-sharp_wave_peaks_shank(2,nevent);
                            elseif skipped_shank(2) == 1 % if shanks [1 3 4] then delay using 3 -> 4
                                peaks_latency(nevent) = sharp_wave_peaks_shank(2,nevent)-sharp_wave_peaks_shank(3,nevent);
                            end
                        else
                            peaks_latency(nevent) = mean(diff(sharp_wave_peaks_shank(:,nevent)),'omitnan');
                        end
                    elseif sum(~isnan(sharp_wave_peaks_shank(:,nevent)))==2 % if delta peaks on two shanks
                        peaks_latency(nevent) = diff(sharp_wave_peaks_shank(1:2,nevent));
                    else % only one peak. Can't calculate latency
                        peaks_latency(nevent) =nan;
                    end

                    % From https://pmc.ncbi.nlm.nih.gov/articles/PMC10125019/
                    % roughly 5 miliseconds to travel 500 micrometer (rough shank spacing)
                    % putatively set the minimum delay threshold to be 3
                    % miliseconds.
                    if session_clusters.probe_hemisphere(nprobe)==1

                        % if direction of mean latency is consistent with the latency more than half
                        % of the shank latency
                        % (i.e. if four shanks, at least 2 jumps between three shanks should be in the same direction as the mean latency)
                        % (if three shanks, then still 2 jumps needed)
                        if peaks_latency(nevent)>0.003 & sum(sharp_wave_peaks_shank(2:end,nevent)-sharp_wave_peaks_shank(1,nevent)>0)>=length(LFP(probe_no).average_V1_shank_id)/2
                            wave_traveling(nevent) = 1; % anterior to posterior
                        elseif peaks_latency(nevent)<-0.003 & sum(sharp_wave_peaks_shank(2:end,nevent)-sharp_wave_peaks_shank(1,nevent)<0)>=length(LFP(probe_no).average_V1_shank_id)/2
                            wave_traveling(nevent) = -1; % posterior to anterior
                        else
                            wave_traveling(nevent) = 0; % noisy or standing wave?
                        end
                        % DOWN_peaks_shank(:,nevent) - DOWN_peaks_shank(1,nevent)
                    elseif session_clusters.probe_hemisphere(nprobe)==2
                        if peaks_latency(nevent)>0.003 & sum(sharp_wave_peaks_shank(2:end,nevent)-sharp_wave_peaks_shank(1,nevent)>0)>=length(LFP(probe_no).average_V1_shank_id)/2
                            wave_traveling(nevent) = -1; % posterior to anterior
                        elseif peaks_latency(nevent)<-0.003 & sum(sharp_wave_peaks_shank(2:end,nevent)-sharp_wave_peaks_shank(1,nevent)<0)>=length(LFP(probe_no).average_V1_shank_id)/2
                            wave_traveling(nevent) = 1; % anterior to posterior
                        else
                            wave_traveling(nevent) = 0; % noisy or standing wave?
                        end
                    end

                    % nexttile
                    % hold on;xline(tvec(idx));
                    % plot(tvec(tidx),zscored_LFP(tidx,:));
                    % % hold on;xline(median(tvec(tidx(1):tidx(end))))
                    % xline(sharp_wave_peaks_shank(:,nevent)')
                end


            end

            ripples(probe_no).sharp_wave_peaktimes = sharp_wave_peaks_shank;
            ripples(probe_no).sharp_wave_latency = peaks_latency;
            ripples(probe_no).sharp_wave_zscore = sharp_wave_zscore_shank;
            ripples(probe_no).SW_traveling = wave_traveling;
            
        end
        disp('HPC sharp wave analysis finished')
        toc

        % histogram(slow_waves(1).DOWN_peaks_latency,100,'Normalization','cdf');
        % hold on; xline(prctile(slow_waves(1).DOWN_peaks_latency,50));
        % histogram(slow_waves(2).DOWN_peaks_latency,100,'Normalization','cdf');
        % hold on; xline(prctile(slow_waves(2).DOWN_peaks_latency,50));
        % 
        % histogram(slow_waves(1).DOWN_travling,3,'Normalization','probability');hold on;
        % histogram(slow_waves(2).DOWN_travling,3,'Normalization','probability');

        % histogram(slow_waves(1).DOWN_peaks_latency,100,'Normalization','probability');
        % hold on; xline(prctile(slow_waves(1).DOWN_peaks_latency,50));
        % histogram(slow_waves(2).DOWN_peaks_latency,100,'Normalization','probability');
        % hold on; xline(prctile(slow_waves(2).DOWN_peaks_latency,50));
        if contains(stimulus_name{n},'Masa2tracks')
            save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'ripples');
            save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_slow_wave_markov_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'slow_waves_markov');
            % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_slow_wave_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'slow_waves');
        else
            % save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'),'slow_waves');
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_markov_events.mat'),'slow_waves_markov');
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'),'ripples');
        end
    end
end





%% UP DOWN state and ripple and spindle analysis (BACKUP)

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

        DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
        DIR1 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));

        if ~isempty(DIR)
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
            session_clusters_RUN=session_clusters;
            clear session_clusters
        end

        if ~isempty(DIR1)
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
            session_clusters_RUN=session_clusters;
            clear session_clusters
        end

        if contains(stimulus_name{n},'Masa2tracks')
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP');
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks4%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            clusters=clusters_ks4;
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
            % load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
%             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
            clusters=clusters_ks4;
        else
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
            
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
            clusters=clusters_ks4;
        end


        % From cell structure back to spike times and spike id
        session_clusters_RUN.spike_id=vertcat(session_clusters_RUN.spike_id{:});
        session_clusters_RUN.spike_times=vertcat(session_clusters_RUN.spike_times{:});
        [session_clusters_RUN.spike_times,index] =sort(session_clusters_RUN.spike_times);
        session_clusters_RUN.spike_id=session_clusters_RUN.spike_id(index);

        clusters_combined= session_clusters;

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

        % HPC frames
        if ~isempty(behavioural_state(probe_no).SWS)
            [HPC_frames] = detect_candidate_frames_masa(LFP(probe_no).tvec,V1_clusters(probe_no).spike_times,behavioural_state(probe_no).SWS,0.01,options);
        end

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




