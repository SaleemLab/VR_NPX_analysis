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
% experiment_info=experiment_info([33 9 10 14]);
% experiment_info=experiment_info([4 5 6 18 19 21 34 35 44 45 58 59 60 71]);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
% experiment_info = experiment_info(4);
% Stimulus_type = 'RUN';
% Stimulus_type = 'SleepChronic';
% all_stimulus_type={'RUN'};
all_stimulus_type={'SleepChronic'};

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
            % LFP_extraction_and_event_detection_pipeline(session_info(n),stimulus_name{n},best_channels)
            % extract_LFP_NPX(session_info(n),stimulus_name{n},best_channels)
            % detect_behavioural_and_brain_states(session_info(n),stimulus_name{n},best_channels)
%             detect_behavioural_and_brain_states_add_on(session_info(n),stimulus_name{n},best_channels)
            extract_LFP_NPX_add_on(session_info(n),stimulus_name{n},best_channels)

        end
    end
end

%% LFP PSD slope and cortical wave direction
% pyversion('C:\Users\masahiro.takigawa\.conda\envs\fooof\python')
% pyversion('C:\Users\masah\anaconda3\envs\fooof\python')

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
pyenv("ExecutionMode","OutOfProcess")

clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar
% experiment_info=experiment_info([4 5 6 ]);
% experiment_info=experiment_info([4 5 6 18 19 21 34 35 44 45 58 59 60 71]);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
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


        for nprobe = 1:length(slow_waves)
            if isfield(slow_waves(nprobe),'ints')
                slow_waves(nprobe).UP_ints = slow_waves(nprobe).ints.UP;
                slow_waves(nprobe).DOWN_ints = slow_waves(nprobe).ints.DOWN;
            end
        end

        if isfield(slow_waves(1),'ints')
            slow_waves = rmfield(slow_waves,'ints');
        end

        params = create_cluster_selection_params('sorting_option','masa');
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

            slow_waves(probe_no).power = [];
            slow_waves(probe_no).frequency = [];
            slow_waves(probe_no).PSD_slope = [];
            slow_waves(probe_no).timebin_edges = [];
            slow_waves(probe_no).DOWN_PSD_slope = [];
            slow_waves(probe_no).UP_PSD_slope = [];

            if ~isempty(behavioural_state_merged.SWS)
                best_channel = find(LFP(probe_no).best_V1_channel==slow_waves(probe_no).best_channel);
                V1_LFP = LFP(probe_no).best_V1(best_channel,:);

            else
                continue
            end

            % if isempty(slow_waves_markov(probe_no).UP_DOWN_index)
            %     continue
            % end
            %
            % if size(slow_waves_markov(probe_no).UP_DOWN_index,1) < 10
            %     continue
            % end

            % slow_waves(probe_no).power=[];
            % slow_waves(probe_no).powerdB=[];
            nClips = floor(length(V1_LFP)/nClipSamps);
            samples_to_pass = 0;

            for clip = 1:nClips
                tidx = [1+samples_to_pass:samples_to_pass+nClipSamps]; % in samples
                % timeWindow = [1+start_samp+samples_to_pass:start_samp+samples_to_pass+nClipSamps]; % in seconds
                % timebin(clip) = round((tvec(tidx(1))+ tvec(tidx(end)))/2);


                [pxx,fxx] = pwelch(V1_LFP(tidx),win,[],nfft,SR);

                slow_waves(probe_no).power(clip,:) = pxx;
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

            slow_waves(probe_no).frequency = fxx;
            slow_waves(probe_no).PSD_slope = PSD_slope;
            slow_waves(probe_no).timebin_edges = timebin_edges;

            delta_power=[];
            delta_power = zscore(mean(10*log10(slow_waves(probe_no).power(:,slow_waves(probe_no).frequency>0.5&slow_waves(probe_no).frequency<4)),2));
            %             delta_power = mean(10*log10(slow_waves(probe_no).power(:,slow_waves(probe_no).frequency>0.5&slow_waves(probe_no).frequency<4)),2);

            %             timebin_centre = timebin_edges(1)+mean(diff(timebin_edges))/2:mean(diff(timebin_edges)):timebin_edges(end)-mean(diff(timebin_edges))/2;
            %             NREM_delta_power = delta_power(find(ismember(timebin_centre,Restrict(timebin_centre,behavioural_state_merged.SWS))));
            %             delta_power = (delta_power-mean(NREM_delta_power))/std(NREM_delta_power);

            event_midpoint = mean(slow_waves(probe_no).UP_ints,2);
            % Find the bin indices for each event time
            [~, ~, bin_indices] = histcounts(event_midpoint, timebin_edges);
            slow_waves(probe_no).UP_PSD_slope= PSD_slope(bin_indices);
            slow_waves(probe_no).UP_delta_power= delta_power(bin_indices);


            event_midpoint = mean(slow_waves(probe_no).DOWN_ints,2);
            % Find the bin indices for each event time
            [~, ~, bin_indices] = histcounts(event_midpoint, timebin_edges);
            slow_waves(probe_no).DOWN_PSD_slope= PSD_slope(bin_indices);
            slow_waves(probe_no).DOWN_delta_power= delta_power(bin_indices);
        end
        toc
        disp('PSD slope for slow waves finished')

        disp('cortical wave analysis started')
        tic
        %%%%%%%%%%%% Cortical wave direction during DOWN state peak
        % -1 is posterior -> anterior, 0 is no delay or noisy delay and 1 is anterior -> posterior
        % https://www.nature.com/articles/s41467-019-10327-5 levenstein et
        % al used 0.5 to 8 Hz for slow wave UP DOWN detection
        filterparms.deltafilter = [0.5 4];%heuristically defined.  room for improvement here.
        filterparms.spindlesfilter = [9 17];%heuristically defined.  room for improvement here.
        filterparms.gammafilter = [100 400];
        filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
        filterparms.gammanormwin = 20; %window for gamma normalization (s)

        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id+1;

            if length(LFP(1).best_HPC) < length(LFP(2).best_HPC)
                time_idx = 1:length(LFP(1).best_HPC);
            elseif length(LFP(1).best_HPC) > length(LFP(2).best_HPC)

                time_idx = 1:length(LFP(2).best_HPC);
            else
                time_idx = 1:length(LFP(1).best_HPC) ;
            end

            lfp.samplingRate = round(1/mean(diff(LFP(probe_no).tvec)));
            lfp.timestamps = LFP(probe_no).tvec(time_idx);
            tvec = LFP(probe_no).tvec(time_idx);
            lfp.data=[];
            DOWN_peaks_shank = [];
            peaks_latency = [];
            DOWN_traveling = [];
            DOWN_peaks_zscore= [];

            slow_waves(probe_no).DOWN_peaks_zscore = [];
            slow_waves(probe_no).DOWN_peaktimes = [];

            % slow_waves(probe_no).DOWN_peaks_latency = [];
            % slow_waves(probe_no).DOWN_traveling = [];
            slow_waves(probe_no).probe_hemisphere = [];
            slow_waves(probe_no).shank_id = [];

            if isempty(slow_waves(probe_no).timestamps)
                continue
            end

            if size(slow_waves(probe_no).timestamps,1) < 10
                continue
            end

            if isfield(LFP(probe_no),'average_V1_xcoord') & ~isempty(behavioural_state_merged.SWS)% if exist best V1 channel for sleep
                probe_id = [];

                if length(LFP)==1
                    lfp.data= [LFP(probe_no).average_V1(:,time_idx)'];
                    probe_hemisphere = probe_no*ones(1,length(LFP(probe_no).average_V1_shank_id));
                    slow_waves(probe_no).shank_id = [LFP(probe_no).average_V1_shank_id];
                else
                    lfp.data= [LFP(1).average_V1(:,time_idx)' LFP(2).average_V1(:,time_idx)'];
                    probe_hemisphere = [ones(1,length(LFP(1).average_V1_shank_id)) 2*ones(1,length(LFP(2).average_V1_shank_id))];
                    if size(LFP(probe_no).average_V1_shank_id,1)==1
                        slow_waves(probe_no).shank_id = [LFP(1).average_V1_shank_id LFP(2).average_V1_shank_id];
                    else
                        slow_waves(probe_no).shank_id = [LFP(1).average_V1_shank_id' LFP(2).average_V1_shank_id'];
                    end
                end

                slow_waves(probe_no).probe_hemisphere = probe_hemisphere;


                if nprobe == 1 % Only need to grab once
                    % grab delta LFP
                    deltaLFP = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);
                    zscored_LFP = zscore(deltaLFP.data);
                    SO_phase_LFP = deltaLFP.phase;
                    SO_amplitude_LFP = zscore(deltaLFP.amp);
                    deltaLFP = [];

                    % grab spindles LFP
                    filter_type  = 'bandpass';
                    passband = [9 17];
                    filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for ripple
                    norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
                    b_spindle = fir1(filter_order, norm_freq_range,filter_type);

                    signal = [];
                    for nShank = 1:length(probe_hemisphere)
                        signal(:,nShank) = filtfilt(b_spindle,1, lfp.data(:,nShank));
                    end
%                     spindle_amplitude_LFP = zscore(envelope(signal,round(lfp.samplingRate/5),'rms'));
%                     spindle_amplitude_LFP = smoothdata(spindle_amplitude_LFP,'gaussian',round(lfp.samplingRate/5)); % envelop amplitude of spindles
                    spindle_amplitude_LFP = zscore(abs(hilbert(signal))); % z scored amplitude
                    spindle_phase_LFP = angle(hilbert(signal)); % spindle phase
                    signal = [];

                    % grab ripples LFP
                    filter_type  = 'bandpass';
                    passband = [125 300];
                    filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for ripple
                    norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
                    b_ripple = fir1(filter_order, norm_freq_range,filter_type);

                    signal = [];
                    for nShank = 1:length(probe_hemisphere)
                        signal(:,nShank) = filtfilt(b_ripple,1, lfp.data(:,nShank));
                    end
                    ripple_amplitude_cortex_LFP = zscore(abs(hilbert(signal))); % z scored amplitude
                    ripple_phase_cortex_LFP = angle(hilbert(signal)); % spindle phase
                    signal = [];



                    % grab gamma LFP
                    filter_type  = 'bandpass';
                    passband = [30 60];
                    filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for theta
                    norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
                    b_theta = fir1(filter_order, norm_freq_range,filter_type);

                    signal = [];
                    for nShank = 1:length(probe_hemisphere)
                        signal(:,nShank) = filtfilt(b_theta,1, lfp.data(:,nShank));
                    end

                    gamma_amplitude_LFP = zscore(abs(hilbert(signal))); % z scored amplitude
                    gamma_phase_LFP = angle(hilbert(signal)); % spindle phase
                end

                %%%%% Calculate Delta peak amplitude and timestamp per event
                for nevent = 1:size(slow_waves(probe_no).DOWN_ints,1)
                    % for nevent = 600:640
                    midpoint = mean(slow_waves(probe_no).timestamps(nevent));
                    % midpoint = mean(slow_waves(probe_no).DOWN_ints(nevent,:));

                    tidx = FindInInterval(tvec,[midpoint-0.15 midpoint+0.15]);
                    % tidx = FindInInterval(tvec,[slow_waves(probe_no).ints.UP(nevent,1)-0.3 slow_waves(probe_no).ints.UP(nevent,1)+0.3]);
                    tidx=tidx(1):tidx(end);

                    for nShank=1:length(probe_hemisphere) % across shanks from both probes if using two probes

                        [~,peak_id] = findpeaks(zscored_LFP(tidx(1):tidx(end),nShank));

                        [~,temp]=min(abs(midpoint-tvec(tidx(peak_id))));
                        if ~isempty(temp)
                            DOWN_peaks_shank(nShank,nevent) = tvec(tidx(peak_id(temp)));
                            DOWN_peaks_zscore(nShank,nevent) = SO_amplitude_LFP(tidx(peak_id(temp)),nShank);
                        else
                            DOWN_peaks_shank(nShank,nevent) = nan;
                            DOWN_peaks_zscore(nShank,nevent) = nan;
                        end
                    end

                    % (diff(ordered_xcoord)/1000000)'./diff(DOWN_peaks_shank(:,nevent))
                    %
                    % diff(DOWN_peaks_shank(:,425))

                    % Putative calculation of average latency across shanks
                    % for each hemisphere
                    % for h = unique(probe_hemisphere)
                    %     DOWN_peaks_shank_temp = DOWN_peaks_shank(probe_hemisphere==h,nevent);
                    %
                    %     if sum(~isnan(DOWN_peaks_shank_temp))==4 % if delta peaks on four shanks
                    %         peaks_latency(h,nevent) = mean(diff(DOWN_peaks_shank_temp),'omitnan');
                    %     elseif sum(~isnan(DOWN_peaks_shank_temp))==3 % if delta peaks on three shanks
                    %         skipped_shank= diff(LFP(h).average_V1_shank_id(~isnan(DOWN_peaks_shank_temp)))>1;
                    %         if sum(skipped_shank)>0 % if delta peak skipped one shank
                    %
                    %             if skipped_shank(1)==1 % if shanks [1 2 4] then delay using 1 -> 2
                    %                 peaks_latency(h,nevent) = DOWN_peaks_shank_temp(1)-DOWN_peaks_shank_temp(2);
                    %             elseif skipped_shank(2) == 1 % if shanks [1 3 4] then delay using 3 -> 4
                    %                 peaks_latency(h,nevent) = DOWN_peaks_shank_temp(2)-DOWN_peaks_shank_temp(3);
                    %             end
                    %         else
                    %             peaks_latency(h,nevent) = mean(diff(DOWN_peaks_shank_temp),'omitnan');
                    %         end
                    %     elseif sum(~isnan(DOWN_peaks_shank_temp))==2 % if delta peaks on two shanks (latency divide by number of shanks skipped to caluclate mean latency per 250 micron)
                    %         peaks_latency(h,nevent) = diff(DOWN_peaks_shank_temp(1:2))/abs(LFP(h).average_V1_shank_id(2)-LFP(h).average_V1_shank_id(1));
                    %     else % only one peak. Can't calculate latency
                    %         peaks_latency(h,nevent) =nan;
                    %     end
                    %
                    %     % From https://www.jneurosci.org/content/34/26/8875#sec-2
                    %     % speed is ~40 milimeter per seconds
                    %     % 6.25 miliseconds to travel 250 micrometer (rough shank spacing)
                    %     % putatively set the minimum delay threshold to be 3
                    %     % miliseconds.
                    %     if h==1
                    %
                    %         % if direction of mean latency is consistent with the latency more than half
                    %         % of the shank latency
                    %         % (i.e. if four shanks, at least 2 jumps between three shanks should be in the same direction as the mean latency)
                    %         % (if three shanks, then still 2 jumps needed)
                    %         if peaks_latency(nevent)>0.003 & sum(DOWN_peaks_shank_temp(2:end)-DOWN_peaks_shank_temp(1)>0)>=length(LFP(h).average_V1_shank_id)/2
                    %             DOWN_traveling(h,nevent) = 1; % anterior to posterior
                    %         elseif peaks_latency(nevent)<-0.003 & sum(DOWN_peaks_shank_temp(2:end)-DOWN_peaks_shank_temp(1)<0)>=length(LFP(h).average_V1_shank_id)/2
                    %             DOWN_traveling(h,nevent) = -1; % posterior to anterior
                    %         else
                    %             DOWN_traveling(h,nevent) = 0; % noisy or standing wave?
                    %         end
                    %         % DOWN_peaks_shank(:,nevent) - DOWN_peaks_shank(1,nevent)
                    %     elseif h==2
                    %         if peaks_latency(h,nevent)>0.003 & sum(DOWN_peaks_shank_temp(2:end)-DOWN_peaks_shank_temp(1)>0)>=length(LFP(h).average_V1_shank_id)/2
                    %             DOWN_traveling(h,nevent) = -1; % posterior to anterior
                    %         elseif peaks_latency(h,nevent)<-0.003 & sum(DOWN_peaks_shank_temp(2:end)-DOWN_peaks_shank_temp(1)<0)>=length(LFP(h).average_V1_shank_id)/2
                    %             DOWN_traveling(h,nevent) = 1; % anterior to posterior
                    %         else
                    %             DOWN_traveling(h,nevent) = 0; % noisy or standing wave?
                    %         end
                    %     end
                    % end
                    % nexttile
                    % hold on;xline(tvec(idx));
                    % plot(tvec(tidx),zscored_LFP(tidx,:));
                    % % hold on;xline(median(tvec(tidx(1):tidx(end))))
                    % xline(DOWN_peaks_shank(:,nevent)','r')
                end


            end


            slow_waves(probe_no).DOWN_peaktimes = DOWN_peaks_shank;
            slow_waves(probe_no).DOWN_peaks_zscore = DOWN_peaks_zscore;
            % slow_waves(probe_no).DOWN_peaks_latency = peaks_latency;
            % slow_waves(probe_no).DOWN_traveling = DOWN_traveling;
        end
        disp('cortical wave analysis finished')
        toc
        %


        disp('Sharp wave ripple analysis started')
        tic
        %%%%%%%%%%%% sharp wave direction during ripple peaktimes
        % -1 is posterior -> anterior, 0 is no delay or noisy delay and 1 is anterior -> posterior
        filterparms.deltafilter = [0.5 4];%heuristically defined.  room for improvement here.

        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id+1;


            if length(LFP(1).best_HPC) < length(LFP(2).best_HPC)
                time_idx = 1:length(LFP(1).best_HPC);
            elseif length(LFP(1).best_HPC) > length(LFP(2).best_HPC)

                time_idx = 1:length(LFP(2).best_HPC);
            else
                time_idx = 1:length(LFP(1).best_HPC) ;
            end

            lfp.samplingRate = round(1/mean(diff(LFP(probe_no).tvec)));
            lfp.timestamps = LFP(probe_no).tvec(time_idx);
            tvec = LFP(probe_no).tvec(time_idx);
            lfp.data=[];
            % DOWN_peaks_shank = [];
            sharp_wave_peaks_shank=[];
            peaks_latency = [];
            wave_traveling = [];
            sharp_wave_zscore_shank=[];


            ripples(probe_no).sharp_wave_peaktimes = [];
            ripples(probe_no).sharp_wave_zscore = [];

            ripples(probe_no).SWR_peaktimes = [];
            ripples(probe_no).SWR_zscore = [];
            % ripples(probe_no).SWR_traveling = [];

            ripples(probe_no).probe_hemisphere = [];
            ripples(probe_no).shank_id = [];

            if isfield(LFP(probe_no),'best_HPC') % if exist best HPC channels
                probe_id = [];

                if length(LFP)==1
                    lfp.data= [LFP(probe_no).best_HPC(:,time_idx)'];
                    probe_hemisphere = probe_no*ones(1,length(LFP(probe_no).best_HPC_shank_id));
                    ripples(probe_no).shank_id = [LFP(1).best_HPC_shank_id LFP(2).best_HPC_shank_id];
                else
                    lfp.data= [LFP(1).best_HPC(:,time_idx)' LFP(2).best_HPC(:,time_idx)'];
                    probe_hemisphere = [ones(1,length(LFP(1).best_HPC_shank_id)) 2*ones(1,length(LFP(2).best_HPC_shank_id))];
                    if size(LFP(probe_no).best_HPC_shank_id,1)==1
                        ripples(probe_no).shank_id = [LFP(1).best_HPC_shank_id LFP(2).best_HPC_shank_id];
                    else
                        ripples(probe_no).shank_id = [LFP(1).best_HPC_shank_id' LFP(2).best_HPC_shank_id'];
                    end
                end

                ripples(probe_no).probe_hemisphere = probe_hemisphere;

                if nprobe == 1 % Only need to grab once
                    % grab delta LFP
                    deltaLFP = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);
                    zscored_LFP = zscore(deltaLFP.data);
                    SO_phase_HPC_LFP = deltaLFP.phase;
                    SO_amplitude_HPC_LFP = zscore(deltaLFP.amp);
                    deltaLFP = [];

                    % grab spindles LFP
                    filter_type  = 'bandpass';
                    passband = [9 17];
                    filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for ripple
                    norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
                    b_spindle = fir1(filter_order, norm_freq_range,filter_type);

                    signal = [];
                    for nShank = 1:length(probe_hemisphere)
                        signal(:,nShank) = filtfilt(b_spindle,1, lfp.data(:,nShank));
                    end
                    %                     spindle_amplitude_HPC_LFP = zscore(envelope(signal,round(lfp.samplingRate/5),'rms'));
                    %                     spindle_amplitude_HPC_LFP = smoothdata(spindle_amplitude_HPC_LFP,'gaussian',round(lfp.samplingRate/5)); % envelop amplitude of spindles
                    spindle_amplitude_HPC_LFP = zscore(abs(hilbert(signal))); % z scored amplitude
                    spindle_phase_HPC_LFP = angle(hilbert(signal)); % spindle phase


                    signal = [];

                    % grab ripples LFP
                    filter_type  = 'bandpass';
                    passband = [125 300];
                    filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for ripple
                    norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
                    b_ripple = fir1(filter_order, norm_freq_range,filter_type);

                    signal = [];
                    for nShank = 1:length(probe_hemisphere)
                        signal(:,nShank) = filtfilt(b_ripple,1, lfp.data(:,nShank));
                    end
                    zscored_LFP = zscore(abs(hilbert(signal)));
                    ripple_amplitude_LFP = zscored_LFP; % z scored amplitude
                    ripple_phase_LFP = angle(hilbert(signal)); % spindle phase
                    signal = [];



                    % grab theta LFP
                    filter_type  = 'bandpass';
                    passband = [4 12];
                    filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for theta
                    norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
                    b_theta = fir1(filter_order, norm_freq_range,filter_type);

                    signal = [];
                    for nShank = 1:length(probe_hemisphere)
                        signal(:,nShank) = filtfilt(b_theta,1, lfp.data(:,nShank));
                    end

                    theta_amplitude_HPC_LFP = zscore(abs(hilbert(signal))); % z scored amplitude
                    theta_phase_HPC_LFP = angle(hilbert(signal)); % spindle phase
                end

                % [ordered_xcoord,~]=sort(LFP(probe_no).best_V1_xcoord);
                zscored_LFP = [];
                zscored_LFP = ripple_amplitude_LFP;
                for nevent = 1:length(ripples(probe_no).onset)

                    tidx = FindInInterval(tvec,[ripples(probe_no).onset(nevent)-0.05 ripples(probe_no).offset(nevent)]);
                    % tidx = FindInInterval(tvec,[slow_waves(probe_no).ints.UP(nevent,1)-0.3 slow_waves(probe_no).ints.UP(nevent,1)+0.3]);
                    tidx=tidx(1):tidx(end);


                    for nShank=1:length(probe_hemisphere)

                        [~,peak_id] = findpeaks(zscored_LFP(tidx(1):tidx(end),nShank));

                        if ~isempty(peak_id)
                            [~,temp]=min(abs(ripples(probe_no).peaktimes(nevent)-tvec(tidx(peak_id))));
                            sharp_wave_peaks_shank(nShank,nevent) = tvec(tidx(peak_id(temp)));
                            sharp_wave_zscore_shank(nShank,nevent) = zscored_LFP(tidx(peak_id(temp)),nShank);
                        else
                            sharp_wave_peaks_shank(nShank,nevent) = nan;
                            sharp_wave_zscore_shank(nShank,nevent) = nan;
                        end
                    end

                    % % (diff(ordered_xcoord)/1000000)'./diff(DOWN_peaks_shank(:,nevent))
                    % %
                    % % diff(DOWN_peaks_shank(:,425))
                    %
                    % % Putatively
                    % for h = unique(probe_hemisphere)
                    %     sharp_wave_peaks_shank_temp = sharp_wave_peaks_shank(probe_hemisphere==h,nevent);
                    %
                    %     if sum(~isnan(sharp_wave_peaks_shank_temp))==4 % if delta peaks on four shanks
                    %         peaks_latency(h,nevent) = mean(diff(sharp_wave_peaks_shank_temp),'omitnan');
                    %     elseif sum(~isnan(sharp_wave_peaks_shank_temp))==3 % if delta peaks on three shanks
                    %         skipped_shank= diff(LFP(probe_no).best_HPC_shank_id(~isnan(sharp_wave_peaks_shank_temp)))>1;
                    %         if sum(skipped_shank)>0 % if delta peak skipped one shank
                    %
                    %             if skipped_shank(1)==1 % if shanks [1 2 4] then delay using 1 -> 2
                    %                 peaks_latency(h,nevent) = sharp_wave_peaks_shank_temp(1)-sharp_wave_peaks_shank_temp(2);
                    %             elseif skipped_shank(2) == 1 % if shanks [1 3 4] then delay using 3 -> 4
                    %                 peaks_latency(h,nevent) = sharp_wave_peaks_shank_temp(2)-sharp_wave_peaks_shank_temp(3);
                    %             end
                    %         else
                    %             peaks_latency(h,nevent) = mean(diff(sharp_wave_peaks_shank_temp),'omitnan');
                    %         end
                    %     elseif sum(~isnan(sharp_wave_peaks_shank_temp))==2 % if delta peaks on two shanks
                    %         peaks_latency(h,nevent) = diff(sharp_wave_peaks_shank_temp(1:2))/abs(LFP(h).best_HPC_shank_id(2)-LFP(h).best_HPC_shank_id(1));
                    %     else % only one peak. Can't calculate latency
                    %         peaks_latency(h,nevent) =nan;
                    %     end
                    %
                    %     % From Patel et al. (2013) https://pmc.ncbi.nlm.nih.gov/articles/PMC3807028/#sec2
                    %     % ripple propogation speed roughly 0.35 m/s or 0.7
                    %     % ms per 250 micron
                    %     % putatively set the minimum delay threshold to be 2
                    %     % miliseconds.
                    %     if h==1
                    %
                    %         % if direction of mean latency is consistent with the latency more than half
                    %         % of the shank latency
                    %         % (i.e. if four shanks, at least 2 jumps between three shanks should be in the same direction as the mean latency)
                    %         % (if three shanks, then still 2 jumps needed)
                    %         if peaks_latency(h,nevent)>0.001 & sum(sharp_wave_peaks_shank_temp(2:end)-sharp_wave_peaks_shank_temp(1)>0)>=length(LFP(h).best_HPC_shank_id)/2
                    %             wave_traveling(h,nevent) = 1; % anterior to posterior
                    %         elseif peaks_latency(h,nevent)<-0.001 & sum(sharp_wave_peaks_shank_temp(2:end)-sharp_wave_peaks_shank_temp(1)<0)>=length(LFP(h).best_HPC_shank_id)/2
                    %             wave_traveling(h,nevent) = -1; % posterior to anterior
                    %         else
                    %             wave_traveling(h,nevent) = 0; % noisy or standing wave?
                    %         end
                    %         % DOWN_peaks_shank(:,nevent) - DOWN_peaks_shank(1,nevent)
                    %     elseif h==2
                    %         if peaks_latency(h,nevent)>0.001 & sum(sharp_wave_peaks_shank_temp(2:end)-sharp_wave_peaks_shank_temp(1)>0)>=length(LFP(h).best_HPC_shank_id)/2
                    %             wave_traveling(h,nevent) = -1; % posterior to anterior
                    %         elseif peaks_latency(h,nevent)<-0.001 & sum(sharp_wave_peaks_shank_temp(2:end)-sharp_wave_peaks_shank_temp(1)<0)>=length(LFP(h).best_HPC_shank_id)/2
                    %             wave_traveling(h,nevent) = 1; % anterior to posterior
                    %         else
                    %             wave_traveling(h,nevent) = 0; % noisy or standing wave?
                    %         end
                    %     end
                    % end
                    % % nexttile
                    % % hold on;xline(tvec(idx));
                    % % plot(tvec(tidx),zscored_LFP(tidx,:));
                    % % % hold on;xline(median(tvec(tidx(1):tidx(end))))
                    % % xline(sharp_wave_peaks_shank(:,nevent)')
                end

                ripples(probe_no).SWR_peaktimes = sharp_wave_peaks_shank;
                % ripples(probe_no).SWR_latency = peaks_latency;
                ripples(probe_no).SWR_zscore = sharp_wave_zscore_shank;
                % ripples(probe_no).SWR_traveling = wave_traveling;


                sharp_wave_zscore_shank = [];
                sharp_wave_peaks_shank = [];


                zscored_LFP = [];
                zscored_LFP = SO_amplitude_HPC_LFP;

                for nevent = 1:length(ripples(probe_no).onset)
                    %                 ripples(probe_no).peaktimes(nevent)

                    tidx = FindInInterval(tvec,[ripples(probe_no).onset(nevent)-0.05 ripples(probe_no).offset(nevent)]);
                    % tidx = FindInInterval(tvec,[slow_waves(probe_no).ints.UP(nevent,1)-0.3 slow_waves(probe_no).ints.UP(nevent,1)+0.3]);
                    tidx=tidx(1):tidx(end);


                    for nShank=1:length(probe_hemisphere)

                        [~,peak_id] = findpeaks(zscored_LFP(tidx(1):tidx(end),nShank));

                        if ~isempty(peak_id)
                            [~,temp]=min(abs(ripples(probe_no).peaktimes(nevent)-tvec(tidx(peak_id))));
                            sharp_wave_peaks_shank(nShank,nevent) = tvec(tidx(peak_id(temp)));
                            sharp_wave_zscore_shank(nShank,nevent) = zscored_LFP(tidx(peak_id(temp)),nShank);
                        else
                            sharp_wave_peaks_shank(nShank,nevent) = nan;
                            sharp_wave_zscore_shank(nShank,nevent) = nan;
                        end
                    end
                end

                ripples(probe_no).sharp_wave_peaktimes = sharp_wave_peaks_shank;
                ripples(probe_no).sharp_wave_zscore = sharp_wave_zscore_shank;
            end
        end

        disp('Sharp wave ripple analysis started')
        toc

        zscored_LFP = [];



        tic
        disp('UP/DOWN and ripples and spindles phase coupling and amplitude correlation')
        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id+1;

            %%%%% Ripples
            ref_shank = find(slow_waves(probe_no).shank_id == slow_waves(probe_no).shank(slow_waves(probe_no).channel == slow_waves(probe_no).best_channel)...
                &slow_waves(probe_no).probe_hemisphere == probe_no);
            cortex_ref_shank = ref_shank;
            spindles(nprobe).best_channel = cortex_ref_shank;


            shank_id = find(ripples(probe_no).probe_hemisphere == probe_no);
            HPC_ref_shank = shank_id(ripples(probe_no).best_channel);
            ref_shank= HPC_ref_shank;

            % ripple amplitude
            step_s = 0.02;
            win_s = step_s*2;
            winSamples = round(win_s * lfp.samplingRate);
            stepSamples = round(step_s * lfp.samplingRate);
            event_tidx = [];
            event_index = [];
            nSteps = floor((size(ripple_amplitude_LFP,1) - winSamples) / stepSamples) + 1;
            tvec_interp1 = tvec(1):step_s:tvec(end);

            meanAmp = [];
            for ntidx = 1:nSteps
                idx = (ntidx-1)*stepSamples + (1:winSamples);
                meanAmp(:,ntidx) = mean(ripple_amplitude_LFP(idx, :));
            end

            % Grab 20ms downsampled idx when events happned
            for nevent = 1:size(ripples(probe_no).peaktimes,1)
                tidx = FindInInterval(tvec_interp1,[ripples(probe_no).onset(nevent) ripples(probe_no).offset(nevent)]);

                event_tidx = [event_tidx tidx(1):tidx(end)];
                event_index = [event_index nevent*ones(size(tidx(1):tidx(end)))];
            end

            amp_corr = [];
            for nchannel = 1:size(meanAmp,1)
                for mchannel = 1:size(meanAmp,1)
                    amp_corr(nchannel,mchannel) = corr(meanAmp(nchannel,event_tidx)',meanAmp(mchannel,event_tidx)');
                end
            end

            ripples(probe_no).amp_corr = amp_corr;


            %%%%%%%%% Phase Locking Value (PLV) during ripples
            %%%%%%%%%% transition and amplitude cross correlation PER EVENT
            plv_ripples = nan(length(ripples(probe_no).shank_id),length(ripples(probe_no).shank_id), length(ripples(probe_no).peaktimes));
            pd_ripples = nan(length(ripples(probe_no).shank_id),length(ripples(probe_no).shank_id), length(ripples(probe_no).peaktimes));
            xcorr_r_ripples = nan(length(ripples(probe_no).shank_id),length(ripples(probe_no).shank_id), length(ripples(probe_no).peaktimes));
            xcorr_lag_ripples = nan(length(ripples(probe_no).shank_id),length(ripples(probe_no).shank_id), length(ripples(probe_no).peaktimes));

            cortex_speed = nan(max(ripples(probe_no).probe_hemisphere), length(ripples(probe_no).peaktimes));
            HPC_speed = nan(max(ripples(probe_no).probe_hemisphere), length(ripples(probe_no).peaktimes));

            for nevent = 1:size(ripples(probe_no).onset,1)
                tidx = FindInInterval(tvec,[ripples(probe_no).onset(nevent) ripples(probe_no).offset(nevent)]);

                for nchannel = 1:length(ripples(probe_no).shank_id)
                    for mchannel = 1:length(ripples(probe_no).shank_id)
                        amp1 = ripple_amplitude_LFP(tidx(1):tidx(end),nchannel); % reference channel
                        amp2 = ripple_amplitude_LFP(tidx(1):tidx(end),mchannel);
                        [r,lag] = xcorr(zscore(amp1),zscore(amp2),'coeff');
                        [~,idx] = max(r);
                        xcorr_lag_ripples(nchannel,mchannel,nevent) = lag(idx)/lfp.samplingRate; % where is max lag
                        xcorr_r_ripples(nchannel,mchannel,nevent) = r(lag==0); % xcorr at zero lag

                        phi1 = ripple_phase_LFP(tidx(1):tidx(end),nchannel); % reference channel
                        phi2 = ripple_phase_LFP(tidx(1):tidx(end),mchannel);
                        dphi = phi1 - phi2;
                        % Circular mean of phase differences
                        pd_ripples(nchannel,mchannel,nevent) = angle(mean(exp(1i * dphi)));
                        % phase locking values
                        plv_ripples(nchannel,mchannel,nevent) = abs(mean(exp(1i * dphi)));
                    end
                end

                for mprobe = 1:max( ripples(probe_no).probe_hemisphere)
                    HPC_phase = angle(mean(exp(1i * ripple_phase_LFP(tidx(1):tidx(end),(ripples(probe_no).probe_hemisphere == mprobe))), 1));
                    HPC_speed(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(ripple_phase_LFP(tidx(1):tidx(end),HPC_ref_shank))))/...
                        mean(diff(unwrap(HPC_phase))./diff(0.25*ripples(probe_no).shank_id(ripples(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)


                    cortex_phase = angle(mean(exp(1i * ripple_phase_cortex_LFP(tidx(1):tidx(end),(slow_waves(probe_no).probe_hemisphere == mprobe))), 1));
                    cortex_speed(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(ripple_phase_cortex_LFP(tidx(1):tidx(end),cortex_ref_shank))))/...
                        mean(diff(unwrap(cortex_phase))./diff(0.25*slow_waves(probe_no).shank_id(slow_waves(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)
                end

            end

            % right hemisphere
            if max(slow_waves(probe_no).probe_hemisphere)==2
                HPC_speed(2,:) = -HPC_speed(2,:);
                cortex_speed(2,:) = -cortex_speed(2,:);
            end

            ripples(probe_no).plv = plv_ripples;
            ripples(probe_no).pd = pd_ripples;
            ripples(probe_no).xcorr_r = xcorr_r_ripples;
            ripples(probe_no).xcorr_lag = xcorr_lag_ripples;
            ripples(probe_no).cortex_speed = cortex_speed;
            ripples(probe_no).HPC_speed = HPC_speed;



            %%%%% Spindles
            ref_shank = cortex_ref_shank;
            % cortex_ref_shank = ref_shank;
            % spindle amplitude
            step_s = 0.02;
            win_s = step_s*2;
            winSamples = round(win_s * lfp.samplingRate);
            stepSamples = round(step_s * lfp.samplingRate);
            event_tidx = [];
            event_index = [];
            nSteps = floor((size(spindle_amplitude_LFP,1) - winSamples) / stepSamples) + 1;
            tvec_interp1 = tvec(1):step_s:tvec(end);

            meanAmp = [];
            for ntidx = 1:nSteps
                idx = (ntidx-1)*stepSamples + (1:winSamples);
                meanAmp(:,ntidx) = mean(spindle_amplitude_LFP(idx, :));
            end

            if ~isempty(spindles(probe_no).peaktimes)
                % Grab 20ms downsampled idx when events happned
                for nevent = 1:size(spindles(probe_no).peaktimes,1)
                    tidx = FindInInterval(tvec_interp1,[spindles(probe_no).onset(nevent) spindles(probe_no).offset(nevent)]);

                    event_tidx = [event_tidx tidx(1):tidx(end)];
                    event_index = [event_index nevent*ones(size(tidx(1):tidx(end)))];
                end

                amp_corr = [];
                for nchannel = 1:size(meanAmp,1)
                    for mchannel = 1:size(meanAmp,1)
                        amp_corr(nchannel,mchannel) = corr(meanAmp(nchannel,event_tidx)',meanAmp(mchannel,event_tidx)');
                    end
                end

                spindles(probe_no).amp_corr = amp_corr;



                %%%%%%%%% Phase Locking Value (PLV) during spindles
                %%%%%%%%%% transition and amplitude cross correlation PER EVENT
                plv_spindles = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(spindles(probe_no).peaktimes));
                pd_spindles = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(spindles(probe_no).peaktimes));
                xcorr_r_spindles = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(spindles(probe_no).peaktimes));
                xcorr_lag_spindles = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(spindles(probe_no).peaktimes));

                cortex_speed = nan(max(slow_waves(probe_no).probe_hemisphere), length(spindles(probe_no).peaktimes));
                HPC_speed = nan(max(slow_waves(probe_no).probe_hemisphere), length(spindles(probe_no).peaktimes));

                for nevent = 1:size(spindles(probe_no).onset,1)
                    tidx = FindInInterval(tvec,[spindles(probe_no).onset(nevent) spindles(probe_no).offset(nevent)]);

                    for nchannel = 1:length(slow_waves(probe_no).shank_id)
                        for mchannel = 1:length(slow_waves(probe_no).shank_id)
                            amp1 = spindle_amplitude_LFP(tidx(1):tidx(end),nchannel); % reference channel
                            amp2 = spindle_amplitude_LFP(tidx(1):tidx(end),mchannel);
                            [r,lag] = xcorr(zscore(amp1),zscore(amp2),'coeff');
                            [~,idx] = max(r);
                            xcorr_lag_spindles(nchannel,mchannel,nevent) = lag(idx)/lfp.samplingRate; % where is max lag
                            xcorr_r_spindles(nchannel,mchannel,nevent) = r(lag==0); % xcorr at zero lag

                            phi1 = spindle_phase_LFP(tidx(1):tidx(end),nchannel); % reference channel
                            phi2 = spindle_phase_LFP(tidx(1):tidx(end),mchannel);
                            dphi = phi1 - phi2;
                            % Circular mean of phase differences
                            pd_spindles(nchannel,mchannel,nevent) = angle(mean(exp(1i * dphi)));
                            % phase locking values
                            plv_spindles(nchannel,mchannel,nevent) = abs(mean(exp(1i * dphi)));
                        end
                    end

                    for mprobe = 1:max( ripples(probe_no).probe_hemisphere)
                        HPC_phase = angle(mean(exp(1i * spindle_phase_HPC_LFP(tidx(1):tidx(end),(ripples(probe_no).probe_hemisphere == mprobe))), 1));
                        HPC_speed(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(spindle_phase_HPC_LFP(tidx(1):tidx(end),HPC_ref_shank))))/...
                            mean(diff(unwrap(HPC_phase))./diff(0.25*ripples(probe_no).shank_id(ripples(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)


                        cortex_phase = angle(mean(exp(1i * spindle_phase_LFP(tidx(1):tidx(end),(slow_waves(probe_no).probe_hemisphere == mprobe))), 1));
                        cortex_speed(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(spindle_phase_LFP(tidx(1):tidx(end),cortex_ref_shank))))/...
                            mean(diff(unwrap(cortex_phase))./diff(0.25*slow_waves(probe_no).shank_id(slow_waves(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)


                    end

                end

                % right hemisphere
                if max(slow_waves(probe_no).probe_hemisphere)==2
                    HPC_speed(2,:) = -HPC_speed(2,:);
                    cortex_speed(2,:) = -cortex_speed(2,:);
                end

                spindles(probe_no).plv = plv_spindles;
                spindles(probe_no).pd = pd_spindles;
                spindles(probe_no).xcorr_r = xcorr_r_spindles;
                spindles(probe_no).xcorr_lag = xcorr_lag_spindles;

                spindles(probe_no).cortex_speed = cortex_speed;
                spindles(probe_no).HPC_speed = HPC_speed;
            else
                spindles(probe_no).amp_corr = [];
                spindles(probe_no).plv = [];
                spindles(probe_no).pd = [];
                spindles(probe_no).xcorr_r = [];
                spindles(probe_no).xcorr_lag = [];

                spindles(probe_no).cortex_speed = [];
                spindles(probe_no).HPC_speed = [];
            end

            %%%%% UP/DOWN phase locking and amplitude correlation
            %%%%% Ripple and spindle coupling with UP DOWN

            if ~isempty(slow_waves(probe_no).timestamps)
                if size(slow_waves(probe_no).timestamps,1) > 10

                    UP_ints = slow_waves(probe_no).UP_ints;
                    DOWN_ints = slow_waves(probe_no).DOWN_ints;

                    ref_shank = cortex_ref_shank;
                    % Amplitude at D-U transition and U-D transition
                    step_s = 0.02;
                    win_s = step_s*2;
                    winSamples = round(win_s * lfp.samplingRate);
                    stepSamples = round(step_s * lfp.samplingRate);
                    event_tidx = [];
                    event_index = [];
                    nSteps = floor((size(SO_amplitude_LFP,1) - winSamples) / stepSamples) + 1;
                    tvec_interp1 = tvec(1):step_s:tvec(end);

                    meanAmp = [];
                    for ntidx = 1:nSteps
                        idx = (ntidx-1)*stepSamples + (1:winSamples);
                        meanAmp(:,ntidx) = mean(SO_amplitude_LFP(idx, :));
                    end

                    %%%% DOWN-UP transition
                    % Grab 20ms downsampled idx when events happned
                    for nevent = 1:size(UP_ints,1)
                        % Default window: 50 ms before UP onset to 100 ms after
                        t_start = UP_ints(nevent,1) - 0.05;
                        t_end   = UP_ints(nevent,1) + 0.1;

                        % Adjust start if overlapping with previous UP offset
                        if nevent > 1 && t_start < UP_ints(nevent-1,2)
                            t_start = UP_ints(nevent-1,2);% Clip to previous UP offset
                        end

                        % Adjust end if overlapping with next UP onset
                        if nevent < size(UP_ints,1) && t_end > UP_ints(nevent+1,1)
                            t_end = UP_ints(nevent,2); % Clip to current UP offset
                        end

                        % Use final t_start and t_end to extract indices
                        tidx = FindInInterval(tvec_interp1, [t_start t_end]);

                        event_tidx = [event_tidx tidx(1):tidx(end)];
                        event_index = [event_index nevent*ones(size(tidx(1):tidx(end)))];
                    end

                    amp_corr = [];
                    for nchannel = 1:size(meanAmp,1)
                        for mchannel = 1:size(meanAmp,1)
                            amp_corr(nchannel,mchannel) = corr(meanAmp(nchannel,event_tidx)',meanAmp(mchannel,event_tidx)');
                        end
                    end

                    slow_waves(probe_no).amp_corr_DU = amp_corr;

                    %%%% UP-DOWN transition
                    event_tidx = [];
                    event_index = [];
                    % Grab 20ms downsampled idx when events happned
                    for nevent = 1:size(DOWN_ints,1)
                        % Default window: 50 ms before DOWN onset to 100 ms after
                        t_start = DOWN_ints(nevent,1) - 0.05;
                        t_end   = DOWN_ints(nevent,1) + 0.1;

                        % Adjust start if overlapping with previous DOWN offset
                        if nevent > 1 && t_start < DOWN_ints(nevent-1,2)
                            t_start = DOWN_ints(nevent-1,2);% Clip to previous DOWN offset
                        end

                        % Adjust end if overlapping with next DOWN onset
                        if nevent < size(DOWN_ints,1) && t_end > DOWN_ints(nevent+1,1)
                            t_end = DOWN_ints(nevent,2); % Clip to current DOWN offset
                        end

                        % Use final t_start and t_end to extract indices
                        tidx = FindInInterval(tvec_interp1, [t_start t_end]);
                        event_tidx = [event_tidx tidx(1):tidx(end)];
                        event_index = [event_index nevent*ones(size(tidx(1):tidx(end)))];
                    end

                    amp_corr = [];
                    for nchannel = 1:size(meanAmp,1)
                        for mchannel = 1:size(meanAmp,1)
                            amp_corr(nchannel,mchannel) = corr(meanAmp(nchannel,event_tidx)',meanAmp(mchannel,event_tidx)');
                        end
                    end

                    slow_waves(probe_no).amp_corr_UD = amp_corr;


                    % UP and DOWN mean Phase
                    phase_UP = nan(length(slow_waves(probe_no).shank_id), length(UP_ints(:,1)));
                    phase_DOWN = nan(length(slow_waves(probe_no).shank_id), length(DOWN_ints(:,1)));

                    % DOWN
                    event_tidx = [];
                    event_index = [];
                    % Grab 20ms downsampled idx when events happned
                    for nevent = 1:size(DOWN_ints,1)
                        % Default window: 50 ms before DOWN onset to 100 ms after
                        t_start = DOWN_ints(nevent,1);
                        t_end   = DOWN_ints(nevent,1) + 0.1;

                        % Adjust end if overlapping with next DOWN onset
                        if nevent < size(DOWN_ints,1) && t_end > DOWN_ints(nevent+1,1)
                            t_end = DOWN_ints(nevent,2); % Clip to current DOWN offset
                        end

                        % Use final t_start and t_end to extract indices
                        tidx = FindInInterval(tvec, [t_start t_end]);
                        for nchannel = 1:length(slow_waves(probe_no).shank_id)
                            phase_DOWN(nchannel,nevent)=angle(mean(exp(1i*SO_phase_LFP(tidx(1):tidx(end),nchannel)))); % phase
                        end

                        % Use final t_start and t_end to extract indices
                        tidx = FindInInterval(tvec_interp1, [t_start t_end]);
                        event_tidx = [event_tidx tidx(1):tidx(end)];
                        event_index = [event_index nevent*ones(size(tidx(1):tidx(end)))];
                    end

                    amp_corr = [];
                    for nchannel = 1:size(meanAmp,1)
                        for mchannel = 1:size(meanAmp,1)
                            amp_corr(nchannel,mchannel) = corr(meanAmp(nchannel,event_tidx)',meanAmp(mchannel,event_tidx)');
                        end
                    end

                    slow_waves(probe_no).amp_corr_DOWN = amp_corr;
                    slow_waves(probe_no).mean_phase_DOWN = phase_DOWN;

                    % UP
                    event_tidx = [];
                    event_index = [];
                    % Grab 20ms downsampled idx when events happned
                    for nevent = 1:size(UP_ints,1)
                        % Default window: UP onset to 100 ms after
                        t_start = UP_ints(nevent,1);
                        t_end   = UP_ints(nevent,1) + 0.1;

                        % Adjust end if overlapping with next UP onset
                        if nevent < size(UP_ints,1) && t_end > UP_ints(nevent+1,1)
                            t_end = UP_ints(nevent,2); % Clip to current UP offset
                        end

                        % Use final t_start and t_end to extract indices
                        tidx = FindInInterval(tvec, [t_start t_end]);
                        for nchannel = 1:length(slow_waves(probe_no).shank_id)
                            phase_UP(nchannel,nevent)=angle(mean(exp(1i*SO_phase_LFP(tidx(1):tidx(end),nchannel)))); % phase
                        end

                        tidx = FindInInterval(tvec_interp1, [t_start t_end]);
                        event_tidx = [event_tidx tidx(1):tidx(end)];
                        event_index = [event_index nevent*ones(size(tidx(1):tidx(end)))];
                    end

                    amp_corr = [];
                    for nchannel = 1:size(meanAmp,1)
                        for mchannel = 1:size(meanAmp,1)
                            amp_corr(nchannel,mchannel) = corr(meanAmp(nchannel,event_tidx)',meanAmp(mchannel,event_tidx)');
                        end
                    end
                    slow_waves(probe_no).amp_corr_UP = amp_corr;
                    slow_waves(probe_no).mean_phase_UP = phase_UP;
                    % polarhistogram(phase_DOWN(1,:));hold on;polarhistogram(phase_UP(1,:))

                    %%%%%%%%% Phase Locking Value (PLV) at D-U transition and U-D
                    %%%%%%%%%% transition and amplitude cross correlation PER EVENT
                    plv_DU = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(UP_ints(:,1)));
                    plv_UD = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(DOWN_ints(:,1)));
                    pd_DU = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(UP_ints(:,1)));
                    pd_UD = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(DOWN_ints(:,1)));
                    xcorr_r_DU = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(UP_ints(:,1)));
                    xcorr_r_UD = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(DOWN_ints(:,1)));
                    xcorr_lag_DU = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(UP_ints(:,1)));
                    xcorr_lag_UD = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(DOWN_ints(:,1)));

                    cortex_speed_UD = nan(max(slow_waves(probe_no).probe_hemisphere), length(DOWN_ints(:,1)));
                    HPC_speed_UD = nan(max(slow_waves(probe_no).probe_hemisphere), length(DOWN_ints(:,1)));
                    cortex_speed_DU = nan(max(slow_waves(probe_no).probe_hemisphere), length(UP_ints(:,1)));
                    HPC_speed_DU = nan(max(slow_waves(probe_no).probe_hemisphere), length(UP_ints(:,1)));

                    for nevent = 1:size(UP_ints,1)
                        % Default window: 50 ms before UP onset to 100 ms after
                        t_start = UP_ints(nevent,1) - 0.05;
                        t_end   = UP_ints(nevent,1) + 0.1;

                        % Adjust start if overlapping with previous UP offset
                        if nevent > 1 && t_start < UP_ints(nevent-1,2)
                            t_start = UP_ints(nevent-1,2);% Clip to previous UP offset
                        end

                        % Adjust end if overlapping with next UP onset
                        if nevent < size(UP_ints,1) && t_end > UP_ints(nevent+1,1)
                            t_end = UP_ints(nevent,2); % Clip to current UP offset
                        end

                        % Use final t_start and t_end to extract indices
                        tidx = FindInInterval(tvec, [t_start t_end]);

                        for nchannel = 1:length(slow_waves(probe_no).shank_id)
                            for mchannel = 1:length(slow_waves(probe_no).shank_id)
                                amp1 = SO_amplitude_LFP(tidx(1):tidx(end),nchannel); % reference channel
                                amp2 = SO_amplitude_LFP(tidx(1):tidx(end),mchannel);
                                [r,lag] = xcorr(zscore(amp1),zscore(amp2),'coeff');
                                [~,idx] = max(r);
                                xcorr_lag_DU(nchannel,mchannel,nevent) = lag(idx)/lfp.samplingRate; % where is max lag
                                xcorr_r_DU(nchannel,mchannel,nevent) = r(lag==0); % xcorr at zero lag

                                phi1 = SO_phase_LFP(tidx(1):tidx(end),nchannel); % reference channel
                                phi2 = SO_phase_LFP(tidx(1):tidx(end),mchannel);
                                dphi = phi1 - phi2;
                                % Circular mean of phase differences
                                pd_DU(nchannel,mchannel,nevent) = angle(mean(exp(1i * dphi)));
                                % phase locking values
                                plv_DU(nchannel,mchannel,nevent) = abs(mean(exp(1i * dphi)));
                            end
                        end


                        for mprobe = 1:max(slow_waves(probe_no).probe_hemisphere)
                            HPC_phase = angle(mean(exp(1i * SO_phase_HPC_LFP(tidx(1):tidx(end),(ripples(probe_no).probe_hemisphere == mprobe))), 1));
                            HPC_speed_DU(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(SO_phase_HPC_LFP(tidx(1):tidx(end),HPC_ref_shank))))/...
                                mean(diff(unwrap(HPC_phase))./diff(0.25*ripples(probe_no).shank_id(ripples(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)

                            cortex_phase = angle(mean(exp(1i * SO_phase_LFP(tidx(1):tidx(end),(slow_waves(probe_no).probe_hemisphere == mprobe))), 1));
                            cortex_speed_DU(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(SO_phase_LFP(tidx(1):tidx(end),cortex_ref_shank))))/...
                                mean(diff(unwrap(cortex_phase))./diff(0.25*slow_waves(probe_no).shank_id(slow_waves(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)
                        end

                    end

                    % right hemisphere
                    if max(slow_waves(probe_no).probe_hemisphere)==2
                        HPC_speed_DU(2,:) = -HPC_speed_DU(2,:);
                        cortex_speed_DU(2,:) = -cortex_speed_DU(2,:);
                    end

                    for nevent = 1:size(DOWN_ints,1)
                        % Default window: 50 ms before DOWN onset to 100 ms after
                        t_start = DOWN_ints(nevent,1) - 0.05;
                        t_end   = DOWN_ints(nevent,1) + 0.1;

                        % Adjust start if overlapping with previous DOWN offset
                        if nevent > 1 && t_start < DOWN_ints(nevent-1,2)
                            t_start = DOWN_ints(nevent-1,2);% Clip to previous DOWN offset
                        end

                        % Adjust end if overlapping with next DOWN onset
                        if nevent < size(DOWN_ints,1) && t_end > DOWN_ints(nevent+1,1)
                            t_end = DOWN_ints(nevent,2); % Clip to current DOWN offset
                        end

                        % Use final t_start and t_end to extract indices
                        tidx = FindInInterval(tvec, [t_start t_end]);

                        for nchannel = 1:length(slow_waves(probe_no).shank_id)
                            for mchannel = 1:length(slow_waves(probe_no).shank_id)
                                amp1 = SO_amplitude_LFP(tidx(1):tidx(end),nchannel); % reference channel
                                amp2 = SO_amplitude_LFP(tidx(1):tidx(end),mchannel);
                                [r,lag] = xcorr(zscore(amp1),zscore(amp2),'coeff');
                                [~,idx] = max(r);
                                xcorr_lag_UD(nchannel,mchannel,nevent) = lag(idx)/lfp.samplingRate; % where is max lag
                                xcorr_r_UD(nchannel,mchannel,nevent) = r(lag==0); % xcorr at zero lag


                                phi1 = SO_phase_LFP(tidx(1):tidx(end),nchannel); % reference channel
                                phi2 = SO_phase_LFP(tidx(1):tidx(end),mchannel);
                                dphi = phi1 - phi2;
                                % Circular mean of phase differences
                                pd_UD(nchannel,mchannel,nevent) = angle(mean(exp(1i * dphi)));
                                % phase locking values
                                plv_UD(nchannel,mchannel,nevent) = abs(mean(exp(1i * dphi)));
                            end
                        end


                        for mprobe = 1:max(slow_waves(probe_no).probe_hemisphere)
                            HPC_phase = angle(mean(exp(1i * SO_phase_HPC_LFP(tidx(1):tidx(end),(ripples(probe_no).probe_hemisphere == mprobe))), 1));
                            HPC_speed_UD(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(SO_phase_HPC_LFP(tidx(1):tidx(end),HPC_ref_shank))))/...
                                mean(diff(unwrap(HPC_phase))./diff(0.25*ripples(probe_no).shank_id(ripples(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)


                            cortex_phase = angle(mean(exp(1i * SO_phase_LFP(tidx(1):tidx(end),(slow_waves(probe_no).probe_hemisphere == mprobe))), 1));
                            cortex_speed_UD(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(SO_phase_LFP(tidx(1):tidx(end),cortex_ref_shank))))/...
                                mean(diff(unwrap(cortex_phase))./diff(0.25*slow_waves(probe_no).shank_id(slow_waves(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)
                        end
                    end

                    % right hemisphere
                    if max(slow_waves(probe_no).probe_hemisphere)==2
                        HPC_speed_UD(2,:) = -HPC_speed_UD(2,:);
                        cortex_speed_UD(2,:) = -cortex_speed_UD(2,:);
                    end

                    slow_waves(probe_no).xcorr_lag_UD = xcorr_lag_UD;
                    slow_waves(probe_no).xcorr_lag_DU = xcorr_lag_DU;

                    slow_waves(probe_no).xcorr_r_UD = xcorr_r_UD;
                    slow_waves(probe_no).xcorr_r_DU = xcorr_r_DU;

                    slow_waves(probe_no).plv_UD = plv_UD;
                    slow_waves(probe_no).pd_UD = pd_UD;

                    slow_waves(probe_no).plv_DU = plv_DU;
                    slow_waves(probe_no).pd_DU = pd_DU;

                    slow_waves(probe_no).cortex_speed_UD = cortex_speed_UD;
                    slow_waves(probe_no).HPC_speed_UD = HPC_speed_UD;

                    slow_waves(probe_no).cortex_speed_DU = cortex_speed_DU;
                    slow_waves(probe_no).HPC_speed_DU = HPC_speed_DU;

                    % phase and amplitude
                    SO_phase_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    SO_phase_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));

                    SO_amplitude_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    SO_amplitude_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));

                    spindle_phase_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    spindle_phase_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).onset));

                    spindle_amplitude_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    spindle_amplitude_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).onset));

                    SO_phase_spindle_onset = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
                    SO_phase_spindle_peaktime = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));

                    SO_amplitude_spindle_onset = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
                    SO_amplitude_spindle_peaktime = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));

                    % SO_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    % spindle_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    % SO_spindles_coupling = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
                    ripple_peak_amplitude = nan(length(ripples(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    ripple_onset_amplitude = nan(length(ripples(probe_no).shank_id), length(ripples(probe_no).onset));
                    spindle_peak_amplitude = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).peaktimes));
                    spindle_onset_amplitude = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));


                    gamma_phase_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    gamma_phase_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    gamma_amplitude_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    gamma_amplitude_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));

                    gamma_phase_spindle_peaktime = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
                    gamma_phase_spindle_onset = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
                    gamma_amplitude_spindle_peaktime = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
                    gamma_amplitude_spindle_onset = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));


                    for nevent = 1:length(ripples(probe_no).onset)
                        [~,tidx]=min(abs(tvec-ripples(probe_no).peaktimes(nevent)));
                        
                        % peak ampltiude and phase
                        SO_phase_ripple_peaktime(:,nevent) = SO_phase_LFP(tidx,:);
                        spindle_phase_ripple_peaktime(:,nevent) = spindle_phase_LFP(tidx,:);
                        SO_amplitude_ripple_peaktime(:,nevent) = SO_amplitude_LFP(tidx,:);
                        spindle_amplitude_ripple_peaktime(:,nevent) = spindle_amplitude_LFP(tidx,:);
                        ripple_peak_amplitude(:,nevent) = ripple_amplitude_LFP(tidx,:);
                        gamma_phase_ripple_peaktime(:, nevent) = gamma_phase_LFP(tidx, :);
                        gamma_amplitude_ripple_peaktime(:, nevent) = gamma_amplitude_LFP(tidx, :);


                        [~,tidx]=min(abs(tvec-ripples(probe_no).onset(nevent)));
                        % onset amplitude and phase
                        SO_phase_ripple_onset(:,nevent) = SO_phase_LFP(tidx,:);
                        spindle_phase_ripple_onset(:,nevent) = spindle_phase_LFP(tidx,:);
                        SO_amplitude_ripple_onset(:,nevent) = SO_amplitude_LFP(tidx,:);
                        spindle_amplitude_ripple_onset(:,nevent) = spindle_amplitude_LFP(tidx,:);
                        ripple_onset_amplitude(:,nevent) = ripple_amplitude_LFP(tidx,:);
                        gamma_phase_ripple_onset(:, nevent) = gamma_phase_LFP(tidx, :);
                        gamma_amplitude_ripple_onset(:, nevent) = gamma_amplitude_LFP(tidx, :);


                        time_ranges = ripples(probe_no).onset(nevent)-(-1:0.1:1);
                        [~,tidx]=min(abs(tvec-ripples(probe_no).onset(nevent)));
                        % onset amplitude and phase
                        SO_phase_ripple_onset(:,nevent) = SO_phase_LFP(tidx,:);
                        spindle_phase_ripple_onset(:,nevent) = spindle_phase_LFP(tidx,:);
                        SO_amplitude_ripple_onset(:,nevent) = SO_amplitude_LFP(tidx,:);
                        spindle_amplitude_ripple_onset(:,nevent) = spindle_amplitude_LFP(tidx,:);
                        ripple_onset_amplitude(:,nevent) = ripple_amplitude_LFP(tidx,:);
                        gamma_phase_ripple_onset(:, nevent) = gamma_phase_LFP(tidx, :);
                        gamma_amplitude_ripple_onset(:, nevent) = gamma_amplitude_LFP(tidx, :);

                    end

                    ripples(probe_no).SO_phase_ripple_peaktime = SO_phase_ripple_peaktime;
                    ripples(probe_no).spindle_phase_ripple_peaktime = spindle_phase_ripple_peaktime;
                    ripples(probe_no).SO_amplitude_ripple_peaktime = SO_amplitude_ripple_peaktime;
                    ripples(probe_no).spindle_amplitude_ripple_peaktime = spindle_amplitude_ripple_peaktime;
                    ripples(probe_no).ripple_peak_amplitude = ripple_peak_amplitude;

                    ripples(probe_no).SO_phase_ripple_onset = SO_phase_ripple_onset;
                    ripples(probe_no).spindle_phase_ripple_onset = spindle_phase_ripple_onset;
                    ripples(probe_no).SO_amplitude_ripple_onset = SO_amplitude_ripple_onset;
                    ripples(probe_no).spindle_amplitude_ripple_onset = spindle_amplitude_ripple_onset;
                    ripples(probe_no).ripple_onset_amplitude = ripple_onset_amplitude;

                    ripples(probe_no).gamma_phase_ripple_peaktime = gamma_phase_ripple_peaktime;
                    ripples(probe_no).gamma_amplitude_ripple_peaktime = gamma_amplitude_ripple_peaktime;
                    ripples(probe_no).gamma_phase_ripple_onset = gamma_phase_ripple_onset;
                    ripples(probe_no).gamma_amplitude_ripple_onset = gamma_amplitude_ripple_onset;

                    if ~isempty(spindles(probe_no).onset)
                        for nevent = 1:length(spindles(probe_no).onset)
                            [~,tidx]=min(abs(tvec-spindles(probe_no).peaktimes(nevent)));
                            % peak ampltiude and phase
                            SO_phase_spindle_peaktime(:,nevent) = SO_phase_LFP(tidx,:);
                            SO_amplitude_spindle_peaktime(:,nevent) = SO_amplitude_LFP(tidx,:);
                            spindle_peak_amplitude(:,nevent) = spindle_amplitude_LFP(tidx,:);
                            gamma_phase_spindle_peaktime(:, nevent) = gamma_phase_LFP(tidx, :);
                            gamma_amplitude_spindle_peaktime(:, nevent) = gamma_amplitude_LFP(tidx, :);


                            [~,tidx]=min(abs(tvec-spindles(probe_no).onset(nevent)));
                            % onset amplitude and phase
                            SO_phase_spindle_onset(:,nevent) = SO_phase_LFP(tidx,:);
                            SO_amplitude_spindle_onset(:,nevent) = SO_amplitude_LFP(tidx,:);
                            spindle_onset_amplitude(:,nevent) = spindle_amplitude_LFP(tidx,:);
                            gamma_phase_spindle_onset(:, nevent) = gamma_phase_LFP(tidx, :);
                            gamma_amplitude_spindle_onset(:, nevent) = gamma_amplitude_LFP(tidx, :);
                        end

                        spindles(probe_no).SO_phase_spindle_peaktime = SO_phase_spindle_peaktime;
                        spindles(probe_no).SO_amplitude_spindle_peaktime = SO_amplitude_spindle_peaktime;
                        spindles(probe_no).spindle_peak_amplitude = spindle_peak_amplitude;

                        spindles(probe_no).SO_phase_spindle_onset = SO_phase_spindle_onset;
                        spindles(probe_no).SO_amplitude_spindle_onset = SO_amplitude_spindle_onset;
                        spindles(probe_no).spindle_onset_amplitude = spindle_onset_amplitude;

                        spindles(probe_no).gamma_phase_spindle_peaktime = gamma_phase_spindle_peaktime;
                        spindles(probe_no).gamma_amplitude_spindle_peaktime = gamma_amplitude_spindle_peaktime;
                        spindles(probe_no).gamma_phase_spindle_onset = gamma_phase_spindle_onset;
                        spindles(probe_no).gamma_amplitude_spindle_onset = gamma_amplitude_spindle_onset;
                    end
                else
                    spindles(probe_no).SO_phase_spindle_peaktime = [];
                    spindles(probe_no).spindle_peak_amplitude = [];
                    spindles(probe_no).SO_amplitude_spindle_peaktime = [];

                    spindles(probe_no).SO_phase_spindle_onset = [];
                    spindles(probe_no).spindle_onset_amplitude = [];
                    spindles(probe_no).SO_amplitude_spindle_onset = [];
                end
            end



        end

        %%%%%%%%%%%%%%%%%%%%%%%%% Grab phase and amplitude per spike

        % From cell structure back to spike times and spike id
        session_clusters.spike_id = vertcat(session_clusters.spike_id{:});
        session_clusters.spike_times = vertcat(session_clusters.spike_times{:});
        [session_clusters.spike_times, index] = sort(session_clusters.spike_times);
        session_clusters.spike_id = session_clusters.spike_id(index);

        % Setup
        nSpikes = length(session_clusters.spike_times);
        n_channels = size(SO_amplitude_LFP, 2);
        ripple_channels = size(ripple_amplitude_LFP, 2);

        chunk_size = 5e5;  % 500,000 spikes per chunk
        n_chunks = ceil(nSpikes / chunk_size);

        % Preallocate cell arrays for chunks
        SO_phase_chunks = cell(n_chunks, 1);
        SO_amplitude_chunks = cell(n_chunks, 1);
        spindle_phase_chunks = cell(n_chunks, 1);
        spindle_amplitude_chunks = cell(n_chunks, 1);
        ripple_amplitude_chunks = cell(n_chunks, 1);
        ripple_phase_chunks = cell(n_chunks, 1);        
        theta_phase_chunks = cell(n_chunks, 1);          
        theta_amplitude_chunks = cell(n_chunks, 1);     

        % Loop over chunks
        for c = 1:n_chunks
            idx_start = (c-1)*chunk_size + 1;
            idx_end = min(c*chunk_size, nSpikes);
            idx_range = idx_start:idx_end;

            % Spike times and nearest LFP index
            spike_times_chunk = session_clusters.spike_times(idx_range);
            spike_idx = interp1(tvec, 1:length(tvec), spike_times_chunk, 'nearest', 'extrap');

            % Extract and convert to single
            SO_phase_chunks{c} = single(SO_phase_LFP(spike_idx, :));
            SO_amplitude_chunks{c} = single(SO_amplitude_LFP(spike_idx, :));
            spindle_phase_chunks{c} = single(spindle_phase_LFP(spike_idx, :));
            spindle_amplitude_chunks{c} = single(spindle_amplitude_LFP(spike_idx, :));
            ripple_amplitude_chunks{c} = single(ripple_amplitude_LFP(spike_idx, :));
            ripple_phase_chunks{c} = single(ripple_phase_LFP(spike_idx, :));           
            theta_phase_chunks{c} = single(theta_phase_HPC_LFP(spike_idx, :));        
            theta_amplitude_chunks{c} = single(theta_amplitude_HPC_LFP(spike_idx, :)); 

            fprintf('Processed chunk %d of %d (%d spikes)\n', c, n_chunks, length(idx_range));
        end

        % Merge chunks
        spike_phase_amplitude.SO_phase = vertcat(SO_phase_chunks{:});
        spike_phase_amplitude.SO_amplitude = vertcat(SO_amplitude_chunks{:});
        spike_phase_amplitude.spindle_phase = vertcat(spindle_phase_chunks{:});
        spike_phase_amplitude.spindle_amplitude = vertcat(spindle_amplitude_chunks{:});
        spike_phase_amplitude.ripple_amplitude = vertcat(ripple_amplitude_chunks{:});
        spike_phase_amplitude.ripple_phase = vertcat(ripple_phase_chunks{:});          
        spike_phase_amplitude.theta_phase = vertcat(theta_phase_chunks{:});             
        spike_phase_amplitude.theta_amplitude = vertcat(theta_amplitude_chunks{:});    

        % Clear chunk variables to free memory
        clear SO_phase_chunks SO_amplitude_chunks spindle_phase_chunks ...
            spindle_amplitude_chunks ripple_amplitude_chunks ripple_phase_chunks ...
            theta_phase_chunks theta_amplitude_chunks




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
            save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_spindle_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'spindles');
            % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_slow_wave_markov_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'slow_waves_markov');
            save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_slow_wave_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'slow_waves');

            save(fullfile(options.ANALYSIS_DATAPATH,sprintf('spike_phase_amplitude%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'spike_phase_amplitude','-v7,3');
        else
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'),'slow_waves');
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'),'spindles');
            % save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_markov_events.mat'),'slow_waves_markov');
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'),'ripples');

            save(fullfile(options.ANALYSIS_DATAPATH,'spike_phase_amplitude.mat'), 'spike_phase_amplitude', '-v7.3');
        end
    end
end




%% LFP time frequency analysis (wavelet method)
% pyversion('C:\Users\masahiro.takigawa\.conda\envs\fooof\python')
% pyversion('C:\Users\masah\anaconda3\envs\fooof\python')

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar
% experiment_info=experiment_info([4 5 6 ]);
% experiment_info=experiment_info([4 5 6 18 19 21 34 35 44 45 58 59 60 71]);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
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
        % PSD slope quantification using fooof ()
        tvec = LFP(1).tvec;
        SR = round(1/mean(diff(tvec)));
        nfft_seconds= 2;
        nfft = 2^(nextpow2(SR*nfft_seconds));
        win  = hanning(nfft);

        %%%%%%%%%%%% Cortical wave direction during DOWN state peak

        filterparms.deltafilter = [0.5 4];%heuristically defined.  room for improvement here.
        filterparms.thetafilter = [4 12];%heuristically defined.  room for improvement here.
        filterparms.spindlesfilter = [9 17];%heuristically defined.  room for improvement here.
        filterparms.ripplesfilter = [125 300];%heuristically defined.  room for improvement here.
        filterparms.gammafilter = [100 400];
        filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
        filterparms.gammanormwin = 20; %window for gamma normalization (s)

        %         for nprobe = 1:length(session_info(n).probe)
        nprobe = 1;
        probe_no = session_info(n).probe(nprobe).probe_id+1;

        if length(LFP(1).best_HPC) < length(LFP(2).best_HPC)
            time_idx = 1:length(LFP(1).best_HPC);
        elseif length(LFP(1).best_HPC) > length(LFP(2).best_HPC)

            time_idx = 1:length(LFP(2).best_HPC);
        else
            time_idx = 1:length(LFP(1).best_HPC) ;
        end

        lfp.samplingRate = round(1/mean(diff(LFP(probe_no).tvec)));
        lfp.timestamps = LFP(probe_no).tvec(time_idx);
        tvec = LFP(probe_no).tvec(time_idx);

        lfp_V1 = lfp;
        lfp_HPC = lfp;
        lfp_V1.data=[];
        lfp_HPC.data=[];

        if isfield(LFP(probe_no),'best_HPC') % if exist best HPC channels
            probe_id = [];

            if length(LFP)==1
                lfp_HPC.data= [LFP(probe_no).best_HPC(:,time_idx)'];
                lfp_V1.data= [LFP(probe_no).best_V1(:,time_idx)'];
            else
                lfp_HPC.data= [LFP(1).best_HPC(:,time_idx)' LFP(2).best_HPC(:,time_idx)'];
                lfp_V1.data= [LFP(1).average_V1(:,time_idx)' LFP(2).average_V1(:,time_idx)'];
            end
        end
        %         end


        %%% get 
        HPC_ref_shank=[];
        ref_shank = [];
        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id+1;

            ref_shank = find(slow_waves(probe_no).shank_id == slow_waves(probe_no).shank(slow_waves(probe_no).channel == slow_waves(probe_no).best_channel)...
                &slow_waves(probe_no).probe_hemisphere == probe_no);
            cortex_ref_shank(probe_no) = ref_shank;

            shank_id = find(ripples(probe_no).probe_hemisphere == probe_no);
            HPC_ref_shank(probe_no) = shank_id(ripples(probe_no).best_channel);
        end

      
        disp('Peri event TF and PPC analysis')

        % Parameters
        win_full = [-2.5 2.5];
        win_save = [-2 2];
        baseline_win = [-2 -1.5];
        downsample_factor = 25;
        fs = lfp_V1.samplingRate/downsample_factor; % 20ms bins

        win_full_samp = round(win_full * fs);
        win_save_samp = round(win_save * fs);
        baseline_samp = round((baseline_win - win_full(1)) * fs);
        timebin = abs(diff(win_save_samp)) + 1;
        tvec = win_save(1):1/fs:win_save(end);

        freq_band = [1 300];
        freqs = logspace(log10(freq_band(1)), log10(freq_band(2)), 100);
        nCycles = logspace(log10(3), log10(20), 100);
        % nCycles = (max(5, freqs * 0.1));
        % nCycles(nCycles > 20) = 20;
%         nCycles = round(min(20, freqs * 0.1));

        event_types = {'ripples','spindles','UP_ints','DOWN_ints'};

        % Output structures
        TF_amp_V1 = struct(); TF_phase_V1 = struct();
        TF_amp_HPC = struct(); TF_phase_HPC = struct();
        PPC_V1_HPC = struct(); PPC_V1 = struct(); PPC_HPC = struct();

        for nprobe = 1:length(session_info(n).probe)
            tic
            disp(['Processing probe ' num2str(nprobe)]);

            % Interpolate event onsets to sample indices
            ripple_sample_idx = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), ripples(nprobe).onset, 'nearest', 'extrap'));
            spindle_sample_idx = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), spindles(nprobe).onset, 'nearest', 'extrap'));
            up_sample_idx      = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), slow_waves(nprobe).UP_ints(:,1), 'nearest', 'extrap'));
            down_sample_idx    = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), slow_waves(nprobe).DOWN_ints(:,1), 'nearest', 'extrap'));

            for eidx = 1:length(event_types)
                event_name = event_types{eidx};
                switch event_name
                    case 'ripples'
                        onsets = ripple_sample_idx;
                    case 'spindles'
                        onsets = spindle_sample_idx;
                    case 'UP_ints'
                        onsets = up_sample_idx;
                    case 'DOWN_ints'
                        onsets = down_sample_idx;
                end

                nevents = length(onsets);

                % Preallocate outputs
                TF_amp_V1(nprobe).(event_name)   = single(nan(2, length(freqs), timebin, nevents));
                TF_phase_V1(nprobe).(event_name) = single(nan(2, length(freqs), timebin, nevents));
                TF_amp_HPC(nprobe).(event_name)  = single(nan(2, length(freqs), timebin, nevents));
                TF_phase_HPC(nprobe).(event_name)= single(nan(2, length(freqs), timebin, nevents));

%                 PPC_V1_HPC(nprobe).(event_name)  = single(nan(2, 2,length(freqs), timebin, nevents));
%                 PPC_V1(nprobe).(event_name)      = single(nan(length(freqs),timebin, nevents));
%                 PPC_HPC(nprobe).(event_name)     = single(nan(length(freqs),timebin, nevents));


                for nevent = 1:nevents
                    t0 = onsets(nevent);
                    idx_full = t0 + downsample_factor*win_full_samp(1):t0 + win_full_samp(2)*downsample_factor;
                    if idx_full(1) < 1 || idx_full(end) > size(lfp_V1.data,1), continue; end

                    for vhemi = 1:2
                        dat = lfp_V1.data(idx_full, cortex_ref_shank(vhemi));
                        wavespec = bz_WaveSpec(dat, 'frange', freq_band, ...
                            'ncyc', nCycles, 'nfreqs', 100,'samplingRate',fs,'downsampleout', downsample_factor);
                        amp_all = abs(wavespec.data);
                        phase_all = angle(wavespec.data);

                        base = mean(amp_all(baseline_samp(1):baseline_samp(2), :), 1);
                        amp_norm = 10 * log10(bsxfun(@rdivide, amp_all, base));

                        TF_amp_V1(nprobe).(event_name)(vhemi, :, :, nevent) = single(amp_norm(win_save_samp(1)-win_full_samp(1)+1:win_save_samp(2)-win_full_samp(1)+1, :)');
                        TF_phase_V1(nprobe).(event_name)(vhemi, :, :, nevent) = single(phase_all(win_save_samp(1)-win_full_samp(1)+1:win_save_samp(2)-win_full_samp(1)+1, :)');
                    end

                    for hhemi = 1:2
                        dat = lfp_HPC.data(idx_full, HPC_ref_shank(hhemi));
                        wavespec = bz_WaveSpec(dat, 'frange', freq_band, ...
                            'ncyc', nCycles, 'nfreqs', 100,'samplingRate',fs,'downsampleout', downsample_factor);
                        amp_all = abs(wavespec.data);
                        phase_all = angle(wavespec.data);

                        base = mean(amp_all(baseline_samp(1):baseline_samp(2), :), 1);
                        amp_norm = 10 * log10(bsxfun(@rdivide, amp_all, base));

                        TF_amp_HPC(nprobe).(event_name)(hhemi, :, :, nevent) = single(amp_norm(win_save_samp(1)-win_full_samp(1)+1:win_save_samp(2)-win_full_samp(1)+1, :)');
                        TF_phase_HPC(nprobe).(event_name)(hhemi, :, :, nevent) = single(phase_all(win_save_samp(1)-win_full_samp(1)+1:win_save_samp(2)-win_full_samp(1)+1, :)');
                    end

%                     % Compute PPC across time and frequency
%                     for vhemi = 1:2
%                         for hhemi = 1:2
%                             v1p = squeeze(TF_phase_V1(nprobe).(event_name)(vhemi, :, :, nevent));   % [freq x time]
%                             hpcp = squeeze(TF_phase_HPC(nprobe).(event_name)(hhemi, :, :, nevent)); % [freq x time]
% 
%                             % Compute phase difference: still [freq x time]
%                             phase_diff = exp(1i * (v1p - hpcp));
% 
%                             % Assign directly
%                             PPC_V1_HPC(nprobe).(event_name)(vhemi, hhemi, :, :, nevent) =phase_diff;  % [time x freq]
%                         end
%                     end
% 
%                     % Cross-hemispheric PPC
%                     v1p_L = squeeze(TF_phase_V1(nprobe).(event_name)(1, :, :, nevent));
%                     v1p_R = squeeze(TF_phase_V1(nprobe).(event_name)(2, :, :, nevent));
%                     hpcp_L = squeeze(TF_phase_HPC(nprobe).(event_name)(1, :, :, nevent));
%                     hpcp_R = squeeze(TF_phase_HPC(nprobe).(event_name)(2, :, :, nevent));
% 
%                     PPC_V1(nprobe).(event_name)(:,:, nevent) = exp(1i * (v1p_L - v1p_R));   % [time x freq]
%                     PPC_HPC(nprobe).(event_name)(:,:, nevent) =exp(1i * (hpcp_L - hpcp_R));
                end
            end
            toc
        end

        if contains(stimulus_name{n},'Masa2tracks')
         
            save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_TF_amplitude%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'TF_amp_HPC','TF_amp_V1','-v7.3');
            save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_TF_phase%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'TF_phase_V1','TF_phase_HPC','-v7.3');

%             save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PPC%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPC_V1','PPC_HPC','PPC_V1_HPC','-v7.3');
        else


            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_TF_amplitude.mat'),'TF_amp_HPC','TF_amp_V1','-v7.3');
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_TF_phase.mat'),'TF_phase_V1','TF_phase_HPC','-v7.3');
            % save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_markov_events.mat'),'slow_waves_markov');
%             save(fullfile(options.ANALYSIS_DATAPATH,'extracted_PPC.mat'),'PPC_V1','PPC_HPC','PPC_V1_HPC','-v7.3');

 
        end
    end

end




%% LFP time frequency analysis (Chronux)
% pyversion('C:\Users\masahiro.takigawa\.conda\envs\fooof\python')
% pyversion('C:\Users\masah\anaconda3\envs\fooof\python')

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar
% experiment_info=experiment_info([4 5 6 ]);
% experiment_info=experiment_info([4 5 6 18 19 21 34 35 44 45 58 59 60 71]);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
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
        % PSD slope quantification using fooof ()
        tvec = LFP(1).tvec;
        SR = round(1/mean(diff(tvec)));
        nfft_seconds= 2;
        nfft = 2^(nextpow2(SR*nfft_seconds));
        win  = hanning(nfft);

        %%%%%%%%%%%% Cortical wave direction during DOWN state peak

        filterparms.deltafilter = [0.5 4];%heuristically defined.  room for improvement here.
        filterparms.thetafilter = [4 12];%heuristically defined.  room for improvement here.
        filterparms.spindlesfilter = [9 17];%heuristically defined.  room for improvement here.
        filterparms.ripplesfilter = [125 300];%heuristically defined.  room for improvement here.
        filterparms.gammafilter = [100 400];
        filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
        filterparms.gammanormwin = 20; %window for gamma normalization (s)

        %         for nprobe = 1:length(session_info(n).probe)
        nprobe = 1;
        probe_no = session_info(n).probe(nprobe).probe_id+1;

        if length(LFP(1).best_HPC) < length(LFP(2).best_HPC)
            time_idx = 1:length(LFP(1).best_HPC);
        elseif length(LFP(1).best_HPC) > length(LFP(2).best_HPC)

            time_idx = 1:length(LFP(2).best_HPC);
        else
            time_idx = 1:length(LFP(1).best_HPC) ;
        end

        lfp.samplingRate = round(1/mean(diff(LFP(probe_no).tvec)));
        lfp.timestamps = LFP(probe_no).tvec(time_idx);
        tvec = LFP(probe_no).tvec(time_idx);

        lfp_V1 = lfp;
        lfp_HPC = lfp;
        lfp_V1.data=[];
        lfp_HPC.data=[];

        if isfield(LFP(probe_no),'best_HPC') % if exist best HPC channels
            probe_id = [];

            if length(LFP)==1
                lfp_HPC.data= [LFP(probe_no).best_HPC(:,time_idx)'];
                lfp_V1.data= [LFP(probe_no).best_V1(:,time_idx)'];
            else
                lfp_HPC.data= [LFP(1).best_HPC(:,time_idx)' LFP(2).best_HPC(:,time_idx)'];
                lfp_V1.data= [LFP(1).average_V1(:,time_idx)' LFP(2).average_V1(:,time_idx)'];
            end
        end
        %         end


        %%% get 
        HPC_ref_shank=[];
        ref_shank = [];
        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id+1;

            ref_shank = find(slow_waves(probe_no).shank_id == slow_waves(probe_no).shank(slow_waves(probe_no).channel == slow_waves(probe_no).best_channel)...
                &slow_waves(probe_no).probe_hemisphere == probe_no);
            cortex_ref_shank(probe_no) = ref_shank;

            shank_id = find(ripples(probe_no).probe_hemisphere == probe_no);
            HPC_ref_shank(probe_no) = shank_id(ripples(probe_no).best_channel);
        end

      
        disp('Peri event TF and PPC analysis')

        % Parameters
        win_full = [-2.5 2.5];% full windows for extraction
        win_save = [-2 2];% windwos for saving
        baseline_win = [-2 -1.5];%windows used for baseline
        % downsample_factor = 25;
        fs = lfp_V1.samplingRate;

        win_full_samp = round(win_full * fs);


        freq_band = [1 300];
        freqs = logspace(log10(freq_band(1)), log10(freq_band(2)), 200);

        params.tapers = [3 5];              % [time-bandwidth product, # of tapers]
        params.pad = 2;                     % pad = 2; pad to 2048 points
        params.Fs = lfp_V1.samplingRate;                   
        params.fpass = [1 300];             % Frequency band of interest
        params.trialave = 0;                % Don't average yet (for spectrogram)
        params.err = 0;                     % No confidence interval
        movingwin = [0.3 0.020];  % 300ms window, 25ms step (Chronux-recommended)

        %%%%%%%%%%% initialise variables
        t0 = win_full_samp(2) + 1;
        idx_full = t0 + win_full_samp(1):t0 + win_full_samp(2);
        dat = lfp_V1.data(idx_full, cortex_ref_shank(1));
        % tic
        [phi,S,t,f] = mtspecgramc(dat, movingwin, params);
        t = t-mean(t); % midpoint aligned to zero

        % find closest frequency
        [~, closest_idx] = min(abs(f(:) - freqs(:)'), [], 1);
        fidx = unique(closest_idx);

        % toc

        % nCycles = (max(5, freqs * 0.1));
        % nCycles(nCycles > 20) = 20;
%         nCycles = round(min(20, freqs * 0.1));

        event_types = {'ripples','spindles','UP','DOWN'};

        % Output structures
        TF_amp_V1 = struct(); TF_phase_V1 = struct();
        TF_amp_HPC = struct(); TF_phase_HPC = struct();
        % PPC_V1_HPC = struct(); PPC_V1 = struct(); PPC_HPC = struct();

        win_save_samp = find(t >=win_save(1) & t <=win_save(2));
        baseline_samp = find(t >=baseline_win(1) & t <=baseline_win(2));

        nprobe = 1;
        TF_amp_V1(nprobe).timebin   = t(win_save_samp); 
        TF_phase_V1(nprobe).timebin =  t(win_save_samp);
        TF_amp_HPC(nprobe).timebin  =  t(win_save_samp);
        TF_phase_HPC(nprobe).timebin=  t(win_save_samp);

        TF_amp_V1(nprobe).freq   = f(fidx);
        TF_phase_V1(nprobe).freq = f(fidx);
        TF_amp_HPC(nprobe).freq  = f(fidx);
        TF_phase_HPC(nprobe).freq= f(fidx);


        for nprobe = 1:length(session_info(n).probe)
            tic
            disp(['Processing probe ' num2str(nprobe)]);
            


            % Interpolate event onsets to sample indices
            ripple_sample_idx = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), ripples(nprobe).onset, 'nearest', 'extrap'));
            spindle_sample_idx = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), spindles(nprobe).onset, 'nearest', 'extrap'));
            up_sample_idx      = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), slow_waves(nprobe).UP_ints(:,1), 'nearest', 'extrap'));
            down_sample_idx    = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), slow_waves(nprobe).DOWN_ints(:,1), 'nearest', 'extrap'));

            for eidx = 1:length(event_types)
                event_name = event_types{eidx};
                switch event_name
                    case 'ripples'
                        onsets = ripple_sample_idx;
                    case 'spindles'
                        onsets = spindle_sample_idx;
                    case 'UP'
                        onsets = up_sample_idx;
                    case 'DOWN'
                        onsets = down_sample_idx;
                end

                nevents = length(onsets);

                % Preallocate outputs
                
                TF_amp_V1(nprobe).(event_name)   = single(nan(2, length(fidx), length(win_save_samp), nevents));
                TF_phase_V1(nprobe).(event_name) = single(nan(2, length(fidx), length(win_save_samp), nevents));
                TF_amp_HPC(nprobe).(event_name)  = single(nan(2, length(fidx), length(win_save_samp), nevents));
                TF_phase_HPC(nprobe).(event_name)= single(nan(2, length(fidx), length(win_save_samp), nevents));

%                 PPC_V1_HPC(nprobe).(event_name)  = single(nan(2, 2,length(freqs), timebin, nevents));
%                 PPC_V1(nprobe).(event_name)      = single(nan(length(freqs),timebin, nevents));
%                 PPC_HPC(nprobe).(event_name)     = single(nan(length(freqs),timebin, nevents));


                for nevent = 1:nevents
                    t0 = onsets(nevent);
                    idx_full = t0 + win_full_samp(1):t0 + win_full_samp(2);
                    if idx_full(1) < 1 || idx_full(end) > size(lfp_V1.data,1), continue; end

                    for hemi = 1:2

                        % V1
                        dat = lfp_V1.data(idx_full, cortex_ref_shank(hemi));
                        [phi,S,t,f] = mtspecgramc(dat, movingwin, params);
                        amp_all = S(win_save_samp,fidx);
                        phase_all = phi(win_save_samp,fidx);

                        % base = mean(amp_all(baseline_samp(1):baseline_samp(2), :), 1);
                        % amp_norm = 10 * log10(bsxfun(@rdivide, amp_all, base));

                        TF_amp_V1(nprobe).(event_name)(hemi, :, :, nevent) = single(amp_all');
                        TF_phase_V1(nprobe).(event_name)(hemi, :, :, nevent) = single(phase_all');


                        % HPC
                        dat = lfp_HPC.data(idx_full, HPC_ref_shank(hemi));
                        [phi,S,t,f] = mtspecgramc(dat, movingwin, params);
                        amp_all = S(win_save_samp,fidx);
                        phase_all = phi(win_save_samp,fidx);

                        % base = mean(amp_all(baseline_samp(1):baseline_samp(2), :), 1);
                        % amp_norm = 10 * log10(bsxfun(@rdivide, amp_all, base));

                        TF_amp_HPC(nprobe).(event_name)(hemi, :, :, nevent) = single(amp_all');
                        TF_phase_HPC(nprobe).(event_name)(hemi, :, :, nevent) = single(phase_all');


                        % [C, phi, S12,S1, S2, t, f] = cohgramc(dat,dat1, movingwin, params);
                    end


                end
            end
            toc
        end

        if contains(stimulus_name{n},'Masa2tracks')
         
            save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_TF_amplitude%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'TF_amp_HPC','TF_amp_V1','-v7.3');
            save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_TF_phase%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'TF_phase_V1','TF_phase_HPC','-v7.3');

%             save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PPC%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPC_V1','PPC_HPC','PPC_V1_HPC','-v7.3');
        else


            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_TF_amplitude.mat'),'TF_amp_HPC','TF_amp_V1','-v7.3');
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_TF_phase.mat'),'TF_phase_V1','TF_phase_HPC','-v7.3');
            % save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_markov_events.mat'),'slow_waves_markov');
%             save(fullfile(options.ANALYSIS_DATAPATH,'extracted_PPC.mat'),'PPC_V1','PPC_HPC','PPC_V1_HPC','-v7.3');
        end
    end
end



%% LFP time frequency analysis (fieldtrip)
% pyversion('C:\Users\masahiro.takigawa\.conda\envs\fooof\python')
% pyversion('C:\Users\masah\anaconda3\envs\fooof\python')

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath('C:\Users\masah\Documents\GitHub\fieldtrip')
addpath('C:\Users\masahiro.takigawa\Documents\GitHub\fieldtrip')


clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar
% experiment_info=experiment_info([4 5 6 ]);
% experiment_info=experiment_info([4 5 6 18 19 21 34 35 44 45 58 59 60 71]);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'Sleep';

ft_defaults
ft_default.interactive   = 'no';     % disables prompts
ft_default.feedback      = 'no';     % disables progress feedback
ft_default.checksize     = inf;      % disables large struct warnings
ft_default.showcallinfo  = 'no';     % disables call stack info
ft_progress('init', 'none');         % disables progress bars

for nsession =1:15
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
        % PSD slope quantification using fooof ()
        tvec = LFP(1).tvec;
        SR = round(1/mean(diff(tvec)));
        nfft_seconds= 2;
        nfft = 2^(nextpow2(SR*nfft_seconds));
        win  = hanning(nfft);

        %%%%%%%%%%%% Cortical wave direction during DOWN state peak

        filterparms.deltafilter = [0.5 4];%heuristically defined.  room for improvement here.
        filterparms.thetafilter = [4 12];%heuristically defined.  room for improvement here.
        filterparms.spindlesfilter = [9 17];%heuristically defined.  room for improvement here.
        filterparms.ripplesfilter = [125 300];%heuristically defined.  room for improvement here.
        filterparms.gammafilter = [100 400];
        filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
        filterparms.gammanormwin = 20; %window for gamma normalization (s)

        %         for nprobe = 1:length(session_info(n).probe)
        nprobe = 1;
        probe_no = session_info(n).probe(nprobe).probe_id+1;

        if length(LFP(1).best_HPC) < length(LFP(2).best_HPC)
            time_idx = 1:length(LFP(1).best_HPC);
        elseif length(LFP(1).best_HPC) > length(LFP(2).best_HPC)

            time_idx = 1:length(LFP(2).best_HPC);
        else
            time_idx = 1:length(LFP(1).best_HPC) ;
        end

        lfp.samplingRate = single(round(1/mean(diff(LFP(probe_no).tvec))));
        lfp.timestamps = single(LFP(probe_no).tvec(time_idx));
        tvec = single(LFP(probe_no).tvec(time_idx));

        lfp_V1 = lfp;
        lfp_HPC = lfp;
        lfp_V1.data=[];
        lfp_HPC.data=[];

        if isfield(LFP(probe_no),'best_HPC') % if exist best HPC channels
            probe_id = [];

            if length(LFP)==1
                lfp_HPC.data= single([LFP(probe_no).best_HPC(:,time_idx)']);
                lfp_V1.data= single([LFP(probe_no).best_V1(:,time_idx)']);
            else
                lfp_HPC.data= single([LFP(1).best_HPC(:,time_idx)' LFP(2).best_HPC(:,time_idx)']);
                lfp_V1.data= single([LFP(1).average_V1(:,time_idx)' LFP(2).average_V1(:,time_idx)']);
            end
        end
        %         end


        %%% get
        HPC_ref_shank=[];
        ref_shank = [];
        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id+1;

            ref_shank = find(slow_waves(probe_no).shank_id == slow_waves(probe_no).shank(slow_waves(probe_no).channel == slow_waves(probe_no).best_channel)...
                &slow_waves(probe_no).probe_hemisphere == probe_no);
            cortex_ref_shank(probe_no) = ref_shank;

            shank_id = find(ripples(probe_no).probe_hemisphere == probe_no);
            HPC_ref_shank(probe_no) = shank_id(ripples(probe_no).best_channel);
        end


        disp('Peri event TF and PPC analysis')

        % Parameters
        win_full = [-3 3];% full windows for extraction
        win_save = [-2 2];% windwos for saving
        baseline_win = [-2 -1.5];%windows used for baseline
        % downsample_factor = 25;
        fs = single(lfp_V1.samplingRate);

        win_full_samp = round(win_full * fs);


        freq_band = [1 300];
        freqs = logspace(log10(freq_band(1)), log10(freq_band(2)), 150);
        nCycles = logspace(log10(3), log10(20), 150);
        nCycles =  max([nCycles;  0.1*(freqs * pi)]);

        % duration = nCycles ./ (freqs * pi);
  


        clear LFP
        %%%%%%%%%%% initialise variables

        event_types = {'ripples','spindles','UP','DOWN'};

        % Output structures
        TF_amp_V1 = struct(); TF_phase_V1 = struct();
        TF_amp_HPC = struct(); TF_phase_HPC = struct();
        % PPC_V1_HPC = struct(); PPC_V1 = struct(); PPC_HPC = struct();

%         win_save_samp = find(t >=win_save(1) & t <=win_save(2));
%         baseline_samp = find(t >=baseline_win(1) & t <=baseline_win(2));

        for nprobe = 1:length(session_info(n).probe)
            tic
            disp(['Processing probe ' num2str(nprobe)]);



            % Interpolate event onsets to sample indices
            ripple_sample_idx = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), ripples(nprobe).onset, 'nearest', 'extrap'));
            spindle_sample_idx = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), spindles(nprobe).onset, 'nearest', 'extrap'));
            up_sample_idx      = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), slow_waves(nprobe).UP_ints(:,1), 'nearest', 'extrap'));
            down_sample_idx    = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), slow_waves(nprobe).DOWN_ints(:,1), 'nearest', 'extrap'));

            for eidx = 1:length(event_types)
                event_name = event_types{eidx};
                switch event_name
                    case 'ripples'
                        onsets = ripple_sample_idx;
                    case 'spindles'
                        onsets = spindle_sample_idx;
                    case 'UP'
                        onsets = up_sample_idx;
                    case 'DOWN'
                        onsets = down_sample_idx;
                end

                nevents = length(onsets);
                
                if nevents == 0
                    TF_amp_V1(nprobe).(event_name) = [];
                    TF_phase_V1(nprobe).(event_name) = [];

                    TF_amp_V1(nprobe).(event_name)= [];
                    TF_phase_V1(nprobe).(event_name) = [];

                    TF_amp_HPC(nprobe).(event_name) = [];
                    TF_phase_HPC(nprobe).(event_name) = [];

                    TF_amp_HPC(nprobe).(event_name) = [];
                    TF_phase_HPC(nprobe).(event_name) = [];
                    continue
                else
  
                end
                % Preallocate outputs

%                 TF_amp_V1(nprobe).(event_name)   = single(nan(2, length(fidx), length(win_save_samp), nevents));
%                 TF_phase_V1(nprobe).(event_name) = single(nan(2, length(fidx), length(win_save_samp), nevents));
%                 TF_amp_HPC(nprobe).(event_name)  = single(nan(2, length(fidx), length(win_save_samp), nevents));
%                 TF_phase_HPC(nprobe).(event_name)= single(nan(2, length(fidx), length(win_save_samp), nevents));

                %                 PPC_V1_HPC(nprobe).(event_name)  = single(nan(2, 2,length(freqs), timebin, nevents));
                %                 PPC_V1(nprobe).(event_name)      = single(nan(length(freqs),timebin, nevents));
                %                 PPC_HPC(nprobe).(event_name)     = single(nan(length(freqs),timebin, nevents));

                cfg = [];
                cfg.method       = 'mtmconvol';
                % cfg.method       = 'wavelet';
                cfg.output       = 'fourier';
                cfg.taper        = 'hanning';
                cfg.foi          = freqs;
                cfg.t_ftimwin    = max(1 ./ freqs * 1, 0.1);
                % cfg.width    = nCycles;
                 
                cfg.toi          = win_full(1)+0.02/2:0.02:win_full(end)-0.02/2;
                cfg.pad          = 'nextpow2';
                cfg.keeptrials   = 'yes';
                cfg.numworkers  = 2;            % or however many CPU cores you want
                % cfg.dpss_keepfv = 'no';

                data = [];
                % data2 = [];

                for nevent = 1:nevents
                    t0 = onsets(nevent);
                    idx_full = t0 + win_full_samp(1):t0 + win_full_samp(2);
                    if idx_full(1) < 1 || idx_full(end) > size(lfp_V1.data,1)
                        data.trial{nevent}  = zeros(4,length(idx_full));
                        data.time{nevent}  = win_full(1):1/fs:win_full(end);
                    else
                        data.trial{nevent}  = [lfp_V1.data(idx_full, cortex_ref_shank)'; lfp_HPC.data(idx_full, HPC_ref_shank)'];
                        data.time{nevent}  = win_full(1):1/fs:win_full(end);
                    end

                    data.sampleinfo(nevent, :) = [t0 + win_full_samp(1), t0 + win_full_samp(2)];
                end

                data.fsample = fs;
                data.label = {'V1_L','V1_R','HPC_L','HPC_R'};  % one channel
                
                % clear TFR
                % TFR = ft_freqanalysis(cfg, data);
                % TFR.fourierspctrm = permute(TFR.fourierspctrm(:, :,:,TFR.time <= 2 &   TFR.time >= -2), [2 3 4 1]);
                % 
                % % Extract [-2 2] seconds
                % FFT=TFR.fourierspctrm;

                [FFT,timebin_out,freqs_out] = run_timefrequency_ft_in_batches(cfg, data, 500);

                TF_amp_V1(nprobe).(event_name)(1, :, :, :) = single(abs(squeeze(FFT(1,:,:,:)).^2));
                TF_phase_V1(nprobe).(event_name)(1, :, :, :) = single(angle(squeeze(FFT(1,:,:,:))));

                TF_amp_V1(nprobe).(event_name)(2, :, :, :) = single(abs(squeeze(FFT(2,:,:,:)).^2));
                TF_phase_V1(nprobe).(event_name)(2, :, :, :) = single(angle(squeeze(FFT(2,:,:,:))));

                TF_amp_HPC(nprobe).(event_name)(1, :, :, :) = single(abs(squeeze(FFT(3,:,:,:)).^2));
                TF_phase_HPC(nprobe).(event_name)(1, :, :, :) = single(angle(squeeze(FFT(3,:,:,:))));

                TF_amp_HPC(nprobe).(event_name)(2, :, :, :) = single(abs(squeeze(FFT(4,:,:,:)).^2));
                TF_phase_HPC(nprobe).(event_name)(2, :, :, :) = single(angle(squeeze(FFT(4,:,:,:))));

            end
            toc
        end

        freqs = freqs_out;
        timebin = timebin_out(timebin_out <= 2 &   timebin_out >= -2);
        clear TFR

        TF_amp_V1(1).timebin   = timebin;
        TF_phase_V1(1).timebin = timebin;
        TF_amp_HPC(1).timebin  =  timebin;
        TF_phase_HPC(1).timebin=  timebin;

        TF_amp_V1(1).freq   = freqs;
        TF_phase_V1(1).freq = freqs;
        TF_amp_HPC(1).freq  = freqs;
        TF_phase_HPC(1).freq= freqs;

        if contains(stimulus_name{n},'Masa2tracks')

            save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_TF_amplitude_ft%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'TF_amp_HPC','TF_amp_V1','-v7.3');
            save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_TF_phase%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'TF_phase_V1','TF_phase_HPC','-v7.3');

            %             save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PPC%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPC_V1','PPC_HPC','PPC_V1_HPC','-v7.3');
        else


            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_TF_amplitude_ft.mat'),'TF_amp_HPC','TF_amp_V1','-v7.3');
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_TF_phase_ft.mat'),'TF_phase_V1','TF_phase_HPC','-v7.3');
            % save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_markov_events.mat'),'slow_waves_markov');
            %             save(fullfile(options.ANALYSIS_DATAPATH,'extracted_PPC.mat'),'PPC_V1','PPC_HPC','PPC_V1_HPC','-v7.3');
        end

        clear TFR data TF_amp_HPC TF_amp_V1 TF_phase_V1 TF_phase_HPC
    end
end


% S_basedline
% timtevec = TFR.time(TFR.time>=-2 & TFR.time<=2);
baseline_win = [-2 -1.5];  % baseline time window
baseline_idx = timebin >= baseline_win(1) & timebin <= baseline_win(2);

S = squeeze(single(squeeze(TF_amp_V1(1).UP(1,:,:,:))));
S = normalize_tf(S, baseline_idx, 'dB');

% S_dB = 10 * log10(S);

S_dB = mean(squeeze(single(squeeze(S(:,:,:)))),3,'omitnan');

S = squeeze(single(squeeze(TF_phase_V1(1).UP(2,:,:,:))));
S_dB = squeeze(S(:,:,500));
% S_dB = TF_amp_V1_ipsi;

% % Create contour plot
% timtevec = linspace(-2, 2, size(amp_all,1));  % time axis


figure;
nexttile
[~, h] = contourf(timebin, log2(freqs), S_dB, 40, 'LineColor', 'none'); % 40 contour levels
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Power Spectrogram (dB) - Contour');
colorbar;
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'});
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
hold on;
xline(0,'r')
yline(log2([9 17]))
yline(log2([1 4]))
yline(log2([125 300]))
clim([-3.2 3.2])
% 

%% Ripple events LFP SO and spindle power and spindle

% pyversion('C:\Users\masahiro.takigawa\.conda\envs\fooof\python')
% pyversion('C:\Users\masah\anaconda3\envs\fooof\python')

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
% addpath(genpath('C:\Users\masah\Documents\GitHub\fieldtrip'))
% addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\fieldtrip'))
clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar
% experiment_info=experiment_info([4 5 6 ]);
% experiment_info=experiment_info([4 5 6 18 19 21 34 35 44 45 58 59 60 71]);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'Sleep';

for nprobe = 1:2
    periripple_LFP_info_V1(nprobe).SO_amplitude = [];
    periripple_LFP_info_V1(nprobe).SO_amplitude{2} = [];

    periripple_LFP_info_HPC(nprobe).SO_amplitude = [];
    periripple_LFP_info_HPC(nprobe).SO_amplitude{2} = [];

    periripple_LFP_info_V1(nprobe).spindle_amplitude = [];
    periripple_LFP_info_HPC(nprobe).spindle_amplitude = [];
    periripple_LFP_info_V1(nprobe).spindle_amplitude{2} = [];
    periripple_LFP_info_HPC(nprobe).spindle_amplitude{2} = [];

    periripple_LFP_info_V1(nprobe).SO_phase = [];
    periripple_LFP_info_HPC(nprobe).SO_phase = [];
    periripple_LFP_info_V1(nprobe).SO_phase{2} = [];
    periripple_LFP_info_HPC(nprobe).SO_phase{2} = [];

    periripple_LFP_info_V1(nprobe).spindle_phase = [];
    periripple_LFP_info_HPC(nprobe).spindle_phase = [];
    periripple_LFP_info_V1(nprobe).spindle_phase{2} = [];
    periripple_LFP_info_HPC(nprobe).spindle_phase{2} = [];

    periripple_LFP_info_HPC(nprobe).tvec = -1:0.02:1;
    periripple_LFP_info_V1(nprobe).tvec = -1:0.02:1;
end

if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end


for nsession =1:22
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT},option);
    % find right date number based on all experiment dates of the subject
    iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    if length(stimulus_name)>1 % Based on if POST or PRE
        if contains(Stimulus_type,'PRE')
            disp('Same stimuli multiple recordings. Will take _2')
            n = find(contains(stimulus_name,'_2')); % Usually because first stimulus is crashed.
        else

            session_info = session_info(~contains(stimulus_name,'PRE'));
            stimulus_name = stimulus_name(~contains(stimulus_name,'PRE'));

            if length(stimulus_name)>1
                disp('Same stimuli multiple recordings. Will take _2')
                n = find(contains(stimulus_name,'_2')); % Usually because first stimulus is crashed.
            else
                n =1;
            end
        end
    else
        n = 1;
    end


    % for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
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

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP_V1_SO.mat'),'LFP');
        LFP1 = LFP;
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');

        Fieldnames = fieldnames(LFP1);
        for nprobe = 1:length(LFP)
            for nfield = 1:length(Fieldnames)
                if contains(Fieldnames{nfield},'V1')==1
                    LFP(nprobe).(Fieldnames{nfield}) = LFP1(nprobe).(Fieldnames{nfield});
                end
            end
        end
        clear LFP1;

        % load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
        %             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters=clusters_ks4;

%         load(fullfile(options.ANALYSIS_DATAPATH,'spike_phase_amplitude.mat'));
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
    % PSD slope quantification using fooof ()
    tvec = LFP(1).tvec;
    SR = round(1/mean(diff(tvec)));
    nfft_seconds= 2;
    nfft = 2^(nextpow2(SR*nfft_seconds));
    win  = hanning(nfft);

    %%%%%%%%%%%% Cortical wave direction during DOWN state peak

    filterparms.deltafilter = [0.5 4];%heuristically defined.  room for improvement here.
    filterparms.thetafilter = [4 12];%heuristically defined.  room for improvement here.
    filterparms.spindlesfilter = [9 17];%heuristically defined.  room for improvement here.
    filterparms.lowgammafilter = [30 60];%heuristically defined.  room for improvement here.
    filterparms.highgammafilter = [60 100];%heuristically defined.  room for improvement here.
    filterparms.ripplesfilter = [125 300];%heuristically defined.  room for improvement here.
    % filterparms.gammafilter = [100 400];
    % filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
    % filterparms.gammanormwin = 20; %window for gamma normalization (s)

    %         for nprobe = 1:length(session_info(n).probe)
    nprobe = 1;
    probe_no = session_info(n).probe(nprobe).probe_id+1;

    if length(LFP(1).best_HPC) < length(LFP(2).best_HPC)
        time_idx = 1:length(LFP(1).best_HPC);
    elseif length(LFP(1).best_HPC) > length(LFP(2).best_HPC)

        time_idx = 1:length(LFP(2).best_HPC);
    else
        time_idx = 1:length(LFP(1).best_HPC) ;
    end
%     if length(LFP(1).best_SO_V1) < length(LFP(2).best_SO_V1)
%         time_idx = 1:length(LFP(1).best_V1_SO);
%     elseif length(LFP(1).best_SO_V1) > length(LFP(2).best_SO_V1)
% 
%         time_idx = 1:length(LFP(2).best_SO_V1);
%     else
%         time_idx = 1:length(LFP(1).best_SO_V1) ;
%     end
    lfp.samplingRate = single(round(1/mean(diff(LFP(probe_no).tvec))));
    lfp.timestamps = single(LFP(probe_no).tvec(time_idx));
    tvec = single(LFP(probe_no).tvec(time_idx));

    lfp_V1 = lfp;
    lfp_HPC = lfp;
    lfp_V1.data=[];
    lfp_HPC.data=[];

    if isfield(LFP(probe_no),'best_HPC') % if exist best HPC channels
        probe_id = [];

        if length(LFP)==1
%             lfp_HPC.data= single([LFP(probe_no).best_HPC(:,time_idx)']);
            %             lfp_V1.data= single([LFP(probe_no).best_V1(:,time_idx)']);
            lfp_V1.data= single([LFP(probe_no).best_SO_V1(:,time_idx)']);
            
        else
            lfp_HPC.data= single([LFP(1).best_HPC(:,time_idx)' LFP(2).best_HPC(:,time_idx)']);
            %             lfp_V1.data= single([LFP(1).average_V1(:,time_idx)' LFP(2).average_V1(:,time_idx)']);
            lfp_V1.data= single([LFP(1).best_SO_V1(:,time_idx)' LFP(2).best_SO_V1(:,time_idx)']);
        end
    end

%     lfp_V1.data= single([LFP(1).best_SO_V1(:,time_idx)' LFP(2).best_SO_V1(:,time_idx)']);
    %         end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% get
    HPC_ref_shank=[];
    ref_shank = [];
    for nprobe = 1:length(session_info(n).probe)
        probe_no = session_info(n).probe(nprobe).probe_id+1;

        [~,ref_shank] = max(LFP(nprobe).best_SO_V1_trough_peak_ratio);
        cortex_ref_shank(probe_no) = ref_shank;
%         ref_shank = find(slow_waves(probe_no).shank_id == slow_waves(probe_no).shank(slow_waves(probe_no).channel == slow_waves(probe_no).best_channel)...
%             &slow_waves(probe_no).probe_hemisphere == probe_no);
%         cortex_ref_shank(probe_no) = ref_shank;

        shank_id = find(ripples(probe_no).probe_hemisphere == probe_no);
        HPC_ref_shank(probe_no) = shank_id(ripples(probe_no).best_channel);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Bandpass filtter %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lfp = [];
    for nregion = 1:2
        if nregion == 1
            lfp = lfp_V1;
            lfp.data = double(lfp.data(:,cortex_ref_shank));
        else
            lfp = lfp_HPC;
            lfp.data = double(lfp.data(:,HPC_ref_shank));
        end

        % grab delta LFP
        deltaLFP = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);
        zscored_LFP = zscore(deltaLFP.data);
        lfp.SO_phase_LFP = deltaLFP.phase;
        lfp.SO_amplitude_LFP = zscore(deltaLFP.amp);
        deltaLFP = [];

        % grab spindles LFP
        filter_type  = 'bandpass';
        passband = [9 17];
        filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for ripple
        norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
        b_spindle = fir1(filter_order, norm_freq_range,filter_type);

        signal = [];
        for nShank = 1:size(lfp.data,2)
            signal(:,nShank) = filtfilt(b_spindle,1, lfp.data(:,nShank));
        end
        %                     spindle_amplitude_LFP = zscore(envelope(signal,round(lfp.samplingRate/5),'rms'));
        %                     spindle_amplitude_LFP = smoothdata(spindle_amplitude_LFP,'gaussian',round(lfp.samplingRate/5)); % envelop amplitude of spindles
        lfp.spindle_amplitude_LFP = zscore(abs(hilbert(signal))); % z scored amplitude
        lfp.spindle_phase_LFP = angle(hilbert(signal)); % phase
        signal = [];


        % grab ripples LFP
        filter_type  = 'bandpass';
        passband = [125 300];
        filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for ripple
        norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
        b_ripple = fir1(filter_order, norm_freq_range,filter_type);

        signal = [];
        for nShank = 1:size(lfp.data,2)
            signal(:,nShank) = filtfilt(b_ripple,1, lfp.data(:,nShank));
        end
        lfp.ripple_amplitude_LFP = zscore(abs(hilbert(signal))); % z scored amplitude
        lfp.ripple_phase_LFP = angle(hilbert(signal)); % phase
        signal = [];



        % grab gamma LFP
        filter_type  = 'bandpass';
        passband = [30 60];
        filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for theta
        norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
        b_theta = fir1(filter_order, norm_freq_range,filter_type);

        signal = [];
        for nShank = 1:size(lfp.data,2)
            signal(:,nShank) = filtfilt(b_theta,1, lfp.data(:,nShank));
        end

        lfp.gamma_amplitude_LFP = zscore(abs(hilbert(signal))); % z scored amplitude
        lfp.gamma_phase_LFP = angle(hilbert(signal)); % phase


        if nregion == 1
            lfp_V1  = lfp;

        else
            lfp_HPC = lfp;
        end

    end



    bin_edges = linspace(-1, 1, 11);  % 10 bins
    nBins = length(bin_edges) - 1;
    nEvents = length(ripples(1).onset);

    % Frequency bands
    band_names = fieldnames(filterparms);
    nBands = numel(band_names);

    % Preallocate: [hemi x band x bin x event]
    V1_amp_all = nan(2, nBands, nBins, nEvents);
    HPC_amp_all = nan(2, nBands, nBins, nEvents);


    fs = single(lfp_V1.samplingRate);
    win_full_samp = round([-1:0.02:1] * fs);

    for nprobe = 1:2
        event_sample_idx = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), ripples(nprobe).onset, 'nearest', 'extrap'));
        % sample_windows = event_sample_idx(nevent) + win_full_samp;

        for nevent = 1:length(event_sample_idx)
            ripple_windows = event_sample_idx(nevent) + win_full_samp; % ripple windows
            if sum(ripple_windows > length(lfp_V1.SO_amplitude_LFP))>0
                ripple_windows(sum(ripple_windows <= length(lfp_V1.SO_amplitude_LFP)):end) = nan;
                % valid_bins = ~isnan(ripple_windows);
            elseif   sum(ripple_windows < 1)>0
                ripple_windows(ripple_windows<1) = nan;

            end
            valid_bins = ~isnan(ripple_windows);

            for mprobe = 1:2
                output = nan(length(ripple_windows),1);
                output(valid_bins) = lfp_V1.SO_amplitude_LFP(ripple_windows(valid_bins),mprobe);
                periripple_LFP_info_V1(nprobe).SO_amplitude{mprobe} = [periripple_LFP_info_V1(nprobe).SO_amplitude{mprobe}    output]; % SO amplitude
                
                output = nan(length(ripple_windows),1);
                output(valid_bins) = lfp_HPC.SO_amplitude_LFP(ripple_windows(valid_bins),mprobe);
                periripple_LFP_info_HPC(nprobe).SO_amplitude{mprobe} = [periripple_LFP_info_HPC(nprobe).SO_amplitude{mprobe}    output]; % SO amplitude

                output = nan(length(ripple_windows),1);
                output(valid_bins) = lfp_V1.spindle_amplitude_LFP(ripple_windows(valid_bins),mprobe);
                periripple_LFP_info_V1(nprobe).spindle_amplitude{mprobe} = [periripple_LFP_info_V1(nprobe).spindle_amplitude{mprobe}    output]; % Spindle amplitude

                output = nan(length(ripple_windows),1);
                output(valid_bins) = lfp_HPC.spindle_amplitude_LFP(ripple_windows(valid_bins),mprobe);
                periripple_LFP_info_HPC(nprobe).spindle_amplitude{mprobe} = [periripple_LFP_info_HPC(nprobe).spindle_amplitude{mprobe}    output]; % Spindle amplitude

                output = nan(length(ripple_windows),1);
                output(valid_bins) = lfp_V1.SO_phase_LFP(ripple_windows(valid_bins),mprobe);
                periripple_LFP_info_V1(nprobe).SO_phase{mprobe} = [periripple_LFP_info_V1(nprobe).SO_phase{mprobe} output];% SO phase

                output = nan(length(ripple_windows),1);
                output(valid_bins) = lfp_HPC.SO_phase_LFP(ripple_windows(valid_bins),mprobe);
                periripple_LFP_info_HPC(nprobe).SO_phase{mprobe} = [periripple_LFP_info_HPC(nprobe).SO_phase{mprobe} output];% SO phase

                output = nan(length(ripple_windows),1);
                output(valid_bins) = lfp_V1.spindle_phase_LFP(ripple_windows(valid_bins),mprobe);
                periripple_LFP_info_V1(nprobe).spindle_phase{mprobe} = [periripple_LFP_info_V1(nprobe).spindle_phase{mprobe} output];% Spindle phase

                output = nan(length(ripple_windows),1);
                output(valid_bins) = lfp_HPC.spindle_phase_LFP(ripple_windows(valid_bins),mprobe);
                periripple_LFP_info_HPC(nprobe).spindle_phase{mprobe} = [periripple_LFP_info_HPC(nprobe).spindle_phase{mprobe} output];% Spindle phase


            end


        end
    end
    % end


    % delete(fullfile(options.ANALYSIS_DATAPATH,'periripple_LFP_info_HPC.mat'));
    % delete(fullfile(options.ANALYSIS_DATAPATH,'periripple_LFP_info_V1.mat'));
end
 periripple_LFP_info_V1_best_SO = periripple_LFP_info_V1;

if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end

% save(fullfile(analysis_folder,'periripple_LFP_info_HPC_RUN.mat'),'periripple_LFP_info_HPC','-v7.3');
save(fullfile(analysis_folder,'periripple_LFP_info_V1_best_SO.mat'),'periripple_LFP_info_V1_best_SO','-v7.3');
% save(fullfile(analysis_folder,'periripple_LFP_info_V1_2Hz.mat'),'periripple_LFP_info_V1_2Hz','-v7.3');

if contains(Stimulus_type,'Sleep') & ~contains(Stimulus_type,'PRE')
    save(fullfile(analysis_folder,'periripple_LFP_info_HPC.mat'),'periripple_LFP_info_HPC','-v7.3');
    save(fullfile(analysis_folder,'periripple_LFP_info_V1.mat'),'periripple_LFP_info_V1','-v7.3');

elseif contains(Stimulus_type,'PRE')
    save(fullfile(analysis_folder,'periripple_LFP_info_HPC_PRE.mat'),'periripple_LFP_info_HPC','-v7.3');
    save(fullfile(analysis_folder,'periripple_LFP_info_V1_PRE.mat'),'periripple_LFP_info_V1','-v7.3');

else
    save(fullfile(analysis_folder,'periripple_LFP_info_HPC_RUN.mat'),'periripple_LFP_info_HPC','-v7.3');
    save(fullfile(analysis_folder,'periripple_LFP_info_V1_RUN.mat'),'periripple_LFP_info_V1','-v7.3');
end



%%
%%


%% Ripple SO phase

% pyversion('C:\Users\masahiro.takigawa\.conda\envs\fooof\python')
% pyversion('C:\Users\masah\anaconda3\envs\fooof\python')

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
% addpath(genpath('C:\Users\masah\Documents\GitHub\fieldtrip'))
% addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\fieldtrip'))
clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar
% experiment_info=experiment_info([4 5 6 ]);
% experiment_info=experiment_info([4 5 6 18 19 21 34 35 44 45 58 59 60 71]);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'Sleep';

SO_phase_ripple_HPC_MUA_spike_rate = [];
SO_phase_HPC_MUA_spike_rate = [];
V1_SO_FR = [];
HPC_SO_FR = [];

if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end


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

    if length(stimulus_name)>1 % Based on if POST or PRE
        if contains(Stimulus_type,'PRE')
            disp('Same stimuli multiple recordings. Will take _2')
            n = find(contains(stimulus_name,'_2')); % Usually because first stimulus is crashed.
        else

            session_info = session_info(~contains(stimulus_name,'PRE'));
            stimulus_name = stimulus_name(~contains(stimulus_name,'PRE'));

            if length(stimulus_name)>1
                disp('Same stimuli multiple recordings. Will take _2')
                n = find(contains(stimulus_name,'_2')); % Usually because first stimulus is crashed.
            else
                n =1;
            end
        end
    else
        n = 1;
    end


    % for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
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

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP_V1_SO.mat'),'LFP');
        LFP1 = LFP;
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');

        Fieldnames = fieldnames(LFP1);
        for nprobe = 1:length(LFP)
            for nfield = 1:length(Fieldnames)
                if contains(Fieldnames{nfield},'V1')==1
                    LFP(nprobe).(Fieldnames{nfield}) = LFP1(nprobe).(Fieldnames{nfield});
                end
            end
        end
        clear LFP1;

        % load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
        %             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters=clusters_ks4;

%         load(fullfile(options.ANALYSIS_DATAPATH,'spike_phase_amplitude.mat'));
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
    % PSD slope quantification using fooof ()
    tvec = LFP(1).tvec;
    SR = round(1/mean(diff(tvec)));
    nfft_seconds= 2;
    nfft = 2^(nextpow2(SR*nfft_seconds));
    win  = hanning(nfft);

    %%%%%%%%%%%% Cortical wave direction during DOWN state peak

    %     filterparms.deltafilter = [0.5 4];%heuristically defined.  room for improvement here.
    filterparms.deltafilter = [0.5 4];%heuristically defined.  room for improvement here.
    filterparms.thetafilter = [4 12];%heuristically defined.  room for improvement here.
    filterparms.spindlesfilter = [9 17];%heuristically defined.  room for improvement here.
    filterparms.lowgammafilter = [30 60];%heuristically defined.  room for improvement here.
    filterparms.highgammafilter = [60 100];%heuristically defined.  room for improvement here.
    filterparms.ripplesfilter = [125 300];%heuristically defined.  room for improvement here.
    % filterparms.gammafilter = [100 400];
    % filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
    % filterparms.gammanormwin = 20; %window for gamma normalization (s)

    %         for nprobe = 1:length(session_info(n).probe)
    nprobe = 1;
    probe_no = session_info(n).probe(nprobe).probe_id+1;

    if length(LFP(1).best_HPC) < length(LFP(2).best_HPC)
        time_idx = 1:length(LFP(1).best_HPC);
    elseif length(LFP(1).best_HPC) > length(LFP(2).best_HPC)

        time_idx = 1:length(LFP(2).best_HPC);
    else
        time_idx = 1:length(LFP(1).best_HPC) ;
    end
    %     if length(LFP(1).best_SO_V1) < length(LFP(2).best_SO_V1)
    %         time_idx = 1:length(LFP(1).best_V1_SO);
    %     elseif length(LFP(1).best_SO_V1) > length(LFP(2).best_SO_V1)
    %
    %         time_idx = 1:length(LFP(2).best_SO_V1);
    %     else
    %         time_idx = 1:length(LFP(1).best_SO_V1) ;
    %     end
    lfp.samplingRate = single(round(1/mean(diff(LFP(probe_no).tvec))));
    lfp.timestamps = single(LFP(probe_no).tvec(time_idx));
    tvec = single(LFP(probe_no).tvec(time_idx));

    lfp_V1 = lfp;
    lfp_HPC = lfp;
    lfp_V1.data=[];
    lfp_HPC.data=[];

    if isfield(LFP(probe_no),'best_HPC') % if exist best HPC channels
        probe_id = [];

        if length(LFP)==1
            %             lfp_HPC.data= single([LFP(probe_no).best_HPC(:,time_idx)']);
            %             lfp_V1.data= single([LFP(probe_no).best_V1(:,time_idx)']);
            lfp_V1.data= single([LFP(probe_no).best_SO_V1(:,time_idx)']);

        else
            lfp_HPC.data= single([LFP(1).best_HPC(:,time_idx)' LFP(2).best_HPC(:,time_idx)']);
            %             lfp_V1.data= single([LFP(1).average_V1(:,time_idx)' LFP(2).average_V1(:,time_idx)']);
            lfp_V1.data= single([LFP(1).best_SO_V1(:,time_idx)' LFP(2).best_SO_V1(:,time_idx)']);
        end
    end

    lfp_V1.data= single([LFP(1).best_SO_V1(:,time_idx)' LFP(2).best_SO_V1(:,time_idx)']);
    %         end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% get
    cortex_ref_shank_probe = [];
    HPC_ref_shank=[];
    ref_shank = [];
    cortex_ref_shank=[];
    for nprobe = 1:length(session_info(n).probe)
        probe_no = session_info(n).probe(nprobe).probe_id+1;

        [~,ref_shank] = max(LFP(nprobe).best_SO_V1_trough_peak_ratio);
        
        cortex_ref_shank(probe_no)=  find(slow_waves(probe_no).probe_hemisphere== probe_no & slow_waves(probe_no).shank_id == ref_shank);
        cortex_ref_shank_probe(probe_no) = ref_shank;
        %         ref_shank = find(slow_waves(probe_no).shank_id == slow_waves(probe_no).shank(slow_waves(probe_no).channel == slow_waves(probe_no).best_channel)...
        %             &slow_waves(probe_no).probe_hemisphere == probe_no);
        %         cortex_ref_shank(probe_no) = ref_shank;

        shank_id = find(ripples(probe_no).probe_hemisphere == probe_no);
        HPC_ref_shank(probe_no) = shank_id(ripples(probe_no).best_channel);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Bandpass filtter %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lfp = [];
    for nregion = 1:2
        if nregion == 1
            lfp = lfp_V1;
            lfp.data = double(lfp.data(:,:));
        else
            lfp = lfp_HPC;
            lfp.data = double(lfp.data(:,:));
        end

        % grab delta LFP
        deltaLFP = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);
        zscored_LFP = zscore(deltaLFP.data);
        lfp.SO_phase_LFP = deltaLFP.phase;
        lfp.SO_amplitude_LFP = zscore(deltaLFP.amp);
        deltaLFP = [];

        % grab spindles LFP
        filter_type  = 'bandpass';
        passband = [9 17];
        filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for ripple
        norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
        b_spindle = fir1(filter_order, norm_freq_range,filter_type);

        signal = [];
        for nShank = 1:size(lfp.data,2)
            signal(:,nShank) = filtfilt(b_spindle,1, lfp.data(:,nShank));
        end
        %                     spindle_amplitude_LFP = zscore(envelope(signal,round(lfp.samplingRate/5),'rms'));
        %                     spindle_amplitude_LFP = smoothdata(spindle_amplitude_LFP,'gaussian',round(lfp.samplingRate/5)); % envelop amplitude of spindles
        lfp.spindle_amplitude_LFP = zscore(abs(hilbert(signal))); % z scored amplitude
        lfp.spindle_phase_LFP = angle(hilbert(signal)); % phase
        signal = [];


        % grab ripples LFP
        filter_type  = 'bandpass';
        passband = [125 300];
        filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for ripple
        norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
        b_ripple = fir1(filter_order, norm_freq_range,filter_type);

        signal = [];
        for nShank = 1:size(lfp.data,2)
            signal(:,nShank) = filtfilt(b_ripple,1, lfp.data(:,nShank));
        end
        lfp.ripple_amplitude_LFP = zscore(abs(hilbert(signal))); % z scored amplitude
        lfp.ripple_phase_LFP = angle(hilbert(signal)); % phase
        signal = [];



        % grab gamma LFP
        filter_type  = 'bandpass';
        passband = [30 60];
        filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for theta
        norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
        b_theta = fir1(filter_order, norm_freq_range,filter_type);

        signal = [];
        for nShank = 1:size(lfp.data,2)
            signal(:,nShank) = filtfilt(b_theta,1, lfp.data(:,nShank));
        end

        lfp.gamma_amplitude_LFP = zscore(abs(hilbert(signal))); % z scored amplitude
        lfp.gamma_phase_LFP = angle(hilbert(signal)); % phase


        if nregion == 1
            lfp_V1  = lfp;

        else
            lfp_HPC = lfp;
        end

    end

    for probe_no = 1:2
        % phase and amplitude
        SO_phase_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
        SO_phase_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));

        SO_amplitude_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
        SO_amplitude_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));

        spindle_phase_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
        spindle_phase_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).onset));

        spindle_amplitude_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
        spindle_amplitude_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).onset));

        SO_phase_spindle_onset = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
        SO_phase_spindle_peaktime = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));

        SO_amplitude_spindle_onset = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
        SO_amplitude_spindle_peaktime = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));

        % SO_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
        % spindle_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
        % SO_spindles_coupling = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
        ripple_peak_amplitude = nan(length(ripples(probe_no).shank_id), length(ripples(probe_no).peaktimes));
        ripple_onset_amplitude = nan(length(ripples(probe_no).shank_id), length(ripples(probe_no).onset));
        spindle_peak_amplitude = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).peaktimes));
        spindle_onset_amplitude = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));


        gamma_phase_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
        gamma_phase_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
        gamma_amplitude_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
        gamma_amplitude_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));

        gamma_phase_spindle_peaktime = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
        gamma_phase_spindle_onset = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
        gamma_amplitude_spindle_peaktime = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
        gamma_amplitude_spindle_onset = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));


        for nevent = 1:length(ripples(probe_no).onset)
            [~,tidx]=min(abs(tvec-ripples(probe_no).peaktimes(nevent)));

            % peak ampltiude and phase
            SO_phase_ripple_peaktime(:,nevent) = lfp_V1.SO_phase_LFP(tidx,:);
            spindle_phase_ripple_peaktime(:,nevent) = lfp_V1.spindle_phase_LFP(tidx,:);
            SO_amplitude_ripple_peaktime(:,nevent) = lfp_V1.SO_amplitude_LFP(tidx,:);
            spindle_amplitude_ripple_peaktime(:,nevent) = lfp_V1.spindle_amplitude_LFP(tidx,:);
            ripple_peak_amplitude(:,nevent) = lfp_HPC.ripple_amplitude_LFP(tidx,:);
            gamma_phase_ripple_peaktime(:, nevent) = lfp_V1.gamma_phase_LFP(tidx, :);
            gamma_amplitude_ripple_peaktime(:, nevent) = lfp_V1.gamma_amplitude_LFP(tidx, :);


            [~,tidx]=min(abs(tvec-ripples(probe_no).onset(nevent)));
            % onset amplitude and phase
            SO_phase_ripple_onset(:,nevent) = lfp_V1.SO_phase_LFP(tidx,:);
            spindle_phase_ripple_onset(:,nevent) = lfp_V1.spindle_phase_LFP(tidx,:);
            SO_amplitude_ripple_onset(:,nevent) = lfp_V1.SO_amplitude_LFP(tidx,:);
            spindle_amplitude_ripple_onset(:,nevent) = lfp_V1.spindle_amplitude_LFP(tidx,:);
            ripple_onset_amplitude(:,nevent) = lfp_HPC.ripple_amplitude_LFP(tidx,:);
            gamma_phase_ripple_onset(:, nevent) = lfp_V1.gamma_phase_LFP(tidx, :);
            gamma_amplitude_ripple_onset(:, nevent) = lfp_V1.gamma_amplitude_LFP(tidx, :);
        end

        ripples(probe_no).SO_phase_ripple_peaktime = SO_phase_ripple_peaktime;
        ripples(probe_no).spindle_phase_ripple_peaktime = spindle_phase_ripple_peaktime;
        ripples(probe_no).SO_amplitude_ripple_peaktime = SO_amplitude_ripple_peaktime;
        ripples(probe_no).spindle_amplitude_ripple_peaktime = spindle_amplitude_ripple_peaktime;
        ripples(probe_no).ripple_peak_amplitude = ripple_peak_amplitude;

        ripples(probe_no).SO_phase_ripple_onset = SO_phase_ripple_onset;
        ripples(probe_no).spindle_phase_ripple_onset = spindle_phase_ripple_onset;
        ripples(probe_no).SO_amplitude_ripple_onset = SO_amplitude_ripple_onset;
        ripples(probe_no).spindle_amplitude_ripple_onset = spindle_amplitude_ripple_onset;
        ripples(probe_no).ripple_onset_amplitude = ripple_onset_amplitude;

        ripples(probe_no).gamma_phase_ripple_peaktime = gamma_phase_ripple_peaktime;
        ripples(probe_no).gamma_amplitude_ripple_peaktime = gamma_amplitude_ripple_peaktime;
        ripples(probe_no).gamma_phase_ripple_onset = gamma_phase_ripple_onset;
        ripples(probe_no).gamma_amplitude_ripple_onset = gamma_amplitude_ripple_onset;

        ripples(probe_no).cortex_SO_ref_shank = cortex_ref_shank_probe;
    end
%     save(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events_V1_best_SO_channel.mat'),'ripples');
%     save(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events_2Hz_SO.mat'),'ripples');
    %   load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'),'ripples');

    % From cell structure back to spike times and spike id
    session_clusters.spike_id = vertcat(session_clusters.spike_id{:});
    session_clusters.spike_times = vertcat(session_clusters.spike_times{:});
    [session_clusters.spike_times, index] = sort(session_clusters.spike_times);
    session_clusters.spike_id = session_clusters.spike_id(index);


    %%%%%%%%%%%%%% Temporary (currently only one ref shank for one probe)
    %     cortex_ref_shank = [1 2];
    %     spike_phase_amplitude.SO_phase(:,cortex_ref_shank(1));
    %
    %
    % Unwrap phase to allow for linear interpolation across the -pi/pi boundary
%     cortex_ref_shank

    phi_unwrapped = unwrap(lfp_V1.SO_phase_LFP(:,cortex_ref_shank));
    % Interpolate LFP phase to the exact spike times
    spike_phases_unwrapped = interp1(tvec, phi_unwrapped, session_clusters.spike_times, 'linear');
    % Re-wrap spike phases back to [-pi, pi]
    spike_phases = mod(spike_phases_unwrapped + pi, 2*pi) - pi;

    % session_clusters.spike_id
    % spike_phase_amplitude
    SO_phase_distribution=[];
    SO_phase_distribution_ripples=[];
    SO_phase_distribution_SO=[];
    SO_phase_distribution_DOWN_UP=[];

    ripple_times = [ripples(1).onset ripples(1).offset; ripples(2).onset ripples(2).offset];
    ripple_times = ConsolidateIntervals(ripple_times);
    sleep_times = [behavioural_state_merged.SWS];
    % sleep_times = ConsolidateIntervals(sleep_times);
    UP_DOWN_times = [slow_waves(1).DOWN_ints(:,1)-0.02 slow_waves(1).DOWN_ints(:,2)+0.02];
    UP_DOWN_times = ConsolidateIntervals(UP_DOWN_times);

    SO_times = SubtractIntervals(sleep_times,UP_DOWN_times); % Sleep without delta

    % SO_times = ConsolidateIntervals(SO_times);

    no_phase_bins = 10;
    phase_bin_edges = -pi:pi/no_phase_bins:pi;
    phase_bins = -pi+pi/no_phase_bins/2:pi/no_phase_bins:pi-pi/no_phase_bins/2;

    for hemi = 1:2
        [~, ~,loc] = histcounts(lfp_V1.timestamps, reshape(sleep_times', [], 1));
        % [~, ~,loc] = histcounts(lfp_V1.timestamps, reshape(behavioural_state_merged.SWS', [], 1));
        tidx = mod(loc, 2) == 1;
        SO_phase_distribution(hemi,:) = histcounts(lfp_V1.SO_phase_LFP(tidx,cortex_ref_shank(hemi)),phase_bin_edges)./lfp_V1.samplingRate;

        [~, ~,loc] = histcounts(lfp_V1.timestamps, reshape(SO_times', [], 1));
        % [~, ~,loc] = histcounts(lfp_V1.timestamps, reshape(behavioural_state_merged.SWS', [], 1));
        tidx = mod(loc, 2) == 1;
        SO_phase_distribution_SO(hemi,:) = histcounts(lfp_V1.SO_phase_LFP(tidx,1),phase_bin_edges)./lfp_V1.samplingRate;

        [~, ~,loc] = histcounts(lfp_V1.timestamps, reshape(UP_DOWN_times', [], 1));
        % [~, ~,loc] = histcounts(lfp_V1.timestamps, reshape(behavioural_state_merged.SWS', [], 1));
        tidx = mod(loc, 2) == 1;
        SO_phase_distribution_DOWN_UP(hemi,:) = histcounts(lfp_V1.SO_phase_LFP(tidx,1),phase_bin_edges)./lfp_V1.samplingRate;

        [~, ~,loc] = histcounts(lfp_V1.timestamps, reshape(ripple_times', [], 1));
        tidx = mod(loc, 2) == 1;
        SO_phase_distribution_ripples(hemi,:) = histcounts(lfp_V1.SO_phase_LFP(tidx,cortex_ref_shank(hemi)),phase_bin_edges)./lfp_V1.samplingRate;
    end

    hemispheres = {'L','R'};
    SO_phase_MUA_spike_rate{nsession} = [];

    for hemi = 1:2
        cell_id = session_clusters.cluster_id(session_clusters.firing_rate<=50 & contains(session_clusters.region,'V1') & contains(session_clusters.region,hemispheres{hemi}));
        spike_index = ismember(session_clusters.spike_id,cell_id);

        [~, ~,loc] = histcounts(session_clusters.spike_times, reshape(behavioural_state_merged.SWS', [], 1));
        spike_index_sleep = mod(loc, 2) == 1;
        spike_index_sleep = (spike_index+spike_index_sleep)>1;

        % ripple_times = [ripples_all(1).onset(ripples_all(1).session_count ==nsession) ripples_all(1).offset(ripples_all(1).session_count ==nsession)
        %     ripples_all(2).onset(ripples_all(2).session_count ==nsession) ripples_all(2).offset(ripples_all(2).session_count ==nsession)];



        [~, ~,loc] = histcounts(session_clusters.spike_times, reshape(ripple_times', [], 1));
        ripple_spike_index = mod(loc, 2) == 1;
        ripple_spike_index = (ripple_spike_index+spike_index)>1;

        [~, ~,loc] = histcounts(session_clusters.spike_times, reshape(SO_times', [], 1));
        SO_spike_index = mod(loc, 2) == 1;
        SO_spike_index = (SO_spike_index+spike_index)>1;

        [~, ~,loc] = histcounts(session_clusters.spike_times, reshape(UP_DOWN_times', [], 1));
        DOWN_spike_index = mod(loc, 2) == 1;
        DOWN_spike_index = (DOWN_spike_index+spike_index)>1;
        %%%%%% SWS
        SO_phase_MUA_spike_rate{nsession}(hemi,1,:) = histcounts(spike_phases(spike_index_sleep,1),phase_bin_edges)./SO_phase_distribution(1,:)/length(cell_id);
        SO_phase_MUA_spike_rate{nsession}(hemi,2,:) = histcounts(spike_phases(spike_index_sleep,2),phase_bin_edges)./SO_phase_distribution(2,:)/length(cell_id);
        SO_phase_SO_MUA_spike_rate{nsession}(hemi,1,:) = histcounts(spike_phases(SO_spike_index,1),phase_bin_edges)./SO_phase_distribution_SO(1,:)/length(cell_id);
        SO_phase_SO_MUA_spike_rate{nsession}(hemi,2,:) = histcounts(spike_phases(SO_spike_index,2),phase_bin_edges)./SO_phase_distribution_SO(2,:)/length(cell_id);
        SO_phase_DOWN_MUA_spike_rate{nsession}(hemi,1,:) = histcounts(spike_phases(DOWN_spike_index,1),phase_bin_edges)./SO_phase_distribution_DOWN_UP(1,:)/length(cell_id);
        SO_phase_DOWN_MUA_spike_rate{nsession}(hemi,2,:) = histcounts(spike_phases(DOWN_spike_index,2),phase_bin_edges)./SO_phase_distribution_DOWN_UP(2,:)/length(cell_id);
        SO_phase_ripple_MUA_spike_rate{nsession}(hemi,1,:) = histcounts(spike_phases(ripple_spike_index,1),phase_bin_edges)./SO_phase_distribution_ripples(1,:)/length(cell_id);
        SO_phase_ripple_MUA_spike_rate{nsession}(hemi,2,:) = histcounts(spike_phases(ripple_spike_index,2),phase_bin_edges)./SO_phase_distribution_ripples(2,:)/length(cell_id);


        %%%% HPC
        cell_id = session_clusters.cluster_id(session_clusters.firing_rate<=50&contains(session_clusters.region,'HPC') & contains(session_clusters.region,hemispheres{hemi}));
        spike_index = ismember(session_clusters.spike_id,cell_id);

        [~, ~,loc] = histcounts(session_clusters.spike_times, reshape(behavioural_state_merged.SWS', [], 1));
        spike_index_sleep = mod(loc, 2) == 1;
        spike_index_sleep = (spike_index+spike_index_sleep)>1;

        % ripple_times = [ripples_all(1).onset(ripples_all(1).session_count ==nsession) ripples_all(1).offset(ripples_all(1).session_count ==nsession)
        %     ripples_all(2).onset(ripples_all(2).session_count ==nsession) ripples_all(2).offset(ripples_all(2).session_count ==nsession)];



        [~, ~,loc] = histcounts(session_clusters.spike_times, reshape(ripple_times', [], 1));
        ripple_spike_index = mod(loc, 2) == 1;
        ripple_spike_index = (ripple_spike_index+spike_index)>1;

        [~, ~,loc] = histcounts(session_clusters.spike_times, reshape(SO_times', [], 1));
        SO_spike_index = mod(loc, 2) == 1;
        SO_spike_index = (SO_spike_index+spike_index)>1;

        [~, ~,loc] = histcounts(session_clusters.spike_times, reshape(UP_DOWN_times', [], 1));
        DOWN_spike_index = mod(loc, 2) == 1;
        DOWN_spike_index = (DOWN_spike_index+spike_index)>1;
        %%%%%% SWS
        SO_phase_HPC_MUA_spike_rate{nsession}(hemi,1,:) = histcounts(spike_phases(spike_index_sleep,1),phase_bin_edges)./SO_phase_distribution(1,:)/length(cell_id);
        SO_phase_HPC_MUA_spike_rate{nsession}(hemi,2,:) = histcounts(spike_phases(spike_index_sleep,2),phase_bin_edges)./SO_phase_distribution(2,:)/length(cell_id);
        SO_phase_SO_HPC_MUA_spike_rate{nsession}(hemi,1,:) = histcounts(spike_phases(SO_spike_index,1),phase_bin_edges)./SO_phase_distribution_SO(1,:)/length(cell_id);
        SO_phase_SO_HPC_MUA_spike_rate{nsession}(hemi,2,:) = histcounts(spike_phases(SO_spike_index,2),phase_bin_edges)./SO_phase_distribution_SO(2,:)/length(cell_id);
        SO_phase_DOWN_HPC_MUA_spike_rate{nsession}(hemi,1,:) = histcounts(spike_phases(DOWN_spike_index,1),phase_bin_edges)./SO_phase_distribution_DOWN_UP(1,:)/length(cell_id);
        SO_phase_DOWN_HPC_MUA_spike_rate{nsession}(hemi,2,:) = histcounts(spike_phases(DOWN_spike_index,2),phase_bin_edges)./SO_phase_distribution_DOWN_UP(2,:)/length(cell_id);
        SO_phase_ripple_HPC_MUA_spike_rate{nsession}(hemi,1,:) = histcounts(spike_phases(ripple_spike_index,1),phase_bin_edges)./SO_phase_distribution_ripples(1,:)/length(cell_id);
        SO_phase_ripple_HPC_MUA_spike_rate{nsession}(hemi,2,:) = histcounts(spike_phases(ripple_spike_index,2),phase_bin_edges)./SO_phase_distribution_ripples(2,:)/length(cell_id);
    end

    fig = figure;
    fig.Name = sprintf('%s %s SO phase MUA firing',options.SUBJECT,options.SESSION)
    sgtitle(fig.Name)
    fig.Position= [844 66 560 906];
    counter = 1;
    
    for hemi = 1:2
        subplot(4,2,counter)
        plot(phase_bins,squeeze(SO_phase_MUA_spike_rate{nsession}(hemi,1,:)),'r');hold on;
        plot(phase_bins,squeeze(SO_phase_MUA_spike_rate{nsession}(hemi,2,:)),'b')
        xlabel('SO Phase');ylabel('FR (Hz)');
        set(gca,'TickDir','out','Box','off','FontSize',12)
        xticks([-pi -pi/2 0 pi/2 pi])
        xticklabels({'-π','-π/2','0','π/2','π'})
        title(sprintf('V1 %s SO Phase',hemispheres{hemi}))
        counter = counter+1;
        
        subplot(4,2,counter)
        plot(phase_bins,squeeze(SO_phase_ripple_MUA_spike_rate{nsession}(hemi,1,:)),'r');hold on;
        plot(phase_bins,squeeze(SO_phase_ripple_MUA_spike_rate{nsession}(hemi,2,:)),'b')
        xticks([-pi -pi/2 0 pi/2 pi])
        xticklabels({'-π','-π/2','0','π/2','π'})
        set(gca,'TickDir','out','Box','off','FontSize',12)
        title(sprintf('V1 %s SO Phase (ripples)',hemispheres{hemi}))
        counter = counter+1;
    end

    % counter = 1;
    for hemi = 1:2
        subplot(4,2,counter)
        plot(phase_bins,squeeze(SO_phase_HPC_MUA_spike_rate{nsession}(hemi,1,:)),'r');hold on;
        plot(phase_bins,squeeze(SO_phase_HPC_MUA_spike_rate{nsession}(hemi,2,:)),'b')
        xlabel('SO Phase');ylabel('FR (Hz)')
        xticks([-pi -pi/2 0 pi/2 pi])
        xticklabels({'-π','-π/2','0','π/2','π'})
        set(gca,'TickDir','out','Box','off','FontSize',12)
        title(sprintf('HPC %s SO Phase',hemispheres{hemi}))
        counter = counter+1;

        subplot(4,2,counter)
        plot(phase_bins,squeeze(SO_phase_ripple_HPC_MUA_spike_rate{nsession}(hemi,1,:)),'r');hold on;
        plot(phase_bins,squeeze(SO_phase_ripple_HPC_MUA_spike_rate{nsession}(hemi,2,:)),'b')
        xlabel('SO Phase');ylabel('FR (Hz)')
        set(gca,'TickDir','out','Box','off','FontSize',12)
        xticks([-pi -pi/2 0 pi/2 pi])
        xticklabels({'-π','-π/2','0','π/2','π'})
        title(sprintf('HPC %s SO Phase (ripples)',hemispheres{hemi}))
        counter = counter+1;
    end


    fig = figure;
    fig.Name = sprintf('%s %s SO phase MUA firing (SO and UP-DOWN)',options.SUBJECT,options.SESSION)
    sgtitle(fig.Name)
    fig.Position= [844 66 560 906];
    counter = 1;

    for hemi = 1:2
        subplot(4,2,counter)
        plot(phase_bins,squeeze(SO_phase_SO_MUA_spike_rate{nsession}(hemi,1,:)),'r');hold on;
        plot(phase_bins,squeeze(SO_phase_SO_MUA_spike_rate{nsession}(hemi,2,:)),'b')
        xlabel('SO Phase');ylabel('FR (Hz)');
        set(gca,'TickDir','out','Box','off','FontSize',12)
        xticks([-pi -pi/2 0 pi/2 pi])
        xticklabels({'-π','-π/2','0','π/2','π'})
        title(sprintf('V1 %s SO Phase (SO)',hemispheres{hemi}))
        counter = counter+1;

        subplot(4,2,counter)
        plot(phase_bins,squeeze(SO_phase_DOWN_MUA_spike_rate{nsession}(hemi,1,:)),'r');hold on;
        plot(phase_bins,squeeze(SO_phase_DOWN_MUA_spike_rate{nsession}(hemi,2,:)),'b')
        xticks([-pi -pi/2 0 pi/2 pi])
        xticklabels({'-π','-π/2','0','π/2','π'})
        set(gca,'TickDir','out','Box','off','FontSize',12)
        title(sprintf('V1 %s SO Phase (DOWN)',hemispheres{hemi}))
        counter = counter+1;
    end

    % counter = 1;
    for hemi = 1:2
        subplot(4,2,counter)
        plot(phase_bins,squeeze(SO_phase_SO_HPC_MUA_spike_rate{nsession}(hemi,1,:)),'r');hold on;
        plot(phase_bins,squeeze(SO_phase_SO_HPC_MUA_spike_rate{nsession}(hemi,2,:)),'b')
        xlabel('SO Phase');ylabel('FR (Hz)')
        xticks([-pi -pi/2 0 pi/2 pi])
        xticklabels({'-π','-π/2','0','π/2','π'})
        set(gca,'TickDir','out','Box','off','FontSize',12)
        title(sprintf('HPC %s SO Phase (SO)',hemispheres{hemi}))
        counter = counter+1;

        subplot(4,2,counter)
        plot(phase_bins,squeeze(SO_phase_DOWN_HPC_MUA_spike_rate{nsession}(hemi,1,:)),'r');hold on;
        plot(phase_bins,squeeze(SO_phase_DOWN_HPC_MUA_spike_rate{nsession}(hemi,2,:)),'b')
        xlabel('SO Phase');ylabel('FR (Hz)')
        set(gca,'TickDir','out','Box','off','FontSize',12)
        xticks([-pi -pi/2 0 pi/2 pi])
        xticklabels({'-π','-π/2','0','π/2','π'})
        title(sprintf('HPC %s SO Phase (DOWN)',hemispheres{hemi}))
        counter = counter+1;
    end


%     for hemi = 1:2
%         temp = squeeze(SO_phase_MUA_spike_rate{nsession}(hemi,1,:));
%         V1_SO_FR(nsession,hemi,1,:) = [mean(temp(phase_bins>-pi/2 & phase_bins<pi/2)) mean(temp(phase_bins<-pi/2 | phase_bins>pi/2))];
%         temp = squeeze(SO_phase_MUA_spike_rate{nsession}(hemi,2,:));
%         V1_SO_FR(nsession,hemi,2,:) = [mean(temp(phase_bins>-pi/2 & phase_bins<pi/2)) mean(temp(phase_bins<-pi/2 | phase_bins>pi/2))];
% 
%         temp = squeeze(SO_phase_HPC_MUA_spike_rate{nsession}(hemi,1,:));
%         HPC_SO_FR(nsession,hemi,1,:) = [mean(temp(phase_bins>-pi/2 & phase_bins<pi/2)) mean(temp(phase_bins<-pi/2 | phase_bins>pi/2))];
%         temp = squeeze(SO_phase_HPC_MUA_spike_rate{nsession}(hemi,2,:));
%         HPC_SO_FR(nsession,hemi,2,:) = [mean(temp(phase_bins>-pi/2 & phase_bins<pi/2)) mean(temp(phase_bins<-pi/2 | phase_bins>pi/2))];
%     end
    
    % save_all_figures('D:\corticohippocampal_replay\V1-HPC sleep reactivation\SO phase plots',[]);
    save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\SO phase plots best'),[])
%     save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\SO phase plots (2Hz)'),[])
end


if contains(Stimulus_type,'Sleep') & ~contains(Stimulus_type,'PRE')
    save(fullfile(analysis_folder,'SO_phase_MUA_best_V1_SO.mat'),...
        'SO_phase_ripple_HPC_MUA_spike_rate','SO_phase_HPC_MUA_spike_rate','SO_phase_SO_HPC_MUA_spike_rate','SO_phase_DOWN_HPC_MUA_spike_rate',...
        'SO_phase_ripple_MUA_spike_rate','SO_phase_MUA_spike_rate','SO_phase_SO_MUA_spike_rate','SO_phase_DOWN_MUA_spike_rate');

%     save(fullfile(analysis_folder,'SO_phase_MUA_SO_2Hz.mat'),...
%         'SO_phase_ripple_HPC_MUA_spike_rate','SO_phase_HPC_MUA_spike_rate','SO_phase_SO_HPC_MUA_spike_rate','SO_phase_DOWN_HPC_MUA_spike_rate',...
%         'SO_phase_ripple_MUA_spike_rate','SO_phase_MUA_spike_rate','SO_phase_SO_MUA_spike_rate','SO_phase_DOWN_MUA_spike_rate');

    % save(fullfile(analysis_folder,'periripple_LFP_info_V1.mat'),'periripple_LFP_info_V1','-v7.3');

elseif contains(Stimulus_type,'PRE')
    % save(fullfile(analysis_folder,'periripple_LFP_info_HPC_PRE.mat'),'periripple_LFP_info_HPC','-v7.3');
    % save(fullfile(analysis_folder,'periripple_LFP_info_V1_PRE.mat'),'periripple_LFP_info_V1','-v7.3');

else
    % save(fullfile(analysis_folder,'periripple_LFP_info_HPC_RUN.mat'),'periripple_LFP_info_HPC','-v7.3');
    % save(fullfile(analysis_folder,'periripple_LFP_info_V1_RUN.mat'),'periripple_LFP_info_V1','-v7.3');
end





%%%%% Plot SO phase relationship
load(fullfile(analysis_folder,'SO_phase_MUA_best_V1_SO.mat'),...
    'SO_phase_ripple_HPC_MUA_spike_rate','SO_phase_HPC_MUA_spike_rate','SO_phase_SO_HPC_MUA_spike_rate','SO_phase_DOWN_HPC_MUA_spike_rate',...
    'SO_phase_ripple_MUA_spike_rate','SO_phase_MUA_spike_rate','SO_phase_SO_MUA_spike_rate','SO_phase_DOWN_MUA_spike_rate');

load(fullfile(analysis_folder,'SO_phase_MUA_SO_2Hz.mat'),...
    'SO_phase_ripple_HPC_MUA_spike_rate','SO_phase_HPC_MUA_spike_rate','SO_phase_SO_HPC_MUA_spike_rate','SO_phase_DOWN_HPC_MUA_spike_rate',...
    'SO_phase_ripple_MUA_spike_rate','SO_phase_MUA_spike_rate','SO_phase_SO_MUA_spike_rate','SO_phase_DOWN_MUA_spike_rate');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
no_phase_bins = 10;
phase_bin_edges = -pi:pi/no_phase_bins:pi;
phase_bins = -pi+pi/no_phase_bins/2:pi/no_phase_bins:pi-pi/no_phase_bins/2;


V1_phase = [];

for hemi = 1:2
    
    for nession = 1:22
        V1_phase{hemi}(1,nession,:) = SO_phase_MUA_spike_rate{nsession}(hemi,1,:);
        V1_phase{hemi}(2,nession,:) = SO_phase_MUA_spike_rate{nsession}(hemi,2,:);
        %         V1_phase{hemi}(1,nession,:) = normalize(SO_phase_MUA_spike_rate{nsession}(hemi,1,:),'range');
        % V1_phase{hemi}(2,nession,:) = normalize(SO_phase_MUA_spike_rate{nsession}(hemi,2,:),'range');
    end
end

mean_phases_ipsi = mean([squeeze((V1_phase{1}(1,:,:))); squeeze((V1_phase{2}(2,:,:)))]);
mean_phases_contra = mean([squeeze((V1_phase{1}(2,:,:))); squeeze((V1_phase{2}(1,:,:)))]);

SE_phases_ipsi = std([squeeze((V1_phase{1}(1,:,:))); squeeze((V1_phase{2}(2,:,:)))])/sqrt(length(SO_phase_MUA_spike_rate{nsession}));
SE_phases_contra = std([squeeze((V1_phase{1}(2,:,:))); squeeze((V1_phase{2}(1,:,:)))])/sqrt(length(SO_phase_MUA_spike_rate{nsession}));


% 1. Setup Colors (Example: Teal for Ipsi, Orange for Contra)
c_ipsi = [35,139,69]/256; 
c_contra = [106,81,163]/256;

% 2. Ensure inputs are row vectors for the fill command
x = phase_bins(:)'; 

% 3. Pre-calculate Upper and Lower bounds for shading
% Ipsi Bounds
ipsi_upper = (mean_phases_ipsi(:)' + SE_phases_ipsi(:)');
ipsi_lower = (mean_phases_ipsi(:)' - SE_phases_ipsi(:)');
% Contra Bounds
contra_upper = (mean_phases_contra(:)' + SE_phases_contra(:)');
contra_lower = (mean_phases_contra(:)' - SE_phases_contra(:)');

%%%% Plot SO V1 phase
fig = figure;
fig.Name = 'V1 SO phase Peak vs Troug (sleep)';
% fig.Name = 'V1 SO phase Peak vs Trough (DOWN 2Hz)';
sgtitle(fig.Name)
fig.Position= [844 66 560 906];

subplot(3,2,1)
hold on;
% --- Plot Ipsi (Shading first, then Line) ---
fill([x fliplr(x)], [ipsi_lower fliplr(ipsi_upper)], ...
    c_ipsi, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h1 = plot(x, mean_phases_ipsi, 'Color', c_ipsi, 'LineWidth', 2);

% --- Plot Contra (Shading first, then Line) ---
fill([x fliplr(x)], [contra_lower fliplr(contra_upper)], ...
    c_contra, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h2 = plot(x, mean_phases_contra, 'Color', c_contra, 'LineWidth', 2);

% --- Formatting ---
yline(0, '--k', 'HandleVisibility', 'off'); % 'HandleVisibility' hides it from legend
xline(-pi/2, '--k', 'HandleVisibility', 'off'); % 'HandleVisibility' hides it from legend
xline(pi/2, '--k', 'HandleVisibility', 'off'); % 'HandleVisibility' hides it from legend

xlabel('Phases');
ylabel('Mean firing rate (Hz)');
title('Ipsi vs Contra Phase Preference');

% Add a legend (using the plot handles h1 and h2)
legend([h1 h2], {'Ipsi', 'Contra'}, 'Location', 'best', 'Box', 'off');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-pi','-pi/2','0','pi/2','pi'})
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-pi pi]);
ylim([1 6])
hold off;



%%%%%% paired t test
for nsession = 1:22
    for hemi = 1:2
        temp = squeeze(SO_phase_SO_MUA_spike_rate{nsession}(hemi,1,:));
        V1_SO_FR(nsession,hemi,1,:) = [mean(temp(phase_bins>-pi/2 & phase_bins<pi/2)) mean(temp(phase_bins<-pi/2 | phase_bins>pi/2))];
        temp = squeeze(SO_phase_SO_MUA_spike_rate{nsession}(hemi,2,:));
        V1_SO_FR(nsession,hemi,2,:) = [mean(temp(phase_bins>-pi/2 & phase_bins<pi/2)) mean(temp(phase_bins<-pi/2 | phase_bins>pi/2))];

        temp = squeeze(SO_phase_SO_HPC_MUA_spike_rate{nsession}(hemi,1,:));
        HPC_SO_FR(nsession,hemi,1,:) = [mean(temp(phase_bins>-pi/2 & phase_bins<pi/2)) mean(temp(phase_bins<-pi/2 | phase_bins>pi/2))];
        temp = squeeze(SO_phase_SO_HPC_MUA_spike_rate{nsession}(hemi,2,:));
        HPC_SO_FR(nsession,hemi,2,:) = [mean(temp(phase_bins>-pi/2 & phase_bins<pi/2)) mean(temp(phase_bins<-pi/2 | phase_bins>pi/2))];
    end
end

peak_vs_trough_ipsi = [];
peak_vs_trough_contra = [];
peak_vs_trough_ipsi(:,1) =  mean([squeeze(V1_SO_FR(:,1,1,1)) squeeze(V1_SO_FR(:,2,2,1))],2);
peak_vs_trough_ipsi(:,2) =  mean([squeeze(V1_SO_FR(:,1,1,2)) squeeze(V1_SO_FR(:,2,2,2))],2);

peak_vs_trough_contra(:,1) =  mean([squeeze(V1_SO_FR(:,1,2,1)) squeeze(V1_SO_FR(:,2,1,1))],2);
peak_vs_trough_contra(:,2) =  mean([squeeze(V1_SO_FR(:,1,2,2)) squeeze(V1_SO_FR(:,2,1,2))],2);


p=[];
p(1) =signrank(peak_vs_trough_ipsi(:,1),peak_vs_trough_ipsi(:,2));
p(2) =signrank(peak_vs_trough_contra(:,1),peak_vs_trough_contra(:,2));

% peak_vs_trough =  squeeze(V1_SO_FR(:,2,1,:));
% peak_vs_trough =  squeeze(V1_SO_FR(:,2,2,:));


%%% Plot peak vs trough phases
% 1. Setup the figure and colors
peak_vs_trough = {peak_vs_trough_ipsi,peak_vs_trough_contra};
conditions = {'V1 ipsi','V1 contra'};

% fig = figure;
% fig.Name = 'V1 SO phase Peak vs Trough';
% sgtitle(fig.Name)
% fig.Position= [844 66 560 906];
% counter = 1;

for iplot = 1:2
    subplot(2,2,iplot+2)
    hold on;

    % Define colors
    color_col1 = [0.1 0.5 0.9];  % Blue for condition 1
    color_col2 = [0.9 0.4 0.1];  % Orange for condition 2
    color_link = [0.5 0.5 0.5];  % Gray linking lines

    % Get the number of pairs (rows)
    num_pairs = size(peak_vs_trough{iplot}, 1);

    jitter_col1 = (rand(num_pairs, 1) - 0.5) * 0.2;
    jitter_col2 = (rand(num_pairs, 1) - 0.5) * 0.2;

    % 2. Draw the linking lines (One line per pair)
    for i = 1:num_pairs
        % X-values are fixed at [1 2] (the categories)
        % Y-values are the data [value_at_1 value_at_2]
        plot([1 + jitter_col1(i), 2 + jitter_col2(i)], peak_vs_trough{iplot}(i, :), ...
            'Color', color_link, ...
            'LineWidth', 0.5);
    end

    % 3. Plot the scatter points on top
    % Scatter for Condition 1 (at X=1)
    scatter(ones(num_pairs, 1)+jitter_col1, peak_vs_trough{iplot}(:, 1), ...
        70, color_col1, 'filled', ...
        'DisplayName', 'Condition 1');

    % Scatter for Condition 2 (at X=2)
    scatter(2 * ones(num_pairs, 1)+jitter_col2, peak_vs_trough{iplot}(:, 2), ...
        70, color_col2, 'filled', ...
        'DisplayName', 'Condition 2');


    % 4. Formatting and Labels
    ylabel('Mean firing rate (Hz)');
    % xlabel('Condition');
    title(conditions{iplot});

    % Set X-axis to show categories 1 and 2 clearly
    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels({'Peak', 'Trough'}); % Use descriptive labels

    text(1.5,6,sprintf('p = %.3e',p(iplot)))
    % Clean up the axes
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
    grid on;

    legend('off'); % Turn off legend since the categories are labeled on the X-axis
    hold off;
end

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\SO phase plots best'),[])
% save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\SO phase plots'),[])
% save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\SO phase plots (2Hz)'),[])








%%%%%% Individual sessions
for nsession = 1:22
    options = experiment_info(nsession).experiment_ID;
    fig = figure;
    fig.Name = sprintf('%s SO phase MUA firing',options(1).experiment_ID)
    sgtitle(fig.Name)
    fig.Position= [844 66 560 906];
    counter = 1;

    for hemi = 1:2
        subplot(4,2,counter)
        plot(phase_bins,squeeze(SO_phase_MUA_spike_rate{nsession}(hemi,1,:)),'r');hold on;
        plot(phase_bins,squeeze(SO_phase_MUA_spike_rate{nsession}(hemi,2,:)),'b')
        xlabel('SO Phase');ylabel('FR (Hz)');
        set(gca,'TickDir','out','Box','off','FontSize',12)
        xticks([-pi -pi/2 0 pi/2 pi])
        xticklabels({'-π','-π/2','0','π/2','π'})
        title(sprintf('V1 %s SO Phase',hemispheres{hemi}))
        counter = counter+1;

        subplot(4,2,counter)
        plot(phase_bins,squeeze(SO_phase_ripple_MUA_spike_rate{nsession}(hemi,1,:)),'r');hold on;
        plot(phase_bins,squeeze(SO_phase_ripple_MUA_spike_rate{nsession}(hemi,2,:)),'b')
        xticks([-pi -pi/2 0 pi/2 pi])
        xticklabels({'-π','-π/2','0','π/2','π'})
        set(gca,'TickDir','out','Box','off','FontSize',12)
        title(sprintf('V1 %s SO Phase (ripples)',hemispheres{hemi}))
        counter = counter+1;
    end

    % counter = 1;
    for hemi = 1:2
        subplot(4,2,counter)
        plot(phase_bins,squeeze(SO_phase_HPC_MUA_spike_rate{nsession}(hemi,1,:)),'r');hold on;
        plot(phase_bins,squeeze(SO_phase_HPC_MUA_spike_rate{nsession}(hemi,2,:)),'b')
        xlabel('SO Phase');ylabel('FR (Hz)')
        xticks([-pi -pi/2 0 pi/2 pi])
        xticklabels({'-π','-π/2','0','π/2','π'})
        set(gca,'TickDir','out','Box','off','FontSize',12)
        title(sprintf('HPC %s SO Phase',hemispheres{hemi}))
        counter = counter+1;

        subplot(4,2,counter)
        plot(phase_bins,squeeze(SO_phase_ripple_HPC_MUA_spike_rate{nsession}(hemi,1,:)),'r');hold on;
        plot(phase_bins,squeeze(SO_phase_ripple_HPC_MUA_spike_rate{nsession}(hemi,2,:)),'b')
        xlabel('SO Phase');ylabel('FR (Hz)')
        set(gca,'TickDir','out','Box','off','FontSize',12)
        xticks([-pi -pi/2 0 pi/2 pi])
        xticklabels({'-π','-π/2','0','π/2','π'})
        title(sprintf('HPC %s SO Phase (ripples)',hemispheres{hemi}))
        counter = counter+1;
    end
end


load((fullfile(analysis_folder,'SO_phase_MUA.mat')));

for nsession = 1:22
    plot(  squeeze(SO_phase_ripple_HPC_MUA_spike_rate{nsession}(hemi,1,:)))

    squeeze(SO_phase_ripple_HPC_MUA_spike_rate{nsession}(hemi,hemi,:))

    squeeze(SO_phase_ripple_HPC_MUA_spike_rate{nsession}(hemi,hemi,:)

    SO_phase_MUA_spike_rate
end

%%
