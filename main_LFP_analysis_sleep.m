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
                    spindle_amplitude_LFP = zscore(envelope(signal,round(lfp.samplingRate/5),'rms'));
                    spindle_amplitude_LFP = smoothdata(spindle_amplitude_LFP,'gaussian',round(lfp.samplingRate/5)); % envelop amplitude of spindles
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
                end

                %%%%% Calculate Delta peak amplitude and timestamp per event
                for nevent = 1:size(slow_waves(probe_no).DOWN_ints,1)
                    % for nevent = 600:640
                    midpoint = mean(slow_waves(probe_no).timestamps(nevent));
                    % midpoint = mean(slow_waves(probe_no).DOWN_ints(nevent,:));

                    tidx = FindInInterval(tvec,[midpoint-0.1 midpoint+0.1]);
                    % tidx = FindInInterval(tvec,[slow_waves(probe_no).ints.UP(nevent,1)-0.3 slow_waves(probe_no).ints.UP(nevent,1)+0.3]);
                    tidx=tidx(1):tidx(end);
                    [~,idx]=min(abs(midpoint-tvec));

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
                    spindle_amplitude_HPC_LFP = zscore(envelope(signal,round(lfp.samplingRate/5),'rms'));
                    spindle_amplitude_HPC_LFP = smoothdata(spindle_amplitude_HPC_LFP,'gaussian',round(lfp.samplingRate/5)); % envelop amplitude of spindles
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
                end

                % [ordered_xcoord,~]=sort(LFP(probe_no).best_V1_xcoord);
                zscored_LFP = [];
                zscored_LFP = ripple_amplitude_LFP;
                for nevent = 1:length(ripples(probe_no).onset)

                    tidx = FindInInterval(tvec,[ripples(probe_no).onset(nevent)-0.05 ripples(probe_no).offset(nevent)]);
                    % tidx = FindInInterval(tvec,[slow_waves(probe_no).ints.UP(nevent,1)-0.3 slow_waves(probe_no).ints.UP(nevent,1)+0.3]);
                    tidx=tidx(1):tidx(end);
                    [~,idx]=min(abs(ripples(probe_no).peaktimes(nevent)-tvec));


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
                    [~,idx]=min(abs(ripples(probe_no).peaktimes(nevent)-tvec));


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
                        % if diff(UP_ints(nevent,:))>0.1
                        if nevent ==1  | nevent ==size(UP_ints,1)
                            tidx = FindInInterval(tvec_interp1,[UP_ints(nevent,1)-0.05 UP_ints(nevent,1)+0.05]);
                        elseif (UP_ints(nevent-1,2)+0.05 < UP_ints(nevent,1) & UP_ints(nevent+1,1)-0.05 > UP_ints(nevent,2))
                            tidx = FindInInterval(tvec_interp1,[UP_ints(nevent,1)-0.05 UP_ints(nevent,1)+0.05]);
                        elseif  UP_ints(nevent-1,2)+0.05 > UP_ints(nevent,1) & UP_ints(nevent+1,1)-0.05 > UP_ints(nevent,2)
                            tidx = FindInInterval(tvec_interp1,[UP_ints(nevent-1,2) UP_ints(nevent,2)+0.05]);
                        elseif  UP_ints(nevent-1,2)+0.05 < UP_ints(nevent,1) & UP_ints(nevent+1,1)-0.05 < UP_ints(nevent,2)
                            tidx = FindInInterval(tvec_interp1,[UP_ints(nevent,1)-0.05 UP_ints(nevent,2)]);
                        end

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
                        if nevent ==1  | nevent ==size(DOWN_ints,1)
                            tidx = FindInInterval(tvec_interp1,[DOWN_ints(nevent,1)-0.05 DOWN_ints(nevent,1)+0.05]);
                        elseif (DOWN_ints(nevent-1,2)+0.05 < DOWN_ints(nevent,1) & DOWN_ints(nevent+1,1)-0.05 > DOWN_ints(nevent,2))
                            tidx = FindInInterval(tvec_interp1,[DOWN_ints(nevent,1)-0.05 DOWN_ints(nevent,1)+0.05]);
                        elseif  DOWN_ints(nevent-1,2)+0.05 > DOWN_ints(nevent,1) & DOWN_ints(nevent+1,1)-0.05 > DOWN_ints(nevent,2)
                            tidx = FindInInterval(tvec_interp1,[DOWN_ints(nevent-1,2) DOWN_ints(nevent,2)+0.05]);
                        elseif  DOWN_ints(nevent-1,2)+0.05 < DOWN_ints(nevent,1) & DOWN_ints(nevent+1,1)-0.05 < DOWN_ints(nevent,2)
                            tidx = FindInInterval(tvec_interp1,[DOWN_ints(nevent,1)-0.05 DOWN_ints(nevent,2)]);
                        end
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
                        if nevent ==1  | nevent ==size(DOWN_ints,1)
                            tidx = FindInInterval(tvec_interp1,[DOWN_ints(nevent,1) DOWN_ints(nevent,1)+0.1]);
                        elseif DOWN_ints(nevent+1,1)-0.1 > DOWN_ints(nevent,2)
                            tidx = FindInInterval(tvec_interp1,[DOWN_ints(nevent,1) DOWN_ints(nevent,1)+0.1]);
                        elseif  DOWN_ints(nevent+1,1)-0.1 < DOWN_ints(nevent,2)
                            tidx = FindInInterval(tvec_interp1,[DOWN_ints(nevent,1) DOWN_ints(nevent,2)]);
                        end

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

                    for nevent = 1:size(DOWN_ints,1)
                        if nevent ==1  | nevent ==size(DOWN_ints,1)
                            tidx = FindInInterval(tvec,[DOWN_ints(nevent,1) DOWN_ints(nevent,1)+0.1]);
                        elseif DOWN_ints(nevent+1,1)-0.1 > DOWN_ints(nevent,2)
                            tidx = FindInInterval(tvec,[DOWN_ints(nevent,1) DOWN_ints(nevent,1)+0.1]);
                        elseif  DOWN_ints(nevent+1,1)-0.1 < DOWN_ints(nevent,2)
                            tidx = FindInInterval(tvec,[DOWN_ints(nevent,1) DOWN_ints(nevent,2)]);
                        end

                        for nchannel = 1:length(slow_waves(probe_no).shank_id)
                            phase_DOWN(nchannel,nevent)=angle(mean(exp(1i*SO_phase_LFP(tidx(1):tidx(end),nchannel)))); % phase
                        end
                    end
                    slow_waves(probe_no).mean_phase_DOWN = phase_DOWN;

                    % UP
                    event_tidx = [];
                    event_index = [];
                    % Grab 20ms downsampled idx when events happned
                    for nevent = 1:size(UP_ints,1)
                        if nevent ==1  | nevent ==size(UP_ints,1)
                            tidx = FindInInterval(tvec_interp1,[UP_ints(nevent,1) UP_ints(nevent,1)+0.1]);
                        elseif UP_ints(nevent+1,1)-0.1 > UP_ints(nevent,2)
                            tidx = FindInInterval(tvec_interp1,[UP_ints(nevent,1) UP_ints(nevent,1)+0.1]);
                        elseif  UP_ints(nevent+1,1)-0.1 < UP_ints(nevent,2)
                            tidx = FindInInterval(tvec_interp1,[UP_ints(nevent,1) UP_ints(nevent,2)]);
                        end
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


                    for nevent = 1:size(UP_ints,1)
                        if nevent ==1  | nevent ==size(UP_ints,1)
                            tidx = FindInInterval(tvec,[UP_ints(nevent,1) UP_ints(nevent,1)+0.1]);
                        elseif UP_ints(nevent+1,1)-0.1 > UP_ints(nevent,2)
                            tidx = FindInInterval(tvec,[UP_ints(nevent,1) UP_ints(nevent,1)+0.1]);
                        elseif  UP_ints(nevent+1,1)-0.1 < UP_ints(nevent,2)
                            tidx = FindInInterval(tvec,[UP_ints(nevent,1) UP_ints(nevent,2)]);
                        end

                        for nchannel = 1:length(slow_waves(probe_no).shank_id)
                            phase_UP(nchannel,nevent)=angle(mean(exp(1i*SO_phase_LFP(tidx(1):tidx(end),nchannel)))); % phase
                        end
                    end
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
                        % if diff(UP_ints(nevent,:))>0.1
                        if nevent ==1  | nevent ==size(UP_ints,1)
                            tidx = FindInInterval(tvec,[UP_ints(nevent,1)-0.05 UP_ints(nevent,1)+0.05]);
                        elseif (UP_ints(nevent-1,2)+0.05 < UP_ints(nevent,1) & UP_ints(nevent+1,1)-0.05 > UP_ints(nevent,2))
                            tidx = FindInInterval(tvec,[UP_ints(nevent,1)-0.05 UP_ints(nevent,1)+0.05]);
                        elseif  UP_ints(nevent-1,2)+0.05 > UP_ints(nevent,1) & UP_ints(nevent+1,1)-0.05 > UP_ints(nevent,2)
                            tidx = FindInInterval(tvec,[UP_ints(nevent-1,2) UP_ints(nevent,2)+0.05]);
                        elseif  UP_ints(nevent-1,2)+0.05 < UP_ints(nevent,1) & UP_ints(nevent+1,1)-0.05 < UP_ints(nevent,2)
                            tidx = FindInInterval(tvec,[UP_ints(nevent,1)-0.05 UP_ints(nevent,2)]);
                        end

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
                        if nevent ==1  | nevent ==size(DOWN_ints,1)
                            tidx = FindInInterval(tvec,[DOWN_ints(nevent,1)-0.05 DOWN_ints(nevent,1)+0.05]);
                        elseif (DOWN_ints(nevent-1,2)+0.05 < DOWN_ints(nevent,1) & DOWN_ints(nevent+1,1)-0.05 > DOWN_ints(nevent,2))
                            tidx = FindInInterval(tvec,[DOWN_ints(nevent,1)-0.05 DOWN_ints(nevent,1)+0.05]);
                        elseif  DOWN_ints(nevent-1,2)+0.05 > DOWN_ints(nevent,1) & DOWN_ints(nevent+1,1)-0.05 > DOWN_ints(nevent,2)
                            tidx = FindInInterval(tvec,[DOWN_ints(nevent-1,2) DOWN_ints(nevent,2)+0.05]);
                        elseif  DOWN_ints(nevent-1,2)+0.05 < DOWN_ints(nevent,1) & DOWN_ints(nevent+1,1)-0.05 < DOWN_ints(nevent,2)
                            tidx = FindInInterval(tvec,[DOWN_ints(nevent,1)-0.05 DOWN_ints(nevent,2)]);
                        end

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
                    slow_waves(probe_no).pd_UD = plv_UD;

                    slow_waves(probe_no).plv_DU = plv_DU;
                    slow_waves(probe_no).pd_DU = plv_DU;

                    slow_waves(probe_no).cortex_speed_UD = cortex_speed_UD;
                    slow_waves(probe_no).HPC_speed_UD = HPC_speed_UD;

                    slow_waves(probe_no).cortex_speed_DU = cortex_speed_DU;
                    slow_waves(probe_no).HPC_speed_DU = HPC_speed_DU;

                    % phase and amplitude
                    SO_phase_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    SO_phase_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));

                    spindle_phase_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    spindle_phase_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).onset));

                    SO_phase_spindle_onset = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
                    SO_phase_spindle_peaktime = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));

                    % SO_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    % spindle_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    % SO_spindles_coupling = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
                    ripple_peak_amplitude = nan(length(ripples(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    ripple_onset_amplitude = nan(length(ripples(probe_no).shank_id), length(ripples(probe_no).onset));
                    spindle_peak_amplitude = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).peaktimes));
                    spindle_onset_amplitude = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));

                    for nevent = 1:length(ripples(probe_no).onset)
                        [~,tidx]=min(abs(ripples(probe_no).peaktimes(nevent)-tvec));
                        % Phase
                        SO_phase_ripple_peaktime(:,nevent) = SO_phase_LFP(tidx,:);
                        spindle_phase_ripple_peaktime(:,nevent) = spindle_phase_LFP(tidx,:);
                        ripple_peak_amplitude(:,nevent) = ripple_amplitude_LFP(tidx,:);

                        [~,tidx]=min(abs(ripples(probe_no).onset(nevent)-tvec));
                        % Phase
                        SO_phase_ripple_onset(:,nevent) = SO_phase_LFP(tidx,:);
                        spindle_phase_ripple_onset(:,nevent) = spindle_phase_LFP(tidx,:);
                        ripple_onset_amplitude(:,nevent) = ripple_amplitude_LFP(tidx,:);
                    end

                    ripples(probe_no).SO_phase_ripple_peaktime = SO_phase_ripple_peaktime;
                    ripples(probe_no).spindle_phase_ripple_peaktime = spindle_phase_ripple_peaktime;
                    ripples(probe_no).ripple_peak_amplitude = ripple_peak_amplitude;

                    ripples(probe_no).SO_phase_ripple_onset = SO_phase_ripple_onset;
                    ripples(probe_no).spindle_phase_ripple_onset = spindle_phase_ripple_onset;
                    ripples(probe_no).ripple_onset_amplitude = ripple_onset_amplitude;

                    if ~isempty(spindles(probe_no).onset)
                        for nevent = 1:length(spindles(probe_no).onset)
                            [~,tidx]=min(abs(spindles(probe_no).peaktimes(nevent)-tvec));
                            % Phase
                            SO_phase_spindle_peaktime(:,nevent) = SO_phase_LFP(tidx,:);
                            spindle_peak_amplitude(:,nevent) = spindle_amplitude_LFP(tidx,:);

                            [~,tidx]=min(abs(spindles(probe_no).onset(nevent)-tvec));
                            % Phase
                            SO_phase_spindle_onset(:,nevent) = SO_phase_LFP(tidx,:);
                            spindle_onset_amplitude(:,nevent) = spindle_amplitude_LFP(tidx,:);
                        end

                        spindles(probe_no).SO_phase_spindle_peaktime = SO_phase_spindle_peaktime;
                        spindles(probe_no).spindle_peak_amplitude = spindle_peak_amplitude;

                        spindles(probe_no).SO_phase_spindle_onset = SO_phase_spindle_onset;
                        spindles(probe_no).spindle_onset_amplitude = spindle_onset_amplitude;
                    end
                else
                    spindles(probe_no).SO_phase_spindle_peaktime = [];
                    spindles(probe_no).spindle_peak_amplitude = [];

                    spindles(probe_no).SO_phase_spindle_onset = [];
                    spindles(probe_no).spindle_onset_amplitude = [];
                end
            end



        end


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
        else
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'),'slow_waves');
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'),'spindles');
            % save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_markov_events.mat'),'slow_waves_markov');
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'),'ripples');
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
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_markov_events.mat'));
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

            if isempty(slow_waves_markov(probe_no).UP_DOWN_index)
                continue
            end

            if size(slow_waves_markov(probe_no).UP_DOWN_index,1) < 10
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

        disp('cortical wave traveling direction analysis started')
        tic
        %%%%%%%%%%%% Cortical wave direction during DOWN state peak
        % -1 is posterior -> anterior, 0 is no delay or noisy delay and 1 is anterior -> posterior
        % https://www.nature.com/articles/s41467-019-10327-5 levenstein et
        % al used 0.5 to 8 Hz for slow wave UP DOWN detection
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
            DOWN_peaks_zscore= [];

            slow_waves_markov(probe_no).DOWN_peaks_zscore = [];
            slow_waves_markov(probe_no).DOWN_peaktimes = [];
            slow_waves_markov(probe_no).DOWN_peaks_latency = [];
            slow_waves_markov(probe_no).DOWN_traveling = [];
            slow_waves_markov(probe_no).probe_hemisphere = [];
            slow_waves_markov(probe_no).shank_id = [];

            if isempty(slow_waves_markov(probe_no).UP_DOWN_index)
                continue
            end

            if size(slow_waves_markov(probe_no).UP_DOWN_index,1) < 10
                continue
            end

            if isfield(LFP(probe_no),'average_V1_xcoord') & ~isempty(behavioural_state_merged.SWS)% if exist best V1 channel for sleep
                probe_id = [];

                if length(LFP)==1
                    lfp.data= [LFP(probe_no).average_V1'];
                    probe_hemisphere = probe_no*ones(1,length(LFP(probe_no).average_V1_shank_id));
                    slow_waves_markov(probe_no).shank_id = [LFP(probe_no).average_V1_shank_id];
                else
                    lfp.data= [LFP(1).average_V1' LFP(2).average_V1'];
                    probe_hemisphere = [ones(1,length(LFP(1).average_V1_shank_id)) 2*ones(1,length(LFP(2).average_V1_shank_id))];
                    slow_waves_markov(probe_no).shank_id = [LFP(1).average_V1_shank_id LFP(2).average_V1_shank_id];
                end

                slow_waves_markov(probe_no).probe_hemisphere = probe_hemisphere;


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

                    for nShank=1:length(probe_hemisphere) % across shanks from both probes if using two probes

                        [~,peak_id] = findpeaks(zscored_LFP(tidx(1):tidx(end),nShank));

                        [~,temp]=min(abs(midpoint-tvec(tidx(peak_id))));
                        if ~isempty(temp)
                            DOWN_peaks_shank(nShank,nevent) = tvec(tidx(peak_id(temp)));
                            DOWN_peaks_zscore(nShank,nevent) = zscored_LFP(tidx(peak_id(temp)),nShank);
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
                    for h = unique(probe_hemisphere)
                        DOWN_peaks_shank_temp = DOWN_peaks_shank(probe_hemisphere==h,nevent);

                        if sum(~isnan(DOWN_peaks_shank_temp))==4 % if delta peaks on four shanks
                            peaks_latency(h,nevent) = mean(diff(DOWN_peaks_shank_temp),'omitnan');
                        elseif sum(~isnan(DOWN_peaks_shank_temp))==3 % if delta peaks on three shanks
                            skipped_shank= diff(LFP(h).average_V1_shank_id(~isnan(DOWN_peaks_shank_temp)))>1;
                            if sum(skipped_shank)>0 % if delta peak skipped one shank

                                if skipped_shank(1)==1 % if shanks [1 2 4] then delay using 1 -> 2
                                    peaks_latency(h,nevent) = DOWN_peaks_shank_temp(1)-DOWN_peaks_shank_temp(2);
                                elseif skipped_shank(2) == 1 % if shanks [1 3 4] then delay using 3 -> 4
                                    peaks_latency(h,nevent) = DOWN_peaks_shank_temp(2)-DOWN_peaks_shank_temp(3);
                                end
                            else
                                peaks_latency(h,nevent) = mean(diff(DOWN_peaks_shank_temp),'omitnan');
                            end
                        elseif sum(~isnan(DOWN_peaks_shank_temp))==2 % if delta peaks on two shanks (latency divide by number of shanks skipped to caluclate mean latency per 250 micron)
                            peaks_latency(h,nevent) = diff(DOWN_peaks_shank_temp(1:2))/abs(LFP(h).average_V1_shank_id(2)-LFP(h).average_V1_shank_id(1));
                        else % only one peak. Can't calculate latency
                            peaks_latency(h,nevent) =nan;
                        end

                        % From https://www.jneurosci.org/content/34/26/8875#sec-2
                        % speed is ~40 milimeter per seconds
                        % 6.25 miliseconds to travel 250 micrometer (rough shank spacing)
                        % putatively set the minimum delay threshold to be 3
                        % miliseconds.
                        if h==1

                            % if direction of mean latency is consistent with the latency more than half
                            % of the shank latency
                            % (i.e. if four shanks, at least 2 jumps between three shanks should be in the same direction as the mean latency)
                            % (if three shanks, then still 2 jumps needed)
                            if peaks_latency(nevent)>0.003 & sum(DOWN_peaks_shank_temp(2:end)-DOWN_peaks_shank_temp(1)>0)>=length(LFP(h).average_V1_shank_id)/2
                                DOWN_traveling(h,nevent) = 1; % anterior to posterior
                            elseif peaks_latency(nevent)<-0.003 & sum(DOWN_peaks_shank_temp(2:end)-DOWN_peaks_shank_temp(1)<0)>=length(LFP(h).average_V1_shank_id)/2
                                DOWN_traveling(h,nevent) = -1; % posterior to anterior
                            else
                                DOWN_traveling(h,nevent) = 0; % noisy or standing wave?
                            end
                            % DOWN_peaks_shank(:,nevent) - DOWN_peaks_shank(1,nevent)
                        elseif h==2
                            if peaks_latency(h,nevent)>0.003 & sum(DOWN_peaks_shank_temp(2:end)-DOWN_peaks_shank_temp(1)>0)>=length(LFP(h).average_V1_shank_id)/2
                                DOWN_traveling(h,nevent) = -1; % posterior to anterior
                            elseif peaks_latency(h,nevent)<-0.003 & sum(DOWN_peaks_shank_temp(2:end)-DOWN_peaks_shank_temp(1)<0)>=length(LFP(h).average_V1_shank_id)/2
                                DOWN_traveling(h,nevent) = 1; % anterior to posterior
                            else
                                DOWN_traveling(h,nevent) = 0; % noisy or standing wave?
                            end
                        end
                    end
                    % nexttile
                    % hold on;xline(tvec(idx));
                    % plot(tvec(tidx),zscored_LFP(tidx,:));
                    % % hold on;xline(median(tvec(tidx(1):tidx(end))))
                    % xline(DOWN_peaks_shank(:,nevent)')
                end


            end


            slow_waves_markov(probe_no).DOWN_peaktimes = DOWN_peaks_shank;
            slow_waves_markov(probe_no).DOWN_peaks_zscore = DOWN_peaks_zscore;
            slow_waves_markov(probe_no).DOWN_peaks_latency = peaks_latency;
            slow_waves_markov(probe_no).DOWN_traveling = DOWN_traveling;
        end
        disp('cortical wave traveling direction analysis finished')
        toc
        %

        disp('Sharp wave ripple traveling direction analysis started')
        tic
        %%%%%%%%%%%% Cortical wave direction during DOWN state peak
        % -1 is posterior -> anterior, 0 is no delay or noisy delay and 1 is anterior -> posterior
        filterparms.deltafilter = [0.5 8];%heuristically defined.  room for improvement here.

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
            ripples(probe_no).sharp_wave_zscore = [];

            ripples(probe_no).SWR_peaktimes = [];
            ripples(probe_no).SWR_zscore = [];
            ripples(probe_no).SWR_traveling = [];

            ripples(probe_no).probe_hemisphere = [];
            ripples(probe_no).shank_id = [];

            if isfield(LFP(probe_no),'best_HPC') % if exist best HPC channels
                probe_id = [];

                if length(LFP)==1
                    lfp.data= [LFP(probe_no).best_HPC'];
                    probe_hemisphere = probe_no*ones(1,length(LFP(probe_no).best_HPC_shank_id));
                    ripples(probe_no).shank_id = [LFP(1).best_HPC_shank_id LFP(2).best_HPC_shank_id];
                else
                    lfp.data= [LFP(1).best_HPC' LFP(2).best_HPC'];
                    probe_hemisphere = [ones(1,length(LFP(1).best_HPC_shank_id)) 2*ones(1,length(LFP(2).best_HPC_shank_id))];
                    ripples(probe_no).shank_id = [LFP(1).best_HPC_shank_id LFP(2).best_HPC_shank_id];
                end

                ripples(probe_no).probe_hemisphere = probe_hemisphere;

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

                % [ordered_xcoord,~]=sort(LFP(probe_no).best_V1_xcoord);

                for nevent = 1:length(ripples(probe_no).onset)

                    tidx = FindInInterval(tvec,[ripples(probe_no).onset(nevent)-0.05 ripples(probe_no).offset(nevent)]);
                    % tidx = FindInInterval(tvec,[slow_waves(probe_no).ints.UP(nevent,1)-0.3 slow_waves(probe_no).ints.UP(nevent,1)+0.3]);
                    tidx=tidx(1):tidx(end);
                    [~,idx]=min(abs(ripples(probe_no).onset(nevent)-tvec));


                    for nShank=1:length(probe_hemisphere)

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
                    for h = unique(probe_hemisphere)
                        sharp_wave_peaks_shank_temp = sharp_wave_peaks_shank(probe_hemisphere==h,nevent);

                        if sum(~isnan(sharp_wave_peaks_shank_temp))==4 % if delta peaks on four shanks
                            peaks_latency(h,nevent) = mean(diff(sharp_wave_peaks_shank_temp),'omitnan');
                        elseif sum(~isnan(sharp_wave_peaks_shank_temp))==3 % if delta peaks on three shanks
                            skipped_shank= diff(LFP(probe_no).best_HPC_shank_id(~isnan(sharp_wave_peaks_shank_temp)))>1;
                            if sum(skipped_shank)>0 % if delta peak skipped one shank

                                if skipped_shank(1)==1 % if shanks [1 2 4] then delay using 1 -> 2
                                    peaks_latency(h,nevent) = sharp_wave_peaks_shank_temp(1)-sharp_wave_peaks_shank_temp(2);
                                elseif skipped_shank(2) == 1 % if shanks [1 3 4] then delay using 3 -> 4
                                    peaks_latency(h,nevent) = sharp_wave_peaks_shank_temp(2)-sharp_wave_peaks_shank_temp(3);
                                end
                            else
                                peaks_latency(h,nevent) = mean(diff(sharp_wave_peaks_shank_temp),'omitnan');
                            end
                        elseif sum(~isnan(sharp_wave_peaks_shank_temp))==2 % if delta peaks on two shanks
                            peaks_latency(h,nevent) = diff(sharp_wave_peaks_shank_temp(1:2))/abs(LFP(h).best_HPC_shank_id(2)-LFP(h).best_HPC_shank_id(1));
                        else % only one peak. Can't calculate latency
                            peaks_latency(h,nevent) =nan;
                        end

                        % From Patel et al. (2013) https://pmc.ncbi.nlm.nih.gov/articles/PMC3807028/#sec2
                        % ripple propogation speed roughly 0.35 m/s or 0.7
                        % ms per 250 micron
                        % putatively set the minimum delay threshold to be 2
                        % miliseconds.
                        if h==1

                            % if direction of mean latency is consistent with the latency more than half
                            % of the shank latency
                            % (i.e. if four shanks, at least 2 jumps between three shanks should be in the same direction as the mean latency)
                            % (if three shanks, then still 2 jumps needed)
                            if peaks_latency(h,nevent)>0.001 & sum(sharp_wave_peaks_shank_temp(2:end)-sharp_wave_peaks_shank_temp(1)>0)>=length(LFP(h).best_HPC_shank_id)/2
                                wave_traveling(h,nevent) = 1; % anterior to posterior
                            elseif peaks_latency(h,nevent)<-0.001 & sum(sharp_wave_peaks_shank_temp(2:end)-sharp_wave_peaks_shank_temp(1)<0)>=length(LFP(h).best_HPC_shank_id)/2
                                wave_traveling(h,nevent) = -1; % posterior to anterior
                            else
                                wave_traveling(h,nevent) = 0; % noisy or standing wave?
                            end
                            % DOWN_peaks_shank(:,nevent) - DOWN_peaks_shank(1,nevent)
                        elseif h==2
                            if peaks_latency(h,nevent)>0.001 & sum(sharp_wave_peaks_shank_temp(2:end)-sharp_wave_peaks_shank_temp(1)>0)>=length(LFP(h).best_HPC_shank_id)/2
                                wave_traveling(h,nevent) = -1; % posterior to anterior
                            elseif peaks_latency(h,nevent)<-0.001 & sum(sharp_wave_peaks_shank_temp(2:end)-sharp_wave_peaks_shank_temp(1)<0)>=length(LFP(h).best_HPC_shank_id)/2
                                wave_traveling(h,nevent) = 1; % anterior to posterior
                            else
                                wave_traveling(h,nevent) = 0; % noisy or standing wave?
                            end
                        end
                    end
                    % nexttile
                    % hold on;xline(tvec(idx));
                    % plot(tvec(tidx),zscored_LFP(tidx,:));
                    % % hold on;xline(median(tvec(tidx(1):tidx(end))))
                    % xline(sharp_wave_peaks_shank(:,nevent)')
                end

                ripples(probe_no).SWR_peaktimes = sharp_wave_peaks_shank;
                ripples(probe_no).SWR_latency = peaks_latency;
                ripples(probe_no).SWR_zscore = sharp_wave_zscore_shank;
                ripples(probe_no).SWR_traveling = wave_traveling;


                sharp_wave_zscore_shank = [];
                sharp_wave_peaks_shank = [];


                deltaLFP = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);

                zscored_LFP = [];
                zscored_LFP = zscore(deltaLFP.amp);

                for nevent = 1:length(ripples(probe_no).onset)
                    %                 ripples(probe_no).peaktimes(nevent)

                    tidx = FindInInterval(tvec,[ripples(probe_no).onset(nevent)-0.05 ripples(probe_no).offset(nevent)]);
                    % tidx = FindInInterval(tvec,[slow_waves(probe_no).ints.UP(nevent,1)-0.3 slow_waves(probe_no).ints.UP(nevent,1)+0.3]);
                    tidx=tidx(1):tidx(end);
                    [~,idx]=min(abs(ripples(probe_no).peaktimes(nevent)-tvec));


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

        disp('Sharp wave ripple traveling direction analysis started')
        toc


        zscored_LFP = [];
        phase_LFP = [];
        amplitude_LFP = [];

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



