
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
            
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            clusters=clusters_ks3;
        elseif contains(stimulus_name{n},'Sleep')
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
%             load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
            % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
            
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state.mat'));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'));
            % load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
%             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
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
        %         
        
        % From cell structure back to spike times and spike id
        session_clusters_RUN.spike_id=vertcat(session_clusters_RUN.spike_id{:});
        session_clusters_RUN.spike_times=vertcat(session_clusters_RUN.spike_times{:});
        [session_clusters_RUN.spike_times,index] =sort(session_clusters_RUN.spike_times);
        session_clusters_RUN.spike_id=session_clusters_RUN.spike_id(index);

        session_clusters.spike_id=vertcat(session_clusters.spike_id{:});
        session_clusters.spike_times=vertcat(session_clusters.spike_times{:});
        [session_clusters.spike_times,index] =sort(session_clusters.spike_times);
        session_clusters.spike_id=session_clusters.spike_id(index);

        clusters_combined= session_clusters;

        if length(clusters) > 1
            MUA_combined = combine_clusters_from_multiple_probes(clusters(1),clusters(2));
        else
            MUA_combined = clusters;
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


        metric_param =[];
        % metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
        metric_param.region = @(x) contains(x,'V1_L');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        V1_spikes{1}=[selected_clusters.spike_id selected_clusters.spike_times];
        metric_param.region = @(x) contains(x,'V1_R');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        V1_spikes{2}=[selected_clusters.spike_id selected_clusters.spike_times];
        metric_param.region = @(x) contains(x,'HPC_L');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        HPC_spikes{1}=[selected_clusters.spike_id selected_clusters.spike_times];
        metric_param.region = @(x) contains(x,'HPC_R');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        HPC_spikes{2}=[selected_clusters.spike_id selected_clusters.spike_times];
        metric_param.region = @(x) contains(x,'HPC');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        HPC_spikes{3}=[selected_clusters.spike_id selected_clusters.spike_times];

        [intervals]= OverlapIntervals(behavioural_state(1).SWS,behavioural_state(2).SWS);


        % HPC frames
        for nprobe = 1:2
            probe_no = session_info(n).probe(nprobe).probe_id+1;
            if ~isempty(behavioural_state(probe_no).SWS)
                [HPC_frames] = detect_candidate_frames_masa(LFP(probe_no).tvec,V1_spikes{1}(:,2),behavioural_state(probe_no).SWS,0.01,options);
            end
        end




        % test detection
        nprobe = 1;
        if ~isempty(behavioural_state(probe_no).SWS)
            best_channel = find(LFP(nprobe).best_V1_channel==slow_waves(nprobe).best_channel);
            temp_SW = DetectSlowWaves_masa('time',tvec,'lfp',LFP(probe_no).best_V1(best_channel,:),'spikes',V1_clusters,'NREMInts',SWS);
            [temp_spindles] = FindSpindles_masa(LFP(probe_no).best_V1(best_channel,:),LFP(probe_no).tvec','behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
                'noise',[],'passband',[9 17],'thresholds',[1 3],'show','off');
        elseif isfield(LFP(nprobe),'best_V1_high_freq')
            [~,best_channel] = max(LFP(nprobe).best_V1_high_freq);
            [temp_spindles] = FindSpindles_masa(LFP(probe_no).best_V1_high_freq(best_channel,:),LFP(probe_no).tvec','behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
                'noise',[],'passband',[9 17],'thresholds',[1 3],'show','off');
        else
            [~,best_channel] = max(LFP(nprobe).L5_power(:,7));
            [temp_spindles] = FindSpindles_masa(LFP(probe_no).L5(best_channel,:),LFP(probe_no).tvec','behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
                'noise',[],'passband',[9 17],'thresholds',[1 3],'show','off');
        end
    end
end
