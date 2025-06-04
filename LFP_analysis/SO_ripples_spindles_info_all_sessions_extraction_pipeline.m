function SO_ripples_spindles_info_all_sessions_extraction_pipeline

%% extract all events for sleep analysis (spiking but also event info)

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))


clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar 
% experiment_info=experiment_info([4 5 6 ]);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'SleepChronic';
% experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
% 1:length(experiment_info)
% [1 2 3 4 6 7 8 9 10 12 14]
session_count = 0;

slow_waves_all = struct();
ripples_all = struct();
spindles_all = struct();
behavioural_state_merged_all = struct();

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


    if length(stimulus_name)>1
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

    session_count = session_count + 1;
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
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
% %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP');
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
% 
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks4%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%         clusters=clusters_ks4;
    elseif contains(stimulus_name{n},'Sleep')
%         load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
        %             load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state_merged.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_markov_events.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
        %             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters=clusters_ks4;
    else
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
%         load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters=clusters_ks4;
    end
    if isfield(slow_waves,'detectorinfo')
        ripples = rmfield(ripples, 'detectorinfo');
        spindles = rmfield(spindles, 'detectorinfo');
        slow_waves = rmfield(slow_waves, 'detectorinfo');
    end
    if isfield(slow_waves,'DOWN_peaks_latency')
        slow_waves = rmfield(slow_waves, 'DOWN_peaks_latency');
        slow_waves = rmfield(slow_waves, 'DOWN_traveling');
        slow_waves = rmfield(slow_waves, 'gammaspikecorr');
        slow_waves = rmfield(slow_waves, 'deltaspikecorr');
        slow_waves = rmfield(slow_waves, 'deltagammacorr');        
    end


    if isfield(ripples,'SWR_traveling')
        ripples = rmfield(ripples, 'SWR_traveling');
        ripples = rmfield(ripples, 'SWR_latency');
    end
    % slow_waves_markov = rmfield(slow_waves_markov, 'NREM_t');
    % slow_waves_markov = rmfield(slow_waves_markov, 'spike_count');
    % slow_waves_markov = rmfield(slow_waves_markov, 'alpha_t');
    % slow_waves_markov = rmfield(slow_waves_markov, 'gamma_t');
    % slow_waves_markov = rmfield(slow_waves_markov, 'beta_t');
    % slow_waves_markov = rmfield(slow_waves_markov, 'xi_t');
    % % 
    % for nprobe = 1:2
    %     slow_waves(nprobe).power = [];
    %     slow_waves(nprobe).frequency = [];
    %     slow_waves(nprobe).PSD_slope = [];
    %     slow_waves(nprobe).timebin_edges = [];
    %     slow_waves(nprobe).DOWN_PSD_slope = [];
    %     slow_waves(nprobe).UP_PSD_slope = [];
    %     slow_waves(nprobe).DOWN_delta_power = [];
    %     slow_waves(nprobe).DOWN_peaks_zscore = [];
    %     slow_waves(nprobe).DOWN_peaktimes = [];
    %     slow_waves(nprobe).DOWN_peaks_latency = [];
    %     slow_waves(nprobe).DOWN_traveling = [];
    %     slow_waves(nprobe).probe_hemisphere = [];
    %     slow_waves(nprobe).shank_id = [];
    %     % slow_waves.power = [];
    % 
    % 
    %     slow_waves(nprobe).power = [];
    %     slow_waves(nprobe).frequency = [];
    %     slow_waves(nprobe).PSD_slope = [];
    %     slow_waves(nprobe).timebin_edges = [];
    %     slow_waves(nprobe).DOWN_PSD_slope = [];
    %     slow_waves(nprobe).UP_PSD_slope = [];
    %     slow_waves(nprobe).DOWN_delta_power = [];
    %     slow_waves(nprobe).DOWN_peaks_zscore = [];
    %     slow_waves(nprobe).DOWN_peaktimes = [];
    %     slow_waves(nprobe).DOWN_peaks_latency = [];
    %     slow_waves(nprobe).DOWN_traveling = [];
    %     slow_waves(nprobe).probe_hemisphere = [];
    %     slow_waves(nprobe).shank_id = [];
    % end
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

    clusters_combined= session_clusters; % SUA from both probes

    if length(clusters) > 1
        MUA_combined = combine_clusters_from_multiple_probes(clusters(1),clusters(2)); % all clusters for MUA
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


    % put all behavioural states in one struct
    field_names = fieldnames(behavioural_state_merged);

    %     for nprobe = 1:length(ripples)
    %         probe_no = session_info(n).probe(nprobe).probe_hemisphere;
    for iField =1:length(field_names)
        behavioural_state_merged_all.(field_names{iField}){session_count} = behavioural_state_merged.(field_names{iField});
    end
    %     end

    % put all ripples in one struct
    field_names = fieldnames(ripples);

    for nprobe = 1:length(ripples)
        probe_no = session_info(n).probe(nprobe).probe_hemisphere;
        ripples_all(probe_no).subject(session_count,:) = options.SUBJECT;
        ripples_all(probe_no).session(session_count,:) = options.SESSION;
        ripples_all(probe_no).session_day(session_count) = iDate;

        for iField =1:length(field_names)
            if session_count == 1
                if  ismember(field_names{iField},{'sharp_wave_zscore','sharp_wave_peaktimes','SWR_zscore','SWR_peaktimes','shank_id','probe_hemisphere',...
                        'amp_corr','xcorr_lag','xcorr_r','pd','plv','SO_phase_ripple_peaktime','spindle_phase_ripple_peaktime','ripple_peak_amplitude',...
                        'SO_phase_ripple_onset','spindle_phase_ripple_onset','ripple_onset_amplitude'})
                    ripples_all(probe_no).(field_names{iField}){session_count} = ripples(nprobe).(field_names{iField});
                else
                    ripples_all(probe_no).(field_names{iField}) = ripples(nprobe).(field_names{iField});
                end


            else
                if  ismember(field_names{iField},{'sharp_wave_zscore','sharp_wave_peaktimes','SWR_zscore','SWR_peaktimes','shank_id','probe_hemisphere',...
                        'amp_corr','xcorr_lag','xcorr_r','pd','plv','SO_phase_ripple_peaktime','spindle_phase_ripple_peaktime','ripple_peak_amplitude',...
                        'SO_phase_ripple_onset','spindle_phase_ripple_onset','ripple_onset_amplitude'})
                    ripples_all(probe_no).(field_names{iField}){session_count} = ripples(nprobe).(field_names{iField});
                else
                    A = ripples_all(probe_no).(field_names{iField});
                    B = ripples(nprobe).(field_names{iField});
                    try
                        ripples_all(probe_no).(field_names{iField}) = [A;B];
                    catch
                        ripples_all(probe_no).(field_names{iField}) = [A B];
                    end
                end
            end
        end

        session_count_events = repmat(session_count,size(ripples(nprobe).onset));
        if session_count == 1
            ripples_all(probe_no).session_count = session_count_events;
        else
            ripples_all(probe_no).session_count = [ripples_all(probe_no).session_count;session_count_events];
        end
    end

    % put all spindles in one struct
    field_names = fieldnames(spindles);

    for nprobe = 1:length(spindles)
        probe_no = session_info(n).probe(nprobe).probe_hemisphere;
        spindles_all(probe_no).subject(session_count,:) = options.SUBJECT;
        spindles_all(probe_no).session(session_count,:) = options.SESSION;
        spindles_all(probe_no).session_day(session_count) = iDate;
        
        for iField =1:length(field_names)
            if session_count == 1
                if  ismember(field_names{iField},{'sharp_wave_zscore','sharp_wave_peaktimes','SWR_zscore','SWR_peaktimes','shank_id','probe_hemisphere',...
                        'amp_corr','xcorr_lag','xcorr_r','pd','plv','SO_phase_spindle_peaktime','SO_phase_spindle_onset','spindle_peak_amplitude','spindle_onset_amplitude'})
                    spindles_all(probe_no).(field_names{iField}){session_count} = spindles(nprobe).(field_names{iField});
                else
                    spindles_all(probe_no).(field_names{iField}) = spindles(nprobe).(field_names{iField});

                end
            else
                if  ismember(field_names{iField},{'sharp_wave_zscore','sharp_wave_peaktimes','SWR_zscore','SWR_peaktimes','shank_id','probe_hemisphere',...
                        'amp_corr','xcorr_lag','xcorr_r','pd','plv','SO_phase_spindle_peaktime','SO_phase_spindle_onset','spindle_peak_amplitude','spindle_onset_amplitude'})
                    spindles_all(probe_no).(field_names{iField}){session_count} = spindles(nprobe).(field_names{iField});
                else
                    A = spindles_all(probe_no).(field_names{iField});
                    B = spindles(nprobe).(field_names{iField});
                    try
                        spindles_all(probe_no).(field_names{iField}) = [A;B];
                    catch
                        spindles_all(probe_no).(field_names{iField}) = [A B];
                    end
                end
            end
        end

        session_count_events = repmat(session_count,size(spindles(nprobe).onset));
        if session_count == 1
            spindles_all(probe_no).session_count = session_count_events;
        else
            spindles_all(probe_no).session_count = [spindles_all(probe_no).session_count;session_count_events];
        end
    end

    % put all slow waves in one struct
    field_names = fieldnames(slow_waves);
%     field_names = fieldnames(slow_waves_markov);
%     slow_waves = slow_waves_markov;
    
    for nprobe = 1:length(slow_waves)
        probe_no = session_info(n).probe(nprobe).probe_hemisphere;
        slow_waves_all(probe_no).subject(session_count,:) = options.SUBJECT;
        slow_waves_all(probe_no).session(session_count,:) = options.SESSION;
        slow_waves_all(probe_no).session_day(session_count) = iDate;

        for iField =1:length(field_names)
            if session_count == 1
                if strcmp(field_names{iField},'ints')
                    slow_waves_all(probe_no).UP_ints = slow_waves(nprobe).ints.UP;
                    slow_waves_all(probe_no).DOWN_ints = slow_waves(nprobe).ints.DOWN;

                elseif  ismember(field_names{iField},{'power','timebin_edges','PSD_slope','frequency',...
                        'deltaspikecorr','gammaspikecorr','deltagammacorr','channel','shank','depth',...
                        'xcoord','gamma_t','viterbi_states','p','DOWN_peaktimes','DOWN_peaks_zscore','shank_id','probe_hemisphere',...
                        'amp_corr_DU','amp_corr_UD','amp_corr_UP','amp_corr_DOWN','mean_phase_UP','mean_phase_DOWN',...
                        'xcorr_lag_UD','xcorr_lag_DU','xcorr_r_UD','xcorr_r_DU','plv_UD','plv_DU','pd_UD','pd_DU'})
                    slow_waves_all(probe_no).(field_names{iField}){session_count} = slow_waves(nprobe).(field_names{iField});
                else
                    slow_waves_all(probe_no).(field_names{iField}) = slow_waves(nprobe).(field_names{iField});
                end
            else
                if strcmp(field_names{iField},'ints')
                    A = slow_waves_all(probe_no).UP_ints;
                    B = slow_waves(nprobe).ints.UP;
                    try
                        slow_waves_all(probe_no).UP_ints = [A;B];
                    catch
                        slow_waves_all(probe_no).UP_ints = [A B];
                    end


                    A = slow_waves_all(probe_no).DOWN_ints;
                    B = slow_waves(nprobe).ints.DOWN;
                    try
                        slow_waves_all(probe_no).DOWN_ints = [A;B];
                    catch
                        slow_waves_all(probe_no).DOWN_ints = [A B];
                    end
                elseif  ismember(field_names{iField},{'power','timebin_edges','PSD_slope','frequency',...
                        'deltaspikecorr','gammaspikecorr','deltagammacorr','channel','shank','depth',...
                        'xcoord','gamma_t','viterbi_states','p','DOWN_peaktimes','DOWN_peaks_zscore','shank_id','probe_hemisphere',...
                        'amp_corr_DU','amp_corr_UD','amp_corr_UP','amp_corr_DOWN','mean_phase_UP','mean_phase_DOWN',...
                        'xcorr_lag_UD','xcorr_lag_DU','xcorr_r_UD','xcorr_r_DU','plv_UD','plv_DU','pd_UD','pd_DU'})

                    slow_waves_all(probe_no).(field_names{iField}){session_count} = slow_waves(nprobe).(field_names{iField});
                else
                    A = slow_waves_all(probe_no).(field_names{iField});
                    B = slow_waves(nprobe).(field_names{iField});


                    try
                        slow_waves_all(probe_no).(field_names{iField}) = [A;B];
                    catch
                        slow_waves_all(probe_no).(field_names{iField}) = [A B];
                    end
                end
            end
        end
        % slow_waves_all(probe_no).UP_PSD_slope
        % session_count_events = repmat(session_count,[size(slow_waves(nprobe).ints.DOWN,1),1]);
        session_count_events = repmat(session_count,[size(slow_waves(nprobe).DOWN_ints,1),1]);
        if session_count == 1
            slow_waves_all(probe_no).DOWN_session_count = session_count_events;
        else
            slow_waves_all(probe_no).DOWN_session_count = [slow_waves_all(probe_no).DOWN_session_count;session_count_events];
        end

        % session_count_events = repmat(session_count,[size(slow_waves(nprobe).ints.UP,1),1]);
        session_count_events = repmat(session_count,[size(slow_waves(nprobe).UP_ints,1),1]);
        if session_count == 1
            slow_waves_all(probe_no).UP_session_count = session_count_events;
        else
            slow_waves_all(probe_no).UP_session_count = [slow_waves_all(probe_no).UP_session_count;session_count_events];
        end

%         session_count_events = repmat(session_count,[size(slow_waves(nprobe).UP_DOWN_index,1),1]);
%         if session_count == 1
%             slow_waves_all(probe_no).UP_DOWN_session_count = session_count_events;
%         else
%             slow_waves_all(probe_no).UP_DOWN_session_count = [slow_waves_all(probe_no).UP_DOWN_session_count;session_count_events];
%         end
% 
%         session_count_events = repmat(session_count,[size(slow_waves(nprobe).DOWN_UP_index,1),1]);
%         if session_count == 1
%             slow_waves_all(probe_no).DOWN_UP_session_count = session_count_events;
%         else
%             slow_waves_all(probe_no).DOWN_UP_session_count = [slow_waves_all(probe_no).DOWN_UP_session_count;session_count_events];
%         end


        if probe_no==1
            metric_param =[];
            % metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
            metric_param.region = @(x) contains(x,'V1_L');
            [selected_clusters,cluster_id] = select_clusters(MUA_combined,metric_param);
            slow_waves_all(probe_no).V1_MUA_spiketimes{session_count} = selected_clusters.spike_times;

            metric_param.region = @(x) contains(x,'HPC_L');
            [selected_clusters,cluster_id] = select_clusters(MUA_combined,metric_param);
            slow_waves_all(probe_no).HPC_MUA_spiketimes{session_count} = selected_clusters.spike_times;

            metric_param.region = @(x) contains(x,'V1_L');
            metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            slow_waves_all(probe_no).V1_spike_times{session_count} = selected_clusters.spike_times;
            slow_waves_all(probe_no).V1_spike_id{session_count} = selected_clusters.spike_times;

            metric_param.region = @(x) contains(x,'HPC_L');
            metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            slow_waves_all(probe_no).HPC_spike_times{session_count} = selected_clusters.spike_times;
            slow_waves_all(probe_no).HPC_spike_id{session_count} = selected_clusters.spike_times;
        elseif probe_no==2
            metric_param =[];
            % metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
            metric_param.region = @(x) contains(x,'V1_R');
            [selected_clusters,cluster_id] = select_clusters(MUA_combined,metric_param);
            slow_waves_all(probe_no).V1_MUA_spiketimes{session_count} = selected_clusters.spike_times;

            metric_param.region = @(x) contains(x,'HPC_R');
            [selected_clusters,cluster_id] = select_clusters(MUA_combined,metric_param);
            slow_waves_all(probe_no).HPC_MUA_spiketimes{session_count} = selected_clusters.spike_times;

            metric_param.region = @(x) contains(x,'V1_L');
            metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            slow_waves_all(probe_no).V1_spike_times{session_count} = selected_clusters.spike_times;
            slow_waves_all(probe_no).V1_spike_id{session_count} = selected_clusters.spike_times;

            metric_param.region = @(x) contains(x,'HPC_L');
            metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            slow_waves_all(probe_no).HPC_spike_times{session_count} = selected_clusters.spike_times;
            slow_waves_all(probe_no).HPC_spike_id{session_count} = selected_clusters.spike_times;
        end
    end
end

% 
% for nsession = 1:max(slow_waves_all(1).DOWN_session_count)
%     UP_index=[];
%     DOWN_index=[];
%     % get events index this session
%     for nprobe = 1:length(slow_waves_all)
%         UP_index{nprobe} = find(slow_waves_all(nprobe).UP_session_count == nsession);
%         DOWN_index{nprobe} = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
% 
%         UP_DOWN_transition = slow_waves_all(nprobe).UP_DOWN_index;
%         DOWN_UP_transition = slow_waves_all(nprobe).DOWN_ints(DOWN_index{nprobe},1);
% 
%         Exclude_index = find(~ismember(DOWN_UP_transition,UP_DOWN_transition));% Only DOWN followed by UP is included for analysis.
% 
%         % remove DOWN events without UP
%         slow_waves_all(nprobe).DOWN_session_count(DOWN_index{nprobe}(Exclude_index))=[];
%         slow_waves_all(nprobe).DOWN_ints(DOWN_index{nprobe}(Exclude_index),:)=[];
% %         slow_waves_all(nprobe).SWpeakmag(DOWN_index{nprobe}(Exclude_index))=[];
%         slow_waves_all(nprobe).timestamps(DOWN_index{nprobe}(Exclude_index))=[];
%         slow_waves_all(nprobe).DOWN_PSD_slope(DOWN_index{nprobe}(Exclude_index))=[];
%         slow_waves_all(nprobe).DOWN_peaks_shank{nsession}(:,DOWN_index{nprobe}(Exclude_index))=[];
%         slow_waves_all(nprobe).DOWN_peaks_latency(DOWN_index{nprobe}(Exclude_index))=[];
%         slow_waves_all(nprobe).DOWN_traveling(DOWN_index{nprobe}(Exclude_index))=[];
% 
%         % DOWN_index{nprobe} = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
%     end
% end

% ripples_all = rmfield(ripples_all, 'detectorinfo');
% spindles_all = rmfield(spindles_all, 'detectorinfo');
if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end

if contains(Stimulus_type,'Sleep') & ~contains(Stimulus_type,'PRE')
    save(fullfile(analysis_folder,'slow_waves_all_POST.mat'),'slow_waves_all','-v7.3')
    save(fullfile(analysis_folder,'ripples_all_POST.mat'),'ripples_all','-v7.3')
    save(fullfile(analysis_folder,'spindles_all_POST.mat'),'spindles_all','-v7.3')
    save(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'),'behavioural_state_merged_all','-v7.3')
elseif contains(Stimulus_type,'PRE')
    save(fullfile(analysis_folder,'slow_waves_all_PRE.mat'),'slow_waves_all','-v7.3')
    save(fullfile(analysis_folder,'ripples_all_PRE.mat'),'ripples_all','-v7.3')
    save(fullfile(analysis_folder,'spindles_all_PRE.mat'),'spindles_all','-v7.3')
    save(fullfile(analysis_folder,'behavioural_state_merged_all_PRE.mat'),'behavioural_state_merged_all','-v7.3')
else
    save(fullfile(analysis_folder,'slow_waves_all.mat'),'slow_waves_all','-v7.3')
    save(fullfile(analysis_folder,'ripples_all.mat'),'ripples_all','-v7.3')
    save(fullfile(analysis_folder,'spindles_all.mat'),'spindles_all','-v7.3')
    save(fullfile(analysis_folder,'behavioural_state_merged_all.mat'),'behavioural_state_merged_all','-v7.3')
end

%% Add on PLS KDE regression

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'SleepChronic';

session_count = 0;

bayesian_reactivation_all = struct();

KDE_reactivation_V1_all= struct();
KDE_reactivation_V1_UP_all= struct();
KDE_reactivation_V1_DOWN_all= struct();

KDE_reactivation_UP_all = struct();
KDE_reactivation_DOWN_all = struct();
KDE_reactivation_all = struct();

for nsession =1:length(experiment_info)

    tic
    disp(sprintf('session %i',nsession))
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT},option);
    iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    if length(stimulus_name)>1
        if contains(Stimulus_type,'PRE')
            disp('Same stimuli multiple recordings. Will take _2')
            n = find(contains(stimulus_name,'_2'));
        else
            session_info = session_info(~contains(stimulus_name,'PRE'));
            stimulus_name = stimulus_name(~contains(stimulus_name,'PRE'));
            if length(stimulus_name)>1
                disp('Same stimuli multiple recordings. Will take _2')
                n = find(contains(stimulus_name,'_2'));
            else
                n =1;
            end
        end
    else
        n = 1;
    end

    session_count = session_count + 1;
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

    if contains(stimulus_name{n},'Sleep')
        load(fullfile(options.ANALYSIS_DATAPATH,'PLS_KDE_reactivation_V1.mat'),'KDE_reactivation_V1');
        load(fullfile(options.ANALYSIS_DATAPATH,'KDE_reactivation_V1_DOWN.mat'),'KDE_reactivation_V1_DOWN');
        load(fullfile(options.ANALYSIS_DATAPATH,'KDE_reactivation_V1_UP.mat'),'KDE_reactivation_V1_UP');
    end

    for nprobe = 1:length(KDE_reactivation_V1)
        bias_distribution_V1 = [KDE_reactivation_V1(nprobe).event_bias(:)' KDE_reactivation_V1_UP(nprobe).event_bias(:)' KDE_reactivation_V1_DOWN(nprobe).event_bias(:)'];

        %%%%% Ripples
        real_bias = KDE_reactivation_V1(nprobe).event_bias(:)';
        shuffled = KDE_reactivation_V1(nprobe).event_T1_probability_shuffled ./ ...
            (KDE_reactivation_V1(nprobe).event_T1_probability_shuffled + KDE_reactivation_V1(nprobe).event_T2_probability_shuffled);

        nbin = size(shuffled, 2);
        percentile = nan(1, nbin);
        zscored_bias = nan(1, nbin);
        for i = 1:nbin
            valid_shuffled = shuffled(:, i);
            valid_shuffled = valid_shuffled(~isnan(valid_shuffled));
            if ~isempty(valid_shuffled) && ~isnan(real_bias(i))
                percentile(i) = sum(valid_shuffled < real_bias(i)) / length(valid_shuffled) * 100;
                zscored_bias(i) = (real_bias(i) - mean(valid_shuffled)) / std(valid_shuffled);
            end
        end

        KDE_reactivation_V1_all(nprobe).event_bins{nsession} = KDE_reactivation_V1(nprobe).event_bins;
        KDE_reactivation_V1_all(nprobe).event_id{nsession} = KDE_reactivation_V1(nprobe).event_id;
        KDE_reactivation_V1_all(nprobe).bias{nsession} = real_bias;
        KDE_reactivation_V1_all(nprobe).zscored_bias{nsession} = (real_bias - mean(bias_distribution_V1,'omitnan'))./std(bias_distribution_V1,'omitnan');
        KDE_reactivation_V1_all(nprobe).zscored_bias_shuffled{nsession} = zscored_bias;

        %%%%% UP
        real_bias_up = KDE_reactivation_V1_UP(nprobe).event_bias(:)';
        shuffled_up = KDE_reactivation_V1_UP(nprobe).event_T1_probability_shuffled ./ ...
            (KDE_reactivation_V1_UP(nprobe).event_T1_probability_shuffled + KDE_reactivation_V1_UP(nprobe).event_T2_probability_shuffled);

        nbin_up = size(shuffled_up, 2);
        zscored_bias_up = nan(1, nbin_up);
        for i = 1:nbin_up
            valid_shuffled = shuffled_up(:, i);
            valid_shuffled = valid_shuffled(~isnan(valid_shuffled));
            if ~isempty(valid_shuffled) && ~isnan(real_bias_up(i))
                zscored_bias_up(i) = (real_bias_up(i) - mean(valid_shuffled)) / std(valid_shuffled);
            end
        end

        KDE_reactivation_V1_UP_all(nprobe).event_bins{nsession} = KDE_reactivation_V1_UP(nprobe).event_bins;
        KDE_reactivation_V1_UP_all(nprobe).event_id{nsession} = KDE_reactivation_V1_UP(nprobe).event_id;
        KDE_reactivation_V1_UP_all(nprobe).bias{nsession} = real_bias_up;
        KDE_reactivation_V1_UP_all(nprobe).zscored_bias{nsession} = (real_bias_up - mean(bias_distribution_V1,'omitnan'))./std(bias_distribution_V1,'omitnan');
        KDE_reactivation_V1_UP_all(nprobe).zscored_bias_shuffled{nsession} = zscored_bias_up;

        %%%%% DOWN
        real_bias_down = KDE_reactivation_V1_DOWN(nprobe).event_bias(:)';
        shuffled_down = KDE_reactivation_V1_DOWN(nprobe).event_T1_probability_shuffled ./ ...
            (KDE_reactivation_V1_DOWN(nprobe).event_T1_probability_shuffled + KDE_reactivation_V1_DOWN(nprobe).event_T2_probability_shuffled);

        nbin_down = size(shuffled_down, 2);
        zscored_bias_down = nan(1, nbin_down);
        for i = 1:nbin_down
            valid_shuffled = shuffled_down(:, i);
            valid_shuffled = valid_shuffled(~isnan(valid_shuffled));
            if ~isempty(valid_shuffled) && ~isnan(real_bias_down(i))
                zscored_bias_down(i) = (real_bias_down(i) - mean(valid_shuffled)) / std(valid_shuffled);
            end
        end

        KDE_reactivation_V1_DOWN_all(nprobe).event_bins{nsession} = KDE_reactivation_V1_DOWN(nprobe).event_bins;
        KDE_reactivation_V1_DOWN_all(nprobe).event_id{nsession} = KDE_reactivation_V1_DOWN(nprobe).event_id;
        KDE_reactivation_V1_DOWN_all(nprobe).bias{nsession} = real_bias_down;
        KDE_reactivation_V1_DOWN_all(nprobe).zscored_bias{nsession} = (real_bias_down - mean(bias_distribution_V1,'omitnan'))./std(bias_distribution_V1,'omitnan');
        KDE_reactivation_V1_DOWN_all(nprobe).zscored_bias_shuffled{nsession} = zscored_bias_down;
    end
    clear KDE_reactivation_V1 KDE_reactivation_V1_DOWN KDE_reactivation_V1_UP

    load(fullfile(options.ANALYSIS_DATAPATH,'PLS_KDE_reactivation.mat'),'KDE_reactivation');
    load(fullfile(options.ANALYSIS_DATAPATH,'KDE_reactivation_DOWN.mat'),'KDE_reactivation_DOWN');
    load(fullfile(options.ANALYSIS_DATAPATH,'KDE_reactivation_UP.mat'),'KDE_reactivation_UP');

    for nprobe = 1:length(KDE_reactivation)
        bias_distribution = [KDE_reactivation(nprobe).event_bias(:)' KDE_reactivation_UP(nprobe).event_bias(:)' KDE_reactivation_DOWN(nprobe).event_bias(:)'];

        %%%%% Ripples
        real_bias = KDE_reactivation(nprobe).event_bias(:)';
        shuffled = KDE_reactivation(nprobe).event_T1_probability_shuffled ./ ...
            (KDE_reactivation(nprobe).event_T1_probability_shuffled + KDE_reactivation(nprobe).event_T2_probability_shuffled);

        nbin = size(shuffled, 2);
        percentile = nan(1, nbin);
        zscored_bias = nan(1, nbin);
        for i = 1:nbin
            valid_shuffled = shuffled(:, i);
            valid_shuffled = valid_shuffled(~isnan(valid_shuffled));
            if ~isempty(valid_shuffled) && ~isnan(real_bias(i))
                percentile(i) = sum(valid_shuffled < real_bias(i)) / length(valid_shuffled) * 100;
                zscored_bias(i) = (real_bias(i) - mean(valid_shuffled)) / std(valid_shuffled);
            end
        end

        KDE_reactivation_all(nprobe).event_bins{nsession} = KDE_reactivation(nprobe).event_bins;
        KDE_reactivation_all(nprobe).event_id{nsession} = KDE_reactivation(nprobe).event_id;
        KDE_reactivation_all(nprobe).bias{nsession} = real_bias;
        KDE_reactivation_all(nprobe).zscored_bias_shuffled{nsession} = zscored_bias;
        KDE_reactivation_all(nprobe).zscored_bias{nsession} = (real_bias - mean(bias_distribution,'omitnan'))./std(bias_distribution,'omitnan');

        %%%%% UP
        real_bias_up = KDE_reactivation_UP(nprobe).event_bias(:)';
        shuffled_up = KDE_reactivation_UP(nprobe).event_T1_probability_shuffled ./ ...
            (KDE_reactivation_UP(nprobe).event_T1_probability_shuffled + KDE_reactivation_UP(nprobe).event_T2_probability_shuffled);

        nbin_up = size(shuffled_up, 2);
        zscored_bias_up = nan(1, nbin_up);
        for i = 1:nbin_up
            valid_shuffled = shuffled_up(:, i);
            valid_shuffled = valid_shuffled(~isnan(valid_shuffled));
            if ~isempty(valid_shuffled) && ~isnan(real_bias_up(i))
                zscored_bias_up(i) = (real_bias_up(i) - mean(valid_shuffled)) / std(valid_shuffled);
            end
        end

        KDE_reactivation_UP_all(nprobe).event_bins{nsession} = KDE_reactivation_UP(nprobe).event_bins;
        KDE_reactivation_UP_all(nprobe).event_id{nsession} = KDE_reactivation_UP(nprobe).event_id;
        KDE_reactivation_UP_all(nprobe).bias{nsession} = real_bias_up;
        KDE_reactivation_UP_all(nprobe).zscored_bias{nsession} = (real_bias_up - mean(bias_distribution,'omitnan'))./std(bias_distribution,'omitnan');
        KDE_reactivation_UP_all(nprobe).zscored_bias_shuffled{nsession} = zscored_bias_up;

        %%%%% DOWN
        real_bias_down = KDE_reactivation_DOWN(nprobe).event_bias(:)';
        shuffled_down = KDE_reactivation_DOWN(nprobe).event_T1_probability_shuffled ./ ...
            (KDE_reactivation_DOWN(nprobe).event_T1_probability_shuffled + KDE_reactivation_DOWN(nprobe).event_T2_probability_shuffled);

        nbin_down = size(shuffled_down, 2);
        zscored_bias_down = nan(1, nbin_down);
        for i = 1:nbin_down
            valid_shuffled = shuffled_down(:, i);
            valid_shuffled = valid_shuffled(~isnan(valid_shuffled));
            if ~isempty(valid_shuffled) && ~isnan(real_bias_down(i))
                zscored_bias_down(i) = (real_bias_down(i) - mean(valid_shuffled)) / std(valid_shuffled);
            end
        end

        KDE_reactivation_DOWN_all(nprobe).event_bins{nsession} = KDE_reactivation_DOWN(nprobe).event_bins;
        KDE_reactivation_DOWN_all(nprobe).event_id{nsession} = KDE_reactivation_DOWN(nprobe).event_id;
        KDE_reactivation_DOWN_all(nprobe).bias{nsession} = real_bias_down;
        KDE_reactivation_DOWN_all(nprobe).zscored_bias{nsession} = (real_bias_down - mean(bias_distribution,'omitnan'))./std(bias_distribution,'omitnan');
        KDE_reactivation_DOWN_all(nprobe).zscored_bias_shuffled{nsession} = zscored_bias_down;
    end
    clear KDE_reactivation KDE_reactivation_DOWN KDE_reactivation_UP
    toc
end

if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end

save(fullfile(analysis_folder,'KDE_reactivation_DOWN_all_POST.mat'),'KDE_reactivation_DOWN_all')
save(fullfile(analysis_folder,'KDE_reactivation_UP_all_POST.mat'),'KDE_reactivation_UP_all')
save(fullfile(analysis_folder,'KDE_reactivation_all_POST.mat'),'KDE_reactivation_all')
save(fullfile(analysis_folder,'KDE_reactivation_V1_DOWN_all_POST.mat'),'KDE_reactivation_V1_DOWN_all')
save(fullfile(analysis_folder,'KDE_reactivation_V1_UP_all_POST.mat'),'KDE_reactivation_V1_UP_all')
save(fullfile(analysis_folder,'KDE_reactivation_V1_all_POST.mat'),'KDE_reactivation_V1_all')


%% Add on Log odds bayesian bias

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'SleepChronic';

session_count = 0;

bayesian_reactivation_all= struct();
bayesian_reactivation_V1_all= struct();


for nsession =1:length(experiment_info)

    tic
    disp(sprintf('session %i',nsession))
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT},option);
    iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    if length(stimulus_name)>1
        if contains(Stimulus_type,'PRE')
            disp('Same stimuli multiple recordings. Will take _2')
            n = find(contains(stimulus_name,'_2'));
        else
            session_info = session_info(~contains(stimulus_name,'PRE'));
            stimulus_name = stimulus_name(~contains(stimulus_name,'PRE'));
            if length(stimulus_name)>1
                disp('Same stimuli multiple recordings. Will take _2')
                n = find(contains(stimulus_name,'_2'));
            else
                n =1;
            end
        end
    else
        n = 1;
    end

    session_count = session_count + 1;
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

    if contains(stimulus_name{n},'Sleep')
        load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'),'decoded_ripple_events');
        load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events_shuffled.mat'),'decoded_ripple_events_shuffled');
        load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events_V1.mat'),'decoded_ripple_events_V1');
        load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events_V1_shuffled.mat'),'decoded_ripple_events_V1_shuffled');
    end

    %%%% Re-enter PV threshold used for decoding
    load(fullfile(options.ANALYSIS_DATAPATH,'..','Masa2tracks','population_vector_corr_HPC_combined_RUN1.mat'),'PPvector');
    PPvector_HPC=PPvector;
    load(fullfile(options.ANALYSIS_DATAPATH,'..','Masa2tracks','population_vector_corr_V1_combined_RUN1.mat'),'PPvector');
    PPvector_V1=PPvector;

    across_tracks_pairs = [3 4 7 8 9 10 13 14];
    within_tracks_pairs = setdiff(1:16,across_tracks_pairs);
    x_bin_width =5;

    
    position_bins = 1:2:140;
    bad_bins = position_bins(sum((PPvector_V1.pval(:,across_tracks_pairs)<0.05).*(PPvector_V1.population_vector(:,across_tracks_pairs)>0),2) +...
        sum(PPvector_V1.pval(:,within_tracks_pairs)>0.05,2)>0);
    good_bins = find(histcounts(bad_bins,0:x_bin_width:140)<1);

    PV_thresholds = [];
    PV_threshold = [];
    if length(good_bins)<5
        PV_thresholds = 0.3:0.05:0.5;
        for i = 1:length(PV_thresholds)
            PV_threshold = PV_thresholds(i);

            % Identify bad bins
            bad_bins = position_bins( ...
                sum((PPvector_V1.population_vector(:,across_tracks_pairs) > PV_threshold) .* ...
                (PPvector_V1.population_vector(:,across_tracks_pairs) > 0), 2) + ...
                sum(PPvector_V1.pval(:,within_tracks_pairs) > 0.05, 2) > 0);

            % Histogram to count which bins are "bad"
            bin_hist = histcounts(bad_bins, 0:x_bin_width:140);

            % Good bins: those not present in bad_bins
            good_bins = find(bin_hist < 1);

            % Exit condition
            if length(good_bins) >= 5
                break;
            end
        end
    end


    for nprobe = 1:length(decoded_ripple_events_V1)

        decoded_ripple_events_V1(nprobe).good_bins = good_bins;
        decoded_ripple_events_V1(nprobe).PV_threshold =   PV_threshold; % empty if based on cell id distribution

        decoded_matrix = [decoded_ripple_events_V1(nprobe).track(1).replay_events(:).replay; decoded_ripple_events_V1(nprobe).track(2).replay_events(:).replay];
        decoded_matrix = reshape(decoded_matrix,56,50,[]);

        %         summed_probability = [good_bins good_bins+28];
        summed_probability_distribution= [decoded_ripple_events_V1(nprobe).track(1).replay_events(:).summed_probability; decoded_ripple_events_V1(nprobe).track(2).replay_events(:).summed_probability];
        summed_probability = reshape(summed_probability_distribution,2,50,[]);

        bayesian_reactivation_V1_all(nprobe).summed_probability = summed_probability;
        bayesian_reactivation_V1_all(nprobe).decoded_matrix = decoded_matrix;

        bayesian_reactivation_V1_all(nprobe).bias = reshape(summed_probability_distribution(1,:)./sum(summed_probability_distribution,1),50,[]);
        bayesian_reactivation_V1_all(nprobe).z_bias = reshape(zscore(summed_probability_distribution(1,:)./sum(summed_probability_distribution,1)),50,[]);

        bayesian_reactivation_V1_all(nprobe).summed_probability = summed_probability;

        for nevent = 1:length(decoded_ripple_events(nprobe).track(1).replay_events)
            data = log(summed_probability(1,:,nevent)./summed_probability(2,:,nevent));

            shuffled_data = log(decoded_ripple_events_shuffled(nprobe).track(1).replay_events(nevent).summed_probability...
                ./decoded_ripple_events_shuffled(nprobe).track(2).replay_events(nevent).summed_probability);

            bayesian_reactivation_V1_all(nprobe).z_log_odds(:,nevent) = (data-mean(shuffled_data))./std(shuffled_data);


            for i = 1:length(data)
                % Get the column distribution
                dist = shuffled_data(:, i);

                % Percentile is the proportion of values less than or equal to data(i)
                bayesian_reactivation_V1_all(nprobe).log_odds_percentile(i,nevent) = sum(dist <= data(i)) / length(dist) * 100;
            end
        end
    end

    save(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events_V1.mat'),'decoded_ripple_events_V1');


    %%%%%%%% HPC
    position_bins = 1:2:140;
    bad_bins = position_bins(sum((PPvector_HPC.pval(:,across_tracks_pairs)<0.05).*(PPvector_HPC.population_vector(:,across_tracks_pairs)>0),2) +...
        sum(PPvector_HPC.pval(:,within_tracks_pairs)>0.05,2)>0);
    good_bins = find(histcounts(bad_bins,0:x_bin_width:140)<1);

    PV_thresholds = [];
    PV_threshold = [];
    if length(good_bins)<5
        PV_thresholds = 0.3:0.05:0.5;
        for i = 1:length(PV_thresholds)
            PV_threshold = PV_thresholds(i);

            % Identify bad bins
            bad_bins = position_bins( ...
                sum((PPvector_HPC.population_vector(:,across_tracks_pairs) > PV_threshold) .* ...
                (PPvector_HPC.population_vector(:,across_tracks_pairs) > 0), 2) + ...
                sum(PPvector_HPC.pval(:,within_tracks_pairs) > 0.05, 2) > 0);

            % Histogram to count which bins are "bad"
            bin_hist = histcounts(bad_bins, 0:x_bin_width:140);

            % Good bins: those not present in bad_bins
            good_bins = find(bin_hist < 1);

            % Exit condition
            if length(good_bins) >= 5
                break;
            end
        end
    end
    save(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'),'decoded_ripple_events');


end

if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end

save(fullfile(analysis_folder,'KDE_reactivation_DOWN_all_POST.mat'),'KDE_reactivation_DOWN_all')
save(fullfile(analysis_folder,'KDE_reactivation_UP_all_POST.mat'),'KDE_reactivation_UP_all')


%% Analyse and plot peri-ripple, peri-spindle and peri-UP activity



%% UP DOWN state and ripple and spindle analysis
clear all
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
% load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
load(fullfile(analysis_folder,'ripples_all_POST.mat'))
load(fullfile(analysis_folder,'spindles_all_POST.mat'))
load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))
% 
% for nsession = 1:max(slow_waves_all(1).DOWN_session_count)
%     UP_index=[];
%     DOWN_index=[];
%     ripples_index=[];
% 
%     % get events index this session
%     for nprobe = 1:length(slow_waves_all)
%         UP_index{nprobe} = find(slow_waves_all(nprobe).UP_session_count == nsession);
%         DOWN_index{nprobe} = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
%         ripples_index{nprobe} = find(ripples_all(nprobe).session_count == nsession& ripples_all(nprobe).SWS_index == 1);
% 
%         UP_DOWN_transition = slow_waves_all(nprobe).UP_ints(UP_index{nprobe},2);
%         DOWN_UP_transition = slow_waves_all(nprobe).DOWN_ints(DOWN_index{nprobe},1);
% 
%         DOWN_index{nprobe} = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
%     end
% end
% 
% 
% 
% for nsession = 1:max(slow_waves_all(1).DOWN_session_count)
%     for nprobe = 1:2
%         shank_id(nprobe,nsession)=ceil(slow_waves_all(nprobe).xcoord{nsession}(slow_waves_all(nprobe).channel{nsession}==slow_waves_all(nprobe).best_channel{nsession})/250);
%         
%     end
% end
% for nprobe = 1:2
%     figure
%     for nsession = 1:all_sessions
% %         UP_DOWN_index = find(slow_waves_all(nprobe).UP_session_count == nsession); % Find UP -> DOWN
%         UP_index = find(slow_waves_all(nprobe).UP_session_count == nsession); % Find UP this session
% %         UP_index = UP_index(slow_waves_all(nprobe).UP_DOWN_index(UP_DOWN_index,1)); % Find UP followed by a DOWN 
%         UP_duration = slow_waves_all(nprobe).UP_ints(UP_index,2)-slow_waves_all(nprobe).UP_ints(UP_index,1);
%         UP_ints = slow_waves_all(nprobe).UP_ints(UP_index(UP_duration<=2),:);
%         UP_duration = UP_ints(:,2)-UP_ints(:,1);
% 
% 
%         DOWN_index = find(slow_waves_all(nprobe).DOWN_session_count == nsession); % Find DOWN this session
%         DOWN_ints = slow_waves_all(nprobe).DOWN_ints(DOWN_index,:);
%         [C,ia,ib] = intersect(UP_ints(:,2),DOWN_ints(:,1));
%         DOWN_index = DOWN_index(ib);
% 
%         DOWN_ints = slow_waves_all(nprobe).DOWN_ints(DOWN_index,:);
%         DOWN_duration = DOWN_ints(:,2)-DOWN_ints(:,1);
% 
%         ripples_index = find(ripples_all(1).session_count == nsession& ripples_all(1).SWS_index == 1);
% 
% 
%         %     UP_ints = slow_waves_all(nprobe).UP_ints(UP_index(UP_duration<10),:);
%         % UP_duration<2
%         subplot(3,5,nsession)
% c
%         hold on
%         histogram(log10(DOWN_duration),-2:0.1:0.5,'EdgeColor','none')
%         histogram(log10(UP_duration),-2:0.1:0.5,'EdgeColor','none')
%     end
% end


%% Caclulating ripple-UP-DOWN temporal probability
time_windows = [-0.5 0.5];
time_bin = 0.01;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
probability=[];
all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;


colour_lines = [215,25,28;44,123,182]/256;

% sessions_to_process = [1 2 3 6 8 9 10 11 12 13];
sessions_to_process = 1:all_sessions;
%%%% P(ripples) during UP DOWN
probability = calculate_UP_DOWN_ripple_probability(slow_waves_all,ripples_all,sessions_to_process,'option','normalised','time_option','peaktimes')
probability_normalised = probability;
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised.mat'),'probability_normalised');
probability = calculate_UP_DOWN_ripple_probability(slow_waves_all,ripples_all,sessions_to_process,'option','absolute','time_option','peaktimes','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability.mat'),'probability');


probability = calculate_UP_DOWN_ripple_probability(slow_waves_all,ripples_all,sessions_to_process,'option','normalised','time_option','whole')
probability_normalised = probability;
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'),'probability_normalised');
probability = calculate_UP_DOWN_ripple_probability(slow_waves_all,ripples_all,sessions_to_process,'option','absolute','time_option','whole','time_windows',[-1 1])
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'),'probability');



%%%% P(UP) and P(DOWN) during ripples
probability = calculate_ripple_UP_DOWN_probability(slow_waves_all,ripples_all,sessions_to_process,'option','absolute','time_option','onset','time_windows',[-1 1])
save(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_SO_probability.mat'),'probability');
probability = calculate_ripple_UP_DOWN_probability(slow_waves_all,ripples_all,sessions_to_process,'option','absolute','time_option','whole','time_windows',[-1 1])
save(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_SO_probability_whole.mat'),'probability');

%%%% P(spindles) during UP DOWN
probability = calculate_UP_DOWN_spindle_probability(slow_waves_all,spindles_all,sessions_to_process,'option','normalised','time_option','onset')
probability_normalised = probability;
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_probability_normalised.mat'),'probability_normalised');
probability = calculate_UP_DOWN_spindle_probability(slow_waves_all,spindles_all,sessions_to_process,'option','absolute','time_option','onset','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_probability.mat'),'probability');


probability = calculate_UP_DOWN_spindle_probability(slow_waves_all,spindles_all,sessions_to_process,'option','normalised','time_option','whole')
probability_normalised = probability;
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_probability_normalised_whole.mat'),'probability_normalised');
probability = calculate_UP_DOWN_spindle_probability(slow_waves_all,spindles_all,sessions_to_process,'option','absolute','time_option','whole','time_windows',[-1 1])
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_probability_whole.mat'),'probability');


%%%% P(spindles) during ripples
probability = calculate_ripple_spindle_probability(slow_waves_all, ripples_all, spindles_all, sessions_to_process,'option','absolute', 'time_option','whole','time_windows',[-1 1])
save(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_spindles_probability_whole.mat'),'probability');

probability = calculate_ripple_spindle_probability(slow_waves_all, ripples_all, spindles_all, sessions_to_process,'option','absolute', 'time_option','onset','time_windows',[-1 1])
save(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_spindles_probability.mat'),'probability');

probability = calculate_ripple_spindle_probability(slow_waves_all, ripples_all, spindles_all, sessions_to_process,'shuffle_option','baseline','option','absolute', 'time_option','whole','time_windows',[-1 1])
save(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_spindles_probability_whole_baseline.mat'),'probability');

probability = calculate_ripple_spindle_probability(slow_waves_all, ripples_all, spindles_all, sessions_to_process,'shuffle_option','baseline','option','absolute', 'time_option','onset','time_windows',[-1 1])
save(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_spindles_probability_baseline.mat'),'probability');


%%%% P(ripples) during spindles
probability = calculate_spindle_ripple_probability(slow_waves_all, ripples_all, spindles_all, sessions_to_process,'option','absolute', 'time_option','whole','time_windows',[-1 1])
save(fullfile(analysis_folder,'V1-HPC sleep interaction','spindles_ripples_probability_whole.mat'),'probability');

probability = calculate_spindle_ripple_probability(slow_waves_all, ripples_all, spindles_all, sessions_to_process,'option','absolute', 'time_option','onset','time_windows',[-1 1])
save(fullfile(analysis_folder,'V1-HPC sleep interaction','spindles_ripples_probability.mat'),'probability');

probability = calculate_spindle_ripple_probability(slow_waves_all, ripples_all, spindles_all, sessions_to_process,'shuffle_option','baseline','option','absolute', 'time_option','whole','time_windows',[-1 1])
save(fullfile(analysis_folder,'V1-HPC sleep interaction','spindles_ripples_probability_whole_baseline.mat'),'probability');

probability = calculate_spindle_ripple_probability(slow_waves_all, ripples_all, spindles_all, sessions_to_process,'shuffle_option','baseline','option','absolute', 'time_option','onset','time_windows',[-1 1])
save(fullfile(analysis_folder,'V1-HPC sleep interaction','spindles_ripples_probability_baseline.mat'),'probability');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% baseline shuffle (-3s before each UP or DOWN events)
%%%% P(ripples) during UP DOWN
probability = calculate_UP_DOWN_ripple_probability(slow_waves_all,ripples_all,sessions_to_process,'option','normalised','time_option','peaktimes','shuffle_option','baseline')
probability_normalised = probability;
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_baseline.mat'),'probability_normalised');
probability = calculate_UP_DOWN_ripple_probability(slow_waves_all,ripples_all,sessions_to_process,'option','absolute','time_option','peaktimes','time_windows',[-1 1],'shuffle_option','baseline');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_baseline.mat'),'probability');


probability = calculate_UP_DOWN_ripple_probability(slow_waves_all,ripples_all,sessions_to_process,'option','normalised','time_option','whole','shuffle_option','baseline')
probability_normalised = probability;
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole_baseline.mat'),'probability_normalised');
probability = calculate_UP_DOWN_ripple_probability(slow_waves_all,ripples_all,sessions_to_process,'option','absolute','time_option','whole','time_windows',[-1 1],'shuffle_option','baseline')
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline.mat'),'probability');

%%%%%% baseline shuffle (-3s before each ripples)
%%%% P(UP) and P(DOWN) during ripples
probability = calculate_ripple_UP_DOWN_probability(slow_waves_all,ripples_all,sessions_to_process,'option','absolute','time_option','onset','shuffle_option','baseline')
save(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_SO_probability_baseline.mat'),'probability');
probability = calculate_ripple_UP_DOWN_probability(slow_waves_all,ripples_all,sessions_to_process,'option','absolute','time_option','whole','shuffle_option','baseline')
save(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_SO_probability_whole_baseline.mat'),'probability');

%%%%%% baseline shuffle (-3s before each ripples)
%%%% P(spindles) during UP DOWN
probability = calculate_UP_DOWN_spindle_probability(slow_waves_all,spindles_all,sessions_to_process,'option','normalised','time_option','onset','shuffle_option','baseline')
probability_normalised = probability;
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_probability_normalised_baseline.mat'),'probability_normalised');
probability = calculate_UP_DOWN_spindle_probability(slow_waves_all,spindles_all,sessions_to_process,'option','absolute','time_option','onset','time_windows',[-1 1],'shuffle_option','baseline');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_probability_baseline.mat'),'probability');


probability = calculate_UP_DOWN_spindle_probability(slow_waves_all,spindles_all,sessions_to_process,'option','normalised','time_option','whole','shuffle_option','baseline')
probability_normalised = probability;
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_probability_normalised_whole_baseline.mat'),'probability_normalised');
probability = calculate_UP_DOWN_spindle_probability(slow_waves_all,spindles_all,sessions_to_process,'option','absolute','time_option','whole','time_windows',[-1 1],'shuffle_option','baseline')
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_probability_whole_baseline.mat'),'probability');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Probability of UP during DOWN and DOWN during UP
time_windows = [-1 1];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
probability=[];
all_sessions = max(slow_waves_all(1).DOWN_session_count);
sessions_to_process = 1:all_sessions;


probability = calculate_SO_SO_probability(slow_waves_all,sessions_to_process,'time_option','whole','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_probability_whole.mat'),'probability');
probability = calculate_SO_SO_probability(slow_waves_all,sessions_to_process,'time_option','whole','shuffle_option','baseline','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_probability_whole_baseline.mat'),'probability');

probability = calculate_SO_SO_probability(slow_waves_all,sessions_to_process,'time_option','onset','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_probability.mat'),'probability');
probability = calculate_SO_SO_probability(slow_waves_all,sessions_to_process,'time_option','onset','shuffle_option','baseline','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_probability_baseline.mat'),'probability');

probability = calculate_SO_SO_contralateral_probability(slow_waves_all,sessions_to_process,'time_option','whole','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_whole_probability.mat'),'probability');
probability = calculate_SO_SO_contralateral_probability(slow_waves_all,sessions_to_process,'time_option','whole','shuffle_option','baseline','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_whole_probability_baseline.mat'),'probability');

probability = calculate_SO_SO_contralateral_probability(slow_waves_all,sessions_to_process,'time_option','onset','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_probability.mat'),'probability');
probability = calculate_SO_SO_contralateral_probability(slow_waves_all,sessions_to_process,'time_option','onset','shuffle_option','baseline','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_probability_baseline.mat'),'probability');

probability = calculate_ripple_ripple_probability(slow_waves_all,ripples_all,sessions_to_process,'time_option','whole','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_ripples_probability_whole.mat'),'probability');
probability = calculate_ripple_ripple_probability(slow_waves_all,ripples_all,sessions_to_process,'time_option','whole','shuffle_option','baseline','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_ripples_probability_whole_baseline.mat'),'probability');

probability = calculate_ripple_ripple_probability(slow_waves_all,ripples_all,sessions_to_process,'time_option','peaktimes','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_ripples_probability.mat'),'probability');
probability = calculate_ripple_ripple_probability(slow_waves_all,ripples_all,sessions_to_process,'time_option','peaktimes','shuffle_option','baseline','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_ripples_probability_baseline.mat'),'probability');


probability = calculate_spindle_spindle_probability(slow_waves_all,spindles_all,sessions_to_process,'time_option','whole','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','spindles_spindles_probability_whole.mat'),'probability');
probability = calculate_spindle_spindle_probability(slow_waves_all,spindles_all,sessions_to_process,'time_option','whole','shuffle_option','baseline','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','spindles_spindles_probability_whole_baseline.mat'),'probability');

probability = calculate_spindle_spindle_probability(slow_waves_all,spindles_all,sessions_to_process,'time_option','onset','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','spindles_spindles_probability.mat'),'probability');
probability = calculate_spindle_spindle_probability(slow_waves_all,spindles_all,sessions_to_process,'time_option','onset','shuffle_option','baseline','time_windows',[-1 1]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','spindles_spindles_probability_baseline.mat'),'probability');


%% Extract key information
temp = [];
sessions_to_process = 1:max(slow_waves_all(1).UP_session_count);
for nprobe = 1:2
    temp{nprobe} = extract_UP_DOWN_ripples_info(slow_waves_all,ripples_all,behavioural_state_merged_all,sessions_to_process,'option','UD','nprobe',nprobe)
    % event_info(2) = extract_UP_DOWN_ripples_info(slow_waves_all,ripples_all,behavioural_state_merged_all,sessions_to_process,'option','DU','nprobe',1)
end
clear event_info
event_info(1) = temp{1};
event_info(2) = temp{2}(2);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');
%% Plotting and calculating MUA relative to UP DOWN and ripple
all_sessions = max(slow_waves_all(1).DOWN_session_count);
sessions_to_process = 1:all_sessions;


PSTH_MUA_1ms = calculate_UP_DOWN_ripple_PSTH_1ms...
    (slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,sessions_to_process,'option','MUA','time_option','absolute','shuffle_option','no');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','PSTH_MUA_1ms.mat'),'PSTH_MUA_1ms','-v7.3');



UP_DOWN_ripple_PSTH_MUA = calculate_UP_DOWN_ripple_PSTH...
    (slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,sessions_to_process,'option','MUA','time_option','absolute','shuffle_option','no');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'),'UP_DOWN_ripple_PSTH_MUA');
UP_DOWN_ripple_PSTH_MUA = calculate_UP_DOWN_ripple_PSTH...
    (slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,sessions_to_process,'option','MUA','time_option','absolute','shuffle_option','baseline');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'),'UP_DOWN_ripple_PSTH_MUA');


UP_DOWN_relative_PSTH_MUA = calculate_UP_DOWN_relative_PSTH...
    (slow_waves_all,ripples_all,behavioural_state_merged_all,sessions_to_process,'option','MUA','shuffle_option','no');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_relative_PSTH_MUA.mat'),'UP_DOWN_relative_PSTH_MUA');
UP_DOWN_relative_PSTH_MUA = calculate_UP_DOWN_relative_PSTH...
    (slow_waves_all,ripples_all,behavioural_state_merged_all,sessions_to_process,'option','MUA','shuffle_option','baseline');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_relative_PSTH_MUA_baseline.mat'),'UP_DOWN_relative_PSTH_MUA');

end
