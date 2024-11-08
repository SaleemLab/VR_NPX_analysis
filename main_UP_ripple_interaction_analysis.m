%% all events for sleep analysis (spiking but also event info)

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
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP');
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        clusters=clusters_ks3;
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
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
        %             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'));
        clusters=clusters_ks3;
    else
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
%         load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'));
        clusters=clusters_ks3;
    end
    ripples = rmfield(ripples, 'detectorinfo');
    spindles = rmfield(spindles, 'detectorinfo');
    slow_waves = rmfield(slow_waves, 'detectorinfo');
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
                ripples_all(probe_no).(field_names{iField}) = ripples(nprobe).(field_names{iField});
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
                spindles_all(probe_no).(field_names{iField}) = spindles(nprobe).(field_names{iField});
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

        session_count_events = repmat(session_count,size(spindles(nprobe).onset));
        if session_count == 1
            spindles_all(probe_no).session_count = session_count_events;
        else
            spindles_all(probe_no).session_count = [spindles_all(probe_no).session_count;session_count_events];
        end
    end

    % put all slow waves in one struct
    field_names = fieldnames(slow_waves);

    for nprobe = 1:length(slow_waves)
        probe_no = session_info(n).probe(nprobe).probe_hemisphere;
        slow_waves_all(probe_no).subject(session_count,:) = options.SUBJECT;
        slow_waves_all(probe_no).session(session_count,:) = options.SESSION;
        slow_waves_all(probe_no).session_day(session_count) = iDate;

        for iField =1:length(field_names)
            if session_count == 1
                if contains(field_names{iField},'ints')
                    slow_waves_all(probe_no).UP_ints = slow_waves(nprobe).ints.UP;
                    slow_waves_all(probe_no).DOWN_ints = slow_waves(nprobe).ints.DOWN;
                    slow_waves_all(probe_no).UP_PSD_slope = slow_waves(nprobe).ints.UP_PSD_slope; 
                elseif  ismember(field_names{iField},{'power','timebin_edges','PSD_slope','frequency',...
                        'deltaspikecorr','gammaspikecorr','deltagammacorr','channel','shank','depth','xcoord','best_channel'})
                    slow_waves_all(probe_no).(field_names{iField}){session_count} = slow_waves(nprobe).(field_names{iField});
                else
                    slow_waves_all(probe_no).(field_names{iField}) = slow_waves(nprobe).(field_names{iField});
                end
            else
                if contains(field_names{iField},'ints')
                    A = slow_waves_all(probe_no).UP_ints;
                    B = slow_waves(nprobe).ints.UP;
                    try
                        slow_waves_all(probe_no).UP_ints = [A;B];
                    catch
                        slow_waves_all(probe_no).UP_ints = [A B];
                    end

                    A = slow_waves_all(probe_no).UP_PSD_slope;
                    B = slow_waves(nprobe).ints.UP_PSD_slope;
                    try
                        slow_waves_all(probe_no).UP_PSD_slope = [A;B];
                    catch
                        slow_waves_all(probe_no).UP_PSD_slope = [A B];
                    end

                    A = slow_waves_all(probe_no).DOWN_ints;
                    B = slow_waves(nprobe).ints.DOWN;
                    try
                        slow_waves_all(probe_no).DOWN_ints = [A;B];
                    catch
                        slow_waves_all(probe_no).DOWN_ints = [A B];
                    end
                elseif  ismember(field_names{iField},{'power','timebin_edges','PSD_slope','frequency',...
                        'deltaspikecorr','gammaspikecorr','deltagammacorr','channel','shank','depth','xcoord','best_channel'})

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

        session_count_events = repmat(session_count,size(slow_waves(nprobe).timestamps));
        if session_count == 1
            slow_waves_all(probe_no).DOWN_session_count = session_count_events;
        else
            slow_waves_all(probe_no).DOWN_session_count = [slow_waves_all(probe_no).DOWN_session_count;session_count_events];
        end

        session_count_events = repmat(session_count,[size(slow_waves(nprobe).ints.UP,1) 1]);
        if session_count == 1
            slow_waves_all(probe_no).UP_session_count = session_count_events;
        else
            slow_waves_all(probe_no).UP_session_count = [slow_waves_all(probe_no).UP_session_count;session_count_events];
        end
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


for nsession = 1:max(slow_waves_all(1).DOWN_session_count)
    UP_index=[];
    DOWN_index=[];
    ripples_index=[];

    % get events index this session
    for nprobe = 1:length(slow_waves_all)
        UP_index{nprobe} = find(slow_waves_all(nprobe).UP_session_count == nsession);
        DOWN_index{nprobe} = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
        ripples_index{nprobe} = find(ripples_all(nprobe).session_count == nsession& ripples_all(nprobe).SWS_index == 1);

        UP_DOWN_transition = slow_waves_all(nprobe).UP_ints(UP_index{nprobe},2);
        DOWN_UP_transition = slow_waves_all(nprobe).DOWN_ints(DOWN_index{nprobe},1);

        Exclude_index = find(~ismember(DOWN_UP_transition,UP_DOWN_transition));% Only DOWN followed by UP is included for analysis.

        % remove DOWN events without UP
        slow_waves_all(nprobe).DOWN_session_count(DOWN_index{nprobe}(Exclude_index))=[];
        slow_waves_all(nprobe).DOWN_ints(DOWN_index{nprobe}(Exclude_index))=[];
        slow_waves_all(nprobe).SWpeakmag(DOWN_index{nprobe}(Exclude_index))=[];
        slow_waves_all(nprobe).timestamps(DOWN_index{nprobe}(Exclude_index))=[];
        slow_waves_all(nprobe).DOWN_PSD_slope(DOWN_index{nprobe}(Exclude_index))=[];
        slow_waves_all(nprobe).DOWN_peaks_shank(:,DOWN_index{nprobe}(Exclude_index))=[];
        slow_waves_all(nprobe).DOWN_peaks_latency(DOWN_index{nprobe}(Exclude_index))=[];
        slow_waves_all(nprobe).DOWN_traveling(DOWN_index{nprobe}(Exclude_index))=[];

        % DOWN_index{nprobe} = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
    end
end

% ripples_all = rmfield(ripples_all, 'detectorinfo');
% spindles_all = rmfield(spindles_all, 'detectorinfo');
if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end

if contains(Stimulus_type,'Sleep') & ~contains(Stimulus_type,'PRE')
    save(fullfile(analysis_folder,'slow_waves_all_POST.mat'),'slow_waves_all')
    save(fullfile(analysis_folder,'ripples_all_POST.mat'),'ripples_all')
    save(fullfile(analysis_folder,'spindles_all_POST.mat'),'spindles_all')
    save(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'),'behavioural_state_merged_all')
elseif contains(Stimulus_type,'PRE')
    save(fullfile(analysis_folder,'slow_waves_all_PRE.mat'),'slow_waves_all')
    save(fullfile(analysis_folder,'ripples_all_PRE.mat'),'ripples_all')
    save(fullfile(analysis_folder,'spindles_all_PRE.mat'),'spindles_all')
    save(fullfile(analysis_folder,'behavioural_state_merged_all_PRE.mat'),'behavioural_state_merged_all')
else
    save(fullfile(analysis_folder,'slow_waves_all.mat'),'slow_waves_all')
    save(fullfile(analysis_folder,'ripples_all.mat'),'ripples_all')
    save(fullfile(analysis_folder,'spindles_all.mat'),'spindles_all')
    save(fullfile(analysis_folder,'behavioural_state_merged_all.mat'),'behavioural_state_merged_all')
end

%% UP DOWN state and ripple and spindle analysis
clear all
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
load(fullfile(analysis_folder,'ripples_all_POST.mat'))
load(fullfile(analysis_folder,'spindles_all_POST.mat'))
load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))

for nsession = 1:max(slow_waves_all(1).DOWN_session_count)
    UP_index=[];
    DOWN_index=[];
    ripples_index=[];

    % get events index this session
    for nprobe = 1:length(slow_waves_all)
        UP_index{nprobe} = find(slow_waves_all(nprobe).UP_session_count == nsession);
        DOWN_index{nprobe} = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
        ripples_index{nprobe} = find(ripples_all(nprobe).session_count == nsession& ripples_all(nprobe).SWS_index == 1);

        UP_DOWN_transition = slow_waves_all(mprobe).UP_ints(UP_index{mprobe},2);
        DOWN_UP_transition = slow_waves_all(mprobe).DOWN_ints(DOWN_index{mprobe},1);

        DOWN_index{nprobe} = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
    end



    %%% Get spike counts
    last_spike = slow_waves_all(nprobe).V1_MUA_spiketimes{nsession}(end);
    tvec = 0:0.01:last_spike;
    % mobility_interp = interp1(selected_clusters.tvec{1},double(selected_clusters.mobility_thresholded{1}),tvec,'previous');
    tvec_edges = [tvec(1)-1/(1/mean(diff(tvec))*2) tvec+1/(1/mean(diff(tvec))*2)];
    if ~isempty(behavioural_state_merged_all.SWS{nsession})
        sleep_tvec = Restrict(tvec,behavioural_state_merged_all.SWS{nsession});
        sleep_index = ismember(tvec,sleep_tvec);
    end
    HPC_spike_counts=[];
    V1_spike_counts=[];
    w = gausswin(0.03*1/mean(diff(tvec)));
    w = w / sum(w);

    HPC_MUA_sleep=[];
    V1_MUA_sleep=[];
    % speed= filtfilt(w,1,zscore(histcounts(spike_times_sleep,tvec_edges))')';
    for nprobe = 1:length(slow_waves_all) % Get distribution of sleep spike counts
        temp_spike_count =[];
        spike_times = slow_waves_all(nprobe).V1_MUA_spiketimes{nsession};
        % temp_spike_count = filtfilt(w,1,histcounts(spike_times,tvec_edges)')';
        temp_spike_count = histcounts(spike_times,tvec_edges);
        V1_MUA_sleep{nprobe} = temp_spike_count(sleep_index==1);

        spike_times = slow_waves_all(nprobe).HPC_MUA_spiketimes{nsession};
        % temp_spike_count = filtfilt(w,1,histcounts(spike_times,tvec_edges)')';
        temp_spike_count = histcounts(spike_times,tvec_edges);
        HPC_MUA_sleep{nprobe} = temp_spike_count(sleep_index==1);
    end



    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(slow_waves_all(nprobe).V1_MUA_spiketimes{nsession}, event_times, [-0.5 0.5], 0.01);
    binnedArray = filtfilt(w,1,binnedArray')';

    zscored_psth= reshape(binnedArray,1,[]);
    zscored_psth = (zscored_psth-mean(V1_MUA_sleep{nprobe}))/(std(V1_MUA_sleep{nprobe}));
    zscored_psth = reshape(zscored_psth,length(event_times),[]);

end



nprobe = 2;
mprobe = 1;


[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(slow_waves_all(nprobe).V1_MUA_spiketimes{nsession}, event_times, [-0.5 0.5], 0.01);
binnedArray = filtfilt(w,1,binnedArray')';

zscored_psth= reshape(binnedArray,1,[]);
zscored_psth = (zscored_psth-mean(V1_MUA_sleep{nprobe}))/(std(V1_MUA_sleep{nprobe}));
zscored_psth = reshape(zscored_psth,length(event_times),[]);
imagesc(zscored_psth)
colorbar
colormap(flipud(gray))

plot_perievent_spiketime_histogram

%% UP DOWN state and ripple and spindle analysis (backup)

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
            % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
%             load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
            % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
            
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state_merged.mat'));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'));
            % load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
%             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'));
            clusters=clusters_ks3;
        else
            % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
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

        metric_param =[];
        % metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
        metric_param.region = @(x) contains(x,'V1_L');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        V1_spikes{1}=[selected_clusters.spike_id selected_clusters.spike_times];
        metric_param.region = @(x) contains(x,'V1_R');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        V1_spikes{2}=[selected_clusters.spike_id selected_clusters.spike_times];

        metric_param.region = @(x) contains(x,'HPC_L');
        % metric_param.cell_type = @(x) x==1;
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        HPC_spikes{1}=[selected_clusters.spike_id selected_clusters.spike_times];
        metric_param.region = @(x) contains(x,'HPC_R');
        % metric_param.cell_type = @(x) x==1;
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        HPC_spikes{2}=[selected_clusters.spike_id selected_clusters.spike_times];
        metric_param.region = @(x) contains(x,'HPC');
        % metric_param.cell_type = @(x) x==1;
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        HPC_spikes{3}=[selected_clusters.spike_id selected_clusters.spike_times];

        %%% Get spike counts
        tvec = LFP(1).tvec;
        mobility_interp = interp1(selected_clusters.tvec{1},double(selected_clusters.mobility_thresholded{1}),tvec,'previous');
        tvec_edges = [tvec(1)-1/(1/mean(diff(tvec))*2) tvec+1/(1/mean(diff(tvec))*2)];
        HPC_spike_counts=[];
        V1_spike_counts=[];
        w = gausswin(0.03*1/mean(diff(tvec)));
        w = w / sum(w);
        % speed= filtfilt(w,1,zscore(histcounts(spike_times_sleep,tvec_edges))')';
        for nprobe = 1:length(V1_spikes)
            spike_times = V1_spikes{nprobe}(:,2);
            spike_speed =  interp1(tvec,mobility_interp,spike_times,'nearest');
            spike_times_sleep = spike_times(spike_speed < 1);
            V1_spike_counts{nprobe} = filtfilt(w,1,zscore(histcounts(spike_times_sleep,tvec_edges))')';

            spike_times = HPC_spikes{nprobe}(:,2);
            spike_speed =  interp1(tvec,mobility_interp,spike_times,'nearest');
            spike_times_sleep = spike_times(spike_speed < 1);
            HPC_spike_counts{nprobe} = filtfilt(w,1,zscore(histcounts(spike_times_sleep,tvec_edges))')';
        end
        spike_times = HPC_spikes{3}(:,2);
        spike_speed =  interp1(tvec,mobility_interp,spike_times,'nearest');
        spike_times_sleep = spike_times(spike_speed < 1);
        HPC_spike_counts{3} = filtfilt(w,1,zscore(histcounts(spike_times_sleep,tvec_edges))')';

        % 


        % Probability of SWR during normalised UP duration
        time_wondows = [-0.5 0.5];
        time_bin = 0.01;
        num_bins=15; % divide one UP event into 20 bins

        for nprobe = 1:2
            figure
            % DOWN state
            probabilities = calculate_relative_event_probability(slow_waves(nprobe).ints.DOWN,ripples(1).SWS_peaktimes,num_bins);
            nexttile
            plot(linspace(0,1,num_bins),cumsum(probabilities));
            xline(0,'r')
            title('cumlative probability of Left ripples during DOWN phases')
            xlabel('Normalised duration of DOWN')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            nexttile
            plot(linspace(0,1,num_bins),probabilities);
            xline(0,'r')
            title('probability of Left ripples during DOWN phases')
            xlabel('Normalised duration of DOWN')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            % UP state
            UP_duration = slow_waves(nprobe).ints.UP(:,2)-slow_waves(nprobe).ints.UP(:,1);

            probabilities = calculate_relative_event_probability(slow_waves(nprobe).ints.UP(find(UP_duration<2),:),ripples(1).SWS_peaktimes,num_bins);
            nexttile
            plot(linspace(0,1,num_bins),cumsum(probabilities));
            xline(0,'r')
            title('cumlative probability of Left ripples during UP phases')
            xlabel('Normalised duration of UP')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            nexttile
            plot(linspace(0,1,num_bins),probabilities);
            ylim([0 max(probabilities)+0.005])
            xline(0,'r')
            title('probability of Left ripples during UP phases')
            xlabel('Normalised duration of UP')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


            % DOWN state
            probabilities = calculate_relative_event_probability(slow_waves(nprobe).ints.DOWN,ripples(2).SWS_peaktimes,num_bins);
            nexttile
            plot(linspace(0,1,num_bins),cumsum(probabilities));
            xline(0,'r')
            title('cumlative probability of Right ripples during DOWN phases')
            xlabel('Normalised duration of DOWN')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            nexttile
            plot(linspace(0,1,num_bins),probabilities);
            xline(0,'r')
            title('probability of Right ripples during DOWN phases')
            xlabel('Normalised duration of DOWN')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            % UP state
            probabilities = calculate_relative_event_probability(slow_waves(nprobe).ints.UP(find(UP_duration<2),:),ripples(2).SWS_peaktimes,num_bins);
            nexttile
            plot(linspace(0,1,num_bins),cumsum(probabilities));
            xline(0,'r')
            title('cumlative probability of Right ripples during UP phases')
            xlabel('Normalised duration of UP')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            nexttile
            plot(linspace(0,1,num_bins),probabilities);
            ylim([0 max(probabilities)+0.005])
            xline(0,'r')
            title('probability of Right ripples during UP phases')
            xlabel('Normalised duration of UP')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)           
        end


        for nprobe = 1:2
            figure
            probabilities = calculate_event_probability(ripples(1).SWS_peaktimes, slow_waves(nprobe).ints.DOWN(:,1), [time_wondows(1):time_bin:time_wondows(2)])

            nexttile
            plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],cumsum(probabilities));
            xline(0,'r')
            title('cumlative probability of Left ripples relative to DOWN phases')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            nexttile
            plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],probabilities);
            xline(0,'r')
            title('probability of Left ripples relative to DOWN phases')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            probabilities = calculate_event_probability(ripples(1).SWS_peaktimes, slow_waves(nprobe).ints.UP(:,1), [time_wondows(1):time_bin:time_wondows(2)])
            nexttile
            plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],cumsum(probabilities));
            xline(0,'r')
            title('cumlative probability of Left ripples relative to DOWN phases')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            nexttile
            plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],probabilities);
            xline(0,'r')
            title('probability of Left ripples relative to DOWN phases')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


            probabilities = calculate_event_probability(ripples(2).SWS_peaktimes, slow_waves(nprobe).ints.DOWN(:,1), [time_wondows(1):time_bin:time_wondows(2)])

            nexttile
            plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],cumsum(probabilities));
            xline(0,'r')
            title('cumlative probability of Left ripples relative to DOWN phases')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            nexttile
            plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],probabilities);
            xline(0,'r')
            title('probability of Left ripples relative to DOWN phases')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            probabilities = calculate_event_probability(ripples(2).SWS_peaktimes, slow_waves(nprobe).ints.UP(:,1), [time_wondows(1):time_bin:time_wondows(2)])
            nexttile
            plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],cumsum(probabilities));
            xline(0,'r')
            title('cumlative probability of Left ripples relative to DOWN phases')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            nexttile
            plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],probabilities);
            xline(0,'r')
            title('probability of Left ripples relative to DOWN phases')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        end


    end
end



%%%%%%%%%%%%%%%% event veiwer
event_times1 = slow_waves(1).ints.UP;
event_times2 = slow_waves(2).ints.UP;

[tvec_epoches,idx] = RestrictInts(tvec',behavioural_state_merged.SWS); %Replace with InInterval
plot(tvec,HPC_spike_counts{3}');hold on
spacer = 1+max(HPC_spike_counts{3}(idx));
plot(tvec,V1_spike_counts{1}'+ spacer,'r');
spacer= 1+max(HPC_spike_counts{3}(idx))+1+max(V1_spike_counts{1}(idx));
plot(tvec,V1_spike_counts{2}' + spacer,'b');
for n = 1:30
    %             xline(min(tvec(tvec>ripples(1).UP_DOWN_transition_onset(n))),'b')
    %             xline(max(tvec(tvec<ripples(1).UP_DOWN_transition_onset(n))),'r')
    %             xline(min(tvec(tvec>event_times(n,2))),'b')
    %             xline(max(tvec(tvec<event_times(n,1))),'r')
    rectangle('Position',[min(tvec(tvec>event_times1(n,1))) 0.5...
        min(tvec(tvec>event_times1(n,2)))-max(tvec(tvec<event_times1(n,1))),...
        20],...
        'FaceColor',[1 0 0],'EdgeColor','none','FaceAlpha',0.2)
end

for n = 1:30
    %             xline(min(tvec(tvec>ripples(1).UP_DOWN_transition_onset(n))),'b')
    %             xline(max(tvec(tvec<ripples(1).UP_DOWN_transition_onset(n))),'r')
    %             xline(min(tvec(tvec>event_times(n,2))),'b')
    %             xline(max(tvec(tvec<event_times(n,1))),'r')
    rectangle('Position',[min(tvec(tvec>event_times2(n,1))) 0.5...
        min(tvec(tvec>event_times2(n,2)))-max(tvec(tvec<event_times2(n,1))),...
        20],...
        'FaceColor',[0 0 1],'EdgeColor','none','FaceAlpha',0.2)
end

% HPC frames
for nprobe = 1:2
    probe_no = session_info(n).probe(nprobe).probe_id+1;
    if ~isempty(behavioural_state_merged.SWS)
        % SWS_no_ripples = SubtractIntervals(behavioural_state_merged.SWS,[ripples(probe_no).onset ripples(probe_no).offset]);
        [HPC_frames] = detect_candidate_frames_masa(tvec,HPC_spikes{probe_no}(:,2),behavioural_state_merged.SWS,0.01,options);
    end
end

%%%%%%%%%%% test detection
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

