function add_on_UP_DOWN_detection(session_info,stimulus_name,best_channels)

options = session_info.probe(1);
% load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));
if isempty(DIR)
    disp('No extracted clusters')
    return
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

% load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_spindle_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'spindles')
% load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'),'ripples')
load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'))
load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'),'spindles')
load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'),'slow_waves')
load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state_merged.mat'),'behavioural_state_merged')
load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'),'PSD','power')% save PSD for the sleep session
% save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP','-v7.3')
load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP')
load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));% clusters for MUA
clusters = clusters_ks4;
load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name,'Chronic'))),'session_clusters'); % Session clusters for SUA

slow_waves2 = slow_waves;
spindles2 = spindles;
% for nprobe = 1:2
%     temp.power = [];
%     temp.frequency = [];
%     temp.PSD_slope = [];
%     temp.timebin_edges = [];
%     temp.DOWN_PSD_slope = [];
%     temp.UP_PSD_slope = [];
%     temp.DOWN_delta_power = [];
%     temp.DOWN_peaks_zscore = [];
%     temp.DOWN_peaktimes = [];
%     temp.DOWN_peaks_latency = [];
%     temp.DOWN_traveling = [];
%     temp.probe_hemisphere = [];
%     temp.shank_id = [];
%     % slow_waves.power = [];
% end

% for nprobe = 1:length(session_info.probe)
for nprobe = 1:2
    probe_no = nprobe;
    
    zscore_min = 0;
    zscore_max = 3;
    %                 metric_param.cluster_id = @(x) ismember(x,session_clusters_RUN.cluster_id(spatial_cell_index));

    %%%%% Detect V1 populational bursting events (Candidate events)
    
    metric_param =[];
    %                  metric_param.cluster_id = @(x) ismember(x,session_clusters_RUN.cluster_id(spatial_cell_index));

    if options.probe_hemisphere==1
        metric_param.region = @(x) contains(x,'V1_L');
        V1_clusters(probe_no) = select_clusters(clusters(nprobe),metric_param);
    elseif options.probe_hemisphere==2
        metric_param.region = @(x) contains(x,'V1_R');
        V1_clusters(probe_no) = select_clusters(clusters(nprobe),metric_param);
    end

    % if length(V1_clusters(probe_no).cluster_id)>200
    %     metric_param =create_cluster_selection_params('sorting_option','masa');
    %     if options.probe_hemisphere==1
    %         metric_param.region = @(x) contains(x,'V1_L');
    %         V1_clusters(probe_no) = select_clusters(clusters(nprobe),metric_param);
    %     elseif options.probe_hemisphere==2
    %         metric_param.region = @(x) contains(x,'V1_R');
    %         V1_clusters(probe_no) = select_clusters(clusters(nprobe),metric_param);
    %     end
    % end
    %%%%%% Detect V1 spindle events and slow waves (combined SWS)
    if ~isempty(behavioural_state_merged.SWS)

        %         best_channel =  find(slow_waves(probe_no).shank(slow_waves(probe_no).channel == slow_waves(probe_no).best_channel));
        best_channel =  2;

        temp_V1_channels(nprobe).deltaspikecorr = [];
        temp_V1_channels(nprobe).gammaspikecorr = [];
        temp_V1_channels(nprobe).deltagammacorr = [];
        temp_V1_channels(nprobe).channel = slow_waves(nprobe).channel;
        temp_V1_channels(nprobe).shank = slow_waves(nprobe).shank;
        temp_V1_channels(nprobe).depth = slow_waves(nprobe).depth;
        temp_V1_channels(nprobe).xcoord = slow_waves(nprobe).xcoord;
        % [~,best_channel] = max(LFP(nprobe).best_V1_high_freq_power(:,7));
        temp_V1_channels(nprobe).best_channel = LFP(probe_no).best_V1_high_freq_channel(best_channel);




        temp_V1_channels(nprobe).best_channel = LFP(probe_no).best_V1_high_freq_channel(best_channel);

        tvec = LFP(probe_no).tvec;
        filterparms.deltafilter = [0.5 4.5];%heuristically defined.  room for improvement here.
        filterparms.gammafilter = [100 400];
        filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
        filterparms.gammanormwin = 20; %window for gamma normalization (s)

        if contains(options.SUBJECT,'62')
            temp = DetectSlowWaves_masa('filterparms',filterparms,'time',tvec,'lfp',LFP(probe_no).best_V1_high_freq(best_channel,:),'spikes',V1_clusters(probe_no),'NREMInts',behavioural_state_merged.SWS,'sensitivity',0.7);
        else
            temp = DetectSlowWaves_masa('filterparms',filterparms,'time',tvec,'lfp',LFP(probe_no).best_V1_high_freq(best_channel,:),'spikes',V1_clusters(probe_no),'NREMInts',behavioural_state_merged.SWS,'sensitivity',0.6);
            %         temp = DetectSlowWaves_masa('time',tvec,'lfp',raw_LFP(15,:),'spikes',V1_clusters(probe_no),'NREMInts',behavioural_state_merged.SWS,'sensitivity',0.6);
        end
        Behaviour.mobility_zscore=session_clusters.mobility_zscore{1};
        spindles_temp(nprobe)= FindSpindles_masa(LFP(probe_no).best_V1_high_freq(best_channel,:),LFP(probe_no).tvec','behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
            'noise',[],'passband',[9 17],'thresholds',[1 3],'show','on','best_channel',best_channel);
% 
%         UP_ints = UPints;
%         for nevent = 1:size(UP_ints,1)
%             if nevent ==1  | nevent ==size(UP_ints,1)
%                 tidx = FindInInterval(deltaLFP.timestamps,[UP_ints(nevent,1) UP_ints(nevent,1)+0.1]);
%             elseif UP_ints(nevent+1,1)-0.1 > UP_ints(nevent,2)
%                 tidx = FindInInterval(deltaLFP.timestamps,[UP_ints(nevent,1) UP_ints(nevent,1)+0.1]);
%             elseif  UP_ints(nevent+1,1)-0.1 < UP_ints(nevent,2)
%                 tidx = FindInInterval(deltaLFP.timestamps,[UP_ints(nevent,1) UP_ints(nevent,2)]);
%             end
% 
%             phase_UP(nevent)=angle(mean(exp(1i*deltaLFP.phase(tidx(1):tidx(end))))); % phase
%         end
%         slow_waves(probe_no).mean_phase_UP = phase_UP;

        temp.deltaspikecorr = [];
        temp.gammaspikecorr = [];
        temp.deltagammacorr = [];
        % temp.channel = [];
        % temp.shank = [];
        % temp.depth = [];
        % temp.xcoord = [];
        % % temp.best_channel = temp_V1_channels(nprobe).best_channel;
        % 
        temp.channel = temp_V1_channels(nprobe).channel;
        temp.shank = temp_V1_channels(nprobe).shank;
        temp.depth = temp_V1_channels(nprobe).depth;
        temp.xcoord = temp_V1_channels(nprobe).xcoord;
        temp.best_channel = LFP(probe_no).best_V1_high_freq_channel(best_channel);

        if ~isempty(temp)
            slow_waves(nprobe) = temp;
        end
    end


    if  contains(stimulus_name,'RUN1')|contains(stimulus_name,'RUN2')
        mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',erase(stimulus_name,'Masa2tracks_')),sprintf('Probe%i',probe_no)))
        save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',erase(stimulus_name,'Masa2tracks_')),sprintf('Probe%i',probe_no)),[])
    else
        mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',stimulus_name),sprintf('Probe%i',probe_no)))
        save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',stimulus_name),sprintf('Probe%i',probe_no)),[])
    end
end


spindles1 = spindles;%backup
clear spindles
spindles = spindles_temp;


for nprobe = 1:2
    if contains(stimulus_name,'Sleep')
        % if ~isempty(CA1_clusters(probe_no).cluster_id)
        %     if length(CA1_clusters(probe_no).cluster_id)>10
        %         [reactivations(nprobe).awake_offset,reactivations(nprobe).awake_index] = RestrictInts(reactivations(nprobe).offset',behavioural_state_merged.quietWake);
        %         reactivations(nprobe).awake_onset = reactivations(nprobe).onset(reactivations(nprobe).awake_index)';
        %     end
        % end
        % 
        % if ~isempty(V1_reactivations(nprobe).onset)
        %     [V1_reactivations(nprobe).awake_offset,V1_reactivations(nprobe).awake_index] = RestrictInts(V1_reactivations(nprobe).offset',behavioural_state_merged.quietWake);
        %     V1_reactivations(nprobe).awake_onset = V1_reactivations(nprobe).onset(V1_reactivations(nprobe).awake_index)';
        % end
        % 
        % [ripples(nprobe).awake_offset,ripples(nprobe).awake_index] = RestrictInts(ripples(nprobe).offset,behavioural_state_merged.quietWake);
        % ripples(nprobe).awake_onset = ripples(nprobe).onset(ripples(nprobe).awake_index);
        % ripples(nprobe).awake_peaktimes = ripples(nprobe).peaktimes(ripples(nprobe).awake_index);

        [spindles(nprobe).awake_offset,spindles(nprobe).awake_index] = RestrictInts(spindles(nprobe).offset,behavioural_state_merged.quietWake);
        spindles(nprobe).awake_onset = spindles(nprobe).onset(spindles(nprobe).awake_index);
        spindles(nprobe).awake_peaktimes = spindles(nprobe).peaktimes(spindles(nprobe).awake_index);
        %                 reactivations(nprobe).awake_peaktimes = reactivations(nprobe).peaktimes(reactivations.awake_index);

        if ~isempty(behavioural_state_merged.SWS)
            % if length(CA1_clusters(probe_no).cluster_id)>10
            %     [reactivations(nprobe).SWS_offset,reactivations(nprobe).SWS_index] = RestrictInts(reactivations(nprobe).offset',behavioural_state_merged.SWS);
            %     reactivations(nprobe).SWS_onset = reactivations(nprobe).onset(reactivations(nprobe).SWS_index)';
            % end
            % 
            % if ~isempty(V1_reactivations(nprobe).onset)
            %     [V1_reactivations(nprobe).SWS_offset,V1_reactivations(nprobe).SWS_index] = RestrictInts(V1_reactivations(nprobe).offset',behavioural_state_merged.SWS);
            %     V1_reactivations(nprobe).SWS_onset = V1_reactivations(nprobe).onset(V1_reactivations(nprobe).SWS_index)';
            % end
            % 
            % [ripples(nprobe).SWS_offset,ripples(nprobe).SWS_index] = RestrictInts(ripples(nprobe).offset,behavioural_state_merged.SWS);
            % ripples(nprobe).SWS_onset = ripples(nprobe).onset(ripples(nprobe).SWS_index);
            % ripples(nprobe).SWS_peaktimes = ripples(nprobe).peaktimes(ripples(nprobe).SWS_index);

            [spindles(nprobe).SWS_offset,spindles(nprobe).SWS_index] = RestrictInts(spindles(nprobe).offset,behavioural_state_merged.SWS);
            spindles(nprobe).SWS_onset = spindles(nprobe).onset(spindles(nprobe).SWS_index);
            spindles(nprobe).SWS_peaktimes = spindles(nprobe).peaktimes(spindles(nprobe).SWS_index);
        end
    else

    end

end

for nprobe =1:2
    probe_no = session_info.probe(nprobe).probe_id + 1;
    options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

    if isfield(Behaviour,'speed')
        % for event = 1:length(ripples(probe_no).onset)
        %     ripples(probe_no).speed(event) = mean(Behaviour.speed(find(Behaviour.sglxTime >= ripples(probe_no).onset(event) & Behaviour.sglxTime <= ripples(probe_no).offset(event))));
        % end
        if ~isempty(spindles(probe_no).onset)
            for event = 1:length(spindles(probe_no).onset)
                spindles(probe_no).speed(event) = mean(Behaviour.speed(find(Behaviour.sglxTime >= spindles(probe_no).onset(event) & Behaviour.sglxTime <= spindles(probe_no).offset(event))));
            end
        end
    else
        % for event = 1:length(ripples(probe_no).onset)
        %     ripples(probe_no).mobility_zscore(event) = mean(Behaviour.mobility_zscore(find(Behaviour.sglxTime >= ripples(probe_no).onset(event) & Behaviour.sglxTime <= ripples(probe_no).offset(event))));
        % end

        if ~isempty(spindles(probe_no).onset)
            for event = 1:length(spindles(probe_no).onset)
                spindles(probe_no).mobility_zscore(event) = mean(Behaviour.mobility_zscore(find(Behaviour.sglxTime >= spindles(probe_no).onset(event) & Behaviour.sglxTime <= spindles(probe_no).offset(event))));
            end
        end
    end

    % if ~contains(stimulus_name,'RUN') % If reactivation events during lap running
    %     continue
    % end
    % 
    % lap_times(1).start = session_clusters_RUN.start_time_all{1}(session_clusters_RUN.track_ID_all{1}==1);
    % lap_times(1).end = session_clusters_RUN.end_time_all{1}(session_clusters_RUN.track_ID_all{1}==1)';
    % 
    % lap_times(2).start = session_clusters_RUN.start_time_all{1}(session_clusters_RUN.track_ID_all{1}==2);
    % lap_times(2).end = session_clusters_RUN.end_time_all{1}(session_clusters_RUN.track_ID_all{1}==2)';
    % 
    % if contains(stimulus_name,'RUN') % If reactivation events during lap running
    %     if length(CA1_clusters(probe_no).cluster_id)>10
    %         if ~isempty(reactivations(probe_no).onset)
    %             [reactivations(probe_no).T1_offset,reactivations(probe_no).T1_index] = RestrictInts(reactivations(probe_no).offset',[lap_times(1).start lap_times(1).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
    %             reactivations(probe_no).T1_onset = reactivations(probe_no).onset(reactivations(probe_no).T1_index);
    %             reactivations(probe_no).T1_midpoint = reactivations(probe_no).midpoint(reactivations(probe_no).T1_index);
    % 
    %             [reactivations(probe_no).T2_offset,reactivations(probe_no).T2_index] = RestrictInts(reactivations(probe_no).offset',[lap_times(2).start lap_times(2).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
    %             reactivations(probe_no).T2_onset = reactivations(probe_no).onset(reactivations(probe_no).T2_index);
    %             reactivations(probe_no).T2_midpoint = reactivations(probe_no).midpoint(reactivations(probe_no).T2_index);
    %         end
    %     end
    % 
    %     if length(V1_clusters(probe_no).cluster_id)>10
    %         if ~isempty(V1_reactivations(probe_no).onset)
    %             [V1_reactivations(probe_no).T1_offset,V1_reactivations(probe_no).T1_index] = RestrictInts(V1_reactivations(probe_no).offset',[lap_times(1).start lap_times(1).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
    %             V1_reactivations(probe_no).T1_onset = V1_reactivations(probe_no).onset(V1_reactivations(probe_no).T1_index);
    %             V1_reactivations(probe_no).T1_midpoint = V1_reactivations(probe_no).midpoint(V1_reactivations(probe_no).T1_index);
    % 
    %             [V1_reactivations(probe_no).T2_offset,V1_reactivations(probe_no).T2_index] = RestrictInts(V1_reactivations(probe_no).offset',[lap_times(2).start lap_times(2).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
    %             V1_reactivations(probe_no).T2_onset = V1_reactivations(probe_no).onset(V1_reactivations(probe_no).T2_index);
    %             V1_reactivations(probe_no).T2_midpoint = V1_reactivations(probe_no).midpoint(V1_reactivations(probe_no).T2_index);
    %         end
    %     end
    % end
    % 
    % if contains(stimulus_name,'RUN') % If reactivation events during lap running
    %     [ripples(probe_no).T1_offset,ripples(probe_no).T1_index] = RestrictInts(ripples(probe_no).offset,[lap_times(1).start lap_times(1).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
    %     ripples(probe_no).T1_onset = ripples(probe_no).onset(ripples(probe_no).T1_index);
    %     ripples(probe_no).T1_peaktimes = ripples(probe_no).peaktimes(ripples(probe_no).T1_index);
    % 
    %     [ripples(probe_no).T2_offset,ripples(probe_no).T2_index] = RestrictInts(ripples(probe_no).offset,[lap_times(2).start lap_times(2).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
    %     ripples(probe_no).T2_onset = ripples(probe_no).onset(ripples(probe_no).T2_index);
    %     ripples(probe_no).T2_peaktimes = ripples(probe_no).peaktimes(ripples(probe_no).T2_index);
    % end
end
% spindles1(1) = spindles(1);%backup
% spindles = spindles1;

clear lap_times
% save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_spindle_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'spindles')
% save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'ripples')
save(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'),'spindles')
save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'),'slow_waves')

save(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events2.mat'),'spindles2')
save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events2.mat'),'slow_waves2')
% save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP','-v7.3')
% save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP','-v7.3')

close all