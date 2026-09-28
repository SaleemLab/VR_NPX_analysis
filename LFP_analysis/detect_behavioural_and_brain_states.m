function detect_behavioural_and_brain_states_add_on(session_info,stimulus_name,best_channels,varargin)
% SETTINGS
p = inputParser;
addParameter(p, 'V1_best_channel', []);
addParameter(p, 'HPC_best_channel', []);
% addParameter(p, 'ripples_info', []);
% addParameter(p, 'ripples_power', []);
parse(p, varargin{:});

V1_best_channel = p.Results.V1_best_channel;
HPC_best_channel = p.Results.HPC_best_channel;




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

if contains(stimulus_name,'Masa2tracks')
    load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters%s.mat',erase(stimulus_name,'Masa2tracks'))),'session_clusters');
    load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name,'Masa2tracks'))));
    load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks4%s.mat',erase(stimulus_name,'Masa2tracks'))));
    load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name,'Masa2tracks'))),'LFP')
    clusters = clusters_ks4;
else
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
    load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name,'Chronic'))),'session_clusters'); % Session clusters for SUA
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));% clusters for MUA
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP')
    clusters = clusters_ks4;
end


% if contains(stimulus_name,'Masa2tracks')
%     load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_candidate_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'replay','reactivations','replay_combined','reactivations_combined')
%     load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_candidate_events_V1%s.mat',erase(stimulus_name,'Masa2tracks'))),'V1_replay','V1_reactivations')
%     load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_spindle_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'spindles')
%     load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'ripples')
%     load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_slow_wave_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'slow_waves')
%     load(fullfile(options.ANALYSIS_DATAPATH,sprintf('behavioural_state_merged%s.mat',erase(stimulus_name,'Masa2tracks'))),'behavioural_state_merged')
% %     save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name,'Masa2tracks'))),'PSD')
% %     save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name,'Masa2tracks'))),'LFP','-v7.3')
% else
%     load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'),'replay','reactivations','replay_combined','reactivations_combined')
%     load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'),'V1_replay','V1_reactivations')
%     load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'),'ripples')
%     load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'),'spindles')
%     load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'),'slow_waves')
%     load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state_merged.mat'),'behavioural_state_merged')
% %     save(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'),'PSD','power')% save PSD for the sleep session
%     % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP','-v7.3')
% %     save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP','-v7.3')
% end

% From cell structure to all cluster spike time and spike id
session_clusters.spike_times=vertcat(session_clusters.spike_times{:});
session_clusters.spike_id=vertcat(session_clusters.spike_id{:});
[session_clusters.spike_times,index] =sort(session_clusters.spike_times);
session_clusters.spike_id=session_clusters.spike_id(index);

% check cluster region label
if ~isfield(clusters(1),'region')
    for nprobe = 1:length(session_info.probe)
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
    end
end

%
if isempty(V1_best_channel)
    V1_best_channel = zeros(1,length(LFP));
end

if isempty(HPC_best_channel)
    HPC_best_channel = zeros(1,length(LFP));
end

for nprobe = 1:length(session_info.probe)
    options = session_info.probe(nprobe);
    tvec = LFP(nprobe).tvec;
    tvec_edges = [tvec(1)-1/(1/mean(diff(tvec))*2) tvec+1/(1/mean(diff(tvec))*2)];

    if exist('slow_waves') == 1 & V1_best_channel(nprobe) == 0;
        V1_best_channel(nprobe) = find(LFP(nprobe).best_V1_high_freq_channel == slow_waves(nprobe).best_channel);

    elseif V1_best_channel(nprobe) == 0
        %         [~,V1_best_channel(nprobe)] = max(LFP(nprobe).best_V1_high_freq_power(:,7));
        %
        %         if nprobe == 1 % Left hemisphere (shank 1 most anterior)
        %             if LFP(nprobe).best_V1_high_freq_power(V1_best_channel(nprobe),7) < 2*LFP(nprobe).best_V1_high_freq_power(1,7)
        %                 V1_best_channel(nprobe) = 1;
        %             end
        %
        %         elseif nprobe ==2 % Right hemisphere (shank 1 most posterior)
        %             if LFP(nprobe).best_V1_high_freq_power(V1_best_channel(nprobe),7) < 2*LFP(nprobe).best_V1_high_freq_power(end,7)
        %                 V1_best_channel(nprobe) = length(LFP(nprobe).best_V1_high_freq_shank_id);
        %             end
        %         end

        [~,V1_best_channel(nprobe)] = max(LFP(nprobe).best_SO_V1_trough_peak_ratio);
    end


    if exist('ripples') == 1 & HPC_best_channel(nprobe) == 0;
        HPC_best_channel(nprobe) = ripples(nprobe).best_channel;

    elseif HPC_best_channel(nprobe) == 0;
        [~,HPC_best_channel(nprobe)] = max(LFP(nprobe).best_HPC_power(:,6));
    end

    speedTreshold = 1;
    %%%%%  Sleep detection
    if contains(stimulus_name,'Sleep')
        speed = double(session_clusters.mobility_thresholded{1});
        Behaviour.mobility_zscore=session_clusters.mobility_zscore{1};
    else
        speed = Behaviour.speed;
        speed(isnan(speed))=0;
        w = gausswin(9);
        w = w / sum(w);
        speed = filtfilt(w,1,speed')';
    end

    speed = interp1(Behaviour.tvec,speed,tvec,'linear');

%     [freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
%         [tvec' LFP(nprobe).best_V1_high_freq(V1_best_channel(nprobe),:)'],[tvec' LFP(nprobe).best_HPC(HPC_best_channel(nprobe),:)'],...
%         [tvec' speed'],speedTreshold);
    [freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
        [tvec' LFP(nprobe).best_SO_V1(V1_best_channel(nprobe),:)'],[tvec' LFP(nprobe).best_HPC(HPC_best_channel(nprobe),:)'],...
        [tvec' speed'],speedTreshold);

    behavioural_state(nprobe).quietWake = quietWake;
    behavioural_state(nprobe).SWS = SWS;
    behavioural_state(nprobe).REM = REM;
    behavioural_state(nprobe).movement = movement;
end


%%%%% Ripple and Spindles and bursting event

% As long as stable firing durng track running. It doesn't have to
% be single-peaked place cell like cells
spatial_cell_index = find(session_clusters_RUN.odd_even_stability(:,1)>0.95 ...
    | session_clusters_RUN.odd_even_stability(:,2)>0.95);

clear replay reactivations ripples spindles CA1_clusters V1_clusters
clear V1_replay V1_reactivations replay_combined replay_combined slow_waves

%

if length(behavioural_state)>1 & isfield(session_info.probe(1),'probe_hemisphere')
    clear behavioural_state_merged
    % [SWS_combined]= OverlapIntervals(behavioural_state(1).SWS,behavioural_state(2).SWS);
    [REM_merged]= MergeIntervals(behavioural_state(1).REM,behavioural_state(2).REM);% Merged REM
    [SWS_merged]= MergeIntervals(behavioural_state(1).SWS,behavioural_state(2).SWS);% Merged SWS

    [SWS_merged] = SubtractIntervals(SWS_merged,REM_merged);% If REM, it should remain REM rather than SWS
    behavioural_state_merged.SWS = SWS_merged;
    behavioural_state_merged.REM = REM_merged;
    behavioural_state_merged.movement = behavioural_state(1).movement;

    [quietWake_merged]= MergeIntervals(behavioural_state(1).quietWake,behavioural_state(2).quietWake);
    [quietWake_merged] = SubtractIntervals(quietWake_merged,REM_merged);% If REM, it should remain REM rather than quiet awake
    [quietWake_merged] = SubtractIntervals(quietWake_merged,SWS_merged);% If SWS, it should remain SWS rather than quiet awake
    behavioural_state_merged.quietWake = quietWake_merged;

else
    % place holder for MEC-V1 recording?
end

for nprobe = 1:length(session_info.probe)
    options = session_info.probe(nprobe);
    probe_no = session_info.probe(nprobe).probe_id + 1;
    options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
    %                 Behavioural state detection

    %%%%%%%%%%%%%%%%%%
    % UP/Down states and ripple and candidate reactivation events detection
    %%%%%%%%%%%%%%%%%%
    tvec = LFP(nprobe).tvec;
    
    zscore_min = 0;
    zscore_max = 3;
    metric_param =[];
    %                 metric_param.cluster_id = @(x) ismember(x,session_clusters_RUN.cluster_id(spatial_cell_index));

    if options.probe_hemisphere==1
        metric_param.region = @(x) contains(x,'HPC_L');
        CA1_clusters(probe_no) = select_clusters(clusters(nprobe),metric_param);
    elseif options.probe_hemisphere==2
        metric_param.region = @(x) contains(x,'HPC_R');
        CA1_clusters(probe_no) = select_clusters(clusters(nprobe),metric_param);
    end

    if ~isempty(CA1_clusters(probe_no).cluster_id)
        if length(CA1_clusters(probe_no).cluster_id)>10
            best_channel = HPC_best_channel(nprobe);
            [replay(probe_no),reactivations(probe_no)] = detect_candidate_events_masa(LFP(probe_no).tvec,LFP(probe_no).best_HPC(best_channel,:),...
                [CA1_clusters(probe_no).spike_id CA1_clusters(probe_no).spike_times],Behaviour,zscore_min,zscore_max,options);
            %                         if ~isempty(behavioural_state(probe_no).SWS)
            %                             [HPC_frames] = detect_candidate_frames_masa(LFP(probe_no).tvec,V1_clusters(probe_no).spike_times,behavioural_state(probe_no).SWS,0.01,options);
            %                         end
        end
    end

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

    if ~isempty(V1_clusters(probe_no).cluster_id)
        if length(V1_clusters(probe_no).cluster_id)>10
            % [~,best_channel] = max(LFP(nprobe).best_HPC_power(:,6));
            [V1_replay(probe_no),V1_reactivations(probe_no)] = detect_candidate_events_masa(LFP(probe_no).tvec,LFP(probe_no).best_HPC(best_channel,:),...
                [V1_clusters(probe_no).spike_id V1_clusters(probe_no).spike_times],Behaviour,zscore_min,zscore_max,options);
        end
    end

    %%%%%% Detect V1 spindle events and slow waves (combined SWS)
    if ~isempty(behavioural_state_merged.SWS)
        best_channel = V1_best_channel(nprobe);
        
        % if ~isfield(LFP(probe_no),'best_V1')|isempty(LFP(nprobe).best_V1) % if best v1 does not exist

        temp = DetectSlowWaves_masa('time',tvec,'lfp',LFP(probe_no).best_SO_V1(best_channel,:),'spikes',V1_clusters(probe_no),'NREMInts',behavioural_state_merged.SWS);
        [spindles(probe_no)] = FindSpindles_masa(LFP(probe_no).best_SO_V1(best_channel,:),LFP(probe_no).tvec','behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
            'noise',[],'passband',[9 17],'thresholds',[1 2.5],'thresholds',[1 2.5],'show','on');
        % else
        %     [~,best_channel] = max(LFP(nprobe).best_V1_power(:,7));
        %     if probe_no == 1 % Left hemisphere (shank 1 most anterior)
        %         if LFP(probe_no).best_V1_power(best_channel,7) < 2*LFP(probe_no).best_V1_power(1,7)
        %             best_channel = 1;
        %             temp_V1_channels(nprobe).best_channel = LFP(probe_no).best_V1_channel(1);
        %         end
        % 
        %     elseif probe_no ==2 % Right hemisphere (shank 1 most posterior)
        %         if LFP(probe_no).best_V1_power(best_channel,7) < 2*LFP(probe_no).best_V1_power(end,7)
        %             best_channel = length(LFP(probe_no).best_V1_shank_id);
        %             temp_V1_channels(nprobe).best_channel = LFP(probe_no).best_V1_channel(end);
        %         end
        %     end
        % 
        %     temp = DetectSlowWaves_masa('time',tvec,'lfp',LFP(probe_no).best_V1(best_channel,:),'spikes',V1_clusters(probe_no),'NREMInts',behavioural_state_merged.SWS);
        %     [spindles(probe_no)] = FindSpindles_masa(LFP(probe_no).best_V1(best_channel,:),LFP(probe_no).tvec','behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
        %         'noise',[],'passband',[9 17],'thresholds',[1 3],'show','on');
        % end

        temp.deltaspikecorr = [];
        temp.gammaspikecorr = [];
        temp.deltagammacorr = [];
        % temp.channel = [];
        % temp.shank = [];
        % temp.depth = [];
        % temp.xcoord = [];
        % % temp.best_channel = temp_V1_channels(nprobe).best_channel;
        % 
        temp.channel = LFP(nprobe).average_V1_channel;
        temp.shank = LFP(nprobe).average_V1_shank_id;
        temp.depth = LFP(nprobe).average_V1_depth;
        temp.xcoord = LFP(nprobe).xcoord;
        temp.best_channel = LFP(probe_no).best_V1_high_freq_channel(best_channel);
        if ~isempty(temp)
            slow_waves(nprobe) = temp;
        end

    elseif isfield(LFP(nprobe),'best_V1_high_freq')
        best_channel = V1_best_channel(nprobe);

        [spindles(probe_no)] = FindSpindles_masa(LFP(probe_no).best_V1_high_freq(best_channel,:),LFP(probe_no).tvec','behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
            'noise',[],'passband',[9 17],'thresholds',[1 3],'show','on');
        slow_waves(nprobe) = struct();
    else

        slow_waves(nprobe) = struct();
    end


    % Detect CA1 ripple events
    best_channel = HPC_best_channel(nprobe);

    [ripples(probe_no)] = FindRipples_masa(LFP(nprobe).best_HPC(best_channel,:),LFP(probe_no).tvec','behaviour',Behaviour,'minDuration',30,'durations',[30 200],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
        'noise',[],'passband',[125 300],'thresholds',[2 5],'show','on','best_channel',best_channel);

    if  contains(stimulus_name,'RUN1')|contains(stimulus_name,'RUN2')
        mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',erase(stimulus_name,'Masa2tracks_')),sprintf('Probe%i',probe_no)))
        save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',erase(stimulus_name,'Masa2tracks_')),sprintf('Probe%i',probe_no)),[])
    else
        mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',stimulus_name),sprintf('Probe%i',probe_no)))
        save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',stimulus_name),sprintf('Probe%i',probe_no)),[])
    end
end

for nprobe = 1:length(session_info.probe)
    if contains(stimulus_name,'Sleep')
        if ~isempty(CA1_clusters(probe_no).cluster_id)
            if length(CA1_clusters(probe_no).cluster_id)>10
                [reactivations(nprobe).awake_offset,reactivations(nprobe).awake_index] = RestrictInts(reactivations(nprobe).offset',behavioural_state_merged.quietWake);
                reactivations(nprobe).awake_onset = reactivations(nprobe).onset(reactivations(nprobe).awake_index)';
            end
        end

        if ~isempty(V1_reactivations(nprobe).onset)
            [V1_reactivations(nprobe).awake_offset,V1_reactivations(nprobe).awake_index] = RestrictInts(V1_reactivations(nprobe).offset',behavioural_state_merged.quietWake);
            V1_reactivations(nprobe).awake_onset = V1_reactivations(nprobe).onset(V1_reactivations(nprobe).awake_index)';
        end

        [ripples(nprobe).awake_offset,ripples(nprobe).awake_index] = RestrictInts(ripples(nprobe).offset,behavioural_state_merged.quietWake);
        ripples(nprobe).awake_onset = ripples(nprobe).onset(ripples(nprobe).awake_index);
        ripples(nprobe).awake_peaktimes = ripples(nprobe).peaktimes(ripples(nprobe).awake_index);

        [spindles(nprobe).awake_offset,spindles(nprobe).awake_index] = RestrictInts(spindles(nprobe).offset,behavioural_state_merged.quietWake);
        spindles(nprobe).awake_onset = spindles(nprobe).onset(spindles(nprobe).awake_index);
        spindles(nprobe).awake_peaktimes = spindles(nprobe).peaktimes(spindles(nprobe).awake_index);
        %                 reactivations(nprobe).awake_peaktimes = reactivations(nprobe).peaktimes(reactivations.awake_index);

        if ~isempty(behavioural_state(nprobe).SWS)
            if length(CA1_clusters(probe_no).cluster_id)>10
                [reactivations(nprobe).SWS_offset,reactivations(nprobe).SWS_index] = RestrictInts(reactivations(nprobe).offset',behavioural_state_merged.SWS);
                reactivations(nprobe).SWS_onset = reactivations(nprobe).onset(reactivations(nprobe).SWS_index)';
            end

            if ~isempty(V1_reactivations(nprobe).onset)
                [V1_reactivations(nprobe).SWS_offset,V1_reactivations(nprobe).SWS_index] = RestrictInts(V1_reactivations(nprobe).offset',behavioural_state_merged.SWS);
                V1_reactivations(nprobe).SWS_onset = V1_reactivations(nprobe).onset(V1_reactivations(nprobe).SWS_index)';
            end

            [ripples(nprobe).SWS_offset,ripples(nprobe).SWS_index] = RestrictInts(ripples(nprobe).offset,behavioural_state_merged.SWS);
            ripples(nprobe).SWS_onset = ripples(nprobe).onset(ripples(nprobe).SWS_index);
            ripples(nprobe).SWS_peaktimes = ripples(nprobe).peaktimes(ripples(nprobe).SWS_index);

            [spindles(nprobe).SWS_offset,spindles(nprobe).SWS_index] = RestrictInts(spindles(nprobe).offset,behavioural_state_merged.SWS);
            spindles(nprobe).SWS_onset = spindles(nprobe).onset(spindles(nprobe).SWS_index);
            spindles(nprobe).SWS_peaktimes = spindles(nprobe).peaktimes(spindles(nprobe).SWS_index);
        end
    else

    end

end


%%%%%%%%%%%%%%%%%%
% Candidate reactivation events detection (probe combined) with
% spatial cells only
%%%%%%%%%%%%%%%%%%
reactivations_combined= [];
replay_combined = [];
clear CA1_clusters_combined

if length(session_info.probe)>1
    zscore_min = 0;
    zscore_max = 3;

    metric_param =[];
    metric_param.cluster_id = @(x) ismember(x,session_clusters_RUN.cluster_id(spatial_cell_index));
    metric_param.region = @(x) contains(x,'HPC');

    CA1_clusters_combined = select_clusters(session_clusters,metric_param);
    best_channel = HPC_best_channel(nprobe);
    [replay_combined,reactivations_combined] = detect_candidate_events_masa(LFP(nprobe).tvec,LFP(nprobe).best_HPC(best_channel,:),...
        [CA1_clusters_combined.spike_id CA1_clusters_combined.spike_times],Behaviour,zscore_min,zscore_max,options);

    if contains(stimulus_name,'Sleep')
        if length(CA1_clusters_combined.cluster_id)>10
            [reactivations_combined.SWS_offset,reactivations_combined.SWS_index] = RestrictInts(reactivations_combined.offset',behavioural_state_merged.SWS);
            reactivations_combined.SWS_onset = reactivations_combined.onset(reactivations_combined.SWS_index)';

            [reactivations_combined.awake_offset,reactivations_combined.awake_index] = RestrictInts(reactivations_combined.offset',behavioural_state_merged.quietWake);
            reactivations_combined.awake_onset = reactivations_combined.onset(reactivations_combined.awake_index)';
        end
    end

end

for nprobe = 1:length(session_info.probe)
    probe_no = session_info.probe(nprobe).probe_id + 1;
    options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

    if isfield(Behaviour,'speed')
        for event = 1:length(ripples(probe_no).onset)
            ripples(probe_no).speed(event) = mean(Behaviour.speed(find(Behaviour.sglxTime >= ripples(probe_no).onset(event) & Behaviour.sglxTime <= ripples(probe_no).offset(event))));
        end
        if ~isempty(spindles(probe_no).onset)
            for event = 1:length(spindles(probe_no).onset)
                spindles(probe_no).speed(event) = mean(Behaviour.speed(find(Behaviour.sglxTime >= spindles(probe_no).onset(event) & Behaviour.sglxTime <= spindles(probe_no).offset(event))));
            end
        end
    else
        for event = 1:length(ripples(probe_no).onset)
            ripples(probe_no).mobility_zscore(event) = mean(Behaviour.mobility_zscore(find(Behaviour.sglxTime >= ripples(probe_no).onset(event) & Behaviour.sglxTime <= ripples(probe_no).offset(event))));
        end

        if ~isempty(spindles(probe_no).onset)
            for event = 1:length(spindles(probe_no).onset)
                spindles(probe_no).mobility_zscore(event) = mean(Behaviour.mobility_zscore(find(Behaviour.sglxTime >= spindles(probe_no).onset(event) & Behaviour.sglxTime <= spindles(probe_no).offset(event))));
            end
        end
    end

    if ~contains(stimulus_name,'RUN') % If reactivation events during lap running
        continue
    end

    lap_times(1).start = session_clusters_RUN.start_time_all{1}(session_clusters_RUN.track_ID_all{1}==1);
    lap_times(1).end = session_clusters_RUN.end_time_all{1}(session_clusters_RUN.track_ID_all{1}==1)';

    lap_times(2).start = session_clusters_RUN.start_time_all{1}(session_clusters_RUN.track_ID_all{1}==2);
    lap_times(2).end = session_clusters_RUN.end_time_all{1}(session_clusters_RUN.track_ID_all{1}==2)';

    if contains(stimulus_name,'RUN') % If reactivation events during lap running
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

    if contains(stimulus_name,'RUN') % If reactivation events during lap running
        [ripples(probe_no).T1_offset,ripples(probe_no).T1_index] = RestrictInts(ripples(probe_no).offset,[lap_times(1).start lap_times(1).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
        ripples(probe_no).T1_onset = ripples(probe_no).onset(ripples(probe_no).T1_index);
        ripples(probe_no).T1_peaktimes = ripples(probe_no).peaktimes(ripples(probe_no).T1_index);

        [ripples(probe_no).T2_offset,ripples(probe_no).T2_index] = RestrictInts(ripples(probe_no).offset,[lap_times(2).start lap_times(2).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
        ripples(probe_no).T2_onset = ripples(probe_no).onset(ripples(probe_no).T2_index);
        ripples(probe_no).T2_peaktimes = ripples(probe_no).peaktimes(ripples(probe_no).T2_index);
    end
end

if contains(stimulus_name,'RUN')
    if length(session_info.probe)>1
        [reactivations_combined.T1_offset,reactivations_combined.T1_index] = RestrictInts(reactivations_combined.offset',[lap_times(1).start lap_times(1).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
        reactivations_combined.T1_onset = reactivations_combined.onset(reactivations_combined.T1_index);
        reactivations_combined.T1_midpoint = reactivations_combined.midpoint(reactivations_combined.T1_index);

        [reactivations_combined.T2_offset,reactivations_combined.T2_index] = RestrictInts(reactivations_combined.offset',[lap_times(2).start lap_times(2).end]); % Including 2 seconds after each lap finishes (it usually takes 3 second before starting next lap)
        reactivations_combined.T2_onset = reactivations_combined.onset(reactivations_combined.T2_index);
        reactivations_combined.T2_midpoint = reactivations_combined.midpoint(reactivations_combined.T2_index);
    end
end

clear lap_times
if contains(stimulus_name,'Masa2tracks')
    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_candidate_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'replay','reactivations','replay_combined','reactivations_combined')
    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_candidate_events_V1%s.mat',erase(stimulus_name,'Masa2tracks'))),'V1_replay','V1_reactivations')
    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_spindle_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'spindles')
    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'ripples')
    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_slow_wave_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'slow_waves')
    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('behavioural_state_merged%s.mat',erase(stimulus_name,'Masa2tracks'))),'behavioural_state_merged')
%     save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name,'Masa2tracks'))),'PSD')
%     save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name,'Masa2tracks'))),'LFP','-v7.3')
else
    save(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'),'replay','reactivations','replay_combined','reactivations_combined')
    save(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'),'V1_replay','V1_reactivations')
    save(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'),'ripples')
    save(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'),'spindles')
    save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'),'slow_waves')
    save(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state_merged.mat'),'behavioural_state_merged')
%     save(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'),'PSD','power')% save PSD for the sleep session
    % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP','-v7.3')
%     save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP','-v7.3')
end

close all