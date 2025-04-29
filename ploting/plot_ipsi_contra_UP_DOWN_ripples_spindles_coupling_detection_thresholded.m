function plot_ipsi_contra_UP_DOWN_spindles_coupling_detection_thresholded(slow_waves_all,ripples_all,spindles_all)

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
% load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
% % load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
% load(fullfile(analysis_folder,'ripples_all_POST.mat'))
% load(fullfile(analysis_folder,'spindles_all_POST.mat'))
% load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov_normalised.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov.mat'));

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole_baseline.mat'));
probability_normalised_whole_baseline = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline.mat'));
probability_psth_whole_baseline = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_baseline.mat'));
probability_normalised_baseline = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_baseline.mat'));
probability_psth_baseline = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
probability_normalised_whole = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised.mat'));
probability_normalised = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability.mat'));
probability_psth = probability;

PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');
PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

% %%%%%%%
% % Phase-amplitude coupling
% nBins = 18;
% edges = linspace(-pi, pi, nBins+1);
% 
% for nchannel = 1:size(SO_phase_ripples,1)
%     for mchannel = 1:size(ripple_peak_amplitude,1)
%         % Phase-amplitude coupling
%         SO_phase_ripples(:,nevent) = ripple_peak_amplitude(tidx,:);
%         spindle_phase_ripples(:,nevent) = ripple_peak_amplitude(tidx,:);
% 
% 
%         [~,~,binIdx] = histcounts(SO_phase_ripples(nchannel,:), edges);
%         % Mean amplitude in each phase bin
%         ampByPhase = accumarray(binIdx(binIdx>0), rippleAmp(binIdx>0), [nBins 1], @mean);
%         p = ampByPhase / sum(ampByPhase); % normalize to get probability distribution
%     end
% end
% if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
%     mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
% end
% save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])


% SO_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
% spindle_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
% SO_spindles_coupling = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));


% Find reference channel/shank
cortex_ref_shank = [];
HPC_ref_shank = [];

for nsession = 1:max(ripples_all(1).session_count)
    for probe_no = 1:2
        cortex_ref_shank(nsession,probe_no) = find(slow_waves_all(probe_no).shank_id{nsession} == slow_waves_all(probe_no).shank{nsession}(slow_waves_all(probe_no).channel{nsession} == slow_waves_all(probe_no).best_channel(nsession))...
            &slow_waves_all(probe_no).probe_hemisphere{nsession} == probe_no);
        % [~,idx] = min(abs(ripples_all(probe_no).SWR_peaktimes{nsession}' - ripples_all(probe_no).peaktimes(ripples_all(probe_no).session_count==nsession))');
        % ripple_counts = histcounts(idx,length(ripples_all(probe_no).shank_id{nsession}));
        % [~,HPC_ref_shank(nsession,probe_no)] = max(ripple_counts);

        shank_id = find(ripples_all(probe_no).probe_hemisphere{nsession} == probe_no);
        HPC_ref_shank(nsession,probe_no) = shank_id(ripples_all(probe_no).best_channel(nsession));

    end
end


%% Grouping of DOWN-UP and UP-DOWN and ripples based on 
load(fullfile(analysis_folder,'V1-HPC sleep interaction','k_cluster_ipsi_contra_events.mat'),'k_cluster');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','MUA_PSTH_merged_thresholded.mat'),'MUA_PSTH_merged_thresholded');


%%%%%%%%%%%%% Overlaps
% Initialize output
L_overlap_idx = [];
R_overlap_idx = [];
overlap_idx = [];
non_overlap_idx = [];
merged_idx = [];
% for nsession = 1:max(slow_waves_all.UP_session_count)

%     L_ints = merged_event_info.UP_ints((merged_event_info.UP_ints(:,1) - nsession * 1000000) > 0 &merged_event_info.UP_hemisphere_id == 1 & (merged_event_info.UP_ints(:,1) - nsession * 1000000) < 1000000,:);
%     R_ints = merged_event_info.UP_ints((merged_event_info.UP_ints(:,1) - nsession * 1000000) > 0 &merged_event_info.UP_hemisphere_id == 2 & (merged_event_info.UP_ints(:,1) - nsession * 1000000) < 1000000,:);
L_idx = find(merged_event_info.UP_hemisphere_id == 1);
R_idx = find(merged_event_info.UP_hemisphere_id == 2);
L_ints = merged_event_info.UP_ints(L_idx,:);
R_ints = merged_event_info.UP_ints(R_idx,:);

windows_threshold = [0.05,0.1,0.2];

for ngroup = 1:4

    % Track indices into merged_event_info.UP_ints
    ipsi_overlap_idx{ngroup} = [];     % leading events involved in overlaps
    contra_overlap_idx{ngroup} = [];     % lagging events involved in overlaps
    unique_lags{ngroup} = [];
    all_lags{ngroup} = [];

    L_overlap_idx{ngroup} = [];
    R_overlap_idx{ngroup} =[];
    L_lags{ngroup} = [];
    R_lags{ngroup} =[];
    lags = [];

    all_overlap_idx{ngroup} = [];
    overlap_idx{ngroup} = [];
    non_overlap_idx{ngroup} = [];

    for iL = 1:size(L_ints,1)
        % Get current L interval

        L_start = L_ints(iL,1);
        if ngroup == 4
            L_end   = L_ints(iL,2);
        else
            L_end   = L_ints(iL,1)+windows_threshold(ngroup);
        end
        % Check overlap with all R intervals
        if ngroup == 4
            is_overlap = (R_ints(:,1) <= L_end) & (R_ints(:,2) >= L_start);
        else
            is_overlap = (R_ints(:,1) <= L_end) & (R_ints(:,1)+windows_threshold(ngroup) >= L_start);
        end

        if any(is_overlap)
            overlapping_R = find(is_overlap);

            for j = 1:length(overlapping_R)
                iR = overlapping_R(j);

                % Store overlapping pair indices
                L_overlap_idx{ngroup}(end+1) = L_idx(iL);
                R_overlap_idx{ngroup}(end+1) = R_idx(iR);
                L_lags{ngroup}(end+1) = L_ints(iL,1) - R_ints(iR,1);
                R_lags{ngroup}(end+1) = R_ints(iR,1) - L_ints(iL,1) ;

            end
        end
    end

    % Unique overlap winners
    lags = [L_lags{ngroup} R_lags{ngroup}];
    all_overlap_idx{ngroup} = [L_overlap_idx{ngroup} R_overlap_idx{ngroup}];
    [overlap_idx{ngroup},~] = unique(all_overlap_idx{ngroup});
    % unique_overlap_idx{ngroup} = lags(index);
    all_lags{ngroup} = lags;

    % Non-overlapping L and R
    nonoverlap_L_idx = setdiff(L_idx, unique(L_overlap_idx{ngroup}));
    nonoverlap_R_idx = setdiff(R_idx, unique(R_overlap_idx{ngroup}));

    non_overlap_idx{ngroup} = unique([nonoverlap_L_idx; nonoverlap_R_idx]);

end
% 
% 
% probability = probability_psth_whole;
%%

%% Ripple probabilities during unilaterally biased and bilaterally synchronised UP
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
probability_normalised_whole = probability;

% ipsi_V1_MUA_baseline = [PSTH_MUA_baseline(1).L_V1_UP; PSTH_MUA_baseline(2).R_V1_UP];
% contra_V1_MUA_baseline = [PSTH_MUA_baseline(1).R_V1_UP; PSTH_MUA_baseline(2).L_V1_UP];

% ipsi_probability = [PSTH_MUA(1).L_V1_UP; PSTH_MUA(2).R_V1_UP];
% contra_probability = [PSTH_MUA(2).L_V1_UP; PSTH_MUA(1).R_V1_UP];



sync_threshold = mean(abs([event_info(1).UP_lag_threshold_low event_info(2).UP_lag_threshold_low event_info(1).UP_lag_threshold_high event_info(2).UP_lag_threshold_high]));
sync_threshold1 = mean(abs([event_info(1).ripples_lag_threshold_low event_info(2).UP_lag_threshold_low event_info(1).UP_lag_threshold_high event_info(2).UP_lag_threshold_high]));

colour_lines=[];
colour_lines = [0,90,50;74,20,134;228,42,168]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
% colour_lines{1} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{2} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 


% Calculate ripples probability according to unilateral vs bilateral events
ripples_times =  merged_event_info.ripples_ints;
ripples_group_info = merged_event_info.ripples_group_id;
ripples_hemisphere_id = merged_event_info.ripples_hemisphere_id;

event_times = merged_event_info.UP_ints;
hemisphere_id = merged_event_info.UP_hemisphere_id;
lag_diff = merged_event_info.UP_lag_diff;
group_id = merged_event_info.UP_group_id;

binnedArray=[];
time_windows = [-1 1];
time_bin = 0.02;
timebin_edge = time_windows(1):time_bin:time_windows(end);
bins_centre = timebin_edge(1)+time_bin/2:time_bin:timebin_edge(end)-time_bin/2;
probability_merged = [];


% Assign ripple hemisphere identity for unilateral and bilateral events
% according to lags
for ngroup = 1
    for nprobe = 1:2
        if ngroup ==1
            index1 = intersect(find(ripples_group_info == ngroup & ripples_hemisphere_id== nprobe),merged_event_info.ripples_index);%ipsi leading
            index2 = intersect(find(ripples_group_info == ngroup & ripples_hemisphere_id~= nprobe),merged_event_info.ripples_index);%contra leading
            index3 = intersect(find(ripples_group_info == ngroup & ripples_hemisphere_id~= nprobe),merged_event_info.ripples_index);%bilateral
        end
        ints = merged_event_info.UP_ints(hemisphere_id,:);


        [~,temp,~] = calculate_event_probability(ripples_times(index1,:), ints(:,1), time_windows,0);
        timebin_edges_all = ints(:,1) + bins_centre;  % Absolute times of peri-event window
        for i = 1:size(ints,1)
            % Previous DOWN (skip if this is the first UP)
            if i > 1
                prev_offset = ints(i-1,2);
                % Find peri-time indices within the previous UP state
                mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                temp(i, mask_prev) = NaN;
            end

            % Next DOWN (skip if this is the last UP)
            if i < size(ints,1)
                next_onset = ints(i+1,1);
                % Find peri-time indices within the next UP state
                mask_next = timebin_edges_all(i,:) >= next_onset;
                temp(i, mask_next) = NaN;
            end
        end
    end
    if ngroup <= max(unique(ripples_group_info)) % first four is ripples by cluster
        probability_merged(nprobe).ripples_UP{ngroup} = temp;
    else
        probability_merged(nprobe).ripples_UP_all = temp;
    end
end

probability_merged = merged_probability_psth_whole;
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_merged.mat'),'probability_merged');

% [probabilities,event_index,normalized_duration,binnedArray] = calculate_relative_event_probability(event_A, event_B,num_bins,shuffle_options);

merged_event_info.ripples_hemisphere_id(merged_event_info.ripples_index)


% DOWN-UP info
lag_diff = merged_event_info.UP_lag_diff;
plv_diff = merged_event_info.UP_plv_diff;
corr_diff = merged_event_info.UP_corr_diff;
ipsi_lag = [ipsi_lag_DU{1} ipsi_lag_DU{2}]';
contra_lag = [contra_lag_DU{1} contra_lag_DU{2}]';
group_id = merged_event_info.UP_group_id;


event_idx = [];
event_idx{1} = {intersect(find((group_id == 3 & lag_diff < -0.1)|(group_id == 3 & lag_diff > 0.1) ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 2 & lag_diff < 0.1 & lag_diff > -0.1 ),merged_event_info.DOWN_index),intersect(find((group_id == 1 & lag_diff < -0.1)|(group_id == 3 & lag_diff > 0.1) ),merged_event_info.DOWN_index)};

event_idx{2} = {intersect(find(group_id == 3 & lag_diff <prctile(lag_diff(group_id == 3& lag_diff <0),50) ),merged_event_info.UP_index),intersect(find(group_id == 3 & lag_diff <0 &lag_diff >prctile(lag_diff(group_id == 3& lag_diff <0),50) ),merged_event_info.UP_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.UP_index),...
    intersect(find(group_id == 1 & lag_diff>0 &lag_diff <prctile(lag_diff(group_id == 1& lag_diff >0),50) ),merged_event_info.UP_index),intersect(find(group_id == 1 &lag_diff >prctile(lag_diff(group_id == 1& lag_diff >0),50) ),merged_event_info.UP_index)};

event_idx{3} = {intersect(find(group_id == 1 & lag_diff <prctile(lag_diff(group_id == 1& lag_diff <0),50) ),merged_event_info.UP_index),intersect(find(group_id == 1 & lag_diff <0 &lag_diff >prctile(lag_diff(group_id == 1& lag_diff <0),50) ),merged_event_info.UP_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.UP_index),...
    intersect(find(group_id == 3 & lag_diff>0 &lag_diff <prctile(lag_diff(group_id == 3& lag_diff >0),50) ),merged_event_info.UP_index),intersect(find(group_id == 3 &lag_diff >prctile(lag_diff(group_id == 3& lag_diff >0),50) ),merged_event_info.UP_index)};

event_idx{4} = {intersect(find(group_id ~= 2 & lag_diff <prctile(lag_diff(group_id ~= 2 & lag_diff <0),50) ),merged_event_info.UP_index),intersect(find(group_id ~= 2 &  lag_diff <0 &lag_diff >prctile(lag_diff(group_id ~= 2 & lag_diff <0),50) ),merged_event_info.UP_index),...
    intersect(find(group_id == 2 & lag_diff <sync_threshold & lag_diff >-sync_threshold),merged_event_info.UP_index),...
    intersect(find(group_id ~= 2 & lag_diff>0 & lag_diff <prctile(lag_diff(group_id ~= 2 & lag_diff >0),50) ),merged_event_info.UP_index),intersect(find(group_id ~= 2 & lag_diff >prctile(lag_diff(group_id ~= 2 & lag_diff >0),50) ),merged_event_info.UP_index)}

event_idx{5} = {intersect(find(group_id == 3 & corr_diff >prctile(corr_diff(group_id == 3),50) ),merged_event_info.UP_index),intersect(find(group_id == 3 & corr_diff <prctile(corr_diff(group_id == 3),50) ),merged_event_info.UP_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.UP_index),...
    intersect(find(group_id == 1 & corr_diff <prctile(corr_diff(group_id == 1),50) ),merged_event_info.UP_index),intersect(find(group_id == 1 & corr_diff >prctile(corr_diff(group_id == 1),50) ),merged_event_info.UP_index)};

event_idx{6} = {intersect(find(group_id == 3 & lag_diff <0 ),merged_event_info.UP_index),intersect(find(group_id == 3 & lag_diff >0 ),merged_event_info.UP_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.UP_index),...
    intersect(find(group_id == 1 & lag_diff <0 ),merged_event_info.UP_index),intersect(find(group_id == 1 & lag_diff >0 ),merged_event_info.UP_index)};


group_name=[];
group_name{1} = {'Ipsi dominant','Bilaterally synchronised','Contra dominant','Shuffled'};
group_name{2} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% dominant clusters
group_name{3} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% non-dominant clusters
group_name{4} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% purely based on lags exlcuding cluster 2
group_name{5} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};

group_name{6} = {'Ipsi dominant ipsi leading','Ipsi dominant contra leading','Bilaterally synchronised','Contra dominant ipsi leading','Contra dominant contra leading','Shuffled'};% dominant clusters

title_names = {'Ipsi-contra DOWN-UP ripples probability by three clusters', 'Ipsi-contra DOWN-UP ripples probability by dominant clusters + 50% lag diff','Ipsi-contra DOWN-UP ripples probability by non-dominant clusters + 50% lag diff',...
    'Ipsi-contra DOWN-UP ripples probability by cluster 1 and 3 merged + 50% lag diff','Ipsi-contra DOWN-UP ripples probability by corr diff','Ipsi-contra DOWN-UP ripples probability by 5 clusters'}
% colour_lines = [0,90,50;74,20,134]/256; % Green Purple

% colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral
% 
% colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 


ipsi_probability = [probability_psth_whole(1).L_ripples_UP; probability_psth_whole(2).R_ripples_UP];
contra_probability = [probability_psth_whole(1).R_ripples_UP; probability_psth_whole(2).L_ripples_UP];

ipsi_probability_baseline = [probability_psth_whole_baseline(1).L_ripples_UP; probability_psth_whole_baseline(2).R_ripples_UP];
contra_probability_baseline = [probability_psth_whole_baseline(1).R_ripples_UP; probability_psth_whole_baseline(2).L_ripples_UP];

%%%%% calculate shuffled MUA baseline
%%% Ipsi
binnedArray = ipsi_probability_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap = temp;

%%% Contra
binnedArray = contra_probability_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

contra_baseline_bootstrap = temp;

%%% Contra
binnedArray = ipsi_probability_baseline-contra_probability_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap = temp;


time_wondows = [-1 1];
time_bin = 0.02;
for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [70 300 1700 500];
    fig.Name =title_names{ngroup};

    if ngroup ==1
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        index =event_idx{ngroup}{i};
        binnedArray = ipsi_probability(index,:);

        % nprobe = 1;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
        temp=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
            temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        y = mean(binnedArray,'omitnan');
        LCI = prctile(temp,2.5);
        UCI = prctile(temp,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline

    y = mean(ipsi_baseline_bootstrap,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(ipsi_baseline_bootstrap,2.5);
    UCI = prctile(ipsi_baseline_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([0 0.1])
    title('Ipsi MUA')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        index =event_idx{ngroup}{i};
        binnedArray = contra_probability(index,:);

        % nprobe = 1;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
        temp=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
            temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        y = mean(binnedArray,'omitnan');
        LCI = prctile(temp,2.5);
        UCI = prctile(temp,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline

    y = mean(contra_baseline_bootstrap,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(contra_baseline_bootstrap,2.5);
    UCI = prctile(contra_baseline_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-1 0.5])

    % xline(0,'r')
    title('Contra MUA')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        index =event_idx{ngroup}{i};
        binnedArray = ipsi_probability(index,:)-contra_probability(index,:);

        % nprobe = 1;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
        temp=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
            temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        y = mean(binnedArray,'omitnan');
        LCI = prctile(temp,2.5);
        UCI = prctile(temp,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline

    y = mean(ipsi_contra_diff_baseline_bootstrap,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(ipsi_contra_diff_baseline_bootstrap,2.5);
    UCI = prctile(ipsi_contra_diff_baseline_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-1 0.5])


    % xline(0,'r')
    title('Ipsi-contra MUA diff')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end


%%%%%%%%%%%%%%%%%%%%%%
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])


%%

ipsi_probability = [probability_normalised_whole(1).L_ripples_UP; probability_normalised_whole(2).R_ripples_UP];
contra_probability = [probability_normalised_whole(1).R_ripples_UP; probability_normalised_whole(2).L_ripples_UP];

ipsi_HPC_MUA_baseline = [PSTH_MUA_baseline(1).L_HPC_UP; PSTH_MUA_baseline(2).R_HPC_UP];
contra_HPC_MUA_baseline = [PSTH_MUA_baseline(1).R_HPC_UP; PSTH_MUA_baseline(2).L_HPC_UP];


%% Visualise and compare first vs final ripples in terms of ipsi-contra PLV diff, corr diff and lag diff
lags = ipsi_lag_ripples{nprobe} - contra_lag_ripples{nprobe};
corrs = ipsi_corr_ripples{nprobe} - contra_corr_ripples{nprobe};
plvs = ipsi_plv_ripples{nprobe} - contra_plv_ripples{nprobe};



% SO_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
% spindle_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
% SO_spindles_coupling = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
