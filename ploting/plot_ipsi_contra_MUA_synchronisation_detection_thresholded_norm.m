function plot_ipsi_contra_MUA_synchronisation_detection_thresholded_norm(slow_waves_all,ripples_all,spindles_all)

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

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_relative_PSTH_MUA.mat'));
PSTH_MUA = UP_DOWN_relative_PSTH_MUA;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_relative_PSTH_MUA_baseline.mat'));
PSTH_MUA_baseline = UP_DOWN_relative_PSTH_MUA;


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


%% Grouping of DOWN-UP and UP-DOWN and ripples based on detetcion lags
load(fullfile(analysis_folder,'V1-HPC sleep interaction','k_cluster_ipsi_contra_events.mat'),'k_cluster');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','MUA_PSTH_merged_thresholded.mat'),'MUA_PSTH_merged_thresholded');


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
event_types = {'UP','DOWN','ripples','spindles'};
for n = 1:4
 
    L_idx = find(merged_event_info.(sprintf('%s_hemisphere_id',event_types{n})) == 1);
    R_idx = find(merged_event_info.(sprintf('%s_hemisphere_id',event_types{n})) == 2);
    L_ints = merged_event_info.(sprintf('%s_ints',event_types{n}))(L_idx,:);
    R_ints = merged_event_info.(sprintf('%s_ints',event_types{n}))(R_idx,:);

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


    merged_event_info.(sprintf('%s_overlap_idx_all',event_types{n})) = all_overlap_idx;
    % (sprintf('%s_lags_all',event_types{n}))
    merged_event_info.(sprintf('%s_non_overlap_idx',event_types{n})) = non_overlap_idx;
    merged_event_info. (sprintf('%s_overlap_idx',event_types{n})) = overlap_idx;
    merged_event_info.(sprintf('%s_lags_all',event_types{n})) = all_lags;
end


%% Ipsi and contra MUA during normalised DOWN UP with different bilateral synchrony


% probability = []
Delta_peaks_zscore_UD = [];
Delta_peaks_zscore_DU = [];
SWpeakmag_UD = [];
SWpeakmag_DU = [];
index = [];
%%%%% UP info
for nprobe = 1:2
    for nsession = 1:max(slow_waves_all(1).UP_session_count)

        % slow_waves_all(1).SWpeakmag(slow_waves_all(1).UP_session_count == nsession)
        % Find DOWN
        [C,ia,ib]  =intersect(find(slow_waves_all(nprobe).DOWN_session_count == nsession),probability(nprobe).DOWN_all_index);

        Delta_peaks_zscore_UD = [Delta_peaks_zscore_UD;...
            slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession,nprobe),ia)'];
        SWpeakmag_UD = [SWpeakmag_UD;...
            slow_waves_all(nprobe).SWpeakmag(C)];

        % Find DOWN before UP
        [C,ia,ib]  = intersect(slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).DOWN_session_count == nsession,2), slow_waves_all(nprobe).UP_ints(intersect(find(slow_waves_all(nprobe).UP_session_count == nsession),probability(nprobe).UP_all_index),1));

        Delta_peaks_zscore_DU = [Delta_peaks_zscore_DU; slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession,nprobe),ia)'];
        index = [index; intersect(find(slow_waves_all(nprobe).DOWN_session_count == nsession),probability(nprobe).DOWN_all_index)];

        [C,ia,ib]  = intersect(slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).DOWN_session_count == nsession,2), slow_waves_all(nprobe).UP_ints(intersect(find(slow_waves_all(nprobe).UP_session_count == nsession),probability(nprobe).UP_all_index),1));

        temp = find(slow_waves_all(nprobe).DOWN_session_count == nsession);

        SWpeakmag_DU = [SWpeakmag_DU;...
            slow_waves_all(nprobe).SWpeakmag(temp(ia))];
    end
end
% scatter( slow_waves_all(1).SWpeakmag(slow_waves_all(1).DOWN_session_count == nsession),slow_waves_all(1).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession,1),:))
% hold on
% scatter(slow_waves_all(1).SWpeakmag(slow_waves_all(1).DOWN_session_count == nsession)',max(slow_waves_all(1).DOWN_peaks_zscore{nsession}(slow_waves_all(1).probe_hemisphere{nsession}==1,:)))


event_times = merged_event_info.UP_ints;
hemisphere_id = merged_event_info.UP_hemisphere_id;
% lag_diff = merged_event_info.UP_lag_diff;
% group_id = merged_event_info.UP_group_id;
all_overlap_idx = merged_event_info.UP_overlap_idx_all{end};
non_overlap_idx = merged_event_info.UP_non_overlap_idx{end};
all_lags =merged_event_info.UP_lags_all{end};


event_idx = [];

lags =all_lags(all_lags>=-0.15 & all_lags<=0.15);
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lags =all_lags;

% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    event_idx{1}{n} =(all_overlap_idx(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end


lags =abs(all_lags(abs(all_lags)>=0 & abs(all_lags)<=0.15));

% lags =all_lags(all_lags>=0 & all_lags<=0.15);
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lags =abs(all_lags);

for n = 1:length(lag_thresholds)-1
    event_idx{2}{n} =(all_overlap_idx(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end


lags =all_lags(all_lags>=-0.05 & all_lags<=0.05);
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lags =all_lags;

% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    event_idx{3}{n} =(all_overlap_idx(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end


% lags =all_lags(all_lags>=-0.05 & all_lags<=0.05);
power_thresholds = prctile(SWpeakmag_DU(all_overlap_idx(all_lags>=-0.15 & all_lags<=0.15)),[0 20 40 60 80 100]);
delta_power =SWpeakmag_DU;
% power_thresholds = prctile(Delta_peaks_zscore_DU(all_lags>=-0.15 & all_lags<=0.15),[0 20 40 60 80 100]);
% delta_power =Delta_peaks_zscore_DU(all_lags>=-0.15 & all_lags<=0.15);

% lag_thresholds = [-0.2]
for n = 1:length(power_thresholds)-1
    event_idx{4}{n} =intersect(all_overlap_idx(abs(all_lags)<=0.15),find(delta_power>power_thresholds(n)&delta_power <power_thresholds(n+1)));
end

power_thresholds = prctile(SWpeakmag_DU,[0 20 40 60 80 100]);
delta_power =SWpeakmag_DU;
% power_thresholds = prctile(Delta_peaks_zscore_DU,[0 20 40 60 80 100]);
% delta_power =Delta_peaks_zscore_DU;

% lag_thresholds = [-0.2]
for n = 1:length(power_thresholds)-1
    event_idx{5}{n} =(delta_power>power_thresholds(n)&delta_power <power_thresholds(n+1));
end


Duration = [event_info(1).previous_DOWN_duration; event_info(2).previous_DOWN_duration];
duration_thresholds = prctile(Duration,[0 20 40 60 80 100]);

for n = 1:length(duration_thresholds)-1
    event_idx{6}{n} =find(Duration>duration_thresholds(n)&Duration <duration_thresholds(n+1));
end


Duration = [event_info(1).previous_DOWN_duration; event_info(2).previous_DOWN_duration];
duration_thresholds = prctile(Duration(all_overlap_idx(all_lags>=-0.15 & all_lags<=0.15)),[0 20 40 60 80 100]);

for n = 1:length(duration_thresholds)-1
    event_idx{7}{n} =intersect(all_overlap_idx(abs(all_lags)<=0.15),find(Duration>duration_thresholds(n)&Duration <duration_thresholds(n+1)));
end

group_name{1} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};
group_name{2} = {'0-20% lagging','20-40% lagging','40-60% lagging','60-80% lagging','80-100% lagging','Shuffled'};
group_name{3} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};
group_name{4} = {'0-20% delta power','20-40% delta power','40-60% delta power','60-80% delta power','80-100% delta power','Shuffled'};
group_name{5} = {'0-20% delta power','20-40% delta power','40-60% delta power','60-80% delta power','80-100% delta power','Shuffled'};
group_name{6} = {'0-20% duration','20-40% duration','40-60% duration','60-80% duration','80-100% duration','Shuffled'};
group_name{7} = {'0-20% duration','20-40% duration','40-60% duration','60-80% duration','80-100% duration','Shuffled'};


title_names{1} = 'Ipsi-contra normalised UP MUA by 5 lags (150ms windows)';
title_names{2} = 'Ipsi-contra normalised UP MUA by 5 lags (150ms windows abs lag)';
title_names{3} = 'Ipsi-contra normalised UP MUA by 5 lags (50ms windows)';
title_names{4} = 'Ipsi-contra normalised UP MUA by 5 powers (150ms windows)';
title_names{5} = 'Ipsi-contra normalised UP MUA by 5 powers (All events)';
title_names{6} = 'Ipsi-contra normalised UP MUA by 5 duration (All events)';
title_names{7} = 'Ipsi-contra normalised UP MUA by 5 duration (150ms events)';


%%%%%%%%%%%%%%%%%%%%%
ipsi_V1_MUA = [PSTH_MUA(1).L_V1_UP; PSTH_MUA(2).R_V1_UP];
contra_V1_MUA = [PSTH_MUA(1).R_V1_UP; PSTH_MUA(2).L_V1_UP];

ipsi_V1_MUA_baseline = [PSTH_MUA_baseline(1).L_V1_UP; PSTH_MUA_baseline(2).R_V1_UP];
contra_V1_MUA_baseline = [PSTH_MUA_baseline(1).R_V1_UP; PSTH_MUA_baseline(2).L_V1_UP];


ipsi_HPC_MUA = [PSTH_MUA(1).L_HPC_UP; PSTH_MUA(2).R_HPC_UP];
contra_HPC_MUA = [PSTH_MUA(1).R_HPC_UP; PSTH_MUA(2).L_HPC_UP];
ipsi_HPC_MUA_baseline = [PSTH_MUA_baseline(1).L_HPC_UP; PSTH_MUA_baseline(2).R_HPC_UP];
contra_HPC_MUA_baseline = [PSTH_MUA_baseline(1).R_HPC_UP; PSTH_MUA_baseline(2).L_HPC_UP];

%%%%%%%%%%%%%%%%%%%%%%%%%%% V1
%%%%% calculate shuffled MUA baseline
%%% Ipsi
binnedArray = ipsi_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap = temp;

%%% Contra
binnedArray = contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

contra_baseline_bootstrap = temp;

%%% Ipsi-Contra baseline
binnedArray = ipsi_V1_MUA_baseline-contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap = temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HPC
%%%%% calculate shuffled MUA baseline
%%% Ipsi
binnedArray = ipsi_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap1 = temp;

%%% Contra
binnedArray = contra_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

contra_baseline_bootstrap1 = temp;



%%% Ipsi-Contra baseline
binnedArray = ipsi_HPC_MUA_baseline-contra_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap1 = temp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Calculate bootstrapped MUA

% clear probability_merged
% time_wondows = [-1 1];
% time_bin = 0.02;
num_bins=20;
x = linspace(0,1,num_bins);
MUA_PSTH_merged_thresholded.x = x;
MUA_PSTH_merged_thresholded.ipsi_UP_V1 = [];
MUA_PSTH_merged_thresholded.contra_UP_V1 = [];
MUA_PSTH_merged_thresholded.ipsi_contra_diff_UP_V1 = [];

MUA_PSTH_merged_thresholded.ipsi_UP_HPC = [];
MUA_PSTH_merged_thresholded.contra_UP_HPC = [];
MUA_PSTH_merged_thresholded.ipsi_contra_diff_UP_HPC = [];


for ngroup = 1:length(event_idx)+1

    group_index = [];
    if ngroup <= length(event_idx)
        group_index = event_idx{ngroup};

    else
        group_index{1} = (1:size(ipsi_V1_MUA,1))';
    end

    for i = 1:length(group_index)
        index =group_index{i};

        binnedArray1 = ipsi_V1_MUA(index,:);
        binnedArray2 = contra_V1_MUA(index,:);
        binnedArray3 = ipsi_V1_MUA(index,:)-contra_V1_MUA(index,:);
        binnedArray4 = ipsi_HPC_MUA(index,:);
        binnedArray5 = contra_HPC_MUA(index,:);
        binnedArray6 = ipsi_HPC_MUA(index,:)-contra_HPC_MUA(index,:);

        temp1=[];
        temp2=[];
        temp3=[];
        temp4=[];
        temp5=[];
        temp6=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray1,1),size(binnedArray1,1));
            temp1(iBoot,:) =  mean(binnedArray1(event_id,:),'omitnan');
            temp2(iBoot,:) =  mean(binnedArray2(event_id,:),'omitnan');
            temp3(iBoot,:) =  mean(binnedArray3(event_id,:),'omitnan');
            temp4(iBoot,:) =  mean(binnedArray4(event_id,:),'omitnan');
            temp5(iBoot,:) =  mean(binnedArray5(event_id,:),'omitnan');
            temp6(iBoot,:) =  mean(binnedArray6(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        MUA_PSTH_merged_thresholded.ipsi_UP_V1{ngroup}{i} = temp1;
        MUA_PSTH_merged_thresholded.contra_UP_V1{ngroup}{i} = temp2;
        MUA_PSTH_merged_thresholded.ipsi_contra_diff_UP_V1{ngroup}{i} = temp3;
        MUA_PSTH_merged_thresholded.ipsi_UP_HPC{ngroup}{i} = temp4;
        MUA_PSTH_merged_thresholded.contra_UP_HPC{ngroup}{i} = temp5;
        MUA_PSTH_merged_thresholded.ipsi_contra_diff_UP_HPC{ngroup}{i} = temp6;
    end
end

MUA_PSTH_merged_thresholded.ipsi_UP_V1_baseline = ipsi_baseline_bootstrap;
MUA_PSTH_merged_thresholded.contra_UP_V1_baseline = contra_baseline_bootstrap;
MUA_PSTH_merged_thresholded.ipsi_contra_diff_UP_V1_baseline = ipsi_contra_diff_baseline_bootstrap;
MUA_PSTH_merged_thresholded.ipsi_UP_HPC_baseline = ipsi_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.contra_UP_HPC_baseline = contra_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.ipsi_contra_diff_UP_HPC_baseline = ipsi_contra_diff_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.UP_groups = [title_names 'All UP'];
MUA_PSTH_merged_thresholded.UP_index = [event_idx {(1:size(ipsi_V1_MUA,1))'}];



% time_wondows = [-1 1];
% time_bin = 0.02;
num_bins=20;
x = linspace(0,1,num_bins);

%%%%%%%%% Plot DU transition
for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 930];
    fig.Name =title_names{ngroup};

    if ngroup ==1  | ngroup ==3
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    elseif ngroup==2
        colour_lines = [228,42,168;222,119,174;184,225,134;127,188,65;0,90,50]/256; % From magenta to dark purple to indicate globally synchrnous to ipsi lagging
    end


    if ngroup >=3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end
    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_UP_V1{ngroup}{i};
        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_UP_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-0.3 0.4])

    title('V1 Ipsi MUA')
    xlabel('Normalised duration of UP')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    if ngroup >=3
%         colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.contra_UP_V1{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.contra_UP_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);
    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.3 0.4])

    % xline(0,'r')
    title('V1 Contra MUA')
    xlabel('Normalised duration of UP')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    if ngroup >=3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_UP_V1{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_UP_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.3 0.4])


    % xline(0,'r')
    title('V1 Ipsi-contra MUA diff')
    xlabel('Normalised duration of UP')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%%%%% HPC MUA
    if ngroup >=3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_UP_HPC{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_UP_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-0.3 0.2])
    title('HPC Ipsi MUA')
    xlabel('Normalised duration of UP')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    if ngroup >=3
%         colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.contra_UP_HPC{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.contra_UP_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);


    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.3 0.2])

    % xline(0,'r')
    title('HPC Contra MUA')
    xlabel('Normalised duration of UP')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    if ngroup >=3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_UP_HPC{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_UP_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);


    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.3 0.2])


    % xline(0,'r')
    title('HPC Ipsi-contra MUA diff')
    xlabel('Normalised duration of UP')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')

%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%% All DOWN UP rpples
event_averaging_scale = 10;
x = linspace(0,1,20);

% for ngroup = 1:length(event_idx)
ngroup = 2;
fig = figure('Color','w');
fig.Position = [350 59 1650/3*2 465];
fig.Name = 'Left-Right combined Ipsi-contra normalised UP MUA'

colour_lines = [0,90,50;74,20,134]/256; % Green Purple


% nexttile
% duration = event_times(:,2) - event_times(:,1);
% [~,sorted_index] = sort(duration);
% imagesc(event_averaging_scale*movmean(ipsi_probability(sorted_index,:),event_averaging_scale,1,'omitnan'))
% % imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
% xticks([0.5 4.5 8.5 12.5 16.5 20.5])
% % xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
% xticklabels([0 0.2 0.4 0.6 0.8 1])
% xline(50.5,'r',LineWidth=1)
% clim([0 1])
% colorbar
% colormap(flipud(gray))
% xlabel('Normalised duration of UP')
% ylabel('Event sorted by UP duration')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% title('ipsi ripples')
% 
% nexttile
% duration = event_times(:,2) - event_times(:,1);
% [~,sorted_index] = sort(duration);
% imagesc(event_averaging_scale*movmean(contra_probability(sorted_index,:),event_averaging_scale,1,'omitnan'))
% % imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
% xticks([0.5 4.5 8.5 12.5 16.5 20.5])
% % xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
% xticklabels([0 0.2 0.4 0.6 0.8 1])
% xline(50.5,'r',LineWidth=1)
% clim([0 1])
% colorbar
% colormap(flipud(gray))
% xlabel('Normalised duration of UP')
% ylabel('Event sorted by UP duration')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% title('contra ripples')


clear ERROR_SHADE
nexttile
binnedArray = MUA_PSTH_merged_thresholded.ipsi_UP_V1{end}{1};
% nprobe = 1;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
% xline(0,'r',LineWidth=1)


binnedArray = MUA_PSTH_merged_thresholded.contra_UP_V1{end}{1};

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
% xline(0,'r',LineWidth=1)

% baseline
binnedArray = MUA_PSTH_merged_thresholded.ipsi_UP_V1_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

ylim([-0.2 0.2])


% xline(0,'r')
legend([ERROR_SHADE(1:end)],{'ipsi','contra','shuffled'},'box','off')
title('V1 MUA')
xlabel('Normalised duration of UP')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
binnedArray = MUA_PSTH_merged_thresholded.ipsi_UP_HPC{end}{1};

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
% xline(0,'r',LineWidth=1)

binnedArray = MUA_PSTH_merged_thresholded.contra_UP_HPC{end}{1};

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
% xline(0,'r',LineWidth=1)

% baseline
binnedArray = MUA_PSTH_merged_thresholded.ipsi_UP_HPC_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

ylim([-0.2 0.2])

% xline(0,'r')
legend([ERROR_SHADE(1:end)],{'ipsi','contra','shuffled'},'box','off')
title('HPC MUA')
xlabel('Normalised duration of UP')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')






%% Ipsi and contra MUA during normalised DOWN with different bilateral synchrony

%%%%% DOWN info
event_times = merged_event_info.DOWN_ints;
hemisphere_id = merged_event_info.DOWN_hemisphere_id;
% lag_diff = merged_event_info.UP_lag_diff;
% group_id = merged_event_info.UP_group_id;
all_overlap_idx = merged_event_info.DOWN_overlap_idx_all{end};
non_overlap_idx = merged_event_info.DOWN_non_overlap_idx{end};
all_lags =merged_event_info.DOWN_lags_all{end};



event_idx = [];

lags =all_lags(all_lags>=-0.15 & all_lags<=0.15);
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lags =all_lags;

% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    event_idx{1}{n} =(all_overlap_idx(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end


lags =abs(all_lags(abs(all_lags)>=0 & abs(all_lags)<=0.15));

% lags =all_lags(all_lags>=0 & all_lags<=0.15);
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lags =abs(all_lags);

for n = 1:length(lag_thresholds)-1
    event_idx{2}{n} =(all_overlap_idx(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end


lags =all_lags(all_lags>=-0.05 & all_lags<=0.05);
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lags =all_lags;

% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    event_idx{3}{n} =(all_overlap_idx(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end


% lags =all_lags(all_lags>=-0.05 & all_lags<=0.05);
power_thresholds = prctile(SWpeakmag_UD(all_overlap_idx(all_lags>=-0.15 & all_lags<=0.15)),[0 20 40 60 80 100]);
delta_power =SWpeakmag_UD;
% power_thresholds = prctile(Delta_peaks_zscore_DU(all_lags>=-0.15 & all_lags<=0.15),[0 20 40 60 80 100]);
% delta_power =Delta_peaks_zscore_DU(all_lags>=-0.15 & all_lags<=0.15);

% lag_thresholds = [-0.2]
for n = 1:length(power_thresholds)-1
    event_idx{4}{n} =intersect(all_overlap_idx(abs(all_lags)<=0.15),find(delta_power>power_thresholds(n)&delta_power <power_thresholds(n+1)));
end

power_thresholds = prctile(SWpeakmag_UD,[0 20 40 60 80 100]);
delta_power =SWpeakmag_UD;
% power_thresholds = prctile(Delta_peaks_zscore_DU,[0 20 40 60 80 100]);
% delta_power =Delta_peaks_zscore_DU;

% lag_thresholds = [-0.2]
for n = 1:length(power_thresholds)-1
    event_idx{5}{n} =(delta_power>power_thresholds(n)&delta_power <power_thresholds(n+1));
end


Duration = [event_info(1).next_DOWN_duration; event_info(2).next_DOWN_duration];
duration_thresholds = prctile(Duration,[0 20 40 60 80 100]);

for n = 1:length(duration_thresholds)-1
    event_idx{6}{n} =find(Duration>duration_thresholds(n)&Duration <duration_thresholds(n+1));
end


Duration = [event_info(1).next_DOWN_duration; event_info(2).next_DOWN_duration];
duration_thresholds = prctile(Duration(all_overlap_idx(all_lags>=-0.15 & all_lags<=0.15)),[0 20 40 60 80 100]);

for n = 1:length(duration_thresholds)-1
    event_idx{7}{n} =intersect(all_overlap_idx(abs(all_lags)<=0.15),find(Duration>duration_thresholds(n)&Duration <duration_thresholds(n+1)));
end

group_name{1} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};
group_name{2} = {'0-20% lagging','20-40% lagging','40-60% lagging','60-80% lagging','80-100% lagging','Shuffled'};
group_name{3} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};
group_name{4} = {'0-20% delta power','20-40% delta power','40-60% delta power','60-80% delta power','80-100% delta power','Shuffled'};
group_name{5} = {'0-20% delta power','20-40% delta power','40-60% delta power','60-80% delta power','80-100% delta power','Shuffled'};
group_name{6} = {'0-20% duration','20-40% duration','40-60% duration','60-80% duration','80-100% duration','Shuffled'};
group_name{7} = {'0-20% duration','20-40% duration','40-60% duration','60-80% duration','80-100% duration','Shuffled'};



title_names{1} = 'Ipsi-contra normalised DOWN MUA by 5 lags (150ms windows)';
title_names{2} = 'Ipsi-contra normalised DOWN MUA by 5 lags (150ms windows abs lag)';
title_names{3} = 'Ipsi-contra normalised DOWN MUA by 5 lags (50ms windows)';
title_names{4} = 'Ipsi-contra normalised DOWN MUA by 5 powers (150ms windows)';
title_names{5} = 'Ipsi-contra normalised DOWN MUA by 5 powers (All events)';
title_names{6} = 'Ipsi-contra normalised DOWN MUA by 5 duration (All events)';
title_names{7} = 'Ipsi-contra normalised DOWN MUA by 5 duration (150ms events)';

%%%%%%%%%%%%%%%%%%%%%
ipsi_V1_MUA = [PSTH_MUA(1).L_V1_DOWN; PSTH_MUA(2).R_V1_DOWN];
contra_V1_MUA = [PSTH_MUA(1).R_V1_DOWN; PSTH_MUA(2).L_V1_DOWN];

ipsi_V1_MUA_baseline = [PSTH_MUA_baseline(1).L_V1_DOWN; PSTH_MUA_baseline(2).R_V1_DOWN];
contra_V1_MUA_baseline = [PSTH_MUA_baseline(1).R_V1_DOWN; PSTH_MUA_baseline(2).L_V1_DOWN];


ipsi_HPC_MUA = [PSTH_MUA(1).L_HPC_DOWN; PSTH_MUA(2).R_HPC_DOWN];
contra_HPC_MUA = [PSTH_MUA(1).R_HPC_DOWN; PSTH_MUA(2).L_HPC_DOWN];
ipsi_HPC_MUA_baseline = [PSTH_MUA_baseline(1).L_HPC_DOWN; PSTH_MUA_baseline(2).R_HPC_DOWN];
contra_HPC_MUA_baseline = [PSTH_MUA_baseline(1).R_HPC_DOWN; PSTH_MUA_baseline(2).L_HPC_DOWN];
%%%%%%%%%%%%%%%%%%%%%%%%%%% V1
%%%%% calculate shuffled MUA baseline
%%% Ipsi
binnedArray = ipsi_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap = temp;

%%% Contra
binnedArray = contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

contra_baseline_bootstrap = temp;

%%% Ipsi-Contra baseline
binnedArray = ipsi_V1_MUA_baseline-contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap = temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HPC
%%%%% calculate shuffled MUA baseline
%%% Ipsi
binnedArray = ipsi_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap1 = temp;

%%% Contra
binnedArray = contra_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

contra_baseline_bootstrap1 = temp;



%%% Ipsi-Contra baseline
binnedArray = ipsi_HPC_MUA_baseline-contra_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap1 = temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Calculate bootstrapped MUA

% clear probability_merged
num_bins=20;
x = linspace(0,1,num_bins);

MUA_PSTH_merged_thresholded.x = x;
MUA_PSTH_merged_thresholded.ipsi_DOWN_V1 = [];
MUA_PSTH_merged_thresholded.contra_DOWN_V1 = [];
MUA_PSTH_merged_thresholded.ipsi_contra_diff_DOWN_V1 = [];

MUA_PSTH_merged_thresholded.ipsi_DOWN_HPC = [];
MUA_PSTH_merged_thresholded.contra_DOWN_HPC = [];
MUA_PSTH_merged_thresholded.ipsi_contra_diff_DOWN_HPC = [];


for ngroup = 1:length(event_idx)+1

    group_index = [];
    if ngroup <= length(event_idx)
        group_index = event_idx{ngroup};

    else
        group_index{1} = (1:size(ipsi_V1_MUA,1))';
    end

    for i = 1:length(group_index)
        index =group_index{i};

        binnedArray1 = ipsi_V1_MUA(index,:);
        binnedArray2 = contra_V1_MUA(index,:);
        binnedArray3 = ipsi_V1_MUA(index,:)-contra_V1_MUA(index,:);
        binnedArray4 = ipsi_HPC_MUA(index,:);
        binnedArray5 = contra_HPC_MUA(index,:);
        binnedArray6 = ipsi_HPC_MUA(index,:)-contra_HPC_MUA(index,:);

        temp1=[];
        temp2=[];
        temp3=[];
        temp4=[];
        temp5=[];
        temp6=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray1,1),size(binnedArray1,1));
            temp1(iBoot,:) =  mean(binnedArray1(event_id,:),'omitnan');
            temp2(iBoot,:) =  mean(binnedArray2(event_id,:),'omitnan');
            temp3(iBoot,:) =  mean(binnedArray3(event_id,:),'omitnan');
            temp4(iBoot,:) =  mean(binnedArray4(event_id,:),'omitnan');
            temp5(iBoot,:) =  mean(binnedArray5(event_id,:),'omitnan');
            temp6(iBoot,:) =  mean(binnedArray6(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        MUA_PSTH_merged_thresholded.ipsi_DOWN_V1{ngroup}{i} = temp1;
        MUA_PSTH_merged_thresholded.contra_DOWN_V1{ngroup}{i} = temp2;
        MUA_PSTH_merged_thresholded.ipsi_contra_diff_DOWN_V1{ngroup}{i} = temp3;
        MUA_PSTH_merged_thresholded.ipsi_DOWN_HPC{ngroup}{i} = temp4;
        MUA_PSTH_merged_thresholded.contra_DOWN_HPC{ngroup}{i} = temp5;
        MUA_PSTH_merged_thresholded.ipsi_contra_diff_DOWN_HPC{ngroup}{i} = temp6;
    end
end

MUA_PSTH_merged_thresholded.ipsi_DOWN_V1_baseline = ipsi_baseline_bootstrap;
MUA_PSTH_merged_thresholded.contra_DOWN_V1_baseline = contra_baseline_bootstrap;
MUA_PSTH_merged_thresholded.ipsi_contra_diff_DOWN_V1_baseline = ipsi_contra_diff_baseline_bootstrap;
MUA_PSTH_merged_thresholded.ipsi_DOWN_HPC_baseline = ipsi_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.contra_DOWN_HPC_baseline = contra_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.ipsi_contra_diff_DOWN_HPC_baseline = ipsi_contra_diff_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.DOWN_groups = [title_names 'All DOWN'];
MUA_PSTH_merged_thresholded.DOWN_index = [event_idx {(1:size(ipsi_V1_MUA,1))'}];


% time_wondows = [-1 1];
% time_bin = 0.02;
num_bins=20;
x = linspace(0,1,num_bins);

%%%%%%%%% Plot DU transition
for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 930];
    fig.Name =title_names{ngroup};

    if ngroup ==1  | ngroup ==3
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    elseif ngroup==2
        colour_lines = [228,42,168;222,119,174;184,225,134;127,188,65;0,90,50]/256; % From magenta to dark purple to indicate globally synchrnous to ipsi lagging
    end


    if ngroup >=3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end
    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_DOWN_V1{ngroup}{i};
        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_DOWN_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-0.3 0.6])

    title('V1 Ipsi MUA')
    xlabel('Normalised duration of DOWN')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    if ngroup >=3
%         colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.contra_DOWN_V1{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.contra_DOWN_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);
    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.3 0.6])

    % xline(0,'r')
    title('V1 Contra MUA')
    xlabel('Normalised duration of DOWN')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    if ngroup >=3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_DOWN_V1{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_DOWN_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.3 0.6])


    % xline(0,'r')
    title('V1 Ipsi-contra MUA diff')
    xlabel('Normalised duration of DOWN')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%%%%% HPC MUA
    if ngroup >=3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_DOWN_HPC{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_DOWN_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-0.3 0.3])
    title('HPC Ipsi MUA')
    xlabel('Normalised duration of DOWN')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    if ngroup >=3
%         colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.contra_DOWN_HPC{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.contra_UP_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);


    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.3 0.3])

    % xline(0,'r')
    title('HPC Contra MUA')
    xlabel('Normalised duration of DOWN')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    if ngroup >=3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_DOWN_HPC{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_DOWN_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);


    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.3 0.3])


    % xline(0,'r')
    title('HPC Ipsi-contra MUA diff')
    xlabel('Normalised duration of DOWN')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

end

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')



%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%% All DOWN UP rpples
event_averaging_scale = 10;
x = linspace(0,1,20);

% for ngroup = 1:length(event_idx)
ngroup = 2;
fig = figure('Color','w');
fig.Position = [350 59 1650/3*2 465];
fig.Name = 'Left-Right combined Ipsi-contra normalised DOWN MUA'

colour_lines = [0,90,50;74,20,134]/256; % Green Purple


% nexttile
% duration = event_times(:,2) - event_times(:,1);
% [~,sorted_index] = sort(duration);
% imagesc(event_averaging_scale*movmean(ipsi_probability(sorted_index,:),event_averaging_scale,1,'omitnan'))
% % imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
% xticks([0.5 4.5 8.5 12.5 16.5 20.5])
% % xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
% xticklabels([0 0.2 0.4 0.6 0.8 1])
% xline(50.5,'r',LineWidth=1)
% clim([0 1])
% colorbar
% colormap(flipud(gray))
% xlabel('Normalised duration of UP')
% ylabel('Event sorted by UP duration')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% title('ipsi ripples')
% 
% nexttile
% duration = event_times(:,2) - event_times(:,1);
% [~,sorted_index] = sort(duration);
% imagesc(event_averaging_scale*movmean(contra_probability(sorted_index,:),event_averaging_scale,1,'omitnan'))
% % imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
% xticks([0.5 4.5 8.5 12.5 16.5 20.5])
% % xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
% xticklabels([0 0.2 0.4 0.6 0.8 1])
% xline(50.5,'r',LineWidth=1)
% clim([0 1])
% colorbar
% colormap(flipud(gray))
% xlabel('Normalised duration of UP')
% ylabel('Event sorted by UP duration')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% title('contra ripples')


clear ERROR_SHADE
nexttile
binnedArray = MUA_PSTH_merged_thresholded.ipsi_DOWN_V1{end}{1};
% nprobe = 1;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
% xline(0,'r',LineWidth=1)


binnedArray = MUA_PSTH_merged_thresholded.contra_DOWN_V1{end}{1};

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
% xline(0,'r',LineWidth=1)

% baseline
binnedArray = MUA_PSTH_merged_thresholded.ipsi_DOWN_V1_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

ylim([-0.2 0.6])


% xline(0,'r')
legend([ERROR_SHADE(1:end)],{'ipsi','contra','shuffled'},'box','off')
title('V1 MUA')
xlabel('Normalised duration of DOWN')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
binnedArray = MUA_PSTH_merged_thresholded.ipsi_DOWN_HPC{end}{1};

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
% xline(0,'r',LineWidth=1)

binnedArray = MUA_PSTH_merged_thresholded.contra_DOWN_HPC{end}{1};

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
% xline(0,'r',LineWidth=1)

% baseline
binnedArray = MUA_PSTH_merged_thresholded.ipsi_DOWN_HPC_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

ylim([-0.1 0.2])

% xline(0,'r')
legend([ERROR_SHADE(1:end)],{'ipsi','contra','shuffled'},'box','off')
title('HPC MUA')
xlabel('Normalised duration of DOWN')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')



%%%%%%%%%% previous UP

title_names{1} = 'Ipsi-contra normalised previous UP MUA by 5 lags (150ms windows)';
title_names{2} = 'Ipsi-contra normalised previous UP MUA by 5 lags (150ms windows abs lag)';
title_names{3} = 'Ipsi-contra normalised previous UP MUA by 5 lags (50ms windows)';
title_names{4} = 'Ipsi-contra normalised previous UP MUA by 5 powers (150ms windows)';
title_names{5} = 'Ipsi-contra normalised previous UP MUA by 5 powers (All events)';
title_names{6} = 'Ipsi-contra normalised previous UP MUA by 5 duration (All events)';
title_names{7} = 'Ipsi-contra normalised previous UP MUA by 5 duration (150ms events)';

ipsi_V1_MUA = [PSTH_MUA(1).L_V1_UP; PSTH_MUA(2).R_V1_UP];
contra_V1_MUA = [PSTH_MUA(1).R_V1_UP; PSTH_MUA(2).L_V1_UP];

ipsi_V1_MUA_baseline = [PSTH_MUA_baseline(1).L_V1_UP; PSTH_MUA_baseline(2).R_V1_UP];
contra_V1_MUA_baseline = [PSTH_MUA_baseline(1).R_V1_UP; PSTH_MUA_baseline(2).L_V1_UP];


ipsi_HPC_MUA = [PSTH_MUA(1).L_HPC_UP; PSTH_MUA(2).R_HPC_UP];
contra_HPC_MUA = [PSTH_MUA(1).R_HPC_UP; PSTH_MUA(2).L_HPC_UP];
ipsi_HPC_MUA_baseline = [PSTH_MUA_baseline(1).L_HPC_UP; PSTH_MUA_baseline(2).R_HPC_UP];
contra_HPC_MUA_baseline = [PSTH_MUA_baseline(1).R_HPC_UP; PSTH_MUA_baseline(2).L_HPC_UP];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Calculate bootstrapped MUA

% clear probability_merged
% time_wondows = [-1 1];
% time_bin = 0.02;
num_bins=20;
x = linspace(0,1,num_bins);
MUA_PSTH_merged_thresholded.x = x;
MUA_PSTH_merged_thresholded.ipsi_previous_UP_V1 = [];
MUA_PSTH_merged_thresholded.contra_previous_UP_V1 = [];
MUA_PSTH_merged_thresholded.ipsi_contra_diff_previous_UP_V1 = [];

MUA_PSTH_merged_thresholded.ipsi_previous_UP_HPC = [];
MUA_PSTH_merged_thresholded.contra_previous_UP_HPC = [];
MUA_PSTH_merged_thresholded.ipsi_contra_diff_previous_UP_HPC = [];


for ngroup = 1:length(event_idx)

    group_index = [];
    if ngroup <= length(event_idx)
        group_index = event_idx{ngroup};

    else
        group_index{1} = (1:size(ipsi_V1_MUA,1))';
    end

    for i = 1:length(group_index)
        index =group_index{i};

        binnedArray1 = ipsi_V1_MUA(index,:);
        binnedArray2 = contra_V1_MUA(index,:);
        binnedArray3 = ipsi_V1_MUA(index,:)-contra_V1_MUA(index,:);
        binnedArray4 = ipsi_HPC_MUA(index,:);
        binnedArray5 = contra_HPC_MUA(index,:);
        binnedArray6 = ipsi_HPC_MUA(index,:)-contra_HPC_MUA(index,:);

        temp1=[];
        temp2=[];
        temp3=[];
        temp4=[];
        temp5=[];
        temp6=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray1,1),size(binnedArray1,1));
            temp1(iBoot,:) =  mean(binnedArray1(event_id,:),'omitnan');
            temp2(iBoot,:) =  mean(binnedArray2(event_id,:),'omitnan');
            temp3(iBoot,:) =  mean(binnedArray3(event_id,:),'omitnan');
            temp4(iBoot,:) =  mean(binnedArray4(event_id,:),'omitnan');
            temp5(iBoot,:) =  mean(binnedArray5(event_id,:),'omitnan');
            temp6(iBoot,:) =  mean(binnedArray6(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        MUA_PSTH_merged_thresholded.ipsi_previous_UP_V1{ngroup}{i} = temp1;
        MUA_PSTH_merged_thresholded.contra_previous_UP_V1{ngroup}{i} = temp2;
        MUA_PSTH_merged_thresholded.ipsi_contra_diff_previous_UP_V1{ngroup}{i} = temp3;
        MUA_PSTH_merged_thresholded.ipsi_previous_UP_HPC{ngroup}{i} = temp4;
        MUA_PSTH_merged_thresholded.contra_previous_UP_HPC{ngroup}{i} = temp5;
        MUA_PSTH_merged_thresholded.ipsi_contra_diff_previous_UP_HPC{ngroup}{i} = temp6;
    end
end

MUA_PSTH_merged_thresholded.previous_UP_groups = [title_names];
MUA_PSTH_merged_thresholded.previous_UP_index = [event_idx];



% time_wondows = [-1 1];
% time_bin = 0.02;
num_bins=20;
x = linspace(0,1,num_bins);

%%%%%%%%% Plot DU transition
for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 930];
    fig.Name =title_names{ngroup};

    if ngroup ==1  | ngroup ==3
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    elseif ngroup==2
        colour_lines = [228,42,168;222,119,174;184,225,134;127,188,65;0,90,50]/256; % From magenta to dark purple to indicate globally synchrnous to ipsi lagging
    end


    if ngroup >=3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end
    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_previous_UP_V1{ngroup}{i};
        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_UP_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-0.3 0.4])

    title('V1 Ipsi MUA')
    xlabel('Normalised duration of UP')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    if ngroup >=3
%         colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.contra_previous_UP_V1{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.contra_UP_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);
    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.3 0.4])

    % xline(0,'r')
    title('V1 Contra MUA')
    xlabel('Normalised duration of UP')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    if ngroup >=3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_previous_UP_V1{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_UP_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.3 0.4])


    % xline(0,'r')
    title('V1 Ipsi-contra MUA diff')
    xlabel('Normalised duration of UP')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%%%%% HPC MUA
    if ngroup >=3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_previous_UP_HPC{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_UP_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-0.3 0.2])
    title('HPC Ipsi MUA')
    xlabel('Normalised duration of UP')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    if ngroup >=3
%         colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.contra_previous_UP_HPC{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.contra_UP_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);


    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.3 0.2])

    % xline(0,'r')
    title('HPC Contra MUA')
    xlabel('Normalised duration of UP')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    if ngroup >=3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_previous_UP_HPC{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_UP_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);


    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.3 0.2])


    % xline(0,'r')
    title('HPC Ipsi-contra MUA diff')
    xlabel('Normalised duration of UP')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')


save(fullfile(analysis_folder,'V1-HPC sleep interaction','MUA_PSTH_merged_thresholded_normalised.mat'),'MUA_PSTH_merged_thresholded');
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
end