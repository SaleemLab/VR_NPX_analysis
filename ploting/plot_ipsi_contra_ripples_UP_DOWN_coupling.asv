function plot_ipsi_contra_ripples_UP_DOWN_coupling(slow_waves_all,ripples_all,spindles_all)

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

load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_SO_probability_whole_baseline.mat'));
probability_psth_whole_baseline = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_SO_probability_whole.mat'));
probability_psth_whole = probability;


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


%%%% lag based on LFP
probability = probability_psth_whole;
nprobe = 1;
event_averaging_scale = 10;
% extract plv, amp corr and lag (latency)
for nprobe=1:2
    mprobe = abs(nprobe-3);

  
    % ripples
    ipsi_lag_ripples{nprobe} = [];
    contra_lag_ripples{nprobe} = [];

    ipsi_corr_ripples{nprobe} = [];
    contra_corr_ripples{nprobe} = [];

    ipsi_plv_ripples{nprobe} = [];
    contra_plv_ripples{nprobe} = [];

    for nsession = 1:max(ripples_all(1).session_count)

  
        %%%%% Ripples
        [C,ia,ib] = intersect(find(ripples_all(nprobe).session_count == sessions_to_process(nsession)),find(ripples_all(nprobe).session_count == sessions_to_process(nsession) & ripples_all(nprobe).SWS_index == 1));
        ipsi_shank = find(ripples_all(nprobe).probe_hemisphere{nsession} == nprobe);
        ipsi_shank(ipsi_shank==HPC_ref_shank(nsession,nprobe))=[];
        contra_shank = find(ripples_all(nprobe).probe_hemisphere{nsession} == mprobe);

        mean_corr_ipsi = [];mean_corr_contra=[];
        for nshank = 1:length(ipsi_shank)
            mean_corr_ipsi(nshank) = mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank(nshank),ia)));
        end

        for nshank = 1:length(contra_shank)
            mean_corr_contra(nshank)= mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank(nshank),ia)));
        end

        [~,id] = max(mean_corr_ipsi);
        ipsi_shank = ipsi_shank(id);
        [~,id] = max(mean_corr_contra);
        contra_shank = contra_shank(id);


        ipsi_lag_ripples{nprobe} = [ipsi_lag_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_lag_ripples{nprobe} = [contra_lag_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_plv_ripples{nprobe} = [ipsi_plv_ripples{nprobe} squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_plv_ripples{nprobe} = [contra_plv_ripples{nprobe} squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_corr_ripples{nprobe} = [ipsi_corr_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_corr_ripples{nprobe} = [contra_corr_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia))'];
    end
end


%% Grouping of DOWN-UP and UP-DOWN and ripples based on detetcion lags
load(fullfile(analysis_folder,'V1-HPC sleep interaction','k_cluster_ipsi_contra_events.mat'),'k_cluster');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');


%%%%%%%%%%%%% Overlaps
% Initialize output
L_overlap_idx = [];
R_overlap_idx = [];
overlap_idx = [];
non_overlap_idx = [];
merged_idx = [];
all_lags = [];
all_overlap_idx = [];
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


%% Ipsi and contra DOWN probabilities during ripples

%%%%% Ripples info
lag_diff = merged_event_info.ripples_lag_diff;
plv_diff = merged_event_info.ripples_plv_diff;
corr_diff = merged_event_info.ripples_corr_diff;
% ipsi_lag = [ipsi_lag_ripples{1} ipsi_lag_ripples{2}]';
contra_lag = [contra_lag_ripples{1} contra_lag_ripples{2}]';
group_id = merged_event_info.ripples_group_id;


%%%%%%% GET events not occuring within 300ms
% Store original indices
[sorted_times, sort_idx] = sort(merged_event_info.ripples_peaktimes);
original_idx = 1:length(merged_event_info.ripples_peaktimes);

dt = diff(sorted_times);

% Find where the difference exceeds 300 ms 
index_sorted = find(diff(sorted_times) > 0.3) + 1;

% Map sorted indices back to original indices
singlet_ripples_index = sort_idx(index_sorted);

%%%%%%% Event grouping based on ipsi-contra difference features
event_idx = [];
% Ipsi-contra ripples HPC MUA by 5 lags (full windows)
lags =all_lags{end};
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    event_idx{1}{n} =(all_overlap_idx{end}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

% Ipsi-contra ripples HPC MUA by 5 abs lags (full windows)
lags =abs(all_lags{end});
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
for n = 1:length(lag_thresholds)-1
    event_idx{2}{n} =(all_overlap_idx{end}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

% Ipsi-contra ripples HPC MUA by 5 lags (full windows singlets)
[~,ia,ib]=intersect(all_overlap_idx{end},singlet_ripples_index);
lags =all_lags{end}(ia);
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lags =all_lags{end}; 

% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    event_idx{3}{n} =intersect(singlet_ripples_index,(all_overlap_idx{end}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1))));
end

% Ipsi-contra ripples HPC MUA by 5 abs lags (full windows singlets)
[~,ia,ib]=intersect(all_overlap_idx{end},singlet_ripples_index);
lags =abs(all_lags{end}(ia));
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lags =abs(all_lags{end}); 

for n = 1:length(lag_thresholds)-1
    event_idx{4}{n} =intersect(singlet_ripples_index,(all_overlap_idx{end}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1))));
end


% Ipsi-contra ripples HPC MUA by 5 contra lags (full windows)
lags =contra_lag;
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    
    event_idx{5}{n} =find(lags>lag_thresholds(n)&lags <=lag_thresholds(n+1));
end

% Ipsi-contra ripples HPC MUA by 5 abs contra lags (full windows)
lags =abs(contra_lag);
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
for n = 1:length(lag_thresholds)-1
    event_idx{6}{n} =find(lags>=lag_thresholds(n)&lags <lag_thresholds(n+1));
end

% Ipsi-contra ripples HPC MUA by 5 contra lags (full windows singlets)

lags =contra_lag(singlet_ripples_index);
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lags =contra_lag;
% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    
    event_idx{7}{n} =intersect(singlet_ripples_index,find(lags>=lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

% Ipsi-contra ripples HPC MUA by 5 abs contra lags (full windows singlets)
lags =abs(contra_lag(singlet_ripples_index));
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lags =abs(contra_lag);
for n = 1:length(lag_thresholds)-1
    event_idx{8}{n} =intersect(singlet_ripples_index,find(lags>=lag_thresholds(n)&lags <lag_thresholds(n+1)));
end


% Ipsi-contra ripples HPC MUA by 5 powers
ripple_powers = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index); ripples_all(2).peak_zscore(ripples_all(2).SWS_index)];
power_thresholds = prctile(ripple_powers,[0 20 40 60 80 100]);
for n = 1:length(power_thresholds)-1
    event_idx{9}{n} =find(ripple_powers>=power_thresholds(n)&ripple_powers <power_thresholds(n+1));
end


% Ipsi-contra ripples HPC MUA by 5 powers (singlets)
ripple_powers = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index); ripples_all(2).peak_zscore(ripples_all(2).SWS_index)];
power_thresholds = prctile(ripple_powers(singlet_ripples_index),[0 20 40 60 80 100]);
for n = 1:length(power_thresholds)-1
    event_idx{10}{n} =intersect(singlet_ripples_index,find(ripple_powers>=power_thresholds(n)&ripple_powers <power_thresholds(n+1)));
end

% Ipsi-contra ripples HPC MUA (singlets)
event_idx{11}{1} = singlet_ripples_index;


group_name{1} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
group_name{2} = {'0-20% lagging','20-40% lagging','40-60% lagging','60-80% lagging','80-100% lagging','Shuffled'};
group_name{3} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
group_name{4} = {'0-20% lagging','20-40% lagging','40-60% lagging','60-80% lagging','80-100% lagging','Shuffled'};
group_name{5} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
group_name{6} = {'0-20% lagging','20-40% lagging','40-60% lagging','60-80% lagging','80-100% lagging','Shuffled'};
group_name{7} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
group_name{8} = {'0-20% lagging','20-40% lagging','40-60% lagging','60-80% lagging','80-100% lagging','Shuffled'};
group_name{9} = {'0-20% ripple power','20-40% ripple power','40-60% ripple power','60-80% ripple power','80-100% ripple power','Shuffled'};
group_name{10} = {'0-20% ripple power','20-40% ripple power','40-60% ripple power','60-80% ripple power','80-100% ripple power','Shuffled'};
group_name{11} = {'Singlets','Shuffled'};

title_names = {
    'Ipsi-contra ripples DOWN probability by 5 lags (full windows)', ...
    'Ipsi-contra ripples DOWN probability by 5 abs lags (full windows)', ...
    'Ipsi-contra ripples DOWN probability by 5 lags (full windows singlets)', ...
    'Ipsi-contra ripples DOWN probability by 5 abs lags (full windows singlets)', ...
    'Ipsi-contra ripples DOWN probability by 5 contra lags (full windows)', ...
    'Ipsi-contra ripples DOWN probability by 5 abs contra lags (full windows)', ...
    'Ipsi-contra ripples DOWN probability by 5 contra lags (full windows singlets)', ...
    'Ipsi-contra ripples DOWN probability by 5 abs contra lags (full windows singlets)', ...
    'Ipsi-contra ripples DOWN probability by 5 powers', ...
    'Ipsi-contra ripples DOWN probability by 5 powers (singlets)', ...
    'Left-Right combined ipsi-contra ripples DOWN probability (singlets)'
};
%%%%%%%%%%%%%%%%%%%%%
% ipsi_probability = [probability_psth_whole(1).L_ripples_UP; probability_psth_whole(2).R_ripples_UP];
% contra_probability = [probability_psth_whole(2).L_ripples_UP; probability_psth_whole(1).R_ripples_UP];

ipsi_probability = [probability_psth_whole(1).L_ripples_DOWN; probability_psth_whole(2).R_ripples_DOWN];
contra_probability = [probability_psth_whole(2).L_ripples_DOWN; probability_psth_whole(1).R_ripples_DOWN];

ipsi_probability_baseline = [probability_psth_whole_baseline(1).L_ripples_DOWN; probability_psth_whole_baseline(2).R_ripples_DOWN];
contra_probability_baseline = [probability_psth_whole_baseline(2).L_ripples_DOWN; probability_psth_whole_baseline(1).R_ripples_DOWN];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Calculate bootstrapped MUA

% clear probability_merged
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

probability_merged.x = x;
probability_merged.ipsi_DOWN_ripples = [];
probability_merged.contra_DOWN_ripples = [];
probability_merged.ipsi_contra_diff_DOWN_ripples = [];

ipsi_baseline_bootstrap=[];
contra_baseline_bootstrap=[];
ipsi_contra_diff_baseline_bootstrap=[];

for ngroup = 1:length(event_idx)+1

    group_index = [];
    if ngroup <= length(event_idx)
        group_index = event_idx{ngroup};

    else
        group_index{1} = (1:size(ipsi_probability,1))';
    end

    for i = 1:length(group_index)
        index =group_index{i};

        binnedArray1 = ipsi_probability(index,:);
        binnedArray2 = contra_probability(index,:);
        binnedArray3 = binnedArray1-binnedArray2;

        temp1=[];
        temp2=[];
        temp3=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray1,1),size(binnedArray1,1));
            temp1(iBoot,:) =  mean(binnedArray1(event_id,:),'omitnan');
            temp2(iBoot,:) =  mean(binnedArray2(event_id,:),'omitnan');
            temp3(iBoot,:) =  mean(binnedArray3(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        probability_merged.ipsi_DOWN_ripples{ngroup}{i} = temp1;
        probability_merged.contra_DOWN_ripples{ngroup}{i} = temp2;
        probability_merged.ipsi_contra_diff_DOWN_ripples{ngroup}{i} = temp3;
    end



    %%%%% Shuffled distribution
    mean_number = [];
    all_index = [];
    for i = 1:length(group_index)
        if size({ngroup},1) == 1
            all_index = [all_index reshape(group_index{i},1,[])];
        else
            all_index = [all_index group_index{i}];
        end

        mean_number(i) = length(group_index{i});
    end
    mean_number = round(mean(mean_number));

    %%%%% V1 Shuffled distribution
    binnedArray = ipsi_probability_baseline;
    temp=[];
    parfor iBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,all_index,mean_number);
        temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
        % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
    end

    ipsi_baseline_bootstrap{ngroup} = temp;

    %%% Contra
    binnedArray = contra_probability_baseline;
    temp=[];
    parfor iBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,all_index,mean_number);
        temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
        % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
    end

    contra_baseline_bootstrap{ngroup} = temp;

    %%% Ipsi-Contra baseline
    binnedArray = ipsi_probability_baseline-contra_probability_baseline;
    temp=[];
    parfor iBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
        temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
        % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
    end

    ipsi_contra_diff_baseline_bootstrap{ngroup} = temp;
end


probability_merged.ipsi_DOWN_baseline_ripples = ipsi_baseline_bootstrap;
probability_merged.contra_DOWN_baseline_ripples = contra_baseline_bootstrap;
probability_merged.ipsi_contra_diff_DOWN_baseline_ripples = ipsi_contra_diff_baseline_bootstrap;
probability_merged.DOWN_ripples_groups = [title_names 'All ripples-DOWN'];
probability_merged.DOWN_ripples_index = [event_idx {(1:size(ipsi_probability,1))'}];


time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

%%%%%%%%% Plot DU transition
for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 465];
    fig.Name =title_names{ngroup};

    % if ngroup ==1
    %     colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
    % elseif ngroup ==5
    %     colour_lines = [0,90,50;228,42,168;74,20,134;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
    % else
    if ismember(ngroup, [1, 3, 5, 7])
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    else
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = probability_merged.ipsi_DOWN_ripples{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = probability_merged.ipsi_DOWN_baseline_ripples{ngroup};
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([0 0.25])
    title('ipsi ripples')
    xlabel('Time relative to ripple peaks (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    if ismember(ngroup, [1, 3, 5, 7])
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    else
        % colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = probability_merged.contra_DOWN_ripples{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = probability_merged.contra_DOWN_baseline_ripples{ngroup};
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);
    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([0 0.25])

    % xline(0,'r')
    title('contra ripples')
    xlabel('Time relative to ripple peaks (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    if ismember(ngroup, [1, 3, 5, 7])
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    else
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = probability_merged.ipsi_contra_diff_DOWN_ripples{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = probability_merged.ipsi_contra_diff_DOWN_baseline_ripples{ngroup};
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.1 0.1])


    % xline(0,'r')
    title('Ipsi-contra diff')
    xlabel('Time relative to ripple peaks (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')

%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%% All DOWN UP rpples
event_averaging_scale = 15;

% for ngroup = 1:length(event_idx)
ngroup = 2;
fig = figure('Color','w');
fig.Position = [350 59 1650 465];
fig.Name = 'Left-Right combined ipsi contra DOWN distribution around ripple peaks'

colour_lines = [0,90,50;74,20,134]/256; % Green Purple


nexttile
event_times = merged_event_info.ripples_ints;
duration = event_times(:,2) - event_times(:,1);
[~,sorted_index] = sort(duration);
imagesc(event_averaging_scale*movmean(ipsi_probability(sorted_index,:),event_averaging_scale,1,'omitnan'))
% imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
xticks([1.5 25.5 50.5 75.5 100.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(50.5,'r',LineWidth=1)
clim([0 1])
colorbar
colormap(flipud(gray))
xlabel('Time relative to ripple peaks (s)')
ylabel('Event sorted by ripple duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('ipsi ripples')

nexttile
duration = event_times(:,2) - event_times(:,1);
[~,sorted_index] = sort(duration);
imagesc(event_averaging_scale*movmean(contra_probability(sorted_index,:),event_averaging_scale,1,'omitnan'))
% imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
xticks([1.5 25.5 50.5 75.5 100.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(50.5,'r',LineWidth=1)
clim([0 1])
colorbar
colormap(flipud(gray))
xlabel('Time relative to ripple peaks (s)')
ylabel('Event sorted by ripple duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('contra ripples')

nexttile
clear ERROR_SHADE

binnedArray = probability_merged.ipsi_DOWN_ripples{end}{1};
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

binnedArray = probability_merged.contra_DOWN_ripples{end}{1};
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline
binnedArray = probability_merged.ipsi_DOWN_baseline_ripples{end};
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

% xline(0,'r')
ylim([0 0.25])
% title('ipsi ripples')
xlabel('Time relative to ripple peaks (s)')
ylabel('Probability')
legend([ERROR_SHADE(1:end)],{'ipsi','contra','shuffled'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')



save(fullfile(analysis_folder,'V1-HPC bilateral interaction','SO_ripples_spindles_probability_merged.mat'),'probability_merged')

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