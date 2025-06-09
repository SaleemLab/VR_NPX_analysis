function plot_ipsi_contra_UP_DOWN_ripples_coupling_detection_thresholded(slow_waves_all,ripples_all,spindles_all)

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


% 
% 
% probability = probability_psth_whole;
%% Compare LFP based and detection based lags
colour_lines = [0,90,50;228,42,168;74,20,134;81 81 81]/256; % Dark Green, Magenta, dark purple
scatter(merged_event_info.UP_lags_all{end},merged_event_info.UP_lag_diff(merged_event_info.UP_overlap_idx_all{end}),5,'MarkerFaceColor',colour_lines(2, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.01)

fig = figure('Color','w');
fig.Position = [30 60 1200 900];
fig.Name = 'All events overlapping vs non-overlapping events lag corr plv distribution';

for n = 1:4
    nexttile
    histogram(merged_event_info.(sprintf('%s_lag_diff',event_types{n}))(merged_event_info.(sprintf('%s_overlap_idx_all',event_types{n})){end}),-0.2:0.005:0.2,'FaceColor',colour_lines(2, :),'EdgeColor','none','Normalization','probability');hold on;
    histogram(merged_event_info.(sprintf('%s_lag_diff',event_types{n}))(merged_event_info.(sprintf('%s_non_overlap_idx',event_types{n})){end}),-0.2:0.005:0.2,'FaceColor',colour_lines(4, :),'EdgeColor','none','Normalization','probability');hold on;
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    title(sprintf('%s lag diff',event_types{n}))

    nexttile
    histogram(merged_event_info.(sprintf('%s_corr_diff',event_types{n}))(merged_event_info.(sprintf('%s_overlap_idx_all',event_types{n})){end}),-2:0.02:2,'FaceColor',colour_lines(2, :),'EdgeColor','none','Normalization','probability');hold on;
    histogram(merged_event_info.(sprintf('%s_corr_diff',event_types{n}))(merged_event_info.(sprintf('%s_non_overlap_idx',event_types{n})){end}),-2:0.02:2,'FaceColor',colour_lines(4, :),'EdgeColor','none','Normalization','probability');hold on;
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    title(sprintf('%s corr diff',event_types{n}))


    nexttile
    histogram(merged_event_info.(sprintf('%s_plv_diff',event_types{n}))(merged_event_info.(sprintf('%s_overlap_idx_all',event_types{n})){end}),-1:0.01:1,'FaceColor',colour_lines(2, :),'EdgeColor','none','Normalization','probability');hold on;
    histogram(merged_event_info.(sprintf('%s_plv_diff',event_types{n}))(merged_event_info.(sprintf('%s_non_overlap_idx',event_types{n})){end}),-1:0.01:1,'FaceColor',colour_lines(4, :),'EdgeColor','none','Normalization','probability');hold on;
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    title(sprintf('%s plv diff',event_types{n}))

end
legend('Overlapping','Non-overlapping','box','off')
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')


% 
% merged_event_info.ripples_lags_all{end}
% merged_event_info.ripples_lag_diff(merged_event_info.ripples_overlap_idx_all{end})


%% Ripple probabilities during unilaterally biased and bilaterally synchronised UP
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
probability_normalised_whole = probability;

% Calculate ripples probability according to unilateral vs bilateral events
%%%%% Ripple info
ripples_times =  merged_event_info.ripples_ints;
ripples_group_info = merged_event_info.ripples_group_id;
ripples_hemisphere_id = merged_event_info.ripples_hemisphere_id;
ripples_lag_diff = merged_event_info.ripples_lag_diff;
ripples_all_overlap_idx = merged_event_info.ripples_overlap_idx_all{end};
ripples_non_overlap_idx = merged_event_info.ripples_non_overlap_idx{end};
ripples_all_lags =merged_event_info.ripples_lags_all{end};


%%%%% UP info
event_times = merged_event_info.UP_ints;
hemisphere_id = merged_event_info.UP_hemisphere_id;
% lag_diff = merged_event_info.UP_lag_diff;
% group_id = merged_event_info.UP_group_id;
all_overlap_idx = merged_event_info.UP_overlap_idx_all{end};
non_overlap_idx = merged_event_info.UP_non_overlap_idx{end};
all_lags =merged_event_info.UP_lags_all{end};


binnedArray=[];
time_windows = [-1 1];
time_bin = 0.02;
timebin_edge = time_windows(1):time_bin:time_windows(end);
bins_centre = timebin_edge(1)+time_bin/2:time_bin:timebin_edge(end)-time_bin/2;
probability_merged = [];

% lag_thresholds = (abs(prctile(all_lags(all_lags>-0.2&all_lags<0.2),[0:10:100])));
% 
% lag_thresholds = (abs(prctile(ripples_all_lags,[0:10:100])));
% 
lag_thresholds = mean(abs(prctile(ripples_lag_diff,[40 60])));
% % 
% lag_thresholds = mean(abs(prctile(ripples_all_lags,[30])));
% 

probability_psth_whole

%%
%%%%%%%%%% Based on detection lag
probability_merged = [];

index=[];

for nprobe = 1:2
    index{1}{nprobe} = intersect(ripples_all_overlap_idx(ripples_all_lags<-0.005& ripples_all_lags>-0.02),find(ripples_hemisphere_id == nprobe))';% leading
    index{2}{nprobe} = intersect(ripples_all_overlap_idx(ripples_all_lags>0.005 & ripples_all_lags<0.02),find(ripples_hemisphere_id == nprobe))';% lagging
    index{3}{nprobe} = ripples_all_overlap_idx(ripples_all_lags<=0&ripples_all_lags>-0.005); % bilaterally synchronised
    index{4}{nprobe} = [intersect(ripples_all_overlap_idx(ripples_all_lags<-0.02 | ripples_all_lags>0.02),find(ripples_hemisphere_id == nprobe))' intersect(ripples_non_overlap_idx,find(ripples_hemisphere_id == nprobe))']; % unilateral
end
group_name = {'ipsi_leading_ripple','contra_leading_ripple','bilateral_ripple','unilateral_ripple'};

for ngroup = 1:length(index)
    tic
    event_times = merged_event_info.UP_ints;
    hemisphere_id = merged_event_info.UP_hemisphere_id;
    temp = [];
    for nprobe = 1:2
        ints = event_times(hemisphere_id==nprobe,:);
        for mprobe = 1:2
            [~,temp{nprobe}{mprobe},~] = calculate_event_probability(ripples_times(index{ngroup}{mprobe},:), ints(:,1), time_windows(1):time_bin:time_windows(end),0);
            timebin_edges_all = ints(:,1) + bins_centre;  % Absolute times of peri-event window
            for i = 1:size(ints,1)
                % Previous DOWN (skip if this is the first UP)
                if i > 1
                    prev_offset = ints(i-1,2);
                    % Find peri-time indices within the previous UP state
                    mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                    temp{nprobe}{mprobe}(i, mask_prev) = NaN;
                end

                % Next DOWN (skip if this is the last UP)
                if i < size(ints,1)
                    next_onset = ints(i+1,1);
                    % Find peri-time indices within the next UP state
                    mask_next = timebin_edges_all(i,:) >= next_onset;
                    temp{nprobe}{mprobe}(i, mask_next) = NaN;
                end
            end


            % probability_merged(nprobe).(sprintf('%s_index',group_name{ngroup})) = temp;
        end
    end
  
    % 1 - ipsi and 2 - contra
    probability_merged(1).(sprintf('%s_UP',group_name{ngroup})) = [temp{1}{1}; temp{2}{2}];
    probability_merged(2).(sprintf('%s_UP',group_name{ngroup})) = [temp{1}{2}; temp{2}{1}];


    % DOWN
    event_times = merged_event_info.DOWN_ints;
    hemisphere_id = merged_event_info.DOWN_hemisphere_id;
    temp = [];
    for nprobe = 1:2
        ints = event_times(hemisphere_id==nprobe,:);
        for mprobe = 1:2
            [~,temp{nprobe}{mprobe},~] = calculate_event_probability(ripples_times(index{ngroup}{mprobe},:), ints(:,1), time_windows(1):time_bin:time_windows(end),0);
            timebin_edges_all = ints(:,1) + bins_centre;  % Absolute times of peri-event window
            for i = 1:size(ints,1)
                % Previous DOWN (skip if this is the first UP)
                if i > 1
                    prev_offset = ints(i-1,2);
                    % Find peri-time indices within the previous UP state
                    mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                    temp{nprobe}{mprobe}(i, mask_prev) = NaN;
                end

                % Next DOWN (skip if this is the last UP)
                if i < size(ints,1)
                    next_onset = ints(i+1,1);
                    % Find peri-time indices within the next UP state
                    mask_next = timebin_edges_all(i,:) >= next_onset;
                    temp{nprobe}{mprobe}(i, mask_next) = NaN;
                end
            end


            % probability_merged(nprobe).(sprintf('%s_index',group_name{ngroup})) = temp;
        end
    end

    % 1 - ipsi and 2 - contra
    probability_merged(1).(sprintf('%s_DOWN',group_name{ngroup})) = [temp{1}{1}; temp{2}{2}];
    probability_merged(2).(sprintf('%s_DOWN',group_name{ngroup})) = [temp{1}{2}; temp{2}{1}];

    % bootstrap distribution
    for nprobe = 1:2
        binnedArrayUP = probability_merged(nprobe).(sprintf('%s_UP',group_name{ngroup}));
        binnedArrayDOWN = probability_merged(nprobe).(sprintf('%s_DOWN',group_name{ngroup}));
        tempUP = [];
        tempDOWN = [];

        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArrayUP,1),size(binnedArrayUP,1));
            tempUP(iBoot,:) = sum(binnedArrayUP(event_id,:),'omitnan')./sum(~isnan(binnedArrayUP(event_id,:)));

            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArrayDOWN,1),size(binnedArrayDOWN,1));
            tempDOWN(iBoot,:) = sum(binnedArrayDOWN(event_id,:),'omitnan')./sum(~isnan(binnedArrayDOWN(event_id,:)));
        end

        probability_merged(nprobe).(sprintf('%s_UP_bootstrap',group_name{ngroup})) = tempUP;
        probability_merged(nprobe).(sprintf('%s_DOWN_bootstrap',group_name{ngroup})) = tempDOWN;
    end
    toc
end



num_bins = 20;
probability_normalised_merged = [];
for ngroup = 1:length(index)
    tic
    event_times = merged_event_info.UP_ints;
    hemisphere_id = merged_event_info.UP_hemisphere_id;
    temp = [];
    for nprobe = 1:2
        ints = event_times(hemisphere_id==nprobe,:);
        for mprobe = 1:2
            [~,temp{nprobe}{mprobe},~] = calculate_relative_event_probability(ints,ripples_times(index{ngroup}{mprobe},:),num_bins,0);
            timebin_edges_all = ints(:,1) + bins_centre;  % Absolute times of peri-event window
            for i = 1:size(ints,1)
                % Previous DOWN (skip if this is the first UP)
                if i > 1
                    prev_offset = ints(i-1,2);
                    % Find peri-time indices within the previous UP state
                    mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                    temp{nprobe}{mprobe}(i, mask_prev) = NaN;
                end

                % Next DOWN (skip if this is the last UP)
                if i < size(ints,1)
                    next_onset = ints(i+1,1);
                    % Find peri-time indices within the next UP state
                    mask_next = timebin_edges_all(i,:) >= next_onset;
                    temp{nprobe}{mprobe}(i, mask_next) = NaN;
                end
            end


            % probability_merged(nprobe).(sprintf('%s_index',group_name{ngroup})) = temp;
        end
    end
  
    % 1 - ipsi and 2 - contra
    probability_normalised_merged(1).(sprintf('%s_UP',group_name{ngroup})) = [temp{1}{1}; temp{2}{2}];
    probability_normalised_merged(2).(sprintf('%s_UP',group_name{ngroup})) = [temp{1}{2}; temp{2}{1}];


    % DOWN
    event_times = merged_event_info.DOWN_ints;
    hemisphere_id = merged_event_info.DOWN_hemisphere_id;
    temp = [];
    for nprobe = 1:2
        ints = event_times(hemisphere_id==nprobe,:);
        for mprobe = 1:2

            [~,temp{nprobe}{mprobe},~] = calculate_relative_event_probability(ints,ripples_times(index{ngroup}{mprobe},:),num_bins,0);
            timebin_edges_all = ints(:,1) + bins_centre;  % Absolute times of peri-event window
            for i = 1:size(ints,1)
                % Previous DOWN (skip if this is the first UP)
                if i > 1
                    prev_offset = ints(i-1,2);
                    % Find peri-time indices within the previous UP state
                    mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                    temp{nprobe}{mprobe}(i, mask_prev) = NaN;
                end

                % Next DOWN (skip if this is the last UP)
                if i < size(ints,1)
                    next_onset = ints(i+1,1);
                    % Find peri-time indices within the next UP state
                    mask_next = timebin_edges_all(i,:) >= next_onset;
                    temp{nprobe}{mprobe}(i, mask_next) = NaN;
                end
            end
            % probability_merged(nprobe).(sprintf('%s_index',group_name{ngroup})) = temp;
        end
    end

    % 1 - ipsi and 2 - contra
    probability_normalised_merged(1).(sprintf('%s_DOWN',group_name{ngroup})) = [temp{1}{1}; temp{2}{2}];
    probability_normalised_merged(2).(sprintf('%s_DOWN',group_name{ngroup})) = [temp{1}{2}; temp{2}{1}];

    % bootstrap distribution
    for nprobe = 1:2
        binnedArrayUP = probability_normalised_merged(nprobe).(sprintf('%s_UP',group_name{ngroup}));
        binnedArrayDOWN = probability_normalised_merged(nprobe).(sprintf('%s_DOWN',group_name{ngroup}));
        tempUP = [];
        tempDOWN = [];

        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArrayUP,1),size(binnedArrayUP,1));
            tempUP(iBoot,:) = sum(binnedArrayUP(event_id,:),'omitnan')./sum(~isnan(binnedArrayUP(event_id,:)));

            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArrayDOWN,1),size(binnedArrayDOWN,1));
            tempDOWN(iBoot,:) = sum(binnedArrayDOWN(event_id,:),'omitnan')./sum(~isnan(binnedArrayDOWN(event_id,:)));
        end

        probability_normalised_merged(nprobe).(sprintf('%s_UP_bootstrap',group_name{ngroup})) = tempUP;
        probability_normalised_merged(nprobe).(sprintf('%s_DOWN_bootstrap',group_name{ngroup})) = tempDOWN;
    end
    toc
end


% end

% probability_merged = merged_probability_psth_whole;
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_merged.mat'),'probability_merged');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole_merged.mat'),'probability_normalised_merged');




%%%%%%%%% LFP based lag

index=[];
probability_merged = [];
for nprobe = 1:2
    lag_thresholds = (prctile(ripples_lag_diff,[12.5 37.5]));
    index{1}{nprobe} = intersect(find(ripples_lag_diff>lag_thresholds(1) & ripples_lag_diff<lag_thresholds(2)),find(ripples_hemisphere_id == nprobe))';% leading
    lag_thresholds = (prctile(ripples_lag_diff,[62.5 87.5]));
    index{2}{nprobe} = intersect(find(ripples_lag_diff>lag_thresholds(1) & ripples_lag_diff<lag_thresholds(2)),find(ripples_hemisphere_id == nprobe))';% leading
    lag_thresholds = (prctile(ripples_lag_diff,[37.5 62.5]));
    index{3}{nprobe} = find(ripples_lag_diff>lag_thresholds(1) & ripples_lag_diff<lag_thresholds(2)); % bilaterally synchronised
    lag_thresholds = (prctile(ripples_lag_diff,[12.5 87.5]));
    index{4}{nprobe} = intersect(find(ripples_lag_diff<lag_thresholds(1) | ripples_lag_diff>lag_thresholds(2)),find(ripples_hemisphere_id == nprobe))';% leading
end
group_name = {'ipsi_leading_ripple','contra_leading_ripple','bilateral_ripple','unilateral_ripple'};

for ngroup = 1:length(index)
    event_times = merged_event_info.UP_ints;
    hemisphere_id = merged_event_info.UP_hemisphere_id;
    temp = [];
    for nprobe = 1:2
        ints = event_times(hemisphere_id==nprobe,:);
        for mprobe = 1:2
            [~,temp{nprobe}{mprobe},~] = calculate_event_probability(ripples_times(index{ngroup}{mprobe},:), ints(:,1), time_windows(1):time_bin:time_windows(end),0);
            timebin_edges_all = ints(:,1) + bins_centre;  % Absolute times of peri-event window
            for i = 1:size(ints,1)
                % Previous DOWN (skip if this is the first UP)
                if i > 1
                    prev_offset = ints(i-1,2);
                    % Find peri-time indices within the previous UP state
                    mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                    temp{nprobe}{mprobe}(i, mask_prev) = NaN;
                end

                % Next DOWN (skip if this is the last UP)
                if i < size(ints,1)
                    next_onset = ints(i+1,1);
                    % Find peri-time indices within the next UP state
                    mask_next = timebin_edges_all(i,:) >= next_onset;
                    temp{nprobe}{mprobe}(i, mask_next) = NaN;
                end
            end


            % probability_merged(nprobe).(sprintf('%s_index',group_name{ngroup})) = temp;
        end
    end
  
    % 1 - ipsi and 2 - contra
    probability_merged(1).(sprintf('%s_UP',group_name{ngroup})) = [temp{1}{1}; temp{2}{2}];
    probability_merged(2).(sprintf('%s_UP',group_name{ngroup})) = [temp{1}{2}; temp{2}{1}];


    % DOWN
    event_times = merged_event_info.DOWN_ints;
    hemisphere_id = merged_event_info.DOWN_hemisphere_id;
    temp = [];
    for nprobe = 1:2
        ints = event_times(hemisphere_id==nprobe,:);
        for mprobe = 1:2
            [~,temp{nprobe}{mprobe},~] = calculate_event_probability(ripples_times(index{ngroup}{mprobe},:), ints(:,1), time_windows(1):time_bin:time_windows(end),0);
            timebin_edges_all = ints(:,1) + bins_centre;  % Absolute times of peri-event window
            for i = 1:size(ints,1)
                % Previous DOWN (skip if this is the first UP)
                if i > 1
                    prev_offset = ints(i-1,2);
                    % Find peri-time indices within the previous UP state
                    mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                    temp{nprobe}{mprobe}(i, mask_prev) = NaN;
                end

                % Next DOWN (skip if this is the last UP)
                if i < size(ints,1)
                    next_onset = ints(i+1,1);
                    % Find peri-time indices within the next UP state
                    mask_next = timebin_edges_all(i,:) >= next_onset;
                    temp{nprobe}{mprobe}(i, mask_next) = NaN;
                end
            end


            % probability_merged(nprobe).(sprintf('%s_index',group_name{ngroup})) = temp;
        end
    end

    % 1 - ipsi and 2 - contra
    probability_merged(1).(sprintf('%s_DOWN',group_name{ngroup})) = [temp{1}{1}; temp{2}{2}];
    probability_merged(2).(sprintf('%s_DOWN',group_name{ngroup})) = [temp{1}{2}; temp{2}{1}];

    % bootstrap distribution
    for nprobe = 1:2
        binnedArrayUP = probability_merged(nprobe).(sprintf('%s_UP',group_name{ngroup}));
        binnedArrayDOWN = probability_merged(nprobe).(sprintf('%s_DOWN',group_name{ngroup}));
        tempUP = [];
        tempDOWN = [];

        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArrayUP,1),size(binnedArrayUP,1));
            tempUP(iBoot,:) = sum(binnedArrayUP(event_id,:),'omitnan')./sum(~isnan(binnedArrayUP(event_id,:)));

            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArrayDOWN,1),size(binnedArrayDOWN,1));
            tempDOWN(iBoot,:) = sum(binnedArrayDOWN(event_id,:),'omitnan')./sum(~isnan(binnedArrayDOWN(event_id,:)));
        end

        probability_merged(nprobe).(sprintf('%s_UP_bootstrap',group_name{ngroup})) = tempUP;
        probability_merged(nprobe).(sprintf('%s_DOWN_bootstrap',group_name{ngroup})) = tempDOWN;
    end

end



num_bins = 20;
probability_normalised_merged = [];
for ngroup = 1:length(index)
    event_times = merged_event_info.UP_ints;
    hemisphere_id = merged_event_info.UP_hemisphere_id;
    temp = [];
    for nprobe = 1:2
        ints = event_times(hemisphere_id==nprobe,:);
        for mprobe = 1:2
            [~,temp{nprobe}{mprobe},~] = calculate_relative_event_probability(ints,ripples_times(index{ngroup}{mprobe},:),num_bins,0);
            timebin_edges_all = ints(:,1) + bins_centre;  % Absolute times of peri-event window
            for i = 1:size(ints,1)
                % Previous DOWN (skip if this is the first UP)
                if i > 1
                    prev_offset = ints(i-1,2);
                    % Find peri-time indices within the previous UP state
                    mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                    temp{nprobe}{mprobe}(i, mask_prev) = NaN;
                end

                % Next DOWN (skip if this is the last UP)
                if i < size(ints,1)
                    next_onset = ints(i+1,1);
                    % Find peri-time indices within the next UP state
                    mask_next = timebin_edges_all(i,:) >= next_onset;
                    temp{nprobe}{mprobe}(i, mask_next) = NaN;
                end
            end


            % probability_merged(nprobe).(sprintf('%s_index',group_name{ngroup})) = temp;
        end
    end
  
    % 1 - ipsi and 2 - contra
    probability_normalised_merged(1).(sprintf('%s_UP',group_name{ngroup})) = [temp{1}{1}; temp{2}{2}];
    probability_normalised_merged(2).(sprintf('%s_UP',group_name{ngroup})) = [temp{1}{2}; temp{2}{1}];


    % DOWN
    event_times = merged_event_info.DOWN_ints;
    hemisphere_id = merged_event_info.DOWN_hemisphere_id;
    temp = [];
    for nprobe = 1:2
        ints = event_times(hemisphere_id==nprobe,:);
        for mprobe = 1:2

            [~,temp{nprobe}{mprobe},~] = calculate_relative_event_probability(ints,ripples_times(index{ngroup}{mprobe},:),num_bins,0);
            timebin_edges_all = ints(:,1) + bins_centre;  % Absolute times of peri-event window
            for i = 1:size(ints,1)
                % Previous DOWN (skip if this is the first UP)
                if i > 1
                    prev_offset = ints(i-1,2);
                    % Find peri-time indices within the previous UP state
                    mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                    temp{nprobe}{mprobe}(i, mask_prev) = NaN;
                end

                % Next DOWN (skip if this is the last UP)
                if i < size(ints,1)
                    next_onset = ints(i+1,1);
                    % Find peri-time indices within the next UP state
                    mask_next = timebin_edges_all(i,:) >= next_onset;
                    temp{nprobe}{mprobe}(i, mask_next) = NaN;
                end
            end
            % probability_merged(nprobe).(sprintf('%s_index',group_name{ngroup})) = temp;
        end
    end

    % 1 - ipsi and 2 - contra
    probability_normalised_merged(1).(sprintf('%s_DOWN',group_name{ngroup})) = [temp{1}{1}; temp{2}{2}];
    probability_normalised_merged(2).(sprintf('%s_DOWN',group_name{ngroup})) = [temp{1}{2}; temp{2}{1}];

    % bootstrap distribution
    for nprobe = 1:2
        binnedArrayUP = probability_normalised_merged(nprobe).(sprintf('%s_UP',group_name{ngroup}));
        binnedArrayDOWN = probability_normalised_merged(nprobe).(sprintf('%s_DOWN',group_name{ngroup}));
        tempUP = [];
        tempDOWN = [];

        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArrayUP,1),size(binnedArrayUP,1));
            tempUP(iBoot,:) = sum(binnedArrayUP(event_id,:),'omitnan')./sum(~isnan(binnedArrayUP(event_id,:)));

            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArrayDOWN,1),size(binnedArrayDOWN,1));
            tempDOWN(iBoot,:) = sum(binnedArrayDOWN(event_id,:),'omitnan')./sum(~isnan(binnedArrayDOWN(event_id,:)));
        end

        probability_normalised_merged(nprobe).(sprintf('%s_UP_bootstrap',group_name{ngroup})) = tempUP;
        probability_normalised_merged(nprobe).(sprintf('%s_DOWN_bootstrap',group_name{ngroup})) = tempDOWN;
    end

end


save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_merged_LFP_lag.mat'),'probability_merged');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole_merged_LFP_lag.mat'),'probability_normalised_merged');


%% Grouping DOWN/UP events based on lags
% UP DOWN based on detection lags 
% Putatively use 0.05s as threshold for bilateral shared
% 0.2s as threshold for bilaterally coordinated but with unilaterally biased
% > 0.2s and/or non overlapping events 

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_merged_LFP_lag.mat'),'probability_merged');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole_merged_LFP_lag.mat'),'probability_normalised_merged');
probability_merged_LFP_lag = probability_merged;
probability_normalised_merged_LFP_lag = probability_normalised_merged;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_merged.mat'),'probability_merged');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole_merged.mat'),'probability_normalised_merged');


%%%%% UP info
event_times = merged_event_info.UP_ints;
hemisphere_id = merged_event_info.UP_hemisphere_id;
% lag_diff = merged_event_info.UP_lag_diff;
% group_id = merged_event_info.UP_group_id;
all_overlap_idx = merged_event_info.UP_overlap_idx_all{end};
non_overlap_idx = merged_event_info.UP_non_overlap_idx{end};
lags =merged_event_info.UP_lags_all{end};

% lags =all_lags{end};



% group_name{4} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% purely based on lags exlcuding cluster 2
% group_name{5} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
% group_name{6} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
% group_name{7} = {'Ipsi dominant ipsi leading','Bilaterally synchronised','Contra dominant contra leading','Shuffled'};


%%%%%%%%%%%%% Grouping events
event_idx = [];
% 1- ipsi leading 2 - contra leading 3 - bilateral 4 - unilateral
lag_thresholds = prctile(lags,[0:5:100]);
event_idx = {all_overlap_idx(lags>-0.2&lags<-0.05),all_overlap_idx(lags<0.05&lags>-0.05),all_overlap_idx(lags>0.05&lags<0.2),non_overlap_idx};

group_name = {'Ipsi leading','Bilaterally synchronised','Contra leading','Non-overlapping','Shuffled'};


%%%%%%%%%% Ripple-UP event pairs
event_pairs_UP = [];
event_pairs_DOWN = [];
group_name = {'Ipsi leading','Bilaterally synchronised','Contra leading','Unilateral'};

for n = 1:length(event_idx)
    for nprobe = 1:2 % 1 is ipsi and 2 is contra
        event_pairs_UP{n}{2}{nprobe} =any(probability_merged(nprobe).bilateral_ripple_UP(event_idx{n},:),2);
        event_pairs_DOWN{n}{2}{nprobe} = any(probability_merged(nprobe).bilateral_ripple_DOWN(event_idx{n},:),2);

        event_pairs_UP{n}{1}{nprobe} =any(probability_merged(nprobe).ipsi_leading_ripple_UP(event_idx{n},:),2);
        event_pairs_DOWN{n}{1}{nprobe} = any(probability_merged(nprobe).ipsi_leading_ripple_DOWN(event_idx{n},:),2);


        event_pairs_UP{n}{3}{nprobe} =any(probability_merged(nprobe).contra_leading_ripple_UP(event_idx{n},:),2);
        event_pairs_DOWN{n}{3}{nprobe} = any(probability_merged(nprobe).contra_leading_ripple_DOWN(event_idx{n},:),2);

        event_pairs_UP{n}{4}{nprobe} =any(probability_merged(nprobe).unilateral_ripple_UP(event_idx{n},:),2);
        event_pairs_DOWN{n}{4}{nprobe} = any(probability_merged(nprobe).unilateral_ripple_DOWN(event_idx{n},:),2);
    end
end

% bilateral UP -bilateral ripples

index=[];

for nprobe = 1:2
    index{1}{nprobe} = intersect(ripples_all_overlap_idx(ripples_all_lags<-0.005& ripples_all_lags>-0.02),find(ripples_hemisphere_id == nprobe))';% leading
    index{2}{nprobe} = intersect(ripples_all_overlap_idx(ripples_all_lags>0.005 & ripples_all_lags<0.02),find(ripples_hemisphere_id == nprobe))';% lagging
    index{3}{nprobe} = ripples_all_overlap_idx(ripples_all_lags<=0&ripples_all_lags>-0.005); % bilaterally synchronised
    index{4}{nprobe} = [intersect(ripples_all_overlap_idx(ripples_all_lags<-0.02 | ripples_all_lags>0.02),find(ripples_hemisphere_id == nprobe))' intersect(ripples_non_overlap_idx,find(ripples_hemisphere_id == nprobe))']; % unilateral
end

% sum(event_pairs_UP{2}{2}{1})/length(index{3}{1})
% % sum(event_pairs_UP{2}{2}{2})
% 
% % ipsi leading UP -ipsi leading ripples
% sum(event_pairs_UP{1}{1}{1})/length(index{1}{1})
% sum(event_pairs_UP{1}{1}{2})
% 
% 
% % unilateral events
% sum(event_pairs_UP{4}{4}{1})/length(index{4}{1})
% sum(event_pairs_UP{4}{4}{2})
% 
% sum(event_pairs_UP{4}{4}{1})
% sum(event_pairs_UP{4}{4}{2})

% figure
% time_wondows = [-1 1];
% time_bin = 0.02;
% x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
% nexttile
% plot(x,mean(probability_merged(1).unilateral_ripple_UP(event_pairs_UP{2}{2}{1},:),'omitnan'));hold on;
% plot(x,mean(probability_merged(2).unilateral_ripple_UP(event_pairs_UP{2}{2}{1},:),'omitnan'));
% title('Unilateral ripples during bilateral UP')
% legend('ipsi','contra')
% 
% nexttile
% plot(x,mean(probability_merged(1).unilateral_ripple_DOWN(event_pairs_DOWN{2}{2}{1},:),'omitnan'));hold on;
% plot(x,mean(probability_merged(2).unilateral_ripple_DOWN(event_pairs_DOWN{2}{2}{1},:),'omitnan'));
% title('Unilateral ripples during bilateral DOWN')
% legend('ipsi','contra')
% 
% nexttile
% plot(x,mean(probability_merged(1).unilateral_ripple_UP(event_pairs_UP{4}{4}{1},:),'omitnan'));hold on;
% plot(x,mean(probability_merged(2).unilateral_ripple_UP(event_pairs_UP{4}{4}{2},:),'omitnan'));
% title('Unilateral ripples during unilateral UP')
% legend('ipsi','contra')
% 
% nexttile
% plot(x,mean(probability_merged(1).unilateral_ripple_DOWN(event_pairs_DOWN{4}{4}{1},:),'omitnan'));hold on;
% plot(x,mean(probability_merged(2).unilateral_ripple_DOWN(event_pairs_DOWN{4}{4}{2},:),'omitnan'));
% title('Unilateral ripples during unilateral DOWN')
% legend('ipsi','contra')
% 
% 
% 
% figure
% time_wondows = [-1 1];
% time_bin = 0.02;
% x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
% nexttile
% plot(x,mean(probability_merged(1).bilateral_ripple_UP(event_pairs_UP{2}{2}{1},:),'omitnan'));hold on;
% plot(x,mean(probability_merged(1).bilateral_ripple_UP(event_pairs_UP{4}{4}{1},:),'omitnan'));
% plot(x,mean(probability_merged(1).bilateral_ripple_UP(event_pairs_UP{4}{4}{2},:),'omitnan'));
% title('Bilateral ripples during UP')
% legend('bilateral','ipsi','contra')
% 
% nexttile
% plot(x,mean(probability_merged(1).bilateral_ripple_DOWN(event_pairs_DOWN{2}{2}{1},:),'omitnan'));hold on;
% plot(x,mean(probability_merged(1).bilateral_ripple_DOWN(event_pairs_DOWN{4}{4}{1},:),'omitnan'));
% plot(x,mean(probability_merged(1).bilateral_ripple_DOWN(event_pairs_DOWN{4}{4}{2},:),'omitnan'));
% title('Bilateral ripples during DOWN')
% legend('bilateral','ipsi','contra')


%% Ipsi and contra ripple probabilities during DOWN UP with different bilateral synchrony

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline.mat'));
probability_psth_whole_baseline = probability;


load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
probability_normalised_whole = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole_baseline.mat'));
probability_normalised_whole_baseline = probability_normalised;

% probability = []
Delta_peaks_zscore_UD = [];
Delta_peaks_zscore_DU = [];
SWpeakmag_UD = [];
SWpeakmag_DU = [];
index = [];
%%%%% UP info
% UP_index_all = [probability(1).UP_index; probability(2).UP_index];
% DOWN_index_all = [probability(1).DOWN_index; probability(2).DOWN_index];

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

%%%%% UP info
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

title_names{1} = 'Ipsi-contra DOWN_UP ripples by 5 lags (150ms windows)';
title_names{2} = 'Ipsi-contra DOWN_UP ripples by 5 lags (150ms windows abs)';
title_names{3} = 'Ipsi-contra DOWN_UP ripples by 5 lags (50ms windows)';
title_names{4} = 'Ipsi-contra DOWN_UP ripples by 5 powers (150ms windows)';
title_names{5} = 'Ipsi-contra DOWN_UP ripples by 5 powers (All events)';
title_names{6} = 'Ipsi-contra DOWN_UP ripples by 5 duration (All windows)';
title_names{7} = 'Ipsi-contra DOWN_UP ripples by 5 duration (150ms events)';


%%%%%%%%%%%%%%%%%%%%%
ipsi_probability = [probability_psth_whole(1).L_ripples_UP; probability_psth_whole(2).R_ripples_UP];
contra_probability = [probability_psth_whole(1).R_ripples_UP; probability_psth_whole(2).L_ripples_UP];

ipsi_probability_baseline = [probability_psth_whole_baseline(1).L_ripples_UP; probability_psth_whole_baseline(2).R_ripples_UP];
contra_probability_baseline = [probability_psth_whole_baseline(1).R_ripples_UP; probability_psth_whole_baseline(2).L_ripples_UP];


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Calculate bootstrapped MUA

% clear probability_merged
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

probability_merged.x = x;
probability_merged.ipsi_ripples_UP = [];
probability_merged.contra_ripples_UP = [];
probability_merged.ipsi_contra_diff_ripples_UP = [];
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

        probability_merged.ipsi_ripples_UP{ngroup}{i} = temp1;
        probability_merged.contra_ripples_UP{ngroup}{i} = temp2;
        probability_merged.ipsi_contra_diff_ripples_UP{ngroup}{i} = temp3;
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


probability_merged.ipsi_ripples_baseline_UP = ipsi_baseline_bootstrap;
probability_merged.contra_ripples_baseline_UP = contra_baseline_bootstrap;
probability_merged.ipsi_contra_diff_ripples_baseline_UP = ipsi_contra_diff_baseline_bootstrap;
probability_merged.ripples_UP_groups = [title_names 'All DOWN-UP'];
probability_merged.ripples_UP_index = [event_idx {(1:size(ipsi_probability,1))'}];


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
    if ngroup ==1  | ngroup ==3
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    elseif ngroup==2
        colour_lines = [228,42,168;222,119,174;184,225,134;127,188,65;0,90,50]/256; % From magenta to dark purple to indicate globally synchrnous to ipsi lagging
    end

    if ngroup >3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = probability_merged.ipsi_ripples_UP{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = probability_merged.ipsi_ripples_baseline_UP{ngroup};
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([0 0.08])
    title('ipsi ripples')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    if ngroup >3
        % colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = probability_merged.contra_ripples_UP{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = probability_merged.contra_ripples_baseline_UP{ngroup};
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);
    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([0 0.08])

    % xline(0,'r')
    title('contra ripples')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    if ngroup >3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = probability_merged.ipsi_contra_diff_ripples_UP{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = probability_merged.ipsi_contra_diff_ripples_baseline_UP{ngroup};
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.04 0.04])


    % xline(0,'r')
    title('Ipsi-contra diff')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')

%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%% All DOWN UP rpples
event_averaging_scale = 30;

% for ngroup = 1:length(event_idx)
ngroup = 2;
fig = figure('Color','w');
fig.Position = [350 59 1650 465];
fig.Name = 'Left-Right combined ipsi contra ripple distribution around DOWN-UP transition'

colour_lines = [0,90,50;74,20,134]/256; % Green Purple


nexttile
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
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Event sorted by UP duration')
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
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Event sorted by UP duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('contra ripples')

nexttile
clear ERROR_SHADE

binnedArray = probability_merged.ipsi_ripples_UP{end}{1};
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

binnedArray = probability_merged.contra_ripples_UP{end}{1};
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline
binnedArray = probability_merged.ipsi_ripples_baseline_UP{ngroup};
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

% xline(0,'r')
ylim([0 0.06])
% title('ipsi ripples')
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Probability')
legend([ERROR_SHADE(1:end)],{'ipsi','contra','shuffled'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')






%% Ipsi and contra ripple probabilities during UP DOWN transition with different bilateral synchrony

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
    event_idx{5}{n} =find(delta_power>power_thresholds(n)&delta_power <power_thresholds(n+1));
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

title_names{1} = 'Ipsi-contra UP_DOWN ripples by 5 lags (150ms windows)';
title_names{2} = 'Ipsi-contra UP_DOWN ripples by 5 lags (150ms windows abs)';
title_names{3} = 'Ipsi-contra UP_DOWN ripples by 5 lags (50ms windows)';
title_names{4} = 'Ipsi-contra UP_DOWN ripples by 5 powers (150ms windows)';
title_names{5} = 'Ipsi-contra UP_DOWN ripples by 5 powers (All events)';
title_names{6} = 'Ipsi-contra UP_DOWN ripples by 5 duration (All windows)';
title_names{7} = 'Ipsi-contra UP_DOWN ripples by 5 duration (150ms events)';

%%%%%%%%%%%%%%%%%%%%%
ipsi_probability = [probability_psth_whole(1).L_ripples_DOWN; probability_psth_whole(2).R_ripples_DOWN];
contra_probability = [probability_psth_whole(1).R_ripples_DOWN; probability_psth_whole(2).L_ripples_DOWN];

ipsi_probability_baseline = [probability_psth_whole_baseline(1).L_ripples_DOWN; probability_psth_whole_baseline(2).R_ripples_DOWN];
contra_probability_baseline = [probability_psth_whole_baseline(1).R_ripples_DOWN; probability_psth_whole_baseline(2).L_ripples_DOWN];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Calculate bootstrapped MUA

% clear probability_merged
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

probability_merged.x = x;
probability_merged.ipsi_ripples_DOWN = [];
probability_merged.contra_ripples_DOWN = [];
probability_merged.ipsi_contra_diff_ripples_DOWN = [];

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

        probability_merged.ipsi_ripples_DOWN{ngroup}{i} = temp1;
        probability_merged.contra_ripples_DOWN{ngroup}{i} = temp2;
        probability_merged.ipsi_contra_diff_ripples_DOWN{ngroup}{i} = temp3;
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


probability_merged.ipsi_ripples_baseline_DOWN = ipsi_baseline_bootstrap;
probability_merged.contra_ripples_baseline_DOWN = contra_baseline_bootstrap;
probability_merged.ipsi_contra_diff_ripples_baseline_DOWN = ipsi_contra_diff_baseline_bootstrap;
probability_merged.ripples_DOWN_groups = [title_names 'All UP-DOWN’'];
probability_merged.ripples_DOWN_index = [event_idx {(1:size(ipsi_probability,1))'}];


time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

%%%%%%%%% Plot DU transition
for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 465];
    fig.Name =title_names{ngroup};

    if ngroup ==1  | ngroup ==3
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    elseif ngroup==2
        colour_lines = [228,42,168;222,119,174;184,225,134;127,188,65;0,90,50]/256; % From magenta to dark purple to indicate globally synchrnous to ipsi lagging
    end


    if ngroup >3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end
    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = probability_merged.ipsi_ripples_DOWN{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = probability_merged.ipsi_ripples_baseline_DOWN{ngroup};
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([0 0.17])
    title('ipsi ripples')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    if ngroup >3
        % colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = probability_merged.contra_ripples_DOWN{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = probability_merged.contra_ripples_baseline_DOWN{ngroup};
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);
    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([0 0.17])

    % xline(0,'r')
    title('contra ripples')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    if ngroup >3
        colour_lines = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for
    end
    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = probability_merged.ipsi_contra_diff_ripples_DOWN{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = probability_merged.ipsi_contra_diff_ripples_baseline_DOWN{ngroup};
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.03 0.03])


    % xline(0,'r')
    title('Ipsi-contra diff')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')




%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%% All UP DOWN rpples
event_averaging_scale = 30;

% for ngroup = 1:length(event_idx)
ngroup = 2;
fig = figure('Color','w');
fig.Position = [350 59 1650 465];
fig.Name = 'Left-Right combined ipsi contra ripple distribution around UP-DOWN transition'

colour_lines = [0,90,50;74,20,134]/256; % Green Purple


nexttile
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
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
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
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('contra ripples')

nexttile
clear ERROR_SHADE

binnedArray = probability_merged.ipsi_ripples_DOWN{end}{1};
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

binnedArray = probability_merged.contra_ripples_DOWN{end}{1};
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline
binnedArray = probability_merged.ipsi_ripples_baseline_DOWN{end};
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

% xline(0,'r')
ylim([0 0.09])
% title('ipsi ripples')
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Probability')
legend([ERROR_SHADE(1:end)],{'ipsi','contra','shuffled'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')


save(fullfile(analysis_folder,'V1-HPC bilateral interaction','SO_ripples_probability_merged.mat'),'probability_merged')

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