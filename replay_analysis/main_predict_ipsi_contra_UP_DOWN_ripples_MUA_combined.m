%%%%%% Main script for mixed effect and cox regression analysis to examine
%%%%%% predictive relationship between ripples and DOWN UP states and
%%%%%% spindles

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('C:\Users\masah\OneDrive\Documents\corticohippocampal_replay')
    analysis_folder = 'C:\Users\masah\OneDrive\Documents\corticohippocampal_replay';
elseif exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
% % load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
load(fullfile(analysis_folder,'ripples_all_POST.mat'))
load(fullfile(analysis_folder,'spindles_all_POST.mat'))
% load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov_normalised.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov.mat'));


load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
% probability_normalised = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_combined.mat'));
probability_psth_whole = probability;


load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'));
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');

load(fullfile(analysis_folder,'periripple_LFP_info_V1.mat'));



all_sessions = max(slow_waves_all(1).DOWN_session_count);
sessions_to_process = 1:all_sessions;

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
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','k_cluster_ipsi_contra_events.mat'),'k_cluster');
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');


% load(fullfile(analysis_folder,'V1-HPC sleep interaction','MUA_PSTH_merged_thresholded.mat'),'MUA_PSTH_merged_thresholded');


%%%%%%%%%%%%% Overlaps
% Initialize output
L_overlap_idx = [];
R_overlap_idx = [];
overlap_idx = [];
non_overlap_idx = [];
merged_idx = [];
all_lags = [];
all_overlap_idx=[];
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


%% Ripples and spindles features and UP DOWN features during UP/DOWN
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_combined.mat'));
probability_psth_whole = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_probability_whole.mat'));
spindle_probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline_combined.mat'));
probability_psth_whole_baseline = probability;

load(fullfile(analysis_folder,'V1-HPC sleep reactivation','ripple_info.mat'),'ripple_info')
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
% probability_normalised_whole = probability;
% contra_lag_ripples

% Get DU and UP delta power
ipsi_Delta_peaks_zscore_UD = [];
ipsi_Delta_peaks_zscore_DU = [];
SWpeakmag_UD = [];
SWpeakmag_DU = [];

contra_Delta_peaks_zscore_UD = [];
contra_Delta_peaks_zscore_DU = [];

index = [];
%%%%% UP info
% UP_index_all = [probability(1).UP_index; probability(2).UP_index];
% DOWN_index_all = [probability(1).DOWN_index; probability(2).DOWN_index];

for nprobe = 1:2
    for nsession = 1:max(slow_waves_all(1).UP_session_count)

        % slow_waves_all(1).SWpeakmag(slow_waves_all(1).UP_session_count == nsession)
        % Find DOWN
        [C,ia,ib]  =intersect(find(slow_waves_all(nprobe).DOWN_session_count == nsession),probability(nprobe).DOWN_all_index);

        ipsi_Delta_peaks_zscore_UD = [ipsi_Delta_peaks_zscore_UD;...
            slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession,nprobe),ia)'];
        SWpeakmag_UD = [SWpeakmag_UD;...
            slow_waves_all(nprobe).SWpeakmag(C)];

        % Find DOWN before UP
        [C,ia,ib]  = intersect(slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).DOWN_session_count == nsession,2), slow_waves_all(nprobe).UP_ints(intersect(find(slow_waves_all(nprobe).UP_session_count == nsession),probability(nprobe).UP_all_index),1));

        ipsi_Delta_peaks_zscore_DU = [ipsi_Delta_peaks_zscore_DU; slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession,nprobe),ia)'];
        index = [index; intersect(find(slow_waves_all(nprobe).DOWN_session_count == nsession),probability(nprobe).DOWN_all_index)];

        [C,ia,ib]  = intersect(slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).DOWN_session_count == nsession,2), slow_waves_all(nprobe).UP_ints(intersect(find(slow_waves_all(nprobe).UP_session_count == nsession),probability(nprobe).UP_all_index),1));

        temp = find(slow_waves_all(nprobe).DOWN_session_count == nsession);

        SWpeakmag_DU = [SWpeakmag_DU;...
            slow_waves_all(nprobe).SWpeakmag(temp(ia))];



        % slow_waves_all(1).SWpeakmag(slow_waves_all(1).UP_session_count == nsession)
        % Find DOWN
        [C,ia,ib]  =intersect(find(slow_waves_all(nprobe).DOWN_session_count == nsession),probability(nprobe).DOWN_all_index);

        contra_Delta_peaks_zscore_UD = [contra_Delta_peaks_zscore_UD;...
            slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession,abs(nprobe-3)),ia)'];

        % Find DOWN before UP
        [C,ia,ib]  = intersect(slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).DOWN_session_count == nsession,2), slow_waves_all(nprobe).UP_ints(intersect(find(slow_waves_all(nprobe).UP_session_count == nsession),probability(nprobe).UP_all_index),1));

        contra_Delta_peaks_zscore_DU = [contra_Delta_peaks_zscore_DU; slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession,abs(nprobe-3)),ia)'];
        index = [index; intersect(find(slow_waves_all(nprobe).DOWN_session_count == nsession),probability(nprobe).DOWN_all_index)];

        [C,ia,ib]  = intersect(slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).DOWN_session_count == nsession,2), slow_waves_all(nprobe).UP_ints(intersect(find(slow_waves_all(nprobe).UP_session_count == nsession),probability(nprobe).UP_all_index),1));

        temp = find(slow_waves_all(nprobe).DOWN_session_count == nsession);

    end
end

%%%%% Ripple info (combined)

UP_ints=[];
DOWN_ints=[];
ripple_peaktimes=[];
ripple_ints=[];
SO_ints=[];
prev_down_idx=[];
next_down_idx=[];

for nprobe = 1:2
    UP_ints{nprobe}=slow_waves_all(nprobe).UP_ints;
    DOWN_ints{nprobe}=slow_waves_all(nprobe).DOWN_ints;
    SO_ints{nprobe} = slow_waves_all(nprobe).DOWN_intervals;
    ripple_peaktimes{nprobe}=ripples_all(nprobe).peaktimes;
    ripple_ints{nprobe}=[ripples_all(nprobe).onset(ripples_all(nprobe).SWS_index == 1) ripples_all(nprobe).offset(ripples_all(nprobe).SWS_index == 1)];
    % spindle_peaktimes{nprobe}=spindles_all(nprobe).peaktimes(spindles_all(nprobe).SWS_index == 1);
    % spindle_ints{nprobe}=[spindles_all(nprobe).onset(spindles_all(nprobe).SWS_index == 1) spindles_all(nprobe).offset(spindles_all(nprobe).SWS_index == 1)];

    for nsession = 1:max(slow_waves_all(1).UP_session_count)
        index = find(slow_waves_all(nprobe).DOWN_intervals_session == sessions_to_process(nsession));
        SO_ints{nprobe}(index,:) = SO_ints{nprobe}(index,:) + nsession * 1000000;

        index = find(slow_waves_all(nprobe).UP_session_count == sessions_to_process(nsession));
        UP_ints{nprobe}(index,:) = UP_ints{nprobe}(index,:) + nsession * 1000000;

        index = find(slow_waves_all(nprobe).DOWN_session_count == sessions_to_process(nsession));
        DOWN_ints{nprobe}(index,:) = DOWN_ints{nprobe}(index,:) + nsession * 1000000;

        index = find(ripples_all(nprobe).session_count(ripples_all(nprobe).SWS_index == 1) == sessions_to_process(nsession));
        ripple_ints{nprobe}(index,:) = ripple_ints{nprobe}(index,:) + nsession * 1000000;
        ripple_peaktimes{nprobe}(index,:) = ripple_peaktimes{nprobe}(index,:) + nsession * 1000000;

        % [C,ia,ib] = intersect(find(spindles_all(nprobe).session_count == sessions_to_process(nsession)),find(spindles_all(nprobe).SWS_index == 1));
        % spindle_ints{nprobe}(ib,:) = spindle_ints{nprobe}(ib,:) + nsession * 1000000;
    end

    nUP = length(UP_ints{nprobe});
    prev_down_idx{nprobe} = nan(nUP, 1);
    next_down_idx{nprobe} = nan(nUP, 1);

    % tol = 1e-6; % Tolerance for floating point comparison

    % 1. Find Previous DOWN: Where DOWN_offset == UP_onset
    % We look for UP(:,1) inside the list of DOWN(:,2)
    [is_prev, prev_down_idx{nprobe}] = ismembertol(UP_ints{nprobe}(:,1), DOWN_ints{nprobe}(:,2), 1e-10);
    prev_down_idx{nprobe}(~is_prev) = nan; % Set zeros to NaN

    % 2. Find Next DOWN: Where DOWN_onset == UP_offset
    % We look for UP(:,2) inside the list of DOWN(:,1)
    [is_next, next_down_idx{nprobe}] = ismembertol(UP_ints{nprobe}(:,2), DOWN_ints{nprobe}(:,1), 1e-10);
    next_down_idx{nprobe}(~is_next) = nan; % Set zeros to NaN
end

ripples_peaktimes =  merged_event_info.ripples_peaktimes;
ripples_times =  merged_event_info.ripples_ints;
% ripples_group_info = merged_event_info.ripples_group_id;
ripples_hemisphere_id = merged_event_info.ripples_hemisphere_id;
ripples_original_index = [find(ripples_all(1).SWS_index); find(ripples_all(2).SWS_index)];
% [event_ids_first,event_ids_second] = merge_bilateral_ripple_events(ripples_hemisphere_id,ripples_times(:,1),0.05);
[event_ids_first,event_ids_second] = merge_bilateral_ripple_events(ripples_hemisphere_id,ripples_peaktimes,0.05);


ripples_hemisphere_id = ripples_hemisphere_id(event_ids_first);
ripples_original_index = ripples_original_index(event_ids_first);
ripples_times = ripples_times(event_ids_first,:);
ripples_peaktimes = ripples_peaktimes(event_ids_first,:);
% ripples_lag_diff = [contra_lag_ripples{1} contra_lag_ripples{2}]';ripples_lag_diff = ripples_lag_diff(event_ids_first);
% ripples_all_overlap_idx = merged_event_info.ripples_overlap_idx_all{end};
% ripples_non_overlap_idx = merged_event_info.ripples_non_overlap_idx{end};ripples_non_overlap_idx = ripples_non_overlap_idx(event_ids_first);
% ripples_all_lags =merged_event_info.ripples_lags_all{end};


%%%% Grab log odds
load(fullfile(analysis_folder,'periripple_LFP_info_V1.mat'));
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_ripples_PSTH.mat'))
spindle_amplitude1 = [periripple_LFP_info_V1(1).spindle_amplitude{1}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).spindle_amplitude{1}(:,ripples_all(2).SWS_index==1)];
spindle_amplitude2 = [periripple_LFP_info_V1(1).spindle_amplitude{2}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).spindle_amplitude{2}(:,ripples_all(2).SWS_index==1)];
spindle_amplitude = nan([size(spindle_amplitude1),2]);
spindle_amplitude(:,:,1) = spindle_amplitude1;
spindle_amplitude(:,:,2) = spindle_amplitude2;

ripple_info.spindle_amplitude_temporal = spindle_amplitude;
ripple_info.spindle_amplitude_temporal = ripple_info.spindle_amplitude_temporal(:,event_ids_first,:);

SO_phase1 = [periripple_LFP_info_V1(1).SO_phase{1}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).SO_phase{1}(:,ripples_all(2).SWS_index==1)];
SO_phase2 = [periripple_LFP_info_V1(1).SO_phase{2}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).SO_phase{2}(:,ripples_all(2).SWS_index==1)];
SO_phase = nan([size(SO_phase1),2]);
SO_phase(:,:,1) = SO_phase1;
SO_phase(:,:,2) = SO_phase2;

ripple_info.SO_phase_temporal = SO_phase;
ripple_info.SO_phase_temporal = ripple_info.SO_phase_temporal(:,event_ids_first,:);
clear SO_phase SO_phase1 SO_phase2 spindle_amplitude spindle_amplitude1 spindle_amplitude2


timebin = 0.01;
time_windows = [-1 1];
% Generate bin edges
bin_edges = time_windows(1):timebin:time_windows(2);
% Generate bin centers
bin_centers = bin_edges(1:end-1) + timebin/2;
% z_bias = KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';
% z_bias_V1 = KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';

z_bias = KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples';
z_bias_V1 = KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples';

z_bias1 = z_bias(isfinite(z_bias));
z_bias(z_bias>=inf) = prctile(z_bias1,99.5);
z_bias(z_bias<=-inf) = prctile(z_bias1,0.5);

z_bias1 = z_bias(isfinite(z_bias_V1));
z_bias_V1(z_bias_V1>=inf) = prctile(z_bias1,99.5);
z_bias_V1(z_bias_V1<=-inf) = prctile(z_bias1,0.5);

%%%% Only grab unique ripple events
z_bias_raw = z_bias;
z_bias_V1_raw = z_bias_V1;
z_bias = z_bias+ KDE_reactivation_ripples_PSTH.nan_mask';
z_bias_V1 = z_bias_V1+ KDE_reactivation_ripples_PSTH.nan_mask';

z_bias = z_bias(:,event_ids_first);
z_bias_V1 = z_bias_V1(:,event_ids_first);

z_bias_raw = z_bias_raw(:,event_ids_first);
z_bias_V1_raw = z_bias_V1_raw(:,event_ids_first);

% mean_z_bias = mean(z_bias(bin_centers>0 & bin_centers<0.1,:),'omitnan');
% mean_z_bias_PRE = mean(z_bias(bin_centers>-0.2 & bin_centers<0,:),'omitnan');
% mean_z_bias_V1= mean(z_bias_V1(bin_centers>0 & bin_centers<0.2,:),'omitnan');
% mean_z_bias_V1_PRE= mean(z_bias_V1(bin_centers>-0.2 & bin_centers<0,:),'omitnan');
% mean_z_bias_V1_PRE_raw= mean(z_bias_V1_raw(bin_centers>-0.2 & bin_centers<0,:),'omitnan');

% merged_event_info.ripples_index_sorted
mean_z_bias = mean(z_bias(bin_centers>0 & bin_centers<0.1,:),'omitnan');
mean_z_bias_PRE = mean(z_bias(bin_centers>-0.1 & bin_centers<0,:),'omitnan');
mean_z_bias_V1= mean(z_bias_V1(bin_centers>0 & bin_centers<0.1,:),'omitnan');
mean_z_bias_V1_PRE= mean(z_bias_V1(bin_centers>-0.1 & bin_centers<0,:),'omitnan');
mean_z_bias_V1_PRE_raw= mean(z_bias_V1_raw(bin_centers>-0.1 & bin_centers<0,:),'omitnan');
% ripple_info.normalised_UP_duration(:,1)

% slow_waves_all
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_PSTH.mat'));

% load(fullfile(analysis_folder,'KDE_reactivation_UP_all_POST.mat'));

UP_HPC_log_odds=[];
UP_V1_log_odds=[];
DOWN_HPC_log_odds=[];
DOWN_V1_log_odds=[];

for nprobe = 1:2
    UP_V1_log_odds{nprobe} = KDE_reactivation_PSTH(nprobe).V1_UP_log_odds;
    temp1 = UP_V1_log_odds{nprobe}(isfinite(UP_V1_log_odds{nprobe}));
    UP_V1_log_odds{nprobe}(UP_V1_log_odds{nprobe}>=inf) = prctile(temp1,99.5);
    UP_V1_log_odds{nprobe}(UP_V1_log_odds{nprobe}<=-inf) = prctile(temp1,0.5);

    UP_HPC_log_odds{nprobe} = KDE_reactivation_PSTH(nprobe).HPC_UP_log_odds;
    temp1 = UP_HPC_log_odds{nprobe}(isfinite(UP_HPC_log_odds{nprobe}));
    UP_HPC_log_odds{nprobe}(UP_HPC_log_odds{nprobe}>=inf) = prctile(temp1,99.5);
    UP_HPC_log_odds{nprobe}(UP_HPC_log_odds{nprobe}<=-inf) = prctile(temp1,0.5);

    DOWN_V1_log_odds{nprobe} = KDE_reactivation_PSTH(nprobe).V1_DOWN_log_odds;
    temp1 = DOWN_V1_log_odds{nprobe}(isfinite(DOWN_V1_log_odds{nprobe}));
    DOWN_V1_log_odds{nprobe}(DOWN_V1_log_odds{nprobe}>=inf) = prctile(temp1,99.5);
    DOWN_V1_log_odds{nprobe}(DOWN_V1_log_odds{nprobe}<=-inf) = prctile(temp1,0.5);

    DOWN_HPC_log_odds{nprobe} = KDE_reactivation_PSTH(nprobe).HPC_DOWN_log_odds;
    temp1 = DOWN_HPC_log_odds{nprobe}(isfinite(DOWN_HPC_log_odds{nprobe}));
    DOWN_HPC_log_odds{nprobe}(DOWN_HPC_log_odds{nprobe}>=inf) = prctile(temp1,99.5);
    DOWN_HPC_log_odds{nprobe}(DOWN_HPC_log_odds{nprobe}<=-inf) = prctile(temp1,0.5);
end




%%%%%%%%%%%%%% Grab ripple info during UP
hemi_labels = {'L', 'R'};
region_labels = {'L_V1','R_V1','HPC'};
varnames = {
    'time_from_last_ripples'
    'first_ripples_index'
    'last_ripples_index'
    'first_ripples_power'
    'last_ripples_power'
    'first_ripples_lag'
    'last_ripples_lag'
    'first_ripples_lag_diff'
    'last_ripples_lag_diff'
    'ripple_counts'
    'ripples_duration'
    'first_ripples_duration'
    'last_ripples_duration'

    'first_ripples_log_odds'
    'last_ripples_log_odds'
    'first_ripples_V1_log_odds'
    'last_ripples_V1_log_odds'
    'first_ripples_PRE_V1_log_odds'
    'last_ripples_PRE_V1_log_odds'

    'second_last_ripples_log_odds'
    'previous_ripples_log_odds'
    'second_last_ripples_PRE_V1_log_odds'
    'previous_ripples_PRE_V1_log_odds'
    'second_last_ripples_V1_log_odds'
    'previous_ripples_V1_log_odds'
    'last_ripple_next_spindle_diff'
    'last_ripple_next_spindle_power'
    'first_ripple_spindle_power'
    'last_ripple_spindle_power'
};
clear data_struct
% Init data_struct
% for h = 1
for v = 1:length(varnames)
    data_struct.(varnames{v}) = cell(1,2);
end

for r = 1:3
    reg = region_labels{r};
    data_struct.(['ripple_' reg '_MUA_cumulative']) = cell(1,2);
    data_struct.(['normalised_ripple_' reg '_MUA_cumulative']) = cell(1,2);

    data_struct.(['first_ripple_' reg '_MUA_peak']) = cell(1,2);
    data_struct.(['last_ripple_' reg '_MUA_peak']) = cell(1,2);
    data_struct.(['first_ripple_' reg '_MUA_cumulative']) = cell(1,2);
    data_struct.(['last_ripple_' reg '_MUA_cumulative']) = cell(1,2);
end
% end

for r = 1:length(region_labels)
    reg = region_labels{r};
    data_struct.([reg '_MUA_cumulative']) = cell(1,2);
    data_struct.(['normalised_' reg '_MUA_cumulative']) = cell(1,2);
    data_struct.(['first_half_' reg '_MUA_cumulative']) = cell(1,2);
    data_struct.(['second_half_' reg '_MUA_cumulative']) = cell(1,2);
end

varnames = fieldnames(data_struct);


% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','ripple_info.mat'),'ripple_info')
% Main loop
for nprobe = 1:2
    UP_indices = probability_psth_whole(nprobe).UP_all_index;
%     event_info
    for h = 1
        ripple_data = ripple_info;
        % Normalised duration of ripples (combined across L and R ripples)
        % for UP less than 2 sec
        index1 = find(ismember(event_info(nprobe).L_ripple_normalised_UP_duration(:,1),ripples_original_index(ripples_hemisphere_id == 1)));
        temp1 = event_info(nprobe).L_ripple_normalised_UP_duration(index1,:);
        [~,ia] = ismember(temp1(:,1),ripples_original_index(ripples_hemisphere_id == 1));

        index2 = find(ismember(event_info(nprobe).R_ripple_normalised_UP_duration(:,1),ripples_original_index(ripples_hemisphere_id == 2)));
        temp2 = event_info(nprobe).R_ripple_normalised_UP_duration(index2,:);
        [~,ib] = ismember(temp2(:,1),ripples_original_index(ripples_hemisphere_id == 2));
 
        ripple_norm_dur = [temp1;temp2];
        ripple_norm_dur(:,1) = [ia; sum(ripples_hemisphere_id == 1)+ib]; % index conversion

        for reg = {'HPC','V1'}
            ripple_data.ripple_HPC_MUA_peak_UP = [event_info(nprobe).(['L_ripple_HPC_MUA_peak_UP'])(index1,:); event_info(nprobe).(['R_ripple_HPC_MUA_peak_UP'])(index2,:)];
            ripple_data.ripple_V1_MUA_peak_UP = [event_info(nprobe).(['L_ripple_V1_MUA_peak_UP'])(index1,:); event_info(nprobe).(['R_ripple_V1_MUA_peak_UP'])(index2,:)];
            
            ripple_data.ripple_HPC_MUA_cumulative_UP = [event_info(nprobe).(['L_ripple_HPC_MUA_cumulative_UP'])(index1,:); event_info(nprobe).(['R_ripple_HPC_MUA_cumulative_UP'])(index2,:)];
            ripple_data.ripple_V1_MUA_cumulative_UP = [event_info(nprobe).(['L_ripple_V1_MUA_cumulative_UP'])(index1,:); event_info(nprobe).(['R_ripple_V1_MUA_cumulative_UP'])(index2,:)];
        end       

        for nevent = 1:length(UP_indices)
            up_idx = UP_indices(nevent);

            % --- Always store UP-wide MUA metrics ---
            up_duration = event_info(nprobe).UP_duration(nevent);
            for r = 1:length(region_labels)
                reg = region_labels{r};
                if contains(reg,'HPC')
                    trace = mean([event_info(nprobe).L_HPC_MUA_UP{nevent}; event_info(nprobe).R_HPC_MUA_UP{nevent}],'omitnan');
                else
                    trace = event_info(nprobe).([reg '_MUA_UP']){nevent};
                end

                npoints = length(trace);
                half_idx = floor(npoints/2);
                cum_sum = sum(trace);
                data_struct.([reg '_MUA_cumulative']){nprobe}(nevent) = cum_sum;
                data_struct.(['normalised_' reg '_MUA_cumulative']){nprobe}(nevent) = cum_sum / up_duration;
                data_struct.(['first_half_' reg '_MUA_cumulative']){nprobe}(nevent) = sum(trace(1:half_idx)) / (up_duration/2);
                data_struct.(['second_half_' reg '_MUA_cumulative']){nprobe}(nevent) = sum(trace(half_idx+1:end)) / (up_duration/2);
            end

            % --- Ripple-based metrics (only if ripples exist in this UP) ---
            event_index = find(up_idx == ripple_norm_dur(:,2));

            if ~isempty(event_index)
                ripples_index = ripple_norm_dur(event_index, 1);
%                 ripples_index = find(ismember(ripples_original_index,ripples_index));

                onset = ripples_times(ripples_index,1);
                [onset,index] = sort(onset);
                ripples_index = ripples_index(index);

                offset = ripples_times(ripples_index,2);
                duration = offset - onset;

                data_struct.first_ripples_index{nprobe}(nevent) = ripples_index(1); % ripple index after combing L and R ripples
                data_struct.last_ripples_index{nprobe}(nevent) = ripples_index(end);

                data_struct.first_ripples_duration{nprobe}(nevent) = offset(1) - onset(1);
                data_struct.last_ripples_duration{nprobe}(nevent) = offset(end) - onset(end);
                data_struct.ripples_duration{nprobe}(nevent) = sum(duration);
                data_struct.ripple_counts{nprobe}(nevent) = length(ripples_index);
                data_struct.time_from_last_ripples{nprobe}(nevent) = ...
                    UP_ints{nprobe}(up_idx,2) - ripples_peaktimes(ripples_index(end));

                data_struct.first_ripples_power{nprobe}(nevent) = ripple_data.ripple_power(ripples_index(1));
                data_struct.last_ripples_power{nprobe}(nevent) = ripple_data.ripple_power(ripples_index(end));

                data_struct.first_ripples_log_odds{nprobe}(nevent) = mean_z_bias(ripples_index(1));
                data_struct.last_ripples_log_odds{nprobe}(nevent) = mean_z_bias(ripples_index(end));

                data_struct.first_ripples_V1_log_odds{nprobe}(nevent) = mean_z_bias_V1(ripples_index(1));
                data_struct.last_ripples_V1_log_odds{nprobe}(nevent) = mean_z_bias_V1(ripples_index(end));

                data_struct.first_ripples_PRE_V1_log_odds{nprobe}(nevent) = mean_z_bias_V1_PRE(ripples_index(1));
                data_struct.last_ripples_PRE_V1_log_odds{nprobe}(nevent) = mean_z_bias_V1_PRE(ripples_index(end));
                if length(ripples_index)>1
                    data_struct.second_last_ripples_log_odds{nprobe}(nevent) = mean_z_bias(ripples_index(end-1));
                    data_struct.second_last_ripples_PRE_V1_log_odds{nprobe}(nevent) = mean_z_bias_V1_PRE(ripples_index(end-1));
                    data_struct.second_last_ripples_V1_log_odds{nprobe}(nevent) = mean_z_bias_V1(ripples_index(end-1));
                else
                     data_struct.second_last_ripples_log_odds{nprobe}(nevent) = nan;
                     data_struct.second_last_ripples_PRE_V1_log_odds{nprobe}(nevent) = nan;
                     data_struct.second_last_ripples_V1_log_odds{nprobe}(nevent) = nan;
                end
                data_struct.previous_ripples_log_odds{nprobe}(nevent) = mean_z_bias(ripples_index(1)-1); % (last ripple before this UP state)
%                 event_index = find(up_idx-1 == ripple_norm_dur(:,2));
                data_struct.previous_ripples_PRE_V1_log_odds{nprobe}(nevent) = mean_z_bias_V1_PRE(ripples_index(1)-1); % (last ripple before this UP state)
                data_struct.previous_ripples_V1_log_odds{nprobe}(nevent) = mean_z_bias_V1(ripples_index(1)-1); % (last ripple before this UP state)

                [~,LFP_bin1] = min(abs(periripple_LFP_info_V1(1).tvec - (UP_ints{nprobe}(up_idx,2) - onset(end)))); % end of UP spindle power
                [~,LFP_bin2] = min(abs(periripple_LFP_info_V1(1).tvec - (UP_ints{nprobe}(up_idx+1,1) - onset(end))));% start of next UP spindle power

                 data_struct.last_ripple_next_spindle_diff{nprobe}(nevent,:) = squeeze(ripple_info.spindle_amplitude_temporal(LFP_bin2,ripples_index(end),:)-ripple_info.spindle_amplitude_temporal(LFP_bin1,ripples_index(end),:))';
                 data_struct.last_ripple_next_spindle_power{nprobe}(nevent,:) = squeeze(ripple_info.spindle_amplitude_temporal(LFP_bin2,ripples_index(end),:))';

                 [~,LFP_bin] = min(abs(periripple_LFP_info_V1(1).tvec - (-0.1)));% spindle power before ripple
                 data_struct.first_ripple_spindle_power{nprobe}(nevent,:) = squeeze(ripple_info.spindle_amplitude_temporal(LFP_bin,ripples_index(1),:));
                 data_struct.last_ripple_spindle_power{nprobe}(nevent,:) = squeeze(ripple_info.spindle_amplitude_temporal(LFP_bin,ripples_index(1),:));

%                 idx1 = find(ripples_original_index == ripples_index(1) & ripples_hemisphere_id == ripple_idx);
%                 idx2 = find(ripples_original_index == ripples_index(end) & ripples_hemisphere_id == ripple_idx);
%                 data_struct.(hemi).first_ripples_lag_diff{nprobe}(nevent) = ripples_lag_diff(idx1);
%                 data_struct.(hemi).last_ripples_lag_diff{nprobe}(nevent) = ripples_lag_diff(idx2);

%                 [C, ia] = intersect(ripples_all_overlap_idx, idx1);
%                 if ~isempty(C)
%                     data_struct.(hemi).first_ripples_lag{nprobe}(nevent) = ripples_all_lags(ia);
%                 end
%                 [C, ia] = intersect(ripples_all_overlap_idx, idx2);
%                 if ~isempty(C)
%                     data_struct.(hemi).last_ripples_lag{nprobe}(nevent) = ripples_all_lags(ia);
%                 end

                % MUA from event_info
                for r = 1:length(region_labels)
                    reg = region_labels{r};  % 'L_V1', 'HPC', etc.
                    region_hemi = strcmp(reg(1), 'R') + 1;  % 1 = L, 2 = R
                    region_type = reg(3:end);              % 'V1' or 'HPC'

                    if contains(region_labels{r},'HPC')
                        data_struct.(['ripple_' reg '_MUA_cumulative']){nprobe}(nevent) = sum(mean(ripple_data.ripple_HPC_MUA_cumulative_UP(event_index, :),2));
                        data_struct.(['normalised_ripple_' reg '_MUA_cumulative']){nprobe}(nevent) = sum(mean(ripple_data.ripple_HPC_MUA_cumulative_UP(event_index, :),2))./ sum(duration);
                        data_struct.(['first_ripple_' reg '_MUA_peak']){nprobe}(nevent) = mean(ripple_data.ripple_HPC_MUA_peak_UP(event_index(1), :));
                        data_struct.(['last_ripple_' reg '_MUA_peak']){nprobe}(nevent) = mean(ripple_data.ripple_HPC_MUA_peak_UP(event_index(end), :));
                    else
                        data_struct.(['ripple_' reg '_MUA_cumulative']){nprobe}(nevent) = sum(ripple_data.ripple_V1_MUA_cumulative_UP(event_index, region_hemi));
                        data_struct.(['normalised_ripple_' reg '_MUA_cumulative']){nprobe}(nevent) = sum(ripple_data.ripple_V1_MUA_cumulative_UP(event_index, region_hemi))./ sum(duration);
                        data_struct.(['first_ripple_' reg '_MUA_peak']){nprobe}(nevent) = (ripple_data.ripple_V1_MUA_peak_UP(event_index(1), region_hemi));
                        data_struct.(['last_ripple_' reg '_MUA_peak']){nprobe}(nevent) = (ripple_data.ripple_V1_MUA_peak_UP(event_index(end), region_hemi));
                    end
                end
            else
                data_struct.time_from_last_ripples{nprobe}(nevent) = 0;
                for v = 1:length(varnames)
                    if contains(varnames{v},'spindle');
                        data_struct.(varnames{v}){nprobe}(nevent,:) = nan(1,2);
                    else contains(varnames{v},'ripple');
                        data_struct.(varnames{v}){nprobe}(nevent) = nan;
                    end
                end

                for r = 1:length(region_labels)
                    reg = region_labels{r};
                    for suffix = {'ripple_', 'normalised_ripple_', 'first_ripple_', 'last_ripple_'}
                        fname = [suffix{1} reg '_MUA_cumulative'];
                        if contains(fname, 'peak')
                            fname = strrep(fname, '_cumulative', '_peak');
                        end
                        data_struct.(fname){nprobe}(nevent) = nan;
                    end
                end
            end
        end
    end
end

% Build UP_DOWN_info
UP_DOWN_info = struct();
varnames = fieldnames(data_struct);
% varnames1 = fieldnames(UP_DOWN_info);
% ~ismember(varnames1,varnames)
for i = 1:length(varnames)
    varname = varnames{i};

    if contains(varname,'spindle')
        ipsi_name = ['ipsi_' varname '_UP'];
        contra_name = ['contra_' varname '_UP'];

        UP_DOWN_info.(ipsi_name) = [data_struct.(varname){1}(:,1); data_struct.(varname){2}(:,2)]';
        UP_DOWN_info.(contra_name) = [data_struct.(varname){1}(:,2); data_struct.(varname){2}(:,1)]';
    elseif (contains(varname,'L_V1') | contains(varname,'R_V1')) 

    else
        UP_DOWN_info.([varname,'_UP']) = [data_struct.(varname){1} data_struct.(varname){2}];
    end
end


% ripple_<hemi>_<region>_MUA_<type>
regions = {'V1'};
suffixes = {
    'ripple_%s_MUA_cumulative'
    'normalised_ripple_%s_MUA_cumulative'
    'first_ripple_%s_MUA_peak'
    'last_ripple_%s_MUA_peak'

    '%s_MUA_cumulative'
    'normalised_%s_MUA_cumulative'
    'first_half_%s_MUA_cumulative'
    'second_half_%s_MUA_cumulative'

    'first_ripple_%s_MUA_cumulative'
    'last_ripple_%s_MUA_cumulative'


    % 'ripple_%s_MUA_cumulative'
};


for i = 1:length(regions)
    region = regions{i};
    for j = 1:length(suffixes)
        suffix_template = suffixes{j};

        % Build field names for L and R
        L_field = sprintf(suffix_template, ['L_' region]);
        R_field = sprintf(suffix_template, ['R_' region]);

        % Build output variable name
        var_base = strrep(suffix_template, '%s', region);
        ipsi_name = ['ipsi_' var_base '_UP'];
        contra_name = ['contra_' var_base '_UP'];

        % Assign
        UP_DOWN_info.(ipsi_name) = [data_struct.(L_field){1}, data_struct.(R_field){2}];
        UP_DOWN_info.(contra_name) = [data_struct.(R_field){1}, data_struct.(L_field){2}];
    end
end



% Cleanup
clear hemi_labels region_labels varnames data_struct
clear h r v f i C ia ripple_idx
clear hemi reg suffix fname field
clear nprobe nevent up_idx up_duration npoints half_idx cum_sum trace
clear UP_indices ripple_norm_dur ripple_data
clear event_index ripples_index onset offset duration vals MUA_cum MUA_peak
clear idx1 idx2
clear ripple_regions ripple_types global_regions global_prefixes
clear region prefix base_name L_field R_field ipsi_name contra_name

% merged_event_info

%%%%% Grab DOWN/UP info
varnames = {
    'ipsi_Delta_peaks_zscore_UD'
    'ipsi_Delta_peaks_zscore_DU'
    'SWpeakmag_UD'
    'SWpeakmag_DU'
    'contra_Delta_peaks_zscore_UD'
    'contra_Delta_peaks_zscore_DU'
    };

% Loop through each variable name
for i = 1:length(varnames)
    varname = varnames{i};
    temp =eval(varname);
    UP_DOWN_info.(varname) = temp;
end
UP_DOWN_info.UP_duration = [event_info(1).UP_duration;event_info(2).UP_duration];
UP_DOWN_info.previous_DOWN_duration = [event_info(1).previous_DOWN_duration;event_info(2).previous_DOWN_duration];
UP_DOWN_info.next_DOWN_duration = [event_info(1).next_DOWN_duration;event_info(2).next_DOWN_duration];

%%%%%%%%%%%%%
%%%%%%%%%%%%% Early UP and DOWN log odds

for nprobe = 1:2
    UP_indices = probability_psth_whole(nprobe).UP_all_index;

    %%%% Late UP Log odds (200ms)
    tidx = KDE_reactivation_PSTH(1).tvec > -0.1 & KDE_reactivation_PSTH(1).tvec < 0; 

    UP_DOWN_info.late_UP_log_odds{nprobe}=mean(DOWN_HPC_log_odds{nprobe}(next_down_idx{nprobe}(UP_indices),tidx),2,'omitnan');
    UP_DOWN_info.late_UP_V1_log_odds{nprobe}=mean(DOWN_V1_log_odds{nprobe}(next_down_idx{nprobe}(UP_indices),tidx),2,'omitnan');

%     prev_down_idx{nprobe}(UP_indices)
    
    %%%% Early UP Log odds (200ms)
    tidx = KDE_reactivation_PSTH(1).tvec > 0 & KDE_reactivation_PSTH(1).tvec < 0.1;
    UP_DOWN_info.early_UP_log_odds{nprobe}=mean(UP_HPC_log_odds{nprobe}(UP_indices,tidx),2,'omitnan');
    UP_DOWN_info.early_UP_V1_log_odds{nprobe}=mean(UP_V1_log_odds{nprobe}(UP_indices,tidx),2,'omitnan');

    %%%% DOWN Log odds (200ms) (following current UP)
    tidx = KDE_reactivation_PSTH(1).tvec > 0 & KDE_reactivation_PSTH(1).tvec < 0.1;
    UP_DOWN_info.DOWN_log_odds{nprobe}=mean(DOWN_HPC_log_odds{nprobe}(next_down_idx{nprobe}(UP_indices),tidx),2,'omitnan');
    UP_DOWN_info.DOWN_V1_log_odds{nprobe}=mean(DOWN_V1_log_odds{nprobe}(next_down_idx{nprobe}(UP_indices),tidx),2,'omitnan');

    %%%% Previous late UP Log odds (200ms)
    tidx = KDE_reactivation_PSTH(1).tvec > -0.1 & KDE_reactivation_PSTH(1).tvec < 0; 
    
    UP_DOWN_info.previous_late_UP_log_odds{nprobe}=mean(DOWN_HPC_log_odds{nprobe}(prev_down_idx{nprobe}(UP_indices),tidx),2,'omitnan');
    UP_DOWN_info.previous_late_UP_V1_log_odds{nprobe}=mean(DOWN_V1_log_odds{nprobe}(prev_down_idx{nprobe}(UP_indices),tidx),2,'omitnan');

    %%%% Next early UP Log odds (200ms)
    tidx = KDE_reactivation_PSTH(1).tvec > 0 & KDE_reactivation_PSTH(1).tvec < 0.1;
    if UP_indices(end)==size(UP_HPC_log_odds{nprobe},1)        
        UP_DOWN_info.next_early_UP_log_odds{nprobe}=[mean(UP_HPC_log_odds{nprobe}(UP_indices(1:end-1)+1,tidx),2,'omitnan'); nan];
        UP_DOWN_info.next_early_UP_V1_log_odds{nprobe}=[mean(UP_V1_log_odds{nprobe}(UP_indices(1:end-1)+1,tidx),2,'omitnan'); nan];
    else
        UP_DOWN_info.next_early_UP_log_odds{nprobe}=[mean(UP_HPC_log_odds{nprobe}(UP_indices+1,tidx),2,'omitnan')];
        UP_DOWN_info.next_early_UP_V1_log_odds{nprobe}=[mean(UP_V1_log_odds{nprobe}(UP_indices+1,tidx),2,'omitnan')];
    end

    UP_indices = probability_psth_whole(nprobe).UP_all_index;
    UP_DOWN_info.edge_event_idx{nprobe} = find(ismember(UP_indices,[find(diff([slow_waves_all(nprobe).UP_session_count])>0); find(diff([slow_waves_all(nprobe).UP_session_count])>0)+1]));
end
varnames =fieldnames(UP_DOWN_info);

for i = 1:length(varnames)
    varname = varnames{i};
    if iscell(UP_DOWN_info.(varname))
       
        UP_DOWN_info.(varname) = vertcat(UP_DOWN_info.(varname){:})';
    end
end
%%%%%%%%%% Spindles
%%%%%%%%%%%%%%%%%%%%
UP_DOWN_info.ipsi_spindles_UP = [nansum(spindle_probability_psth_whole(1).L_spindles_UP(:,50:60)')'>0; nansum(spindle_probability_psth_whole(2).R_spindles_UP(:,50:60)')'>0];
UP_DOWN_info.contra_spindles_UP = [nansum(spindle_probability_psth_whole(1).R_spindles_UP(:,50:60)')'>0; nansum(spindle_probability_psth_whole(2).L_spindles_UP(:,50:60)')'>0];

UP_DOWN_info.ipsi_spindles_DOWN = [nansum(spindle_probability_psth_whole(1).L_spindles_DOWN(:,50:60)')'>0; nansum(spindle_probability_psth_whole(2).R_spindles_DOWN(:,50:60)')'>0];
UP_DOWN_info.contra_spindles_DOWN = [nansum(spindle_probability_psth_whole(1).R_spindles_DOWN(:,50:60)')'>0; nansum(spindle_probability_psth_whole(2).L_spindles_DOWN(:,50:60)')'>0];

save(fullfile(analysis_folder,'V1-HPC sleep reactivation','UP_DOWN_info.mat'),'UP_DOWN_info');

save(fullfile(analysis_folder,'V1-HPC sleep reactivation','UP_DOWN_info_all.mat'),'UP_DOWN_info');
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','UP_DOWN_info.mat'),'UP_DOWN_info');
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction\survival analysis','UP_DOWN_info.mat'),'UP_DOWN_info');

%% Predict Ripples from DOWN-UP info

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'));
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_combined.mat'));
probability_psth_whole = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_probability_whole.mat'));
spindle_probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline_combined.mat'));
probability_psth_whole_baseline = probability;

load(fullfile(analysis_folder,'V1-HPC sleep reactivation','UP_DOWN_info.mat'),'UP_DOWN_info');

%%%%% UP info
event_times = merged_event_info.UP_ints;
hemisphere_id = merged_event_info.UP_hemisphere_id;
% lag_diff = merged_event_info.UP_lag_diff;
% group_id = merged_event_info.UP_group_id;
all_overlap_idx = merged_event_info.UP_overlap_idx_all{end};
non_overlap_idx = merged_event_info.UP_non_overlap_idx{end};
lags =merged_event_info.UP_lags_all{end};

UP_session_count = [slow_waves_all(1).UP_session_count(probability(1).UP_all_index); slow_waves_all(2).UP_session_count(probability(2).UP_all_index)];
subject_id = str2double(cellstr(slow_waves_all(1).subject(UP_session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

%%%%%%%%%%%
%%%%% V1 MUA
ipsi_V1_MUA = [PSTH_MUA(1).L_V1_UP; PSTH_MUA(2).R_V1_UP];
contra_V1_MUA = [PSTH_MUA(1).R_V1_UP; PSTH_MUA(2).L_V1_UP];

% ipsi_V1_MUA_baseline = [PSTH_MUA_baseline(1).L_V1_UP; PSTH_MUA_baseline(2).R_V1_UP];
% contra_V1_MUA_baseline = [PSTH_MUA_baseline(1).R_V1_UP; PSTH_MUA_baseline(2).L_V1_UP];

%%%%% HPC MUA
ipsi_HPC_MUA = [PSTH_MUA(1).L_HPC_UP; PSTH_MUA(2).R_HPC_UP];
contra_HPC_MUA = [PSTH_MUA(1).R_HPC_UP; PSTH_MUA(2).L_HPC_UP];
ipsi_HPC_MUA = (ipsi_HPC_MUA+contra_HPC_MUA)/2;
% ipsi_HPC_MUA_baseline = [PSTH_MUA_baseline(1).L_HPC_UP; PSTH_MUA_baseline(2).R_HPC_UP];
% contra_HPC_MUA_baseline = [PSTH_MUA_baseline(1).R_HPC_UP; PSTH_MUA_baseline(2).L_HPC_UP];

%%%%%%%%%%% ripples
ipsi_probability = [probability_psth_whole(1).ripples_UP; probability_psth_whole(2).ripples_UP];
% contra_probability = [probability_psth_whole(1).R_ripples_UP; probability_psth_whole(2).L_ripples_UP];

%%%%%%%%%% Predict HPC MUA and ripples during DOWN-UP transition based on
%%%%%%%%%% V1 DOWN UP bilateral synchrony and magnitude
%%%%%%%%% 0 to 100ms ripples
% 
index = all_overlap_idx(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);
output = predict_bilateral_ripples_by_DOWN_UP_synchrony(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:),ipsi_probability(index,:),UP_DOWN_info,...
    'DOWN_UP_index',index,'DOWN_UP_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),'session_id',UP_session_count(index));
% output = predict_ripples_by_DOWN_UP_synchrony(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,'DOWN_UP_index',index,'DOWN_UP_lag',abs(lags(lag_index)'),'subject_id',subject_id(index));
save(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','DU_synchrony_Bilateral_Ripples_output.mat'),'output');


% 
% output = predict_ripples_by_DOWN_UP(ipsi_V1_MUA, contra_V1_MUA, ipsi_HPC_MUA, contra_HPC_MUA,ipsi_probability, contra_probability,UP_DOWN_info,'subject_id',subject_id);

output = predict_bilateral_ripples_by_DOWN_UP(ipsi_V1_MUA, contra_V1_MUA, ipsi_HPC_MUA,ipsi_probability,UP_DOWN_info,'subject_id',subject_id,'session_id',UP_session_count);
save(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','DU_Bilateral_Ripples_output.mat'),'output');


%%%%%%%%% 0 to 200ms ripples
index = all_overlap_idx(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);
output = predict_bilateral_ripples_by_DOWN_UP_synchrony(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:),ipsi_probability(index,:),UP_DOWN_info,...
    'DOWN_UP_index',index,'DOWN_UP_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),'session_id',UP_session_count(index),'time_window', 0.2);
% output = predict_ripples_by_DOWN_UP_synchrony(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,'DOWN_UP_index',index,'DOWN_UP_lag',abs(lags(lag_index)'),'subject_id',subject_id(index));
save(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','DU_synchrony_Bilateral_Ripples_output_200ms.mat'),'output');


output = predict_bilateral_ripples_by_DOWN_UP(ipsi_V1_MUA, contra_V1_MUA, ipsi_HPC_MUA,ipsi_probability,UP_DOWN_info,'subject_id',subject_id,'session_id',UP_session_count,'time_window', 0.2);
save(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','DU_Bilateral_Ripples_output_200ms.mat'),'output');

% 
% %%%%%%%%%%%%%%%%%%% Plotting
index = all_overlap_idx(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);

load(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','DU_synchrony_Bilateral_Ripples_output.mat'),'output');
% 

output = predict_bilateral_ripples_by_DOWN_UP_synchrony(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:),ipsi_probability(index,:),UP_DOWN_info,...
    'DOWN_UP_index',index,'DOWN_UP_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),'session_id',UP_session_count(index),'output',output,'plot_option',1);
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (150ms windows)'),[],'ContentType','image')
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','mixed effect regression (150ms lag windows)'),[],'SVG_option',1)


load(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','DU_synchrony_Bilateral_Ripples_output_200ms.mat'),'output');
output = predict_bilateral_ripples_by_DOWN_UP_synchrony(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:),ipsi_probability(index,:),UP_DOWN_info,...
    'DOWN_UP_index',index,'DOWN_UP_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),'session_id',UP_session_count(index),'output',output,'plot_option',1,'time_window', 0.2);
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','mixed effect regression (150ms lag windows 200ms)'),[],'SVG_option',1)



load(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','DU_Bilateral_Ripples_output.mat'),'output');
output = predict_bilateral_ripples_by_DOWN_UP(ipsi_V1_MUA, contra_V1_MUA, ipsi_HPC_MUA,ipsi_probability,UP_DOWN_info,'subject_id',subject_id,'session_id',UP_session_count,'output',output,'plot_option',1);
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (full windows)'),[],'ContentType','image')
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','mixed effect regression (full windows)'),[],'SVG_option',1)


load(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','DU_Bilateral_Ripples_output_200ms.mat'),'output');
output = predict_bilateral_ripples_by_DOWN_UP(ipsi_V1_MUA, contra_V1_MUA, ipsi_HPC_MUA,ipsi_probability,UP_DOWN_info,'subject_id',subject_id,'session_id',UP_session_count,'output',output,'plot_option',1,'time_window', 0.2);
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (full windows)'),[],'ContentType','image')
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','mixed effect regression (full windows 200ms)'),[],'SVG_option',1)


%% UP DOWN


%%%%% DOWN info
event_times = merged_event_info.DOWN_ints;
hemisphere_id = merged_event_info.DOWN_hemisphere_id;
% lag_diff = merged_event_info.UP_lag_diff;
% group_id = merged_event_info.UP_group_id;
all_overlap_idx = merged_event_info.DOWN_overlap_idx_all{end};
non_overlap_idx = merged_event_info.DOWN_non_overlap_idx{end};
lags =merged_event_info.DOWN_lags_all{end};

DOWN_session_count = [slow_waves_all(1).DOWN_session_count(probability(1).DOWN_all_index); slow_waves_all(2).DOWN_session_count(probability(2).DOWN_all_index)];
subject_id = str2double(cellstr(slow_waves_all(1).subject(DOWN_session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

%%%%%%%%%%%%%% Grab ripple info during UP DOWN
%%%%%%%%%%%
%%%%% V1 MUA
ipsi_V1_MUA = [PSTH_MUA(1).L_V1_DOWN; PSTH_MUA(2).R_V1_DOWN];
contra_V1_MUA = [PSTH_MUA(1).R_V1_DOWN; PSTH_MUA(2).L_V1_DOWN];

ipsi_V1_MUA_baseline = [PSTH_MUA_baseline(1).L_V1_DOWN; PSTH_MUA_baseline(2).R_V1_DOWN];
contra_V1_MUA_baseline = [PSTH_MUA_baseline(1).R_V1_DOWN; PSTH_MUA_baseline(2).L_V1_DOWN];

%%%%% HPC MUA
ipsi_HPC_MUA = [PSTH_MUA(1).L_HPC_DOWN; PSTH_MUA(2).R_HPC_DOWN];
contra_HPC_MUA = [PSTH_MUA(1).R_HPC_DOWN; PSTH_MUA(2).L_HPC_DOWN];
ipsi_HPC_MUA = (ipsi_HPC_MUA+contra_HPC_MUA)/2;
% ipsi_HPC_MUA_baseline = [PSTH_MUA_baseline(1).L_HPC_DOWN; PSTH_MUA_baseline(2).R_HPC_DOWN];
% contra_HPC_MUA_baseline = [PSTH_MUA_baseline(1).R_HPC_DOWN; PSTH_MUA_baseline(2).L_HPC_DOWN];

%%%%%%%%%%% ripples
ipsi_probability = [probability_psth_whole(1).ripples_DOWN; probability_psth_whole(2).ripples_DOWN];
% contra_probability = [probability_psth_whole(1).R_ripples_DOWN; probability_psth_whole(2).L_ripples_DOWN];

% ipsi_probability_baseline = [probability_psth_whole_baseline(1).L_ripples_DOWN; probability_psth_whole_baseline(2).R_ripples_DOWN];
% contra_probability_baseline = [probability_psth_whole_baseline(1).R_ripples_DOWN; probability_psth_whole_baseline(2).L_ripples_DOWN];

% 
% %%%%%%%%%%%%%%%% -0.1 to 0 ripples
index = 1:size(ipsi_V1_MUA,1);
% lag_index = (abs(lags)<=2);
output = predict_UP_DOWN_V1_MUA_by_bilateral_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:),ipsi_probability(index,:),UP_DOWN_info,...
    'subject_id',subject_id(index),'session_id',DOWN_session_count(index),'time_window', -0.1);
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output.mat'),'output');
save(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','Ripples_V1depression_output.mat'),'output');

% predict_UP_DOWN_synchrony_by_bilateral_ripples
index = all_overlap_idx(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);
output = predict_UP_DOWN_synchrony_by_bilateral_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:),ipsi_probability(index,:),UP_DOWN_info,...
    'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),'session_id',DOWN_session_count(index),'time_window', -0.1);
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output.mat'),'output');
save(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','Ripples_V1synchrony_output.mat'),'output');


% %%%%%%%%%%%%%%%% -0.2 to 0 ripples
index = 1:size(ipsi_V1_MUA,1);
% lag_index = (abs(lags)<=2);
output = predict_UP_DOWN_V1_MUA_by_bilateral_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:),ipsi_probability(index,:),UP_DOWN_info,...
    'subject_id',subject_id(index),'session_id',DOWN_session_count(index),'time_window', -0.2);
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output.mat'),'output');
save(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','Ripples_V1depression_output_200ms.mat'),'output');

% predict_UP_DOWN_synchrony_by_bilateral_ripples
index = all_overlap_idx(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);
output = predict_UP_DOWN_synchrony_by_bilateral_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:),ipsi_probability(index,:),UP_DOWN_info,...
    'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),'session_id',DOWN_session_count(index),'time_window', -0.2);
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output.mat'),'output');
save(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','Ripples_V1synchrony_output_200ms.mat'),'output');



% 
% 
% %%%%%%%%%%%% Plot
index = 1:size(ipsi_V1_MUA,1);
% lag_index = (abs(lags)<=2);
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output.mat'),'output');
load(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','Ripples_V1depression_output.mat'),'output');

% output = predict_UP_DOWN_V1_MUA_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,...
%     'output',output,'plot_option',1,'subject_id',subject_id(index));
output = predict_UP_DOWN_V1_MUA_by_bilateral_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:),ipsi_probability(index,:),UP_DOWN_info,...
    'output',output,'plot_option',1,'subject_id',subject_id(index),'session_id',DOWN_session_count(index));

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','mixed effect regression (full windows)'),[],'SVG_option',1)


index = 1:size(ipsi_V1_MUA,1);
% lag_index = (abs(lags)<=2);
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output.mat'),'output');
load(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','Ripples_V1depression_output_200ms.mat'),'output');

% output = predict_UP_DOWN_V1_MUA_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,...
%     'output',output,'plot_option',1,'subject_id',subject_id(index));
output = predict_UP_DOWN_V1_MUA_by_bilateral_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:),ipsi_probability(index,:),UP_DOWN_info,...
    'output',output,'plot_option',1,'subject_id',subject_id(index),'session_id',DOWN_session_count(index),'time_window',-0.2);

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','mixed effect regression (full windows 200ms)'),[],'SVG_option',1)





index = all_overlap_idx(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);
load(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','Ripples_V1synchrony_output.mat'),'output');
% output = predict_UP_DOWN_synchrony_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,...
%     'output',output,'plot_option',1,'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index));
output = predict_UP_DOWN_synchrony_by_bilateral_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:),ipsi_probability(index,:),UP_DOWN_info,...
    'output',output,'plot_option',1,'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),'session_id',DOWN_session_count(index));

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','mixed effect regression (150ms lag windows)'),[],'SVG_option',1)


load(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','Ripples_V1synchrony_output_200ms.mat'),'output');
% output = predict_UP_DOWN_synchrony_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,...
%     'output',output,'plot_option',1,'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index));
output = predict_UP_DOWN_synchrony_by_bilateral_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:),ipsi_probability(index,:),UP_DOWN_info,...
    'output',output,'plot_option',1,'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),'session_id',DOWN_session_count(index),'time_window',-0.2);

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','mixed effect regression (150ms lag windows 200ms)'),[],'SVG_option',1)


% %%%%%%%%%%%%%%%% -0.1 to 0 ripples (control)

UP_DOWN_info.SWpeakmag_UD1 = UP_DOWN_info.SWpeakmag_UD ;
UP_DOWN_info.SWpeakmag_UD = UP_DOWN_info.SWpeakmag_DU;

all_overlap_idx_UP = merged_event_info.UP_overlap_idx_all{end};
lags =merged_event_info.UP_lags_all{end};
index = all_overlap_idx_UP(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);
output = predict_UP_DOWN_synchrony_by_bilateral_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:),ipsi_probability(index,:),UP_DOWN_info,...
    'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),'session_id',DOWN_session_count(index),'time_window', -0.1);
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output.mat'),'output');
save(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','Ripples_V1synchrony_output_previous_DU.mat'),'output');

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','mixed effect regression (150ms lag windows)'),[],'SVG_option',1)



index = 1:size(ipsi_V1_MUA,1);
% lag_index = (abs(lags)<=2);
output = predict_UP_DOWN_V1_MUA_by_bilateral_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:),ipsi_probability(index,:),UP_DOWN_info,...
    'subject_id',subject_id(index),'session_id',DOWN_session_count(index),'time_window', -0.1);
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output.mat'),'output');
save(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','Ripples_V1depression_output_previous_DU.mat'),'output');
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction\UP_DOWN_ripples_lme','mixed effect regression (full windows)'),[],'SVG_option',1)

%% Generate CSV files

model_indices = [1 2 7 8 17 18 5 6 11 12 21 22 3 4 9 10 19 20];
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression 250ms ripples (full windows)'));
load(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_Ripples_250ms_output.mat'),'output');
% load('Ripples_V1synchrony_output.mat');  % Load your 'output' variable
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'DU_ripples_Summary.csv');


model_indices = [1 2 7 8 17 18 5 6 11 12 21 22 3 4 9 10 19 20];
load(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_Ripples_output.mat'),'output');
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (full windows)'));
% load('Ripples_V1synchrony_output.mat');  % Load your 'output' variable
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'DU_ripples_Summary.csv');


model_indices = [1 2 9 10 25 26 3 4 11 12 27 28 7 8 15 16 31 32 5 6 13 14 29 30];
load(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_synchrony_Ripples_output.mat'),'output');
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (150ms windows)'));
% load('Ripples_V1synchrony_output.mat');  % Load your 'output' variable
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'DU_synchrony Ripples_Summary.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output_based_on_UP.mat'),'output');
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (150ms windows UP lag control)'));
model_indices = [13 14 29 30 5 6 2 3 25 26 18 19 21 22 32 33 7 23 15];
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'Ripples_V1synchrony_Summary.csv');


model_indices = [13 14 29 30 5 6 2 3 25 26 18 19 21 22 32 33 7 23 15];
load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output.mat'),'output');
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (150ms windows)'));
% load('Ripples_V1synchrony_output.mat');  % Load your 'output' variable
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'Ripples_V1synchrony_Summary.csv');

% [10 11 13 14 21 22 5 6 2 3 17 18 7 15]
model_indices = [10 11 13 14 21 22 5 6 2 3 17 18 7 15];
load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output_based_on_UP.mat'),'output');
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (full windows UP lag control)'));
% load('Ripples_V1synchrony_output.mat');  % Load your 'output' variable
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'Ripples_V1depression_Summary.csv');


model_indices = [10 11 13 14 21 22 5 6 2 3 17 18 7 15];
load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output.mat'),'output');
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (full windows)'));
% load('Ripples_V1synchrony_output.mat');  % Load your 'output' variable
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'Ripples_V1depression_Summary.csv');



load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_amp_phase_Ripples_output.mat'),'output');
model_indices = 1:length(output(1).model);
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression SO spindle phase power'));
% load('Ripples_V1synchrony_output.mat');  % Load your 'output' variable
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'SO_spindles_amp_phase_Ripples_Summary.csv');


% output = predict_UP_DOWN_synchrony_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,...
%     plot_option',1,'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index));
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output_based_on_UP.mat'),'output');

%% Effect of ripples and cumulative HPC activities  on UP survival probability

load(fullfile(analysis_folder,'V1-HPC sleep reactivation','UP_DOWN_info.mat'),'UP_DOWN_info');

UP_session_count = [slow_waves_all(1).UP_session_count(probability(1).UP_all_index); slow_waves_all(2).UP_session_count(probability(2).UP_all_index)];
subject_id = str2double(cellstr(slow_waves_all(1).subject(UP_session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);
UP_DOWN_info.subject_id = subject_id;

%%%%%%%%%%%%% Ripple power predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.last_ripples_power_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple peak power and last ripple to UP-DOWN transition','subject_id',subject_id...
);
save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_power_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.last_ripples_power_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple peak power and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'subject_id',subject_id...
);

output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_power_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])



%%%%%%%%%%%%% Ripple duration predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.last_ripples_duration_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple duration and last ripple to UP-DOWN transition','subject_id',subject_id...
);
save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_duration_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option')



% scatter([UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP], ...
%     [UP_DOWN_info.ipsi_last_ripples_power_UP, UP_DOWN_info.contra_last_ripples_power_UP],...
%     'filled','MarkerFaceAlpha','0.02')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% ylim([0 25])
% xlim([0 0.3])


% histogram(UP_DOWN_info.ipsi_last_ripples_power_UP)
% plot_UP_survival_probability( ...
%     {UP_DOWN_info.ipsi_last_ripples_power_UP, UP_DOWN_info.contra_last_ripples_power_UP}, ...
%     {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
%     {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',1,'subject_id',subject_id, ...
%     'event_option',[],'title_name', 'Last ripple peak power and last ripple to UP-DOWN transition (1 ripple)'...
% );

%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ripple MUA peak predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.last_ripple_HPC_MUA_peak_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id,'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple peak MUA and last ripple to UP-DOWN transition'...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_HPC_MUA_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.last_ripple_HPC_MUA_peak_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple peak MUA and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'subject_id',subject_id...
);

output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_HPC_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])

%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Last ripple V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_last_ripple_V1_MUA_peak_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Last ripple ipsi V1 peak MUA and last ripple to UP-DOWN transition'...
);

output = plot_UP_survival_probability( ...
    {UP_DOWN_info.contra_last_ripple_V1_MUA_peak_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Last ripple contra V1 peak MUA and last ripple to UP-DOWN transition'...
);

output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_last_ripple_V1_MUA_peak_UP - UP_DOWN_info.contra_last_ripple_V1_MUA_peak_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Last ripple ipsi contra V1 peak MUA diff and last ripple to UP-DOWN transition'...
);
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])



%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalised cumulative ripple HPC MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.normalised_ripple_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative ripple activity and last ripple to UP-DOWN transition','timebin', 0.015...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_normalised_HPC_MUA_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.normalised_ripple_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Normalised cumulative ripple activity and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'subject_id',subject_id...
);

output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','normalised_cumulative_ripple_HPC_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalised cumulative ripple V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_ripple_V1_MUA_cumulative_UP,UP_DOWN_info.contra_normalised_ripple_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 ripple activity and last ripple to UP-DOWN transition','timebin', 0.015...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_normalised_HPC_MUA_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_ripple_V1_MUA_cumulative_UP,UP_DOWN_info.contra_normalised_ripple_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 ripple activity and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);

for i = 1:length(output)
output(i).p_shuffled = output_shuffled(i).p;
output(i).b_shuffled = output_shuffled(i).b;
end

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','normalised_cumulative_ripple_V1_MUA_survival.mat'),'output');
% load(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','normalised_cumulative_ripple_V1_MUA_survival.mat'),'output');

% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalised cumulative ripple V1 MUA ipsi contra diff predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_ripple_V1_MUA_cumulative_UP - UP_DOWN_info.contra_normalised_ripple_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative ipsi contra V1 ripple activity diff and last ripple to UP-DOWN transition','timebin', 0.015...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_normalised_HPC_MUA_survival.mat'),'output');

output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_ripple_V1_MUA_cumulative_UP - UP_DOWN_info.contra_normalised_ripple_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative ipsi contra V1 ripple activity diff and last ripple to UP-DOWN transition (Shuffled)','timebin', 0.015...
);

output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','normalised_cumulative_ripple_V1_MUA_diff_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cumulative ripple HPC MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ripple_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Cumulative ripple activity and last ripple to UP-DOWN transition','timebin', 0.015...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_normalised_HPC_MUA_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ripple_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Cumulative ripple activity and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'subject_id',subject_id...
);

output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','cumulative_ripple_HPC_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cumulative ripple V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_ripple_V1_MUA_cumulative_UP,UP_DOWN_info.contra_ripple_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Cumulative ripple V1 activity and last ripple to UP-DOWN transition','timebin', 0.015...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_normalised_HPC_MUA_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_ripple_V1_MUA_cumulative_UP,UP_DOWN_info.contra_ripple_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Cumulative ripple V1 activity and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);

for i = 1:length(output)
output(i).p_shuffled = output_shuffled(i).p;
output(i).b_shuffled = output_shuffled(i).b;
end


save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','cumulative_ripple_V1_MUA_survival.mat'),'output');
load(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','cumulative_ripple_V1_MUA_survival.mat'),'output');

% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])






%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%% Normalised cumulative V1 MUA predicts UP probability
% Normalised cumulative V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_V1_MUA_cumulative_UP,UP_DOWN_info.contra_normalised_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative UP V1 activity and last ripple to UP-DOWN transition','timebin', 0.015...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_normalised_HPC_MUA_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_V1_MUA_cumulative_UP,UP_DOWN_info.contra_normalised_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative UP V1 activity and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);

for i = 1:length(output)
output(i).p_shuffled = output_shuffled(i).p;
output(i).b_shuffled = output_shuffled(i).b;
end


save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','normalised_cumulative_UP_V1_MUA_survival.mat'),'output');
% load(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','normalised_cumulative_UP_V1_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])



%%%%%%%%%%%%%%%% Normalised cumulative HPC MUA predicts UP probability
% Normalised cumulative V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.normalised_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative UP HPC activity and last ripple to UP-DOWN transition','timebin', 0.015...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_normalised_HPC_MUA_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.normalised_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative UP HPC activity and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);

output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','normalised_cumulative_UP_HPC_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])



%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%% 1st half Normalised cumulative HPC MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.first_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 1st half and last ripple to UP-DOWN transition','timebin', 0.015...
);

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.first_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 1st half and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);


output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','1st_half_normalised_cumulative_UP_HPC_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])

%%%%%%%%%%%%%%%% 1st half Normalised cumulative V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP, UP_DOWN_info.contra_first_half_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 1st half and last ripple to UP-DOWN transition','timebin', 0.015...
);

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP, UP_DOWN_info.contra_first_half_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 1st half and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);


for i = 1:length(output)
output(i).p_shuffled = output_shuffled(i).p;
output(i).b_shuffled = output_shuffled(i).b;
end


save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','1st_half_normalised_cumulative_UP_V1_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%% 2nd half Normalised cumulative HPC MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.second_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 2nd half and last ripple to UP-DOWN transition','timebin', 0.015...
);

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.second_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 2nd half and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);


output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','2nd_half_normalised_cumulative_UP_HPC_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


%%%%%%%%%%%%%%%% 2nd half Normalised cumulative V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP, UP_DOWN_info.contra_second_half_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 2nd half and last ripple to UP-DOWN transition','timebin', 0.015...
);

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP, UP_DOWN_info.contra_second_half_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 2nd half and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);


for i = 1:length(output)
output(i).p_shuffled = output_shuffled(i).p;
output(i).b_shuffled = output_shuffled(i).b;
end


save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','2nd_half_normalised_cumulative_UP_V1_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])




%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%% 2nd half - 1st half Normalised cumulative HPC MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.second_half_HPC_MUA_cumulative_UP-UP_DOWN_info.first_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 2nd - 1st half and last ripple to UP-DOWN transition','timebin', 0.015...
);

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.second_half_HPC_MUA_cumulative_UP-UP_DOWN_info.first_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 2nd - 1st half and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);


output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','2nd_1st_half_normalised_cumulative_UP_HPC_MUA_diff_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


%%%%%%%%%%%%%%%% 2nd half Normalised cumulative V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP - UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP, UP_DOWN_info.contra_second_half_V1_MUA_cumulative_UP - UP_DOWN_info.contra_first_half_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 2nd - 1st half and last ripple to UP-DOWN transition','timebin', 0.015...
);

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP - UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP, UP_DOWN_info.contra_second_half_V1_MUA_cumulative_UP - UP_DOWN_info.contra_first_half_V1_MUA_cumulative_UP}, ....
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 2nd - 1st half and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);


for i = 1:length(output)
output(i).p_shuffled = output_shuffled(i).p;
output(i).b_shuffled = output_shuffled(i).b;
end


save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','2nd_1st_half_normalised_cumulative_UP_V1_MUA_diff_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


%%%%%%%%%%%%%%%%%%%%%%
%%%%%
% 
% UP_DOWN_info.ipsi_last_ripple_spindle_power_UP
%% Log odds and UP termination
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_log_odds_UP.mat'),'KDE_log_odds_UP')
UP_session_count = [slow_waves_all(1).UP_session_count(probability(1).UP_all_index); slow_waves_all(2).UP_session_count(probability(2).UP_all_index)];
subject_id = str2double(cellstr(slow_waves_all(1).subject(UP_session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);
UP_DOWN_info.subject_id = subject_id;


V1_direction = [-1*ones(1,length(probability_psth_whole(1).UP_all_index)) 1*ones(1,length(probability_psth_whole(2).UP_all_index))];
log_odds_threshold = prctile(UP_DOWN_info.last_ripples_log_odds_UP,0:10:100);
log_odds_prod = UP_DOWN_info.last_ripples_log_odds_UP .* UP_DOWN_info.last_ripples_V1_log_odds_UP;
log_odds_prod = UP_DOWN_info.first_ripples_log_odds_UP .* UP_DOWN_info.first_ripples_V1_log_odds_UP;


matching_index = zeros(1,length(UP_session_count));
non_matching_index = zeros(1,length(UP_session_count));

HPC_reactivation_index=[];HPC_non_reactivation_index=[];
for n = 1:max(UP_session_count)
    %     matching_index = [matching_index ];

    session_index = find(UP_session_count == n);

    HPC_threhold = prctile(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n),[25 75]) ;

    HPC_reactivation_index = [HPC_reactivation_index;...
        session_index(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) > HPC_threhold(2) |UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) < HPC_threhold(1))];

    HPC_non_reactivation_index = [HPC_non_reactivation_index;...
        session_index(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) > HPC_threhold(1) &UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) < HPC_threhold(2))];

    V1_threshold = prctile(UP_DOWN_info.last_ripples_V1_log_odds_UP( UP_session_count == n),[50 50]) ;
    %     V1_threshold = median(UP_DOWN_info.last_ripples_V1_log_odds_UP( UP_session_count == n),'omitnan') ;

%     coherent_index = [find(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) > HPC_threhold(2) & UP_DOWN_info.last_ripples_V1_log_odds_UP( UP_session_count == n) > V1_threshold(2))'
%         find(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) < HPC_threhold(1) & UP_DOWN_info.last_ripples_V1_log_odds_UP( UP_session_count == n) < V1_threshold(1))'];
% 
    coherent_index = [find(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) > HPC_threhold(2) & V1_direction( UP_session_count == n) > 0)'
        find(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) < HPC_threhold(1) & V1_direction( UP_session_count == n) <0)'];
    matching_index(session_index(coherent_index))=1;
    
% 
%     incoherent_index = [find(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) > HPC_threhold(2) & UP_DOWN_info.last_ripples_V1_log_odds_UP( UP_session_count == n) < V1_threshold(2))'
%         find(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) < HPC_threhold(1) & UP_DOWN_info.last_ripples_V1_log_odds_UP( UP_session_count == n) > V1_threshold(1))'];

    incoherent_index = [find(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) > HPC_threhold(2) & V1_direction( UP_session_count == n) < 0)'
        find(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) < HPC_threhold(1) & V1_direction( UP_session_count == n) >0)'];
    non_matching_index(session_index(incoherent_index))=1;



% %     HPC_threhold = prctile(UP_DOWN_info.late_UP_log_odds( UP_session_count == n),[25 75]) ;
%     V1_threshold = prctile(UP_DOWN_info.late_UP_V1_log_odds( UP_session_count == n),[50 50]) ;
% %     coherent_index = [find(UP_DOWN_info.late_UP_log_odds( UP_session_count == n) > HPC_threhold(2) & V1_direction( UP_session_count == n) > 0)'
% %         find(UP_DOWN_info.late_UP_log_odds( UP_session_count == n) < HPC_threhold(1) & V1_direction( UP_session_count == n) <0)'];
%         coherent_index = [find(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) > HPC_threhold(2) & UP_DOWN_info.late_UP_V1_log_odds( UP_session_count == n) > V1_threshold(2))'
%             find(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) < HPC_threhold(1) & UP_DOWN_info.late_UP_V1_log_odds( UP_session_count == n) < V1_threshold(1))'];
%     matching_index(session_index(coherent_index))=1;
% %     incoherent_index = [find(UP_DOWN_info.late_UP_log_odds( UP_session_count == n) > HPC_threhold(2) & V1_direction( UP_session_count == n) < 0)'
% %         find(UP_DOWN_info.late_UP_log_odds( UP_session_count == n) < HPC_threhold(1) & V1_direction( UP_session_count == n) >0)'];
%         incoherent_index = [find(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) > HPC_threhold(2) & UP_DOWN_info.late_UP_V1_log_odds( UP_session_count == n) < V1_threshold(2))'
%             find(UP_DOWN_info.last_ripples_log_odds_UP( UP_session_count == n) < HPC_threhold(1) & UP_DOWN_info.late_UP_V1_log_odds( UP_session_count == n) > V1_threshold(1))'];
% 
%     non_matching_index(session_index(incoherent_index))=1;
end

index = find((non_matching_index + matching_index)>0);
output = plot_UP_survival_probability( ...
    {matching_index(index)}, ...
    {UP_DOWN_info.time_from_last_ripples_UP(index)}, ...
    {UP_DOWN_info.ripple_counts_UP(index)},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple coherence and last ripple to UP-DOWN transition','event_option','binary','subject_id',subject_id(index)...
    );



output = plot_UP_survival_probability( ...
    {abs(UP_DOWN_info.last_ripples_log_odds_UP(HPC_reactivation_index))}, ...
    {UP_DOWN_info.time_from_last_ripples_UP(HPC_reactivation_index)}, ...
    {UP_DOWN_info.ripple_counts_UP(HPC_reactivation_index)},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple coherence (reactivation) and last ripple to UP-DOWN transition','event_option',[],'subject_id',subject_id(HPC_reactivation_index)...
    );

output = plot_UP_survival_probability( ...
    {UP_DOWN_info.last_ripples_log_odds_UP(HPC_non_reactivation_index).*V1_direction(HPC_non_reactivation_index)}, ...
    {UP_DOWN_info.time_from_last_ripples_UP(HPC_non_reactivation_index)}, ...
    {UP_DOWN_info.ripple_counts_UP(HPC_non_reactivation_index)},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple coherence (non reactivation) and last ripple to UP-DOWN transition','event_option',[],'subject_id',subject_id(HPC_non_reactivation_index)...
    );


output = plot_UP_survival_probability( ...
    {UP_DOWN_info.late_UP_V1_log_odds.*V1_direction}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple coherence and last ripple to UP-DOWN transition','event_option',[],'subject_id',subject_id...
    );

output = plot_UP_survival_probability( ...
    {UP_DOWN_info.last_ripples_log_odds_UP.*V1_direction}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple V1 log odds and last ripple to UP-DOWN transition','event_option',[],'subject_id',subject_id...
    );

output = plot_UP_survival_probability( ...
    {UP_DOWN_info.last_ripples_PRE_V1_log_odds_UP.*V1_direction}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple V1 log odds PRE and last ripple to UP-DOWN transition','event_option',[],'subject_id',subject_id...
    );

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_duration_survival.mat'),'output');
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])





for nprobe = 1:2
    UP_indices = probability_psth_whole(nprobe).UP_all_index;
    for nevent = 1:length(UP_indices)

        KDE_log_odds_UP(nprobe).HPC{UP_indices(nevent)}

        sum(KDE_log_odds_UP(nprobe).V1{UP_indices(nevent)})

        KDE_log_odds_UP(nprobe).V1{UP_indices(nevent)+1}

        KDE_log_odds_UP(nprobe).V1{UP_indices(nevent)-1}

    end
end







%% Relationship between neighbouring UP states
UP_session_count = [slow_waves_all(1).UP_session_count(probability(1).UP_all_index); slow_waves_all(2).UP_session_count(probability(2).UP_all_index)];
subject_id = str2double(cellstr(slow_waves_all(1).subject(UP_session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);
UP_DOWN_info.subject_id = subject_id;
% UP_DOWN_info.previous_late_UP_V1_log_odds
% UP_DOWN_info.early_UP_V1_log_odds

scatter(UP_DOWN_info.previous_late_UP_V1_log_odds,UP_DOWN_info.early_UP_V1_log_odds)


scatter(UP_DOWN_info.late_UP_V1_log_odds,UP_DOWN_info.early_UP_V1_log_odds)

% V1_direction = [zeros(1,length(probability_psth_whole(1).UP_all_index)) 1*ones(1,length(probability_psth_whole(2).UP_all_index))];
V1_direction = [1*ones(1,length(probability_psth_whole(1).UP_all_index)) 2*ones(1,length(probability_psth_whole(2).UP_all_index))];

low_index = UP_DOWN_info.last_ripples_power_UP < prctile(UP_DOWN_info.last_ripples_power_UP,[25]);
high_index = UP_DOWN_info.last_ripples_power_UP > prctile(UP_DOWN_info.last_ripples_power_UP,[75]);
tbl = table(V1_direction',normalize(UP_DOWN_info.last_ripples_power_UP)',normalize(UP_DOWN_info.next_early_UP_V1_log_odds)',normalize(UP_DOWN_info.late_UP_V1_log_odds)',normalize(UP_DOWN_info.last_ripples_V1_log_odds_UP)',...
    normalize(UP_DOWN_info.late_UP_log_odds)',normalize(UP_DOWN_info.last_ripples_log_odds_UP)',...
    normalize(UP_DOWN_info.early_UP_V1_log_odds)',normalize(UP_DOWN_info.previous_late_UP_V1_log_odds)',normalize(UP_DOWN_info.previous_late_UP_log_odds)',...
    subject_id,'VariableNames',{'V1Track','ripplePower','nextUP','lateUP','lastRippleV1','lateUPHPC','lastRippleHPC','earlyUP','previousUP','previousUPHPC','subject_id'});
tbl = tbl(UP_DOWN_info.ripple_counts_UP>0,:);

subplot(2,2,1)
scatter(UP_DOWN_info.early_UP_V1_log_odds(low_index),UP_DOWN_info.next_early_UP_log_odds(low_index));
% glme = fitglme(tbl, 'nextUP ~ lateUP*lateUPHPC*ripplePower + (1|subject_id)');
% glme.Coefficients

glme = fitglme(tbl, 'nextUP ~ lastRippleV1*lastRippleHPC*ripplePower + (1|subject_id)');
glme.Coefficients
glme.Rsquared.Adjusted

% glme = fitglme(tbl, 'nextUP ~ lastRippleV1*lastRippleHPC + ripplePower*lastRippleHPC + ripplePower*lastRippleV1 + (1|subject_id)');
% glme.Rsquared.Adjusted

glme = fitglme(tbl, 'nextUP ~ lateUP + earlyUP+ (1|subject_id)');
glme.Coefficients

glme = fitglme(tbl, 'nextUP ~ earlyUP*lastRippleHPC*ripplePower+ (1|subject_id)');
glme.Coefficients


glme = fitglme(tbl, 'nextUP ~ lateUP*lastRippleHPC*ripplePower+ (1|subject_id)');
glme.Coefficients
% 
% glme = fitglme(tbl, 'nextUP ~ V1Track*lastRippleHPC*ripplePower+ (1|subject_id)');
% glme.Coefficients


glme = fitglme(tbl, 'earlyUP ~ previousUP*previousUPHPC + (1|subject_id)');
glme.Coefficients
% tbl.V1Track

subplot(2,2,2)
scatter(UP_DOWN_info.late_UP_V1_log_odds(high_index),UP_DOWN_info.next_early_UP_log_odds(high_index));

%% UP DOWN ripple log odds distribution
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);
subject_id = subject_id(event_ids_first);
session_count = session_count(event_ids_first);

duration_thresholds = [0:0.1:1];
% duration_thresholds = [0:1/5:1];
xbins = duration_thresholds(1) + mean(diff(duration_thresholds))/2:mean(diff(duration_thresholds)):duration_thresholds(end)-mean(diff(duration_thresholds))/2;

mean_z_bias = mean(z_bias(bin_centers>0 & bin_centers<0.1,:),'omitnan');
mean_z_bias_PRE = mean(z_bias(bin_centers>-0.2 & bin_centers<0,:),'omitnan');
mean_z_bias_V1= mean(z_bias_V1(bin_centers>0 & bin_centers<0.2,:),'omitnan');

mean_z_bias_V1_PRE= mean(z_bias_V1(bin_centers>-0.2 & bin_centers<0,:),'omitnan');
mean_z_bias_V1_PRE_raw= mean(z_bias_V1_raw(bin_centers>-0.2 & bin_centers<0,:),'omitnan');

mean_z_bias_200ms = mean(z_bias(bin_centers>0 & bin_centers<0.2,:),'omitnan');

subject_id(subject_id>0)=1;
hemi_id_shuffled = ripple_info.normalised_UP_duration;

log_odds_normalised_duration_session = [];
V1_log_odds_normalised_duration_session = [];

% for n = 1:22
    V1_log_odds_normalised_duration = [];
    V1_log_odds_normalised_duration_shuffled = [];
    log_odds_normalised_duration = [];
    log_odds_normalised_duration_shuffled = [];
    for hemi = 1:2
        for i = 1:length(duration_thresholds)-1
            index = find(ripple_info.UP_duration(:,hemi) <=2 & ripple_info.normalised_UP_duration(:,hemi) > duration_thresholds(i) & ripple_info.normalised_UP_duration(:,hemi) <= duration_thresholds(i+1));
%             index = find(ripple_info.UP_duration(:,hemi) <=2 &subject_id == n & ripple_info.normalised_UP_duration(:,hemi) > duration_thresholds(i) & ripple_info.normalised_UP_duration(:,hemi) <= duration_thresholds(i+1));
            
            for iboot = 1:1000
                s = RandStream('philox4x32_10', 'Seed', iboot);
                temp_index = randsample(s, index, size(index,1), true);
                %                 V1_log_odds_normalised_duration{hemi}(i,iboot) = mean(mean_z_bias_V1_PRE(temp_index),'omitnan');
                V1_log_odds_normalised_duration{hemi}(i,iboot) = mean(mean_z_bias_V1(temp_index),'omitnan');
                log_odds_normalised_duration{hemi}(i,iboot) = mean(mean_z_bias_200ms(temp_index),'omitnan');
                %                 log_odds_normalised_duration{hemi}(i,iboot) = mean(mean_z_bias(temp_index),'omitnan');

                s = RandStream('philox4x32_10', 'Seed', iboot);
                hemi_id_shuffled = ripple_info.normalised_UP_duration;
                swap_mask = rand(s,size(ripple_info.normalised_UP_duration, 1), 1) > 0.5;
                hemi_id_shuffled(swap_mask, :) = hemi_id_shuffled(swap_mask, [2, 1]);

                temp_index = find(ripple_info.UP_duration(:,hemi) <=2 &subject_id == n & hemi_id_shuffled(:,hemi) > duration_thresholds(i) & hemi_id_shuffled(:,hemi) <= duration_thresholds(i+1));
%                 V1_log_odds_normalised_duration_shuffled{hemi}(i,iboot) = mean(mean_z_bias_V1_PRE(temp_index),'omitnan');
                V1_log_odds_normalised_duration_shuffled{hemi}(i,iboot) = mean(mean_z_bias_V1(temp_index),'omitnan');
                log_odds_normalised_duration_shuffled{hemi}(i,iboot) = mean(mean_z_bias_200ms(temp_index),'omitnan');
%                 log_odds_normalised_duration_shuffled{hemi}(i,iboot) = mean(mean_z_bias(temp_index),'omitnan');

            end
        end
    end


%     for hemi = 1:2
%         for i = 1:length(duration_thresholds)-1
%             index = find(ripple_info.UP_duration(:,hemi) <=2 &session_count == n & ripple_info.normalised_UP_duration(:,hemi) > duration_thresholds(i) & ripple_info.normalised_UP_duration(:,hemi) <= duration_thresholds(i+1));
%             V1_log_odds_normalised_duration_session(hemi,n,i) = mean(mean_z_bias_V1(index),'omitnan');
%             log_odds_normalised_duration_session(hemi,n,i) = mean(mean_z_bias_200ms(index),'omitnan');
%         end
%     end
% end
 
tvec = duration_thresholds(2:end);

% 
% fig = figure('Name','Normalised UP duration ripple log odds distribution session');
% subplot(2,2,1)
% hold on
% dist1 = squeeze(V1_log_odds_normalised_duration_session(1,:,:));
% dist2 = squeeze(V1_log_odds_normalised_duration_session(2,:,:));
% dist3 = dist2-dist1;
% 
% plot(tvec,dist1,'r')
% plot(tvec,dist2,'b')
% ylim([-0.4 0.4])
% % plot(tvec,dist3,'m')
% % plot(tvec,nanmean(dist1),'r')
% % ci = [nanmean(dist1) - std(dist1)/sqrt(max(subject_id))'; nanmean(dist1) + std(dist1)/sqrt(max(subject_id))']';
% % 
% % fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
% %     'r', 'EdgeColor','none','FaceAlpha',0.3);
% xlabel('Normalised UP duration')
% ylabel('V1 ripple Log Odds (z)')
% set(gca,'TickDir','out','Box','off','FontSize',12);
% 
% subplot(2,2,2)
% hold on
% hold on
% dist1 = squeeze(log_odds_normalised_duration_session(1,:,:));
% dist2 = squeeze(log_odds_normalised_duration_session(2,:,:));
% dist3 = dist2-dist1;
% 
% plot(tvec,dist1,'r')
% plot(tvec,dist2,'b')
% ylim([-0.3 0.3])
% % plot(tvec,dist3,'m')
% set(gca,'TickDir','out','Box','off','FontSize',12);
% xlabel('Normalised UP duration')
% ylabel('HPC ripple Log Odds (z)')
% 
% 
% subplot(2,2,3)
% hold on
% dist1 = squeeze(V1_log_odds_normalised_duration_session(1,:,:));
% dist2 = squeeze(V1_log_odds_normalised_duration_session(2,:,:));
% dist3 = dist2-dist1;
% plot(tvec,dist3,'m')
% ylim([-0.3 0.3])
% xlabel('Normalised UP duration')
% ylabel('V1 ripple Log Odds diff(z)')
% set(gca,'TickDir','out','Box','off','FontSize',12);
% 
% 
% subplot(2,2,4)
% dist1 = squeeze(log_odds_normalised_duration_session(1,:,:));
% dist2 = squeeze(log_odds_normalised_duration_session(2,:,:));
% dist3 = dist2-dist1;
% plot(tvec,dist3,'m')
% ylim([-0.4 0.4])
% xlabel('Normalised UP duration')
% ylabel('ripple Log Odds diff (z)')
% set(gca,'TickDir','out','Box','off','FontSize',12);


tvec = duration_thresholds(2:end);

%     fig = figure('Name','Normalised UP duration ripple log odds distribution PRE');
%     fig = figure('Name','Normalised UP duration ripple log odds distribution (3s)');
fig = figure('Name','Normalised UP duration ripple log odds distribution (3s)');
subplot(2,2,1)
hold on
plot(tvec,nanmean(V1_log_odds_normalised_duration{1}'),'r')
ci = prctile(V1_log_odds_normalised_duration{1}',[2.5 97.5])';

fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
    'r', 'EdgeColor','none','FaceAlpha',0.3);


plot(tvec,nanmean(V1_log_odds_normalised_duration{2}'),'b')
ci = prctile(V1_log_odds_normalised_duration{2}',[2.5 97.5])';
fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
    'b', 'EdgeColor','none','FaceAlpha',0.3);
xlabel('Normalised UP duration')
ylabel('V1 ripple Log Odds (z)')
set(gca,'TickDir','out','Box','off','FontSize',12);

subplot(2,2,2)
hold on
plot(tvec,nanmean(log_odds_normalised_duration{1}'),'r')
ci = prctile(log_odds_normalised_duration{1}',[2.5 97.5])';
fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
    'r', 'EdgeColor','none','FaceAlpha',0.3);


plot(tvec,nanmean(log_odds_normalised_duration{2}'),'b')
ci = prctile(log_odds_normalised_duration{2}',[2.5 97.5])';

fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
    'b', 'EdgeColor','none','FaceAlpha',0.3);
set(gca,'TickDir','out','Box','off','FontSize',12);
xlabel('Normalised UP duration')
ylabel('HPC ripple Log Odds (z)')


subplot(2,2,3)
hold on
plot(tvec,nanmean(V1_log_odds_normalised_duration{2}'-V1_log_odds_normalised_duration{1}'),'m')
ci = prctile(V1_log_odds_normalised_duration{2}'-V1_log_odds_normalised_duration{1}',[2.5 97.5])';

fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
    'm', 'EdgeColor','none','FaceAlpha',0.3);

hold on
plot(tvec,nanmean(V1_log_odds_normalised_duration_shuffled{2}'-V1_log_odds_normalised_duration_shuffled{1}'),'k')
ci = prctile(V1_log_odds_normalised_duration_shuffled{2}'-V1_log_odds_normalised_duration_shuffled{1}',[2.5 97.5])';

fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
    'k', 'EdgeColor','none','FaceAlpha',0.3);


xlabel('Normalised UP duration')
ylabel('V1 ripple Log Odds diff(z)')
set(gca,'TickDir','out','Box','off','FontSize',12);


subplot(2,2,4)
hold on
plot(tvec,nanmean(log_odds_normalised_duration{2}'-log_odds_normalised_duration{1}'),'m')
ci = prctile(log_odds_normalised_duration{2}'-log_odds_normalised_duration{1}',[2.5 97.5])';

fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
    'm', 'EdgeColor','none','FaceAlpha',0.3);
hold on
plot(tvec,nanmean(log_odds_normalised_duration_shuffled{2}'-log_odds_normalised_duration_shuffled{1}'),'k')
ci = prctile(log_odds_normalised_duration_shuffled{2}'-log_odds_normalised_duration_shuffled{1}',[2.5 97.5])';

fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
    'k', 'EdgeColor','none','FaceAlpha',0.3);
xlabel('Normalised UP duration')
ylabel('ripple Log Odds diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12);


    save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','UP-DOWN log odds'),[])

% D:\corticohippocampal_replay\V1-HPC bilateral interaction\UP-DOWN log odds

%% Ripple predicted by spindle and SO phase and amplitude


%%%%%%%%%%%%%%%%%% Ripple info (based in ipsi or contra)
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_spindles_probability_whole.mat'));
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

%%% ripple power
ripple_info.ripple_power = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).SWS_index==1)];

%%% spindle co-occurance
[~,spindle_index,~,index] =RestrictInts(merged_event_info.ripples_ints,merged_event_info.spindles_ints);
ripple_info.spindle_presence = spindle_index;
ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.spindle_presence_hemi(find(spindle_index)) = merged_event_info.spindles_hemisphere_id(index);

%%% spindle phase
spindle_phase=[];
spindle_phase1=[];
for probe_no = 1:2
    spindle_phase1{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        spindle_phase1{probe_no} = [spindle_phase1{probe_no} ripples_all(probe_no).spindle_phase_ripple_onset{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end


spindle_phase(1,:) = [spindle_phase1{1}(1,ripples_all(1).SWS_index==1) spindle_phase1{2}(2,ripples_all(2).SWS_index==1)];
spindle_phase(2,:) = [spindle_phase1{1}(2,ripples_all(1).SWS_index==1) spindle_phase1{2}(1,ripples_all(2).SWS_index==1)];


ripple_info.spindle_phase = spindle_phase';

%%% spindle power
spindle_amplitude=[];
spindle_amplitude1=[];
for probe_no = 1:2
    spindle_amplitude1{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        spindle_amplitude1{probe_no} = [spindle_amplitude1{probe_no} ripples_all(probe_no).spindle_amplitude_ripple_onset{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

spindle_amplitude(1,:) = [spindle_amplitude1{1}(1,ripples_all(1).SWS_index==1) spindle_amplitude1{2}(2,ripples_all(2).SWS_index==1)];
spindle_amplitude(2,:) = [spindle_amplitude1{1}(2,ripples_all(1).SWS_index==1) spindle_amplitude1{2}(1,ripples_all(2).SWS_index==1)];

ripple_info.spindle_amplitude = spindle_amplitude';

%%% SO phase
SO_phase=[];
SO_phase1=[];
for probe_no = 1:2
    SO_phase1{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        SO_phase1{probe_no} = [SO_phase1{probe_no} ripples_all(probe_no).SO_phase_ripple_onset{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

SO_phase(1,:) = [SO_phase1{1}(1,ripples_all(1).SWS_index==1) SO_phase1{2}(2,ripples_all(2).SWS_index==1)];
SO_phase(2,:) = [SO_phase1{1}(2,ripples_all(1).SWS_index==1) SO_phase1{2}(1,ripples_all(2).SWS_index==1)];
ripple_info.SO_phase = SO_phase';

%%% SO power
SO_amplitude=[];
SO_amplitude1=[];
for probe_no = 1:2
    SO_amplitude1{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        SO_amplitude1{probe_no} = [SO_amplitude1{probe_no} ripples_all(probe_no).SO_amplitude_ripple_onset{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

SO_amplitude(1,:) = [SO_amplitude1{1}(1,ripples_all(1).SWS_index==1) SO_amplitude1{2}(2,ripples_all(2).SWS_index==1)];
SO_amplitude(2,:) = [SO_amplitude1{1}(2,ripples_all(1).SWS_index==1) SO_amplitude1{2}(1,ripples_all(2).SWS_index==1)];

ripple_info.SO_amplitude = SO_amplitude';

%%%%%% Ripple session count and subject_id
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);


% singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));
singlet_index = logical(ones(length(merged_event_info.ripples_peaktimes),1));


%%%%%%%%%%%%% Mixed effect model (Does spindle and SO predict ripple power?)

output = predict_ripples_by_SO_spindles(ripple_info,'subject_id',subject_id);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_amp_phase_Ripples_output.mat'),'output');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_amp_phase_Ripples_output.mat'),'output')
predict_ripples_by_SO_spindles(ripple_info,'subject_id',subject_id,'output',output,'plot_option',1);

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')



%%%%%%%%%%%%% UP transition
output = predict_ripples_by_SO_spindles(ripple_info,'subject_id',subject_id);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_amp_phase_Ripples_UP_transition_output.mat'),'output');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_amp_phase_Ripples_DOWN_transition_output.mat'),'output')
predict_ripples_by_SO_spindles(ripple_info,'subject_id',subject_id,'output',output,'plot_option',1);

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')


predict_ripples_by_SO_spindles_TF(ripple)




% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis'),[])
