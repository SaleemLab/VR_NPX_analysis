%% Loading extracted SO ripples spindles info for plotting and analysis 
clear all
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
% load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
load(fullfile(analysis_folder,'ripples_all_POST.mat'))
load(fullfile(analysis_folder,'spindles_all_POST.mat'))
load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'),'UP_DOWN_ripple_PSTH_MUA');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));

load(fullfile(analysis_folder,'KDE_reactivation_DOWN_all_POST.mat'),'KDE_reactivation_DOWN_all')
load(fullfile(analysis_folder,'KDE_reactivation_UP_all_POST.mat'),'KDE_reactivation_UP_all')
load(fullfile(analysis_folder,'KDE_reactivation_all_POST.mat'),'KDE_reactivation_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_DOWN_all_POST.mat'),'KDE_reactivation_V1_DOWN_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_UP_all_POST.mat'),'KDE_reactivation_V1_UP_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_all_POST.mat'),'KDE_reactivation_V1_all')


%%%% Combined V1 and HPC PSTH extraction for both bias and log odds

all_sessions = max(slow_waves_all(1).DOWN_session_count);
sessions_to_process = 1:all_sessions;

% Parameters
psth_window = [-1 1];
psth_step = 0.01;
nbins = round(diff(psth_window) / psth_step);
half_bins = nbins / 2;

% Initialize PSTH structure
KDE_reactivation_PSTH = struct();

for region = ["V1", "HPC"]
    up_bias_all = [];
    down_bias_all = [];
    up_bias_all_z = [];
    down_bias_all_z = [];
    up_bias_all_zshuff = [];
    down_bias_all_zshuff = [];

    up_logodds_all = [];
    down_logodds_all = [];
    up_logodds_all_z = [];
    down_logodds_all_z = [];
    up_logodds_all_zshuff = [];
    down_logodds_all_zshuff = [];
    for nprobe = 1:2
        if region == "V1"
            up_struct = KDE_reactivation_V1_UP_all;
            down_struct = KDE_reactivation_V1_DOWN_all;
        else
            up_struct = KDE_reactivation_UP_all;
            down_struct = KDE_reactivation_DOWN_all;
        end

        for nsession = 1:length(up_struct(nprobe).bias)
            bias_up = up_struct(nprobe).bias{nsession};
            bias_down = down_struct(nprobe).bias{nsession};
            z_up = up_struct(nprobe).zscored_bias{nsession};
            z_down = down_struct(nprobe).zscored_bias{nsession};
            zshuff_up = up_struct(nprobe).zscored_bias_shuffled{nsession};
            zshuff_down = down_struct(nprobe).zscored_bias_shuffled{nsession};

            logodds_up = up_struct(nprobe).log_odds{nsession};
            logodds_down = down_struct(nprobe).log_odds{nsession};
            zlog_up = up_struct(nprobe).zscored_log_odds{nsession};
            zlog_down = down_struct(nprobe).zscored_log_odds{nsession};
            zshuff_log_up = up_struct(nprobe).zscored_log_odds_shuffled{nsession};
            zshuff_log_down = down_struct(nprobe).zscored_log_odds_shuffled{nsession};

            time_up = up_struct(nprobe).event_bins{nsession}(:,1)';
            time_down = down_struct(nprobe).event_bins{nsession}(:,1)';
            event_id_up = up_struct(nprobe).event_id{nsession};
            event_id_down = down_struct(nprobe).event_id{nsession};

            up_all = find(slow_waves_all(nprobe).UP_session_count == nsession);
            UP_event_index = intersect(up_all, probability(nprobe).UP_all_index);

            temp_index = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
            [~, ia, ib] = intersect(slow_waves_all(nprobe).DOWN_ints(temp_index, 2), ...
                                    slow_waves_all(nprobe).UP_ints(UP_event_index, 1));

            up_ints = slow_waves_all(nprobe).UP_ints(up_all, :);
            previous_DOWN_event_index = ia;
            UP_event_index = find(ismember(up_all, UP_event_index));

            down_all = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
            DOWN_event_index = intersect(down_all, probability(nprobe).DOWN_all_index);
            DOWN_event_index = find(ismember(down_all, DOWN_event_index));
            down_ints = slow_waves_all(nprobe).DOWN_ints(down_all, :);

            for i = 1:length(UP_event_index)
                up_bias = nan(1, nbins);
                up_z = nan(1, nbins);
                up_zshuff = nan(1, nbins);

                up_log = nan(1, nbins);
                up_log_z = nan(1, nbins);
                up_log_zshuff = nan(1, nbins);

                if i <= length(previous_DOWN_event_index)
                    eid = previous_DOWN_event_index(i);
                    event_end = down_ints(eid, 2);
                    rel_time = time_down - event_end;
                    idx = find(event_id_down == eid & rel_time >= -1 & rel_time < 0);
                    rel_times = rel_time(idx);
                    bin_idx = round((rel_times + 1) / psth_step) + 1;
                    valid = bin_idx > 0 & bin_idx <= half_bins;

                    up_bias(bin_idx(valid)) = bias_down(idx(valid));
                    up_z(bin_idx(valid)) = z_down(idx(valid));
                    up_zshuff(bin_idx(valid)) = zshuff_down(idx(valid));

                    up_log(bin_idx(valid)) = logodds_down(idx(valid));
                    up_log_z(bin_idx(valid)) = zlog_down(idx(valid));
                    up_log_zshuff(bin_idx(valid)) = zshuff_log_down(idx(valid));
                end

                eid = UP_event_index(i);
                event_start = up_ints(eid, 1);
                rel_time = time_up - event_start;
                idx = find(event_id_up == eid & rel_time >= 0 & rel_time < 1);
                rel_times = rel_time(idx);
                bin_idx = round(rel_times / psth_step) + half_bins;
                valid = bin_idx > half_bins & bin_idx <= nbins;

                up_bias(bin_idx(valid)) = bias_up(idx(valid));
                up_z(bin_idx(valid)) = z_up(idx(valid));
                up_zshuff(bin_idx(valid)) = zshuff_up(idx(valid));

                up_log(bin_idx(valid)) = logodds_up(idx(valid));
                up_log_z(bin_idx(valid)) = zlog_up(idx(valid));
                up_log_zshuff(bin_idx(valid)) = zshuff_log_up(idx(valid));

                up_bias_all = [up_bias_all; up_bias];
                up_bias_all_z = [up_bias_all_z; up_z];
                up_bias_all_zshuff = [up_bias_all_zshuff; up_zshuff];

                up_logodds_all = [up_logodds_all; up_log];
                up_logodds_all_z = [up_logodds_all_z; up_log_z];
                up_logodds_all_zshuff = [up_logodds_all_zshuff; up_log_zshuff];
            end

            for i = 1:length(DOWN_event_index)
                down_bias = nan(1, nbins);
                down_z = nan(1, nbins);
                down_zshuff = nan(1, nbins);

                down_log = nan(1, nbins);
                down_log_z = nan(1, nbins);
                down_log_zshuff = nan(1, nbins);

                if i <= length(UP_event_index)
                    eid = UP_event_index(i);
                    event_end = up_ints(eid, 2);
                    rel_time = time_up - event_end;
                    idx = find(event_id_up == eid & rel_time >= -1 & rel_time < 0);
                    rel_times = rel_time(idx);
                    bin_idx = round((rel_times + 1) / psth_step) + 1;
                    valid = bin_idx > 0 & bin_idx <= half_bins;

                    down_bias(bin_idx(valid)) = bias_up(idx(valid));
                    down_z(bin_idx(valid)) = z_up(idx(valid));
                    down_zshuff(bin_idx(valid)) = zshuff_up(idx(valid));

                    down_log(bin_idx(valid)) = logodds_up(idx(valid));
                    down_log_z(bin_idx(valid)) = zlog_up(idx(valid));
                    down_log_zshuff(bin_idx(valid)) = zshuff_log_up(idx(valid));
                end

                eid = DOWN_event_index(i);
                event_start = down_ints(eid, 1);
                rel_time = time_down - event_start;
                idx = find(event_id_down == eid & rel_time >= 0 & rel_time < 1);
                rel_times = rel_time(idx);
                bin_idx = round(rel_times / psth_step) + half_bins;
                valid = bin_idx > half_bins & bin_idx <= nbins;

                down_bias(bin_idx(valid)) = bias_down(idx(valid));
                down_z(bin_idx(valid)) = z_down(idx(valid));
                down_zshuff(bin_idx(valid)) = zshuff_down(idx(valid));

                down_log(bin_idx(valid)) = logodds_down(idx(valid));
                down_log_z(bin_idx(valid)) = zlog_down(idx(valid));
                down_log_zshuff(bin_idx(valid)) = zshuff_log_down(idx(valid));

                down_bias_all = [down_bias_all; down_bias];
                down_bias_all_z = [down_bias_all_z; down_z];
                down_bias_all_zshuff = [down_bias_all_zshuff; down_zshuff];

                down_logodds_all = [down_logodds_all; down_log];
                down_logodds_all_z = [down_logodds_all_z; down_log_z];
                down_logodds_all_zshuff = [down_logodds_all_zshuff; down_log_zshuff];
            end
        end

        prefix = sprintf('%s_', region);
        KDE_reactivation_PSTH.([prefix 'UP']) = up_bias_all;
        KDE_reactivation_PSTH.([prefix 'DOWN']) = down_bias_all;
        KDE_reactivation_PSTH.([prefix 'UP_z']) = up_bias_all_z;
        KDE_reactivation_PSTH.([prefix 'DOWN_z']) = down_bias_all_z;
        KDE_reactivation_PSTH.([prefix 'UP_zshuff']) = up_bias_all_zshuff;
        KDE_reactivation_PSTH.([prefix 'DOWN_zshuff']) = down_bias_all_zshuff;

        KDE_reactivation_PSTH.([prefix 'UP_log_odds']) = up_logodds_all;
        KDE_reactivation_PSTH.([prefix 'DOWN_log_odds']) = down_logodds_all;
        KDE_reactivation_PSTH.([prefix 'UP_log_odds_z']) = up_logodds_all_z;
        KDE_reactivation_PSTH.([prefix 'DOWN_log_odds_z']) = down_logodds_all_z;
        KDE_reactivation_PSTH.([prefix 'UP_log_odds_zshuff']) = up_logodds_all_zshuff;
        KDE_reactivation_PSTH.([prefix 'DOWN_log_odds_zshuff']) = down_logodds_all_zshuff;
    end
end

save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_PSTH.mat'),'KDE_reactivation_PSTH');


%%%%% HPC ripples
% Parameters
ripple_window = [0, 0.2];  % seconds
psth_step = 0.01;          % 10 ms
nbins = round(diff(ripple_window) / psth_step);

% Metrics to extract
% metric_names = {'bias', 'zscored_bias', 'zscored_bias_shuffled'};
metric_names = {'bias', 'zscored_bias', 'log_odds', 'zscored_log_odds','zscored_log_odds_shuffled','log_odds_percentile'};
% Output containers
ripple_bias_masked_SWS = struct();
ripple_bias_masked_nonSWS = struct();
for region = ["HPC", "V1"]
    for metric = metric_names
        ripple_bias_masked_SWS.(region).(metric{1})     = cell(1, 2);
        ripple_bias_masked_nonSWS.(region).(metric{1})  = cell(1, 2);
    end
end

% Main loop
for isSWS = [1, 0]
    for nprobe = 1:2
        for region = ["HPC", "V1"]

            if contains(region,'V1')
                temp = KDE_reactivation_V1_all;
            else
                temp = KDE_reactivation_all;
            end

            for metric = metric_names
                metric_str = metric{1};
                masked_all = [];

                for nsession = 1:length(temp(nprobe).(metric_str))
                    % Get all ripple indices for this session
                    session_mask = ripples_all(nprobe).session_count == nsession;
                    ripple_idx_all = find(session_mask);

                    % Get SWS or non-SWS ripple subset
                    state_mask = ripples_all(nprobe).SWS_index == isSWS;
                    [ripple_idx_all,ripple_idx_session,~] = intersect(ripple_idx_all,find(state_mask));

                    if isempty(ripple_idx_session), continue; end

                    ripple_times = [ ...
                        ripples_all(nprobe).onset(ripple_idx_all), ...
                        ripples_all(nprobe).offset(ripple_idx_all)];

                    % Get data
                    bias_vals = temp(nprobe).(metric_str){nsession};
                    time_bins = temp(nprobe).event_bins{nsession}(:,1)';
                    event_ids = temp(nprobe).event_id{nsession};

                    for i = 1:length(ripple_idx_session)
                        ripple_id = ripple_idx_session(i);
                        t_start = ripple_times(i, 1);

                        rel_time = time_bins - t_start;
                        idx = find(event_ids == ripple_id & rel_time >= 0 & rel_time <= 0.2);
                   
                        ripple_bias = nan(1, nbins);

                        ripple_bias(1:length(idx)) = bias_vals(idx);
                        % Mask if another ripple starts within this window
                        all_starts = ripples_all(nprobe).onset(ripple_idx_all);
                        next_ripple = all_starts > t_start & all_starts < (t_start + 0.2);
                        if any(next_ripple)
                            t_next = min(all_starts(next_ripple));
                            mask_start = round((t_next - t_start) / psth_step) + 1;
                            ripple_bias(mask_start:end) = nan;
                        end

                        masked_all = [masked_all; ripple_bias];
                    end
                end

                % Save result
                if isSWS
                    ripple_bias_masked_SWS.(region).(metric_str){nprobe} = masked_all;
                else
                    ripple_bias_masked_nonSWS.(region).(metric_str){nprobe} = masked_all;
                end
            end
        end
    end
end


% SWS ripples
% SWS ripples
% bias
KDE_reactivation_content.HPC_ripples                  = [ripple_bias_masked_SWS.HPC.bias{1}; ripple_bias_masked_SWS.HPC.bias{2}];
KDE_reactivation_content.V1_ripples                   = [ripple_bias_masked_SWS.V1.bias{1}; ripple_bias_masked_SWS.V1.bias{2}];
KDE_reactivation_content.HPC_awake_ripples            = [ripple_bias_masked_nonSWS.HPC.bias{1}; ripple_bias_masked_nonSWS.HPC.bias{2}];
KDE_reactivation_content.V1_awake_ripples             = [ripple_bias_masked_nonSWS.V1.bias{1}; ripple_bias_masked_nonSWS.V1.bias{2}];

% zscored_bias
KDE_reactivation_content.HPC_z_ripples                = [ripple_bias_masked_SWS.HPC.zscored_bias{1}; ripple_bias_masked_SWS.HPC.zscored_bias{2}];
KDE_reactivation_content.V1_z_ripples                 = [ripple_bias_masked_SWS.V1.zscored_bias{1}; ripple_bias_masked_SWS.V1.zscored_bias{2}];
KDE_reactivation_content.HPC_z_awake_ripples          = [ripple_bias_masked_nonSWS.HPC.zscored_bias{1}; ripple_bias_masked_nonSWS.HPC.zscored_bias{2}];
KDE_reactivation_content.V1_z_awake_ripples           = [ripple_bias_masked_nonSWS.V1.zscored_bias{1}; ripple_bias_masked_nonSWS.V1.zscored_bias{2}];

% log_odds
KDE_reactivation_content.HPC_logodds_ripples         = [ripple_bias_masked_SWS.HPC.log_odds{1}; ripple_bias_masked_SWS.HPC.log_odds{2}];
KDE_reactivation_content.V1_logodds_ripples          = [ripple_bias_masked_SWS.V1.log_odds{1}; ripple_bias_masked_SWS.V1.log_odds{2}];
KDE_reactivation_content.HPC_logodds_awake_ripples   = [ripple_bias_masked_nonSWS.HPC.log_odds{1}; ripple_bias_masked_nonSWS.HPC.log_odds{2}];
KDE_reactivation_content.V1_logodds_awake_ripples    = [ripple_bias_masked_nonSWS.V1.log_odds{1}; ripple_bias_masked_nonSWS.V1.log_odds{2}];

% zscored_log_odds
KDE_reactivation_content.HPC_z_logodds_ripples       = [ripple_bias_masked_SWS.HPC.zscored_log_odds{1}; ripple_bias_masked_SWS.HPC.zscored_log_odds{2}];
KDE_reactivation_content.V1_z_logodds_ripples        = [ripple_bias_masked_SWS.V1.zscored_log_odds{1}; ripple_bias_masked_SWS.V1.zscored_log_odds{2}];
KDE_reactivation_content.HPC_z_logodds_awake_ripples = [ripple_bias_masked_nonSWS.HPC.zscored_log_odds{1}; ripple_bias_masked_nonSWS.HPC.zscored_log_odds{2}];
KDE_reactivation_content.V1_z_logodds_awake_ripples  = [ripple_bias_masked_nonSWS.V1.zscored_log_odds{1}; ripple_bias_masked_nonSWS.V1.zscored_log_odds{2}];

% zscored_log_odds_shuffled
KDE_reactivation_content.HPC_zshuff_logodds_ripples       = [ripple_bias_masked_SWS.HPC.zscored_log_odds_shuffled{1}; ripple_bias_masked_SWS.HPC.zscored_log_odds_shuffled{2}];
KDE_reactivation_content.V1_zshuff_logodds_ripples        = [ripple_bias_masked_SWS.V1.zscored_log_odds_shuffled{1}; ripple_bias_masked_SWS.V1.zscored_log_odds_shuffled{2}];
KDE_reactivation_content.HPC_zshuff_logodds_awake_ripples = [ripple_bias_masked_nonSWS.HPC.zscored_log_odds_shuffled{1}; ripple_bias_masked_nonSWS.HPC.zscored_log_odds_shuffled{2}];
KDE_reactivation_content.V1_zshuff_logodds_awake_ripples  = [ripple_bias_masked_nonSWS.V1.zscored_log_odds_shuffled{1}; ripple_bias_masked_nonSWS.V1.zscored_log_odds_shuffled{2}];

% log_odds_percentile
KDE_reactivation_content.HPC_logodds_percentile_ripples         = [ripple_bias_masked_SWS.HPC.log_odds_percentile{1}; ripple_bias_masked_SWS.HPC.log_odds_percentile{2}];
KDE_reactivation_content.V1_logodds_percentile_ripples          = [ripple_bias_masked_SWS.V1.log_odds_percentile{1}; ripple_bias_masked_SWS.V1.log_odds_percentile{2}];
KDE_reactivation_content.HPC_logodds_percentile_awake_ripples   = [ripple_bias_masked_nonSWS.HPC.log_odds_percentile{1}; ripple_bias_masked_nonSWS.HPC.log_odds_percentile{2}];
KDE_reactivation_content.V1_logodds_percentile_awake_ripples    = [ripple_bias_masked_nonSWS.V1.log_odds_percentile{1}; ripple_bias_masked_nonSWS.V1.log_odds_percentile{2}];

save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_content.mat'),'KDE_reactivation_content');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data
load(fullfile(analysis_folder,'KDE_reactivation_ripples_all_POST.mat'),'KDE_reactivation_ripples_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_ripples_all_POST.mat'),'KDE_reactivation_V1_ripples_all')

ripple_window = [-1, 1];
psth_step = 0.01;
nbins = round(diff(ripple_window) / psth_step);
metric_names = {'bias', 'zscored_bias', 'log_odds', 'zscored_log_odds'};
regions = ["HPC", "V1"];

% Initialize
ripple_bias_masked_SWS = struct();
ripple_bias_masked_nonSWS = struct();
KDE_reactivation_ripples_PSTH = struct();
nan_mask = cell(1,2);          % for probe 1,2 SWS
nan_mask_awake = cell(1,2);    % for probe 1,2 awake

% Loop over states
for isSWS = [1, 0]
    for nprobe = 1:2
        for region = regions
            for metric = metric_names
                metric_str = metric{1};
                if isSWS
                    ripple_bias_masked_SWS.(region).(metric_str){nprobe} = [];
                else
                    ripple_bias_masked_nonSWS.(region).(metric_str){nprobe} = [];
                end
            end
        end

        nan_mask_probe = [];  % accumulate here

        for region = regions
            if region == "V1"
                temp = KDE_reactivation_V1_ripples_all;
            else
                temp = KDE_reactivation_ripples_all;
            end

            for metric = metric_names
                metric_str = metric{1};
                metric_all = [];

                for nsession = 1:length(temp(nprobe).(metric_str))
                    session_mask = ripples_all(nprobe).session_count == nsession;
                    ripple_idx_all = find(session_mask);
                    state_mask = ripples_all(nprobe).SWS_index == isSWS;
                    [ripple_idx_all, ripple_idx_session] = intersect(ripple_idx_all, find(state_mask));
                    if isempty(ripple_idx_session), continue; end

                    ripple_times = [ripples_all(nprobe).onset(ripple_idx_all), ...
                                    ripples_all(nprobe).offset(ripple_idx_all)];
                    metric_vals = temp(nprobe).(metric_str){nsession};
                    time_bins = temp(nprobe).event_bins{nsession}(:,1)';
                    event_ids = temp(nprobe).event_id{nsession};

                    for i = 1:length(ripple_idx_session)
                        ripple_id = ripple_idx_session(i);
                        t_start = ripple_times(i, 1);
                        rel_time = time_bins - t_start;
                        idx = find(event_ids == ripple_id & rel_time >= -1 & rel_time < 1);

                        ripple_vec = nan(1, nbins);
                        ripple_vec(1:length(idx)) = metric_vals(idx);
                        metric_all = [metric_all; ripple_vec];

                        % Only compute mask once (on HPC region, one metric)
                        if strcmp(region, 'HPC') && strcmp(metric_str, 'bias')
                            nan_row = zeros(1, nbins);  % 0 = keep, NaN = mask

                            all_onsets = ripple_times(:,1);
                            all_offsets = ripple_times(:,2);

                            % Mask next ripple
                            next_ripples = all_onsets > t_start & all_onsets < (t_start + 1);
                            if any(next_ripples)
                                t_next = min(all_onsets(next_ripples));
                                mask_start = round((t_next - t_start) / psth_step) + 101;
                                if mask_start <= nbins
                                    nan_row(mask_start:end) = NaN;
                                end
                            end

                            % Mask previous ripple
                            prev_ripples = all_offsets < t_start & all_offsets > (t_start - 1);
                            if any(prev_ripples)
                                t_prev = max(all_offsets(prev_ripples));
                                mask_end = round((t_prev - t_start) / psth_step) + 101;
                                if mask_end >= 1
                                    nan_row(1:mask_end) = NaN;
                                end
                            end

                            nan_mask_probe = [nan_mask_probe; nan_row];
                        end
                    end
                end

                % Save metric data per region/state
                if isSWS
                    ripple_bias_masked_SWS.(region).(metric_str){nprobe} = metric_all;
                else
                    ripple_bias_masked_nonSWS.(region).(metric_str){nprobe} = metric_all;
                end
            end
        end

        % Store mask for this probe
        if isSWS
            nan_mask{nprobe} = nan_mask_probe;
        else
            nan_mask_awake{nprobe} = nan_mask_probe;
        end
    end
end

% Combine and store
for region = regions
    region = char(region)
    KDE_reactivation_ripples_PSTH.([region '_ripples'])             = [ripple_bias_masked_SWS.(region).bias{1};              ripple_bias_masked_SWS.(region).bias{2}];
    KDE_reactivation_ripples_PSTH.([region '_z_ripples'])           = [ripple_bias_masked_SWS.(region).zscored_bias{1};      ripple_bias_masked_SWS.(region).zscored_bias{2}];
    KDE_reactivation_ripples_PSTH.([region '_logodds_ripples'])     = [ripple_bias_masked_SWS.(region).log_odds{1};          ripple_bias_masked_SWS.(region).log_odds{2}];
    KDE_reactivation_ripples_PSTH.([region '_z_logodds_ripples'])   = [ripple_bias_masked_SWS.(region).zscored_log_odds{1};  ripple_bias_masked_SWS.(region).zscored_log_odds{2}];

    KDE_reactivation_ripples_PSTH.([region '_awake_ripples'])             = [ripple_bias_masked_nonSWS.(region).bias{1};              ripple_bias_masked_nonSWS.(region).bias{2}];
    KDE_reactivation_ripples_PSTH.([region '_z_awake_ripples'])           = [ripple_bias_masked_nonSWS.(region).zscored_bias{1};      ripple_bias_masked_nonSWS.(region).zscored_bias{2}];
    KDE_reactivation_ripples_PSTH.([region '_logodds_awake_ripples'])     = [ripple_bias_masked_nonSWS.(region).log_odds{1};          ripple_bias_masked_nonSWS.(region).log_odds{2}];
    KDE_reactivation_ripples_PSTH.([region '_z_logodds_awake_ripples'])   = [ripple_bias_masked_nonSWS.(region).zscored_log_odds{1};  ripple_bias_masked_nonSWS.(region).zscored_log_odds{2}];
end

% Add nan masks
KDE_reactivation_ripples_PSTH.nan_mask        = [nan_mask{1}; nan_mask{2}];
KDE_reactivation_ripples_PSTH.nan_mask_awake  = [nan_mask_awake{1}; nan_mask_awake{2}];

% Save
save(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'KDE_reactivation_ripples_PSTH.mat'), 'KDE_reactivation_ripples_PSTH');


