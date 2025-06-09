%% Loading extracted SO ripples spindles info for plotting and analysis 
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




%%%%% V1: Grabbing bias at DOWN-UP and UP-DOWN transition

all_sessions = max(slow_waves_all(1).DOWN_session_count);
sessions_to_process = 1:all_sessions;

% Initialize output variables
UP_bias_V1_all = cell(1, 2);
DOWN_bias_V1_all = cell(1, 2);
UP_bias_V1_all_z = cell(1, 2);
DOWN_bias_V1_all_z = cell(1, 2);
UP_bias_V1_all_zshuff = cell(1, 2);
DOWN_bias_V1_all_zshuff = cell(1, 2);

% Parameters
psth_window = [-1 1];
psth_step = 0.01;
nbins = round(diff(psth_window) / psth_step);
half_bins = nbins / 2;

for nprobe = 1:2
    tic
    up_bias_all = [];
    down_bias_all = [];
    up_bias_all_z = [];
    down_bias_all_z = [];
    up_bias_all_zshuff = [];
    down_bias_all_zshuff = [];

    for nsession = 1:length(KDE_reactivation_V1_UP_all(nprobe).bias)
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

        % Bias data
        bias_up = KDE_reactivation_V1_UP_all(nprobe).bias{nsession};
        bias_down = KDE_reactivation_V1_DOWN_all(nprobe).bias{nsession};
        z_up = KDE_reactivation_V1_UP_all(nprobe).zscored_bias{nsession};
        z_down = KDE_reactivation_V1_DOWN_all(nprobe).zscored_bias{nsession};
        zshuff_up = KDE_reactivation_V1_UP_all(nprobe).zscored_bias_shuffled{nsession};
        zshuff_down = KDE_reactivation_V1_DOWN_all(nprobe).zscored_bias_shuffled{nsession};
        time_up = KDE_reactivation_V1_UP_all(nprobe).event_bins{nsession}(:,1)';
        time_down = KDE_reactivation_V1_DOWN_all(nprobe).event_bins{nsession}(:,1)';
        event_id_up = KDE_reactivation_V1_UP_all(nprobe).event_id{nsession};
        event_id_down = KDE_reactivation_V1_DOWN_all(nprobe).event_id{nsession};

        for i = 1:length(UP_event_index)
            up_bias = nan(1, nbins);
            up_z = nan(1, nbins);
            up_zshuff = nan(1, nbins);

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

            up_bias_all = [up_bias_all; up_bias];
            up_bias_all_z = [up_bias_all_z; up_z];
            up_bias_all_zshuff = [up_bias_all_zshuff; up_zshuff];
        end

        for i = 1:length(DOWN_event_index)
            down_bias = nan(1, nbins);
            down_z = nan(1, nbins);
            down_zshuff = nan(1, nbins);

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

            down_bias_all = [down_bias_all; down_bias];
            down_bias_all_z = [down_bias_all_z; down_z];
            down_bias_all_zshuff = [down_bias_all_zshuff; down_zshuff];
        end
    end

    UP_bias_V1_all{nprobe} = up_bias_all;
    DOWN_bias_V1_all{nprobe} = down_bias_all;
    UP_bias_V1_all_z{nprobe} = up_bias_all_z;
    DOWN_bias_V1_all_z{nprobe} = down_bias_all_z;
    UP_bias_V1_all_zshuff{nprobe} = up_bias_all_zshuff;
    DOWN_bias_V1_all_zshuff{nprobe} = down_bias_all_zshuff;
    toc
end

KDE_reactivation_PSTH.V1_UP = [UP_bias_V1_all{1}; UP_bias_V1_all{2}];
KDE_reactivation_PSTH.V1_DOWN = [DOWN_bias_V1_all{1}; DOWN_bias_V1_all{2}];
KDE_reactivation_PSTH.V1_UP_z = [UP_bias_V1_all_z{1}; UP_bias_V1_all_z{2}];
KDE_reactivation_PSTH.V1_DOWN_z = [DOWN_bias_V1_all_z{1}; DOWN_bias_V1_all_z{2}];
KDE_reactivation_PSTH.V1_UP_zshuff = [UP_bias_V1_all_zshuff{1}; UP_bias_V1_all_zshuff{2}];
KDE_reactivation_PSTH.V1_DOWN_zshuff = [DOWN_bias_V1_all_zshuff{1}; DOWN_bias_V1_all_zshuff{2}];


%%%%% HPC: Same logic as above (just change variable names)
% Initialize output
UP_bias_all = cell(1, 2);
DOWN_bias_all = cell(1, 2);
UP_bias_all_z = cell(1, 2);
DOWN_bias_all_z = cell(1, 2);
UP_bias_all_zshuff = cell(1, 2);
DOWN_bias_all_zshuff = cell(1, 2);

for nprobe = 1:2
    tic
    up_bias_all = [];
    down_bias_all = [];
    up_bias_all_z = [];
    down_bias_all_z = [];
    up_bias_all_zshuff = [];
    down_bias_all_zshuff = [];

    for nsession = 1:length(KDE_reactivation_UP_all(nprobe).bias)
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

        bias_up = KDE_reactivation_UP_all(nprobe).bias{nsession};
        bias_down = KDE_reactivation_DOWN_all(nprobe).bias{nsession};
        z_up = KDE_reactivation_UP_all(nprobe).zscored_bias{nsession};
        z_down = KDE_reactivation_DOWN_all(nprobe).zscored_bias{nsession};
        zshuff_up = KDE_reactivation_UP_all(nprobe).zscored_bias_shuffled{nsession};
        zshuff_down = KDE_reactivation_DOWN_all(nprobe).zscored_bias_shuffled{nsession};
        time_up = KDE_reactivation_UP_all(nprobe).event_bins{nsession}(:,1)';
        time_down = KDE_reactivation_DOWN_all(nprobe).event_bins{nsession}(:,1)';
        event_id_up = KDE_reactivation_UP_all(nprobe).event_id{nsession};
        event_id_down = KDE_reactivation_DOWN_all(nprobe).event_id{nsession};

        for i = 1:length(UP_event_index)
            up_bias = nan(1, nbins);
            up_z = nan(1, nbins);
            up_zshuff = nan(1, nbins);

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

            up_bias_all = [up_bias_all; up_bias];
            up_bias_all_z = [up_bias_all_z; up_z];
            up_bias_all_zshuff = [up_bias_all_zshuff; up_zshuff];
        end

        for i = 1:length(DOWN_event_index)
            down_bias = nan(1, nbins);
            down_z = nan(1, nbins);
            down_zshuff = nan(1, nbins);

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

            down_bias_all = [down_bias_all; down_bias];
            down_bias_all_z = [down_bias_all_z; down_z];
            down_bias_all_zshuff = [down_bias_all_zshuff; down_zshuff];
        end
    end

    UP_bias_all{nprobe} = up_bias_all;
    DOWN_bias_all{nprobe} = down_bias_all;
    UP_bias_all_z{nprobe} = up_bias_all_z;
    DOWN_bias_all_z{nprobe} = down_bias_all_z;
    UP_bias_all_zshuff{nprobe} = up_bias_all_zshuff;
    DOWN_bias_all_zshuff{nprobe} = down_bias_all_zshuff;
    toc
end

KDE_reactivation_PSTH.HPC_UP = [UP_bias_all{1}; UP_bias_all{2}];
KDE_reactivation_PSTH.HPC_DOWN = [DOWN_bias_all{1}; DOWN_bias_all{2}];
KDE_reactivation_PSTH.HPC_UP_z = [UP_bias_all_z{1}; UP_bias_all_z{2}];
KDE_reactivation_PSTH.HPC_DOWN_z = [DOWN_bias_all_z{1}; DOWN_bias_all_z{2}];
KDE_reactivation_PSTH.HPC_UP_zshuff = [UP_bias_all_zshuff{1}; UP_bias_all_zshuff{2}];
KDE_reactivation_PSTH.HPC_DOWN_zshuff = [DOWN_bias_all_zshuff{1}; DOWN_bias_all_zshuff{2}];


save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_PSTH.mat'),'KDE_reactivation_PSTH');


%%%%% HPC ripples
% Parameters
ripple_window = [0, 0.2];  % seconds
psth_step = 0.01;          % 10 ms
nbins = round(diff(ripple_window) / psth_step);

% Metrics to extract
metric_names = {'bias', 'zscored_bias', 'zscored_bias_shuffled'};

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
KDE_reactivation_content.HPC_ripples             = [ripple_bias_masked_SWS.HPC.bias{1};             ripple_bias_masked_SWS.HPC.bias{2}];
KDE_reactivation_content.V1_ripples              = [ripple_bias_masked_SWS.V1.bias{1};              ripple_bias_masked_SWS.V1.bias{2}];
KDE_reactivation_content.HPC_z_ripples     = [ripple_bias_masked_SWS.HPC.zscored_bias{1};     ripple_bias_masked_SWS.HPC.zscored_bias{2}];
KDE_reactivation_content.V1_z_ripples      = [ripple_bias_masked_SWS.V1.zscored_bias{1};      ripple_bias_masked_SWS.V1.zscored_bias{2}];
KDE_reactivation_content.HPC_zshuff_ripples    = [ripple_bias_masked_SWS.HPC.zscored_bias_shuffled{1}; ripple_bias_masked_SWS.HPC.zscored_bias_shuffled{2}];
KDE_reactivation_content.V1_zshuff_ripples     = [ripple_bias_masked_SWS.V1.zscored_bias_shuffled{1}; ripple_bias_masked_SWS.V1.zscored_bias_shuffled{2}];

% Non-SWS (awake/REM)
KDE_reactivation_content.HPC_awake_ripples             = [ripple_bias_masked_nonSWS.HPC.bias{1};             ripple_bias_masked_nonSWS.HPC.bias{2}];
KDE_reactivation_content.V1_awake_ripples              = [ripple_bias_masked_nonSWS.V1.bias{1};              ripple_bias_masked_nonSWS.V1.bias{2}];
KDE_reactivation_content.HPC_z_awake_ripples     = [ripple_bias_masked_nonSWS.HPC.zscored_bias{1};     ripple_bias_masked_nonSWS.HPC.zscored_bias{2}];
KDE_reactivation_content.V1_z_awake_ripples      = [ripple_bias_masked_nonSWS.V1.zscored_bias{1};      ripple_bias_masked_nonSWS.V1.zscored_bias{2}];
KDE_reactivation_content.HPC_zshuff_awake_ripples    = [ripple_bias_masked_nonSWS.HPC.zscored_bias_shuffled{1}; ripple_bias_masked_nonSWS.HPC.zscored_bias_shuffled{2}];
KDE_reactivation_content.V1_zshuff_awake_ripples     = [ripple_bias_masked_nonSWS.V1.zscored_bias_shuffled{1}; ripple_bias_masked_nonSWS.V1.zscored_bias_shuffled{2}];

save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_content.mat'),'KDE_reactivation_content');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(analysis_folder,'KDE_reactivation_ripples_all_POST.mat'),'KDE_reactivation_ripples_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_ripples_all_POST.mat'),'KDE_reactivation_V1_ripples_all')


%%%%% HPC ripples
% Parameters
ripple_window = [-1, 1];  % seconds
psth_step = 0.01;          % 10 ms
nbins = round(diff(ripple_window) / psth_step);

% Metrics to extract
metric_names = {'bias', 'zscored_bias'};

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
                temp = KDE_reactivation_V1_ripples_all;
            else
                temp = KDE_reactivation_ripples_all;
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
                        idx = find(event_ids == ripple_id & rel_time >= -1 & rel_time < 1);
                   
                        ripple_bias = nan(1, nbins);

                        ripple_bias(1:length(idx)) = bias_vals(idx);
                        % Mask if another ripple starts within this window
                        all_starts = ripples_all(nprobe).onset(ripple_idx_all);
                        next_ripple = all_starts > t_start & all_starts < (t_start + 1);
                        if any(next_ripple)
                            t_next = min(all_starts(next_ripple));
                            mask_start = round((t_next - t_start) / psth_step) + 101;
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
KDE_reactivation_ripples_PSTH.HPC_ripples             = [ripple_bias_masked_SWS.HPC.bias{1};             ripple_bias_masked_SWS.HPC.bias{2}];
KDE_reactivation_ripples_PSTH.V1_ripples              = [ripple_bias_masked_SWS.V1.bias{1};              ripple_bias_masked_SWS.V1.bias{2}];
KDE_reactivation_ripples_PSTH.HPC_z_ripples     = [ripple_bias_masked_SWS.HPC.zscored_bias{1};     ripple_bias_masked_SWS.HPC.zscored_bias{2}];
KDE_reactivation_ripples_PSTH.V1_z_ripples      = [ripple_bias_masked_SWS.V1.zscored_bias{1};      ripple_bias_masked_SWS.V1.zscored_bias{2}];
% KDE_reactivation_ripples_PSTH.HPC_zshuff_ripples    = [ripple_bias_masked_SWS.HPC.zscored_bias_shuffled{1}; ripple_bias_masked_SWS.HPC.zscored_bias_shuffled{2}];
% KDE_reactivation_ripples_PSTH.V1_zshuff_ripples     = [ripple_bias_masked_SWS.V1.zscored_bias_shuffled{1}; ripple_bias_masked_SWS.V1.zscored_bias_shuffled{2}];

% Non-SWS (awake/REM)
KDE_reactivation_ripples_PSTH.HPC_awake_ripples             = [ripple_bias_masked_nonSWS.HPC.bias{1};             ripple_bias_masked_nonSWS.HPC.bias{2}];
KDE_reactivation_ripples_PSTH.V1_awake_ripples              = [ripple_bias_masked_nonSWS.V1.bias{1};              ripple_bias_masked_nonSWS.V1.bias{2}];
KDE_reactivation_ripples_PSTH.HPC_z_awake_ripples     = [ripple_bias_masked_nonSWS.HPC.zscored_bias{1};     ripple_bias_masked_nonSWS.HPC.zscored_bias{2}];
KDE_reactivation_ripples_PSTH.V1_z_awake_ripples      = [ripple_bias_masked_nonSWS.V1.zscored_bias{1};      ripple_bias_masked_nonSWS.V1.zscored_bias{2}];
% KDE_reactivation_ripples_PSTH.HPC_zshuff_awake_ripples    = [ripple_bias_masked_nonSWS.HPC.zscored_bias_shuffled{1}; ripple_bias_masked_nonSWS.HPC.zscored_bias_shuffled{2}];
% KDE_reactivation_ripples_PSTH.V1_zshuff_awake_ripples     = [ripple_bias_masked_nonSWS.V1.zscored_bias_shuffled{1}; ripple_bias_masked_nonSWS.V1.zscored_bias_shuffled{2}];

save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_ripples_PSTH.mat'),'KDE_reactivation_ripples_PSTH');


