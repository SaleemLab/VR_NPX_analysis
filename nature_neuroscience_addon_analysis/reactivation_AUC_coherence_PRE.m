%% V1-HC reactivation coherence AUC analysis for PRE data (or POST, same sessions)
% Uses non-duplicated ripple events (make_unique_ripple_data) with capped log odds.
% HPC->V1 (HPC 0-0.1 s selects Track 1/2, V1 bias difference measured):
%   all ripples, within-session ripple power quartiles
% V1->HPC (V1 bias selects, hemisphere-matched condition, HPC 0-0.1 s measured):
%   within-session spindle percentile quartiles, SO peak vs trough, and
%   dominant SO trough with non-dominant SO peak (unilateral) vs trough (bilateral)
% Static AUC uses a pre (-0.2-0 s) and post (0-0.2 s) V1 window; temporal AUC
% slides the V1 window, as in POST.
% sleep: 'PRE' (default) or 'POST_matched' (POST events from the 12 PRE
%   sessions only; spindle percentile of spindle_amplitude with the bottom 10%
%   excluded, as in the original paper), 'POST' (all 22 POST
%   sessions, same processing) or 'POST_rate_lt10' (as 'POST' with the PLS-KDE
%   log odds of decoders built from neurons < 10 Hz, assemble_PLS_KDE_rate_lt10_
%   PSTH_POST), 'POST_Bayes' / 'POST_Bayes_rate_lt10' (as 'POST' with the
%   multinomial Bayesian log odds of run_Bayesian_logodds_PSTH_POST, all cells /
%   cells < 10 Hz; windows without spikes = log odds 0); 'POST_Bayes_raw' /
%   'POST_Bayes_raw_rate_lt10' use the raw instead of the z-scored Bayesian log odds.
%   'awake': non-SWS POST ripples (make_unique_awake_ripple_data), all 22 sessions;
%   spindle percentile is within session among awake ripples (0-100). Adds a
%   temporal AUC of the ripple power quartiles.
%   Output folder = sleep.
% equal_n: if true, each session is randomly subsampled (fixed seed) to
%   min(PRE, POST) events of that session, so PRE and POST have the same number
%   of events per session. Ripple power quartiles and spindle percentiles are
%   assigned on all events of the session before subsampling. Output folder =
%   [sleep '_equal_n'].
% only_spindle: if true, run only the spindle sections (static and temporal)
%   and leave the other outputs in the folder unchanged (default false).
function reactivation_AUC_coherence_PRE(sleep, equal_n, only_spindle)
    if nargin < 1
        sleep = 'PRE';
    end
    if nargin < 2
        equal_n = false;
    end
    if nargin < 3
        only_spindle = false;
    end
    is_awake = strcmp(sleep, 'awake');
    if is_awake && (equal_n || only_spindle)
        error('equal_n and only_spindle are not available for awake ripples');
    end
    run_all = ~only_spindle;
    tag = sleep;
    if equal_n
        tag = [sleep '_equal_n'];
    end
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    out_dir = fullfile(rev_dir, tag);
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % --- 1. Load Data ---
    fprintf('Loading %s data...\n', sleep);
    load(fullfile(rev_dir, 'ripple_info_PRE_unique.mat'), 'ripple_info_PRE');
    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    pre_sessions = unique(ripple_info_PRE.session_id);
    info = struct();
    switch sleep
        case 'PRE'
            load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_PRE_unique.mat'), 'KDE_reactivation_ripples_PSTH_PRE');
            PSTH = KDE_reactivation_ripples_PSTH_PRE;
            keep = true(numel(ripple_info_PRE.session_id), 1);
            info.session_id = ripple_info_PRE.session_id(:);
            info.ripple_power = ripple_info_PRE.ripple_power(:);
            info.SO_phase = ripple_info_PRE.SO_phase;
            info.spindle_percentile = ripple_info_PRE.spindle_percentile;
        case {'POST_matched', 'POST', 'POST_rate_lt10', 'POST_Bayes', 'POST_Bayes_rate_lt10', 'POST_Bayes_raw', 'POST_Bayes_raw_rate_lt10'}
            switch erase(sleep, '_raw')
                case 'POST_rate_lt10'
                    load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique_rate_lt10.mat'), 'KDE_reactivation_ripples_PSTH');
                case 'POST_Bayes'
                    load(fullfile(rev_dir, 'Bayes_reactivation_ripples_PSTH_POST_unique.mat'), 'Bayes_reactivation_ripples_PSTH');
                    KDE_reactivation_ripples_PSTH = Bayes_reactivation_ripples_PSTH;
                case 'POST_Bayes_rate_lt10'
                    load(fullfile(rev_dir, 'Bayes_reactivation_ripples_PSTH_POST_unique_rate_lt10.mat'), 'Bayes_reactivation_ripples_PSTH');
                    KDE_reactivation_ripples_PSTH = Bayes_reactivation_ripples_PSTH;
                otherwise
                    load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
            end
            PSTH = KDE_reactivation_ripples_PSTH;
            if strcmp(sleep, 'POST_matched')
                keep = ismember(ripple_info.session_id(:), pre_sessions);
            else
                keep = true(numel(ripple_info.session_id), 1);   % all 22 POST sessions
            end
            info.session_id = ripple_info.session_id(keep);
            info.session_id = info.session_id(:);
            info.ripple_power = ripple_info.ripple_power(keep);
            info.ripple_power = info.ripple_power(:);
            info.SO_phase = ripple_info.SO_phase(keep, :);
            % POST spindle percentile of spindle_amplitude as in the original paper:
            % events below the 10th within-session percentile are NaN (excluded),
            % so Q1 covers the 10-25th percentile
            info.spindle_percentile = ripple_info.spindle_percentile(keep, :);
        case 'awake'
            load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_awake_unique.mat'), 'KDE_reactivation_ripples_PSTH_awake');
            load(fullfile(rev_dir, 'ripple_info_awake_unique.mat'), 'ripple_info_awake');
            PSTH = KDE_reactivation_ripples_PSTH_awake;
            clear KDE_reactivation_ripples_PSTH_awake
            keep = true(numel(ripple_info_awake.session_id), 1);
            info.session_id = ripple_info_awake.session_id(:);
            info.ripple_power = ripple_info_awake.ripple_power(:);
            info.SO_phase = ripple_info_awake.SO_phase;
            info.spindle_percentile = ripple_info_awake.spindle_percentile;
        otherwise
            error('Unknown sleep %s', sleep);
    end

    % Within-session ripple power quartiles on all events of the session
    q_masks = within_session_quartiles(info.ripple_power, info.session_id);

    % Equal number of PRE and POST events per session
    sub = true(numel(info.session_id), 1);
    if equal_n
        rng(1);
        sub = false(numel(info.session_id), 1);
        for s = pre_sessions(:)'
            n = min(sum(ripple_info_PRE.session_id == s), sum(ripple_info.session_id == s));
            idx = find(info.session_id == s);
            sub(idx(randperm(numel(idx), n))) = true;
        end
        info = structfun(@(x) x(sub, :), info, 'UniformOutput', false);
        q_masks = cellfun(@(m) m(sub'), q_masks, 'UniformOutput', false);
    end
    keep = find(keep);
    keep = keep(sub);
    fprintf('%s: %d events from %d sessions\n', tag, numel(keep), numel(unique(info.session_id)));

    % Track bias matrices (nBins x nEvents)
    % '_raw' modes: raw (not z-scored) log odds, so windows without spikes stay at 0
    if contains(sleep, '_raw')
        lo_field = '_logodds_ripples';
    else
        lo_field = '_z_logodds_ripples';
    end
    track_bias_HC = (PSTH.(['HPC' lo_field])(keep, :) + PSTH.nan_mask(keep, :))';
    track_bias_V1 = (PSTH.(['V1' lo_field])(keep, :) + PSTH.nan_mask(keep, :))';
    clear PSTH KDE_reactivation_ripples_PSTH KDE_reactivation_ripples_PSTH_PRE

    % Time parameters
    timebin = 0.01;
    bin_edges = -1:timebin:1;
    bin_centers = bin_edges(1:end-1) + timebin/2;

    % HPC window (fixed) and V1 windows (pre / post ripple)
    bins_HPC = bin_centers > 0 & bin_centers < 0.1;
    v1_windows = struct( ...
        'name',  {'pre', 'post'}, ...
        'label', {'V1 -0.2-0 s (pre)', 'V1 0-0.2 s (post)'}, ...
        'bins',  {bin_centers > -0.2 & bin_centers < 0, bin_centers > 0 & bin_centers < 0.2});

    nBoot = 1000;
    main_color = [231, 41, 138] / 256;
    mean_HPC = mean(track_bias_HC(bins_HPC, :), 'omitnan');
    nEv = numel(mean_HPC);

    %% SECTION 1: ALL RIPPLES (HPC->V1)
    all_set = struct('label', 'All Ripples', 'mask', true(1, nEv), 'c1', [], 'c2', []);
    if run_all
        fprintf('Processing All Ripples...\n');
        run_condition_set(tag, all_set, 'HPC->V1', main_color, sprintf('All %s Ripples (HPC->V1)', tag), 'AUC_all_events', ...
            out_dir, mean_HPC, track_bias_V1, v1_windows, nBoot);

        %% SECTION 2: RIPPLE POWER QUARTILES, within session (HPC->V1)
        fprintf('Processing Ripple Percentiles...\n');
        q_labels = {'Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)'};
        rp_set = struct('label', q_labels, 'mask', q_masks, 'c1', [], 'c2', []);
        run_condition_set(tag, rp_set, 'HPC->V1', cool(4), 'Ripple Power Percentiles, within session (HPC->V1)', 'AUC_ripple_percentile', ...
            out_dir, mean_HPC, track_bias_V1, v1_windows, nBoot);
    end

    %% SECTION 3: SPINDLE PERCENTILE QUARTILES, hemisphere-matched (V1->HPC)
    fprintf('Processing Spindle Percentiles...\n');
    sp_set = add_all_mask(hemisphere_conditions('spindle', info), nEv);
    run_condition_set(tag, sp_set, 'V1->HPC', cool(4), 'Spindle Percentiles, within session (V1->HPC)', 'AUC_spindle_percentile', ...
        out_dir, mean_HPC, track_bias_V1, v1_windows, nBoot);

    so_set = add_all_mask(hemisphere_conditions('SO', info), nEv);
    sob_set = add_all_mask(hemisphere_conditions('SO_trough_bilateral', info), nEv);
    sob_colors = [0.95 0.65 0.82; 0.8 0.1 0.5];
    if run_all
        %% SECTION 4: SO PHASE, hemisphere-matched (V1->HPC)
        fprintf('Processing SO Phase...\n');
        run_condition_set(tag, so_set, 'V1->HPC', [0.2 0.6 0.8; 0.8 0.2 0.2], 'SO Phase (V1->HPC)', 'AUC_SO_phase', ...
            out_dir, mean_HPC, track_bias_V1, v1_windows, nBoot);

        %% SECTION 4b: DOMINANT SO TROUGH, non-dominant hemisphere peak vs trough (V1->HPC)
        fprintf('Processing SO Trough (non-dominant peak vs trough)...\n');
        run_condition_set(tag, sob_set, 'V1->HPC', sob_colors, 'Dominant SO Trough, non-dominant Peak vs Trough (V1->HPC)', ...
            'AUC_SO_trough_bilateral', out_dir, mean_HPC, track_bias_V1, v1_windows, nBoot);
    end

    %% SECTION 5: TEMPORAL AUC
    % V1 window centred on each time bin (t +/- win_size/2) for all sets, as in
    %   corticohippocampal_reactivation_temporal_KDE_diff.m and
    %   reactivation_coherence_temporal_KDE_diff_V1_detection_SO_phase.m
    tp.win_size = 0.1;
    tp.step_size = 0.02;
    tp.time_bins = -1:tp.step_size:1;
    tp.bin_centers = bin_centers;

    if run_all
        fprintf('Processing Temporal AUC (all ripples)...\n');
        run_temporal_set(tag, all_set, 'HPC->V1', 'centred', main_color, sprintf('All %s Ripples (HPC->V1)', tag), ...
            ['KDE_temporal_AUC_all_ripples_' tag], out_dir, mean_HPC, track_bias_V1, tp, nBoot);
    end
    if is_awake
        fprintf('Processing Temporal AUC (ripple power)...\n');
        run_temporal_set(tag, rp_set, 'HPC->V1', 'centred', cool(4), 'Ripple Power Percentiles (HPC->V1)', ...
            ['KDE_temporal_AUC_ripple_percentile_' tag], out_dir, mean_HPC, track_bias_V1, tp, nBoot);
    end
    fprintf('Processing Temporal AUC (spindle)...\n');
    run_temporal_set(tag, sp_set, 'V1->HPC', 'centred', cool(4), 'Spindle Percentiles (V1->HPC)', ...
        ['KDE_temporal_AUC_spindle_percentile_' tag], out_dir, mean_HPC, track_bias_V1, tp, nBoot);
    if run_all
        fprintf('Processing Temporal AUC (SO)...\n');
        run_temporal_set(tag, so_set, 'V1->HPC', 'centred', [0.2 0.6 0.8; 0.8 0.2 0.2], 'SO Phase (V1->HPC)', ...
            ['KDE_temporal_AUC_SO_phase_' tag], out_dir, mean_HPC, track_bias_V1, tp, nBoot);
        fprintf('Processing Temporal AUC (SO trough, non-dominant peak vs trough)...\n');
        run_temporal_set(tag, sob_set, 'V1->HPC', 'centred', sob_colors, 'Dominant SO Trough (V1->HPC)', ...
            ['KDE_temporal_AUC_SO_trough_bilateral_' tag], out_dir, mean_HPC, track_bias_V1, tp, nBoot);
    end

    fprintf('All analyses finished successfully!\n');
end

function S = add_all_mask(C, nEv)
    % V1->HPC conditions use all events; the condition enters via c1/c2
    S = struct('label', {C.label}, 'mask', true(1, nEv), 'c1', {C.c1}, 'c2', {C.c2});
end

function q_masks = within_session_quartiles(metric, session_ids)
    % Logical mask per quartile (1-4) of the metric WITHIN each session
    nBins = 4;
    assigned_quantile = NaN(size(metric));
    unique_sessions = unique(session_ids);
    for i = 1:length(unique_sessions)
        idx = find(session_ids == unique_sessions(i));
        sess_metric = metric(idx);
        pcts = prctile(sess_metric, 0:25:100);
        for iBin = 1:nBins
            % First bin is closed on the left so the session minimum is included
            if iBin == 1
                q_mask = sess_metric >= pcts(iBin) & sess_metric <= pcts(iBin+1);
            else
                q_mask = sess_metric > pcts(iBin) & sess_metric <= pcts(iBin+1);
            end
            assigned_quantile(idx(q_mask)) = iBin;
        end
    end
    q_masks = arrayfun(@(q) (assigned_quantile(:) == q)', 1:nBins, 'UniformOutput', false);
end

function r = condition_auc(cond, direction, mean_HPC, mean_V1, nBoot)
    % Run reactivation_coherence_auc for one condition in the given direction
    m = cond.mask;
    c1 = cond.c1; c2 = cond.c2;
    if ~isempty(c1), c1 = c1(m); c2 = c2(m); end
    if strcmp(direction, 'HPC->V1')
        r = reactivation_coherence_auc(mean_HPC(m), mean_V1(m), nBoot, c1, c2);
    else
        r = reactivation_coherence_auc(mean_V1(m), mean_HPC(m), nBoot, c1, c2);
    end
end

function run_condition_set(sleep, conds, direction, colors, title_str, save_name, out_dir, mean_HPC, track_bias_V1, v1_windows, nBoot)
    % Static AUC for each condition and V1 window; save and plot paired bars
    nCond = numel(conds);
    res = struct('label', {conds.label});
    for iCond = 1:nCond
        for w = 1:numel(v1_windows)
            mean_V1 = mean(track_bias_V1(v1_windows(w).bins, :), 'omitnan');
            res(iCond).(v1_windows(w).name) = condition_auc(conds(iCond), direction, mean_HPC, mean_V1, nBoot);
        end
    end
    window_labels = {v1_windows.label};
    save(fullfile(out_dir, [save_name '.mat']), 'res', 'window_labels', 'direction');

    fig = figure('Name', [sleep ' ' title_str], 'Position', [640 100 300 + 200*nCond 400]);
    tiledlayout(1, numel(v1_windows), 'TileSpacing', 'compact');
    bar_width = 0.35; offset = 0.2;
    for w = 1:numel(v1_windows)
        ax = nexttile; hold on;
        for iCond = 1:nCond
            r = res(iCond).(v1_windows(w).name);
            if ~isfinite(r.AUC_mean), continue; end
            % Shuffled (grey) on the left, real (coloured) on the right
            bar(iCond - offset, r.AUC_mean_shuffled, bar_width, 'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
            errorbar(iCond - offset, r.AUC_mean_shuffled, r.AUC_mean_shuffled - r.AUC_CI_shuffled(1), r.AUC_CI_shuffled(2) - r.AUC_mean_shuffled, ...
                'o', 'Color', [0.4 0.4 0.4], 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerSize', 4, 'LineWidth', 1.2, 'CapSize', 6);
            bar(iCond + offset, r.AUC_mean, bar_width, 'FaceColor', colors(iCond,:), 'FaceAlpha', 0.7, 'EdgeColor', 'none');
            errorbar(iCond + offset, r.AUC_mean, r.AUC_mean - r.AUC_CI(1), r.AUC_CI(2) - r.AUC_mean, 'k', 'LineWidth', 1.5, 'CapSize', 8);
        end
        yline(0, 'k-');
        set(ax, 'XTick', 1:nCond, 'XTickLabel', {conds.label});
        xlim([0.4, nCond + 0.6]); ylabel('Reactivation coherence AUC');
        title(v1_windows(w).label);
        set(ax, 'TickDir', 'out', 'box', 'off', 'Color', 'none', 'FontSize', 12);
    end
    linkaxes(findobj(fig, 'Type', 'axes'), 'y');
    sgtitle([title_str '  (grey: shuffled)'], 'FontSize', 13);
    savefig(fig, fullfile(out_dir, [save_name '.fig']));
    exportgraphics(fig, fullfile(out_dir, [save_name '.png']), 'Resolution', 300);
    exportgraphics(fig, fullfile(out_dir, [save_name '.pdf']), 'ContentType', 'vector');
end

function run_temporal_set(sleep, conds, direction, alignment, colors, title_str, save_name, out_dir, mean_HPC, track_bias_V1, tp, nBoot)
    % Temporal AUC with a sliding V1 window (HPC fixed 0-0.1 s); one tile per condition
    nCond = numel(conds);
    nTime = numel(tp.time_bins);
    AUC.mean = nan(nTime, nCond);
    AUC.ci = nan(nTime, nCond, 2);
    AUC.shifted_mean = nan(nTime, nCond);
    AUC.shifted_ci = nan(nTime, nCond, 2);

    for t = 1:nTime
        if strcmp(alignment, 'centred')
            t0 = tp.time_bins(t) - tp.win_size/2;
        else
            t0 = tp.time_bins(t);
        end
        bins_V1 = tp.bin_centers >= t0 & tp.bin_centers < t0 + tp.win_size;
        if ~any(bins_V1), continue; end
        mean_V1 = mean(track_bias_V1(bins_V1, :), 'omitnan');
        for iCond = 1:nCond
            r = condition_auc(conds(iCond), direction, mean_HPC, mean_V1, nBoot);
            AUC.mean(t, iCond) = r.AUC_mean;
            AUC.ci(t, iCond, :) = r.AUC_CI;
            AUC.shifted_mean(t, iCond) = r.AUC_mean_shuffled;
            AUC.shifted_ci(t, iCond, :) = r.AUC_CI_shuffled;
        end
    end
    AUC.labels = {conds.label};
    AUC.direction = direction;
    AUC.window_alignment = alignment;
    AUC.time_bins = tp.time_bins;
    AUC.win_size = tp.win_size;
    AUC.step_size = tp.step_size;
    save(fullfile(out_dir, [save_name '.mat']), 'AUC');

    if strcmp(alignment, 'centred')
        xlab = 'V1 window centre (s rel. ripple onset)';
    else
        xlab = 'V1 window start (s rel. ripple onset)';
    end
    fig = figure('Name', [sleep ' temporal ' title_str], 'Position', [640 100 400 60 + 200*nCond]);
    tiledlayout(nCond, 1, 'TileSpacing', 'compact');
    for iCond = 1:nCond
        nexttile; hold on;
        ok = ~isnan(AUC.mean(:, iCond));
        tvec = tp.time_bins(ok);
        ci = squeeze(AUC.ci(ok, iCond, :));
        ci_s = squeeze(AUC.shifted_ci(ok, iCond, :));
        fill([tvec fliplr(tvec)], [ci_s(:,1)' fliplr(ci_s(:,2)')], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
        plot(tvec, AUC.shifted_mean(ok, iCond), 'k', 'LineWidth', 1.2);
        fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], colors(iCond,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        plot(tvec, AUC.mean(ok, iCond), 'Color', colors(iCond,:), 'LineWidth', 2);
        yline(0, '--r'); xline(0, '--k');
        ylabel('Reactivation coherence AUC');
        title(conds(iCond).label);
        set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 11);
        xlim([-0.5 0.5]);
    end
    xlabel(xlab);
    linkaxes(findobj(fig, 'Type', 'axes'), 'y');
    sgtitle([sleep ' temporal: ' title_str '  (black: shuffled)'], 'FontSize', 12);
    savefig(fig, fullfile(out_dir, [save_name '.fig']));
    exportgraphics(fig, fullfile(out_dir, [save_name '.png']), 'Resolution', 300);
    exportgraphics(fig, fullfile(out_dir, [save_name '.pdf']), 'ContentType', 'vector');
end
