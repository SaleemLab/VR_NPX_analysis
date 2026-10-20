%% MUA-equalised event control: pooled reactivation coherence AUC (POST)
% Same conditions and AUC as reactivation_AUC_coherence_PRE (POST data), with
% Low and High events equalised for HC MUA and dominant-hemisphere V1 MUA within
% each session (as MUA_equalised_event_control), then pooled across sessions.
%   ripple_power : HPC->V1. Events selected by within-session ripple power
%                  percentile <= 25 (Low) vs > 75 (High); thresholds from |HPC| of
%                  each condition's events. Temporal V1 window centred on each bin.
%   spindle      : V1->HPC. Events selected by V1 bias sign and the spindle
%                  percentile (re-ranked within session) of the hemisphere matched
%                  to it (V1 > 0 -> hemisphere 2, V1 < 0 -> hemisphere 1), <= 25 vs > 75;
%                  thresholds from |V1| of all events. Temporal V1 window centred on each bin.
%   SO_phase     : V1->HPC as spindle, SO peak [-pi/2, pi/2] vs trough.
%   SO_trough_uni_bi : V1->HPC as spindle, dominant hemisphere in SO trough with
%                  the non-dominant hemisphere in peak (unilateral, 'Peak') vs
%                  also in trough (bilateral, 'Trough'), as hemisphere_conditions
%                  'SO_trough_bilateral'.
% HPC =mean z log odds 0-0.1 s; V1 = mean z log odds over the V1 window.
% MUA (z, 10 ms bins): HC 0-0.1 s; V1 of the dominant hemisphere over the V1 window
% (dominance from the V1 bias sign in that window, so matching is redone per window).
% Each iteration (seeded): one MUA-matched draw of Low and High events in every
% session (nQ x nQ joint quantile bins, min(nA, nB) kept per bin; sessions with
% < min_events matched events dropped), pooled, then one bootstrap resample of the
% pooled events. Shuffled control keeps the selector un-resampled while resampling
% the measure (as reactivation_coherence_auc). Unmatched AUC uses the same bootstrap
% without matching. Matched second - first condition difference (High - Low,
% Trough - Peak; bilateral - unilateral for SO_trough_uni_bi) is taken within each
% iteration; diff_p is one-sided (difference <= 0).
% Output: nature_neuroscience_revision_analysis\MUA_equalised_event_control\AUC,
% plus the low vs high bar plot (plot_AUC_low_vs_high_controls).
%   which : analyses to run (default all: 'ripple_power', 'spindle', 'SO_phase',
%           'SO_trough_uni_bi'); static results of the others are kept from the
%           saved MUA_equalised_AUC_static.mat
function MUA_equalised_AUC(which)
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    out_dir = fullfile(rev_dir, 'MUA_equalised_event_control', 'AUC');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    nIter = 1000;
    nQ = 5;
    min_events = 10;

    % --- Load non-duplicated POST data (rows aligned to ripple_info) ---
    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
    HC = (KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples + KDE_reactivation_ripples_PSTH.nan_mask)';
    V1 = (KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples + KDE_reactivation_ripples_PSTH.nan_mask)';
    clear KDE_reactivation_ripples_PSTH
    load(fullfile(rev_dir, 'ripple_MUA_HPC.mat'), 'ripple_MUA_HPC');
    load(fullfile(rev_dir, 'ripple_MUA_V1_L.mat'), 'ripple_MUA_V1_L');
    load(fullfile(rev_dir, 'ripple_MUA_V1_R.mat'), 'ripple_MUA_V1_R');
    n_ev = numel(ripple_info.session_id);
    assert(size(ripple_MUA_HPC, 1) == n_ev && size(HC, 2) == n_ev, 'MUA / PSTH rows do not match ripple_info');

    timebin = 0.01;
    bin_edges = -1:timebin:1;
    bin_centers = bin_edges(1:end-1) + timebin/2;
    mua_centers = -2 + timebin/2:timebin:2 - timebin/2;

    E.sess = ripple_info.session_id(:);
    E.HPC_b = mean(HC(bin_centers > 0 & bin_centers < 0.1, :), 'omitnan')';
    E.HC_MUA = double(mean(ripple_MUA_HPC(:, mua_centers > 0 & mua_centers < 0.1), 2));
    E.ripple_pct = within_session_percentile(ripple_info.ripple_power(:), E.sess);
    E.spindle_pct = within_session_percentile(ripple_info.spindle_percentile, E.sess);
    E.SO_phase = ripple_info.SO_phase;
    % Cumulative sums for fast window means (bins x events)
    cV1 = [zeros(1, n_ev); cumsum(fillmissing_zero(V1), 1)];
    nV1 = [zeros(1, n_ev); cumsum(isfinite(V1), 1)];
    cL = [zeros(n_ev, 1) cumsum(double(ripple_MUA_V1_L), 2)];
    cR = [zeros(n_ev, 1) cumsum(double(ripple_MUA_V1_R), 2)];
    clear HC V1 ripple_MUA_HPC ripple_MUA_V1_L ripple_MUA_V1_R

    style = AUC_low_high_style();
    an = struct( ...
        'name',      {'ripple_power', 'spindle', 'SO_phase', 'SO_trough_uni_bi'}, ...
        'title',     {'Ripple power, within session (HPC->V1)', 'Spindle-band power, matched hemisphere (V1->HPC)', ...
                      'SO phase, matched hemisphere (V1->HPC)', 'Dominant SO trough, non-dominant peak vs trough (V1->HPC)'}, ...
        'direction', {'HPC->V1', 'V1->HPC', 'V1->HPC', 'V1->HPC'}, ...
        'alignment', {'centred', 'centred', 'centred', 'centred'}, ...
        'colors',    {style.ripple_power.colors, style.spindle.colors, style.SO_phase.colors, style.SO_trough_uni_bi.colors});
    an_all = an;
    if nargin >= 1, an = an(ismember({an.name}, which)); end
    cfg = struct('nIter', nIter, 'nQ', nQ, 'min_events', min_events);

    %% Static AUC: pre and post V1 windows
    v1_windows = struct('name', {'pre', 'post'}, 'label', {'V1 -0.2-0 s (pre)', 'V1 0-0.2 s (post)'}, ...
        'range', {[-0.2 0], [0 0.2]});
    Tstat = table();
    res = struct();
    static_file = fullfile(out_dir, 'MUA_equalised_AUC_static.mat');
    if nargin >= 1 && exist(static_file, 'file')
        % keep the saved results of the analyses not rerun
        S = load(static_file, 'res', 'Tstat');
        res = S.res;
        rerun = intersect(fieldnames(res), {an.name});
        if ~isempty(rerun), res = rmfield(res, rerun); end
        Tstat = S.Tstat(~ismember(S.Tstat.analysis, {an.name}), :);
    end
    for iA = 1:numel(an)
        for w = 1:numel(v1_windows)
            W = window_data(E, cV1, nV1, cL, cR, bin_centers, mua_centers, v1_windows(w).range);
            r = matched_auc(an(iA), E, W, cfg);
            res.(an(iA).name).(v1_windows(w).name) = r;
            Tstat = [Tstat; summary_row(an(iA).name, v1_windows(w).name, r)]; %#ok<AGROW>
            fprintf('%s %s: matched %s %.3f vs %s %.3f, diff %.3f [%.3f %.3f], p = %.3f\n', an(iA).name, ...
                v1_windows(w).name, r.labels{1}, r.matched_mean(1), r.labels{2}, r.matched_mean(2), ...
                r.diff_mean, r.diff_CI(1), r.diff_CI(2), r.diff_p);
        end
    end
    disp(Tstat);
    save(static_file, 'res', 'Tstat', 'v1_windows', 'cfg');
    writetable(Tstat, fullfile(out_dir, 'MUA_equalised_AUC_static.csv'));
    plot_static(res, an_all(isfield(res, {an_all.name})), v1_windows, fullfile(out_dir, 'MUA_equalised_AUC_static'));
    plot_AUC_low_vs_high_controls('MUA_equalised');

    %% Temporal AUC: sliding V1 window (HPC fixed 0-0.1 s)
    tp.win_size = 0.1;
    tp.step_size = 0.02;
    tp.time_bins = -0.6:tp.step_size:0.6;
    for iA = 1:numel(an)
        fprintf('Temporal AUC: %s\n', an(iA).name);
        AUC = temporal_auc(an(iA), E, cV1, nV1, cL, cR, bin_centers, mua_centers, tp, cfg);
        save(fullfile(out_dir, ['MUA_equalised_AUC_temporal_' an(iA).name '.mat']), 'AUC');
        plot_temporal(AUC, an(iA), fullfile(out_dir, ['MUA_equalised_AUC_temporal_' an(iA).name]));
    end
    fprintf('Saved to %s\n', out_dir);
end

function X = fillmissing_zero(X)
    X(~isfinite(X)) = 0;
end

function W = window_data(E, cV1, nV1, cL, cR, bin_centers, mua_centers, range)
    % V1 bias, dominant hemisphere and its V1 MUA for one V1 window
    b = find(bin_centers > range(1) & bin_centers < range(2));
    m = find(mua_centers > range(1) & mua_centers < range(2));
    W.V1_b = ((cV1(b(end) + 1, :) - cV1(b(1), :)) ./ (nV1(b(end) + 1, :) - nV1(b(1), :)))';
    W.hemi = NaN(size(W.V1_b)); W.hemi(W.V1_b < 0) = 1; W.hemi(W.V1_b > 0) = 2;
    V1_L = (cL(:, m(end) + 1) - cL(:, m(1))) / numel(m);
    V1_R = (cR(:, m(end) + 1) - cR(:, m(1))) / numel(m);
    W.V1_MUA = NaN(size(W.V1_b));
    W.V1_MUA(W.hemi == 1) = V1_L(W.hemi == 1);
    W.V1_MUA(W.hemi == 2) = V1_R(W.hemi == 2);
end

function r = matched_auc(a, E, W, cfg)
    % Bootstrapped AUC per condition: MUA-matched, shuffled (matched) and unmatched
    [masks, r.labels] = condition_masks(a.name, E, W.hemi);
    valid = isfinite(W.V1_b) & isfinite(E.HPC_b) & isfinite(E.HC_MUA) & isfinite(W.V1_MUA);
    if strcmp(a.direction, 'HPC->V1')
        sel = E.HPC_b; meas = W.V1_b;
    else
        sel = W.V1_b; meas = E.HPC_b;
    end
    idxU = cell(1, 2); th = cell(1, 2);
    for c = 1:2
        idxU{c} = find(masks{c} & valid);
        if strcmp(a.direction, 'HPC->V1')
            ref = abs(sel(idxU{c}));            % thresholds from the condition's events
        else
            ref = abs(sel(isfinite(sel)));      % thresholds from all events
        end
        t = prctile(ref, 0:10:100); th{c} = t(1:end-1);
    end
    P = prepare_matching(idxU{1}, idxU{2}, E.sess, E.HC_MUA, W.V1_MUA, cfg.nQ, cfg.min_events);

    nIter = cfg.nIter;
    au_m = NaN(nIter, 2); au_s = NaN(nIter, 2); au_u = NaN(nIter, 2);
    parfor i = 1:nIter
        s = RandStream('philox4x32_10', 'Seed', i);
        idxM = draw_matched(P, s);
        am = NaN(1, 2); as = NaN(1, 2); au = NaN(1, 2);
        for c = 1:2
            m = idxM{c};
            if isempty(m), continue; end
            b = m(randi(s, numel(m), numel(m), 1));
            am(c) = auc_from_diff(bias_diff(sel(b), meas(b), th{c}), th{c});
            as(c) = auc_from_diff(bias_diff(sel(m), meas(b), th{c}), th{c});
            u = idxU{c};
            bu = u(randi(s, numel(u), numel(u), 1));
            au(c) = auc_from_diff(bias_diff(sel(bu), meas(bu), th{c}), th{c});
        end
        au_m(i, :) = am; au_s(i, :) = as; au_u(i, :) = au;
    end
    d = au_m(:, 2) - au_m(:, 1);
    du = au_u(:, 2) - au_u(:, 1);
    r.matched_mean = mean(au_m, 'omitnan');   r.matched_CI = prctile(au_m, [2.5 97.5]);
    r.shuffled_mean = mean(au_s, 'omitnan');  r.shuffled_CI = prctile(au_s, [2.5 97.5]);
    r.unmatched_mean = mean(au_u, 'omitnan'); r.unmatched_CI = prctile(au_u, [2.5 97.5]);
    r.diff_mean = mean(d, 'omitnan'); r.diff_CI = prctile(d, [2.5 97.5]);
    r.diff_p = mean(d <= 0);                  % one-sided bootstrap p (High/Trough > Low/Peak)
    r.diff_unmatched_mean = mean(du, 'omitnan'); r.diff_unmatched_CI = prctile(du, [2.5 97.5]);
    r.diff_unmatched_p = mean(du <= 0);
    r.p_vs_shuffled = mean(au_m <= au_s);     % per condition, matched real vs shuffled
    r.n_unmatched = cellfun(@numel, idxU);
    r.n_matched = sum(P.k);
    r.n_sessions_matched = numel(unique(P.session));
end

function [masks, labels] = condition_masks(name, E, hemi)
    switch name
        case 'ripple_power'
            masks = {E.ripple_pct <= 25, E.ripple_pct > 75};
            labels = {'Low', 'High'};
        case 'spindle'
            f = dominant_feature(E.spindle_pct, hemi);
            masks = {f <= 25, f > 75};
            labels = {'Low', 'High'};
        case 'SO_phase'
            f = dominant_feature(E.SO_phase, hemi);
            masks = {peak(f), trough(f)};
            labels = {'Peak', 'Trough'};
        case 'SO_trough_uni_bi'
            % dominant (f) in trough; non-dominant (g) in peak (unilateral) or trough (bilateral)
            f = dominant_feature(E.SO_phase, hemi);
            g = dominant_feature(E.SO_phase, 3 - hemi);
            masks = {trough(f) & peak(g), trough(f) & trough(g)};
            labels = {'Peak', 'Trough'};   % non-dominant hemisphere's phase
    end
end

function m = peak(x)
    m = x >= -pi/2 & x <= pi/2;
end

function m = trough(x)
    m = (x >= -pi & x < -pi/2) | (x > pi/2 & x <= pi);
end

function f = dominant_feature(F, hemi)
    f = NaN(size(hemi));
    f(hemi == 1) = F(hemi == 1, 1);
    f(hemi == 2) = F(hemi == 2, 2);
end

function pct = within_session_percentile(X, session_ids)
    % Percentile (0-100) of each column of X within each session
    pct = NaN(size(X));
    for s = unique(session_ids)'
        idx = session_ids == s;
        for c = 1:size(X, 2)
            x = X(idx, c);
            ok = isfinite(x);
            p = NaN(size(x));
            p(ok) = (tiedrank(x(ok)) - 0.5) ./ sum(ok) * 100;
            pct(idx, c) = p;
        end
    end
end

function P = prepare_matching(ia, ib, sess, x1, x2, nQ, min_events)
    % Per session joint (x1, x2) quantile bins; keep bins with events in both sets
    P.ea = {}; P.eb = {}; P.k = []; P.session = [];
    for s = unique(sess([ia; ib]))'
        a = ia(sess(ia) == s); b = ib(sess(ib) == s);
        if isempty(a) || isempty(b), continue; end
        pooled = [a; b];
        ca = joint_bin(x1, x2, pooled, a, nQ);
        cb = joint_bin(x1, x2, pooled, b, nQ);
        ea = {}; eb = {}; k = [];
        for c = unique([ca; cb])'
            n = min(sum(ca == c), sum(cb == c));
            if n > 0
                ea{end+1} = a(ca == c); eb{end+1} = b(cb == c); k(end+1) = n; %#ok<AGROW>
            end
        end
        if sum(k) >= min_events
            P.ea = [P.ea ea]; P.eb = [P.eb eb]; P.k = [P.k k];
            P.session = [P.session s * ones(1, numel(k))];
        end
    end
end

function idx = draw_matched(P, s)
    sa = cell(numel(P.k), 1); sb = sa;
    for j = 1:numel(P.k)
        sa{j} = P.ea{j}(randperm(s, numel(P.ea{j}), P.k(j)));
        sb{j} = P.eb{j}(randperm(s, numel(P.eb{j}), P.k(j)));
    end
    idx = {vertcat(sa{:}), vertcat(sb{:})};
end

function cid = joint_bin(x1, x2, pooled, idx, nQ)
    cid = (quantile_bin(x1(pooled), x1(idx), nQ) - 1) * nQ + quantile_bin(x2(pooled), x2(idx), nQ);
end

function q = quantile_bin(ref, x, nQ)
    edges = unique(prctile(ref, linspace(0, 100, nQ + 1)));
    edges(1) = -Inf; edges(end) = Inf;
    q = discretize(x, edges);
end

function d = bias_diff(sel, meas, th)
    % mean(meas | sel >= th) - mean(meas | sel <= -th), as reactivation_coherence_auc
    d = NaN(1, numel(th));
    for j = 1:numel(th)
        t1 = sel >= th(j); t2 = sel <= -th(j);
        if any(t1) && any(t2)
            d(j) = mean(meas(t1), 'omitnan') - mean(meas(t2), 'omitnan');
        end
    end
end

function a = auc_from_diff(d, th)
    a = trapz(th, d) / (max(th) - min(th));
end

function row = summary_row(name, window, r)
    row = table(string(name), string(window), string(r.labels{1}), string(r.labels{2}), ...
        r.unmatched_mean(1), r.unmatched_mean(2), r.matched_mean(1), r.matched_CI(1, 1), r.matched_CI(2, 1), ...
        r.matched_mean(2), r.matched_CI(1, 2), r.matched_CI(2, 2), r.shuffled_mean(1), r.shuffled_mean(2), ...
        r.p_vs_shuffled(1), r.p_vs_shuffled(2), r.diff_unmatched_mean, r.diff_unmatched_p, ...
        r.diff_mean, r.diff_CI(1), r.diff_CI(2), r.diff_p, r.n_unmatched(1), r.n_unmatched(2), r.n_matched, r.n_sessions_matched, ...
        'VariableNames', {'analysis', 'V1_window', 'cond_A', 'cond_B', 'AUC_A_unmatched', 'AUC_B_unmatched', ...
        'AUC_A_matched', 'AUC_A_CI_lo', 'AUC_A_CI_hi', 'AUC_B_matched', 'AUC_B_CI_lo', 'AUC_B_CI_hi', ...
        'AUC_A_shuffled', 'AUC_B_shuffled', 'p_A_vs_shuffled', 'p_B_vs_shuffled', 'diff_unmatched', 'p_diff_unmatched', ...
        'diff_matched', 'diff_CI_lo', 'diff_CI_hi', 'p_diff_matched', 'n_A_unmatched', 'n_B_unmatched', ...
        'n_per_cond_matched', 'n_sessions_matched'});
end

function AUC = temporal_auc(a, E, cV1, nV1, cL, cR, bin_centers, mua_centers, tp, cfg)
    nT = numel(tp.time_bins);
    f = {'matched_mean', 'shuffled_mean', 'unmatched_mean'};
    for k = 1:numel(f), AUC.(f{k}) = NaN(nT, 2); end
    fci = {'matched_CI', 'shuffled_CI', 'unmatched_CI'};
    for k = 1:numel(fci), AUC.(fci{k}) = NaN(nT, 2, 2); end
    AUC.diff_mean = NaN(nT, 1); AUC.diff_CI = NaN(nT, 2); AUC.diff_p = NaN(nT, 1);
    AUC.diff_unmatched_mean = NaN(nT, 1); AUC.n_matched = NaN(nT, 1);
    for t = 1:nT
        if strcmp(a.alignment, 'centred')
            t0 = tp.time_bins(t) - tp.win_size/2;
        else
            t0 = tp.time_bins(t);
        end
        W = window_data(E, cV1, nV1, cL, cR, bin_centers, mua_centers, [t0 - 1e-9, t0 + tp.win_size - 1e-9]);
        r = matched_auc(a, E, W, cfg);
        for k = 1:numel(f), AUC.(f{k})(t, :) = r.(f{k}); end
        for k = 1:numel(fci), AUC.(fci{k})(t, :, :) = permute(r.(fci{k}), [3 2 1]); end
        AUC.diff_mean(t) = r.diff_mean; AUC.diff_CI(t, :) = r.diff_CI; AUC.diff_p(t) = r.diff_p;
        AUC.diff_unmatched_mean(t) = r.diff_unmatched_mean; AUC.n_matched(t) = r.n_matched;
    end
    AUC.labels = r.labels;
    AUC.direction = a.direction;
    AUC.window_alignment = a.alignment;
    AUC.time_bins = tp.time_bins;
    AUC.win_size = tp.win_size;
    AUC.step_size = tp.step_size;
    AUC.cfg = cfg;
end

function plot_static(res, an, v1_windows, save_path)
    % Per condition: shuffled (grey), unmatched (outline), MUA-matched (filled) with 95% CI
    fig = figure('Name', 'MUA-equalised AUC', 'Position', [100 50 860 330 * numel(an)]);
    tiledlayout(numel(an), numel(v1_windows), 'TileSpacing', 'compact');
    bw = 0.22;
    for iA = 1:numel(an)
        cond_colors = num2cell(an(iA).colors, 2);
        for w = 1:numel(v1_windows)
            ax = nexttile; hold on;
            r = res.(an(iA).name).(v1_windows(w).name);
            for c = 1:2
                bar(c - bw, r.shuffled_mean(c), bw, 'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
                errorbar(c - bw, r.shuffled_mean(c), r.shuffled_mean(c) - r.shuffled_CI(1, c), r.shuffled_CI(2, c) - r.shuffled_mean(c), ...
                    'Color', [0.4 0.4 0.4], 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 5);
                bar(c, r.unmatched_mean(c), bw, 'FaceColor', 'none', 'EdgeColor', cond_colors{c}, 'LineWidth', 1.5);
                errorbar(c, r.unmatched_mean(c), r.unmatched_mean(c) - r.unmatched_CI(1, c), r.unmatched_CI(2, c) - r.unmatched_mean(c), ...
                    'Color', cond_colors{c}, 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 5);
                bar(c + bw, r.matched_mean(c), bw, 'FaceColor', cond_colors{c}, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
                errorbar(c + bw, r.matched_mean(c), r.matched_mean(c) - r.matched_CI(1, c), r.matched_CI(2, c) - r.matched_mean(c), ...
                    'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 6);
            end
            yline(0, 'k-');
            set(ax, 'XTick', [1 2], 'XTickLabel', r.labels, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 11);
            xlim([0.5 2.5]); ylabel('Reactivation coherence AUC');
            title({an(iA).title, v1_windows(w).label, ...
                sprintf('%s-%s matched %.3f [%.3f %.3f], p = %.3f', r.labels{2}, r.labels{1}, ...
                r.diff_mean, r.diff_CI(1), r.diff_CI(2), r.diff_p), ...
                sprintf('(unmatched %.3f, p = %.3f)', r.diff_unmatched_mean, r.diff_unmatched_p)}, ...
                'FontSize', 8, 'FontWeight', 'normal');
        end
    end
    sgtitle('Pooled AUC, within-session MUA-equalised (grey: shuffled, outline: unmatched, filled: matched)', 'FontSize', 11);
    savefig(fig, [save_path '.fig']);
    export_png(fig, [save_path '.png']);
end

function plot_temporal(AUC, a, save_path)
    % Tiles: each condition (matched + shuffled, unmatched dashed) and matched difference
    cond_colors = num2cell(a.colors, 2);
    tvec = AUC.time_bins;
    if strcmp(AUC.window_alignment, 'centred')
        xlab = 'V1 window centre (s rel. ripple onset)';
    else
        xlab = 'V1 window start (s rel. ripple onset)';
    end
    fig = figure('Name', ['MUA-equalised temporal ' a.title], 'Position', [100 100 420 760]);
    tiledlayout(3, 1, 'TileSpacing', 'compact');
    for c = 1:2
        ax(c) = nexttile; hold on; %#ok<AGROW>
        ok = isfinite(AUC.matched_mean(:, c))';
        ci = squeeze(AUC.matched_CI(ok, c, :)); ci_s = squeeze(AUC.shuffled_CI(ok, c, :));
        fill([tvec(ok) fliplr(tvec(ok))], [ci_s(:, 1)' fliplr(ci_s(:, 2)')], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
        plot(tvec(ok), AUC.shuffled_mean(ok, c), 'k', 'LineWidth', 1.2);
        fill([tvec(ok) fliplr(tvec(ok))], [ci(:, 1)' fliplr(ci(:, 2)')], cond_colors{c}, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        h1 = plot(tvec(ok), AUC.matched_mean(ok, c), 'Color', cond_colors{c}, 'LineWidth', 2);
        h2 = plot(tvec(ok), AUC.unmatched_mean(ok, c), '--', 'Color', cond_colors{c}, 'LineWidth', 1);
        yline(0, '--r'); xline(0, '--k');
        ylabel('Reactivation coherence AUC'); title(AUC.labels{c});
        if c == 1, legend([h1 h2], {'MUA-matched', 'unmatched'}, 'Box', 'off', 'Location', 'best'); end
        set(ax(c), 'TickDir', 'out', 'Box', 'off', 'FontSize', 11);
    end
    linkaxes(ax, 'y');
    ax3 = nexttile; hold on;
    ok = isfinite(AUC.diff_mean)';
    fill([tvec(ok) fliplr(tvec(ok))], [AUC.diff_CI(ok, 1)' fliplr(AUC.diff_CI(ok, 2)')], cond_colors{2}, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(tvec(ok), AUC.diff_mean(ok), 'Color', cond_colors{2}, 'LineWidth', 2);
    plot(tvec(ok), AUC.diff_unmatched_mean(ok), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
    sig = ok & (AUC.diff_CI(:, 1)' > 0);
    if any(sig)
        yl = ylim; plot(tvec(sig), yl(2) * ones(1, sum(sig)), 'k.', 'MarkerSize', 8);
    end
    yline(0, 'k-'); xline(0, '--k');
    ylabel(sprintf('AUC %s - %s', AUC.labels{2}, AUC.labels{1}));
    title('MUA-matched difference (dashed: unmatched; dots: 95% CI > 0)', 'FontWeight', 'normal', 'FontSize', 9);
    set(ax3, 'TickDir', 'out', 'Box', 'off', 'FontSize', 11);
    xlabel(xlab);
    set([ax ax3], 'XLim', [-0.5 0.5]);
    sgtitle(['MUA-equalised temporal: ' a.title], 'FontSize', 11);
    savefig(fig, [save_path '.fig']);
    export_png(fig, [save_path '.png']);
end

function export_png(fig, png_path)
    % PNG plus vector PDF (same name). Retry: files can be briefly locked by a
    % viewer / sync client
    [p, f] = fileparts(png_path);
    for k = 1:5
        try
            exportgraphics(fig, png_path, 'Resolution', 300);
            exportgraphics(fig, fullfile(p, [f '.pdf']), 'ContentType', 'vector');
            return
        catch err
            fprintf('PNG export failed (%s), retry %d/5\n', err.message, k);
            pause(5);
        end
    end
    rethrow(err);
end
