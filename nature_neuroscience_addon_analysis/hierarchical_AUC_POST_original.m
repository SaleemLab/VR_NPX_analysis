%% Reactivation coherence AUC with hierarchical (session) bootstrap, original POST decoder
% Original POST V1 and HPC z log odds (KDE_reactivation_ripples_PSTH_POST_unique.mat,
% non-duplicated events), same AUC as reactivation_AUC_coherence_PRE /
% reactivation_coherence_auc, but resampling respects events nested in sessions:
%   each iteration draws the 22 sessions with replacement, then events with
%   replacement within each drawn session; all conditions of an analysis use the
%   same resample (paired differences). Shuffled control: HPC->V1 / V1->HPC measure
%   permuted across events within each drawn session (selector and conditions kept).
% Analyses:
%   all_ripples  (HPC->V1) : all events; selector HPC (0-0.1 s), measure V1 window
%   ripple_power (HPC->V1) : within-session ripple power quartiles; thresholds
%                            from |HPC| of each quartile's events (as original)
%   spindle      (V1->HPC) : dominant-hemisphere (V1 bias sign) spindle percentile
%                            quartiles (spindle_percentile, 10-100: bottom
%                            10% excluded as in the original paper)
%   SO_phase     (V1->HPC) : dominant-hemisphere SO peak [-pi/2, pi/2] vs trough
%   SO_trough_uni_bi (V1->HPC): trough bilateral (both trough) vs unilateral
%                            (dominant trough, non-dominant peak)
%   SO_peak_uni_bi   (V1->HPC): peak bilateral vs unilateral (dominant peak,
%                            non-dominant trough)
%   V1->HPC thresholds from |V1| of all events (as original).
% p: one-sided fraction of iterations with (last - first condition) <= 0;
% two-sided = 2 * min(p, 1 - p).
% Static: V1 -0.2-0 s and 0-0.2 s. Temporal: 0.1 s V1 window, centred on each bin
% for HPC->V1 and starting at each bin for V1->HPC (as original), -0.5:0.02:0.5 s.
% Output: nature_neuroscience_revision_analysis\POST_original_AUC_hierarchical
%   which : analyses to run (default all)
function hierarchical_AUC_POST_original(which)
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    out_dir = fullfile(rev_dir, 'POST_original_AUC_hierarchical');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    nBoot = 1000;

    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
    P = KDE_reactivation_ripples_PSTH; clear KDE_reactivation_ripples_PSTH
    HC = (P.HPC_z_logodds_ripples + P.nan_mask)';
    V1 = (P.V1_z_logodds_ripples + P.nan_mask)';
    clear P

    bc = (-1:0.01:1); bc = bc(1:end-1) + 0.005;
    HPC_b = mean(HC(bc > 0 & bc < 0.1, :), 'omitnan');
    sessions = unique(ripple_info.session_id(:))';
    sess_events = arrayfun(@(s) find(ripple_info.session_id(:)' == s), sessions, 'UniformOutput', false);
    rp_q = within_session_quartile(ripple_info.ripple_power(:), ripple_info.session_id(:))';
    fprintf('%d sessions, %d events\n', numel(sessions), numel(HPC_b));

    an = struct('name', {}, 'title', {}, 'direction', {}, 'labels', {}, 'colors', {});
    an(end+1) = struct('name', 'all_ripples', 'title', 'All ripples (HPC->V1)', 'direction', 'HPC->V1', ...
        'labels', {{'All'}}, 'colors', [231 41 138] / 256);
    an(end+1) = struct('name', 'ripple_power', 'title', 'Ripple power, within-session quartiles (HPC->V1)', 'direction', 'HPC->V1', ...
        'labels', {{'Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)'}}, 'colors', cool(4));
    an(end+1) = struct('name', 'spindle', 'title', 'Spindle-band power, dominant hemisphere (V1->HPC)', 'direction', 'V1->HPC', ...
        'labels', {{'Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)'}}, 'colors', cool(4));
    an(end+1) = struct('name', 'SO_phase', 'title', 'SO phase, dominant hemisphere (V1->HPC)', 'direction', 'V1->HPC', ...
        'labels', {{'SO Peak', 'SO Trough'}}, 'colors', [0.2 0.6 0.8; 0.8 0.2 0.2]);
    an(end+1) = struct('name', 'SO_trough_uni_bi', 'title', 'SO trough (dominant): bilateral vs unilateral (V1->HPC)', 'direction', 'V1->HPC', ...
        'labels', {{'Trough bilateral', 'Trough unilateral'}}, 'colors', [0.95 0.6 0.6; 0.8 0.2 0.2]);
    an(end+1) = struct('name', 'SO_peak_uni_bi', 'title', 'SO peak (dominant): bilateral vs unilateral (V1->HPC)', 'direction', 'V1->HPC', ...
        'labels', {{'Peak bilateral', 'Peak unilateral'}}, 'colors', [0.6 0.8 0.95; 0.2 0.6 0.8]);
    if nargin >= 1, an = an(ismember({an.name}, which)); end

    %% Static AUC
    v1_windows = struct('name', {'pre', 'post'}, 'label', {'V1 -0.2-0 s (pre)', 'V1 0-0.2 s (post)'}, ...
        'range', {[-0.2 0], [0 0.2]});
    for iA = 1:numel(an)
        res = struct(); T = table();
        for w = 1:numel(v1_windows)
            V1_b = mean(V1(bc > v1_windows(w).range(1) & bc < v1_windows(w).range(2), :), 'omitnan');
            [sel, meas, C, th] = analysis_inputs(an(iA), V1_b, HPC_b, ripple_info, rp_q);
            r = condition_auc(sel, meas, C, th, sess_events, nBoot);
            res.(v1_windows(w).name) = r;
            for c = 1:numel(C)
                T = [T; table(string(an(iA).name), string(an(iA).direction), string(v1_windows(w).name), string(an(iA).labels{c}), ...
                    r.n(c), r.AUC_mean(c), r.AUC_CI(1, c), r.AUC_CI(2, c), r.AUC_shuffled_mean(c), r.AUC_shuffled_CI(1, c), ...
                    r.AUC_shuffled_CI(2, c), r.p_vs_shuffled(c), r.diff_mean, r.diff_CI(1), r.diff_CI(2), r.diff_p, r.diff_p_two, ...
                    'VariableNames', {'analysis', 'direction', 'V1_window', 'condition', 'n_events', 'AUC_mean', 'AUC_CI_lo', ...
                    'AUC_CI_hi', 'AUC_shuffled_mean', 'AUC_shuffled_CI_lo', 'AUC_shuffled_CI_hi', 'p_vs_shuffled', ...
                    'diff_last_minus_first', 'diff_CI_lo', 'diff_CI_hi', 'diff_p_one_sided', 'diff_p_two_sided'})]; %#ok<AGROW>
            end
            fprintf('%s %s: AUC %s (p vs shuffle %s) | %s - %s = %.3f [%.3f %.3f], p = %.3f\n', an(iA).name, v1_windows(w).name, ...
                sprintf('%.3f ', r.AUC_mean), sprintf('%.3f ', r.p_vs_shuffled), an(iA).labels{end}, an(iA).labels{1}, ...
                r.diff_mean, r.diff_CI(1), r.diff_CI(2), r.diff_p);
        end
        save(fullfile(out_dir, ['hier_AUC_static_' an(iA).name '.mat']), 'res', 'v1_windows');
        writetable(T, fullfile(out_dir, ['hier_AUC_static_' an(iA).name '.csv']));
        plot_static(res, an(iA), v1_windows, fullfile(out_dir, ['hier_AUC_static_' an(iA).name]));
    end

    %% Temporal AUC
    time_bins = -0.5:0.02:0.5; win = 0.1;
    for iA = 1:numel(an)
        fprintf('Temporal %s\n', an(iA).name);
        nC = numel(an(iA).labels);
        A.mean = NaN(numel(time_bins), nC); A.ci = NaN(numel(time_bins), nC, 2);
        A.shuffled_mean = A.mean; A.shuffled_ci = A.ci; A.p_vs_shuffled = A.mean;
        A.diff_mean = NaN(numel(time_bins), 1); A.diff_ci = NaN(numel(time_bins), 2); A.diff_p = A.diff_mean;
        for t = 1:numel(time_bins)
            if strcmp(an(iA).direction, 'HPC->V1')
                t0 = time_bins(t) - win / 2;   % centred, as original HPC->V1
            else
                t0 = time_bins(t);             % starting, as original V1->HPC
            end
            V1_b = mean(V1(bc >= t0 & bc < t0 + win, :), 'omitnan');
            [sel, meas, C, th] = analysis_inputs(an(iA), V1_b, HPC_b, ripple_info, rp_q);
            r = condition_auc(sel, meas, C, th, sess_events, nBoot);
            A.mean(t, :) = r.AUC_mean; A.ci(t, :, :) = permute(r.AUC_CI, [3 2 1]);
            A.shuffled_mean(t, :) = r.AUC_shuffled_mean; A.shuffled_ci(t, :, :) = permute(r.AUC_shuffled_CI, [3 2 1]);
            A.p_vs_shuffled(t, :) = r.p_vs_shuffled;
            A.diff_mean(t) = r.diff_mean; A.diff_ci(t, :) = r.diff_CI; A.diff_p(t) = r.diff_p;
        end
        A.labels = an(iA).labels; A.time_bins = time_bins; A.win_size = win; A.direction = an(iA).direction;
        AUC = A;
        save(fullfile(out_dir, ['hier_AUC_temporal_' an(iA).name '.mat']), 'AUC');
        plot_temporal(AUC, an(iA), fullfile(out_dir, ['hier_AUC_temporal_' an(iA).name]));
    end
    fprintf('Saved to %s\n', out_dir);
end

function [sel, meas, C, th] = analysis_inputs(a, V1_b, HPC_b, info, rp_q)
% Selector, measure, condition masks and per-condition thresholds
if strcmp(a.direction, 'HPC->V1')
    sel = HPC_b; meas = V1_b;
    if strcmp(a.name, 'all_ripples')
        C = {true(size(sel))};
    else
        C = arrayfun(@(q) rp_q == q, 1:4, 'UniformOutput', false);
    end
    th = cellfun(@(c) thresholds(abs(sel(c))), C, 'UniformOutput', false);   % per condition
else
    sel = V1_b; meas = HPC_b;
    dom = NaN(size(V1_b)); dom(V1_b < 0) = 1; dom(V1_b > 0) = 2;
    switch a.name
        case 'spindle', Fv = info.spindle_percentile;   % bottom 10% NaN, as original
        otherwise, Fv = info.SO_phase;
    end
    f = NaN(size(V1_b)); g = f;
    f(dom == 1) = Fv(dom == 1, 1); g(dom == 1) = Fv(dom == 1, 2);
    f(dom == 2) = Fv(dom == 2, 2); g(dom == 2) = Fv(dom == 2, 1);
    peak = @(x) x >= -pi/2 & x <= pi/2;
    trough = @(x) (x >= -pi & x < -pi/2) | (x > pi/2 & x <= pi);
    switch a.name
        case 'spindle', C = {f <= 25, f > 25 & f <= 50, f > 50 & f <= 75, f > 75};
        case 'SO_phase', C = {peak(f), trough(f)};
        case 'SO_trough_uni_bi', C = {trough(f) & trough(g), trough(f) & peak(g)};
        case 'SO_peak_uni_bi', C = {peak(f) & peak(g), peak(f) & trough(g)};
    end
    th = repmat({thresholds(abs(sel))}, 1, numel(C));                           % all events
end
end

function th = thresholds(x)
th = prctile(x, 0:10:100); th = th(1:end-1);
end

function q = within_session_quartile(metric, session_ids)
% Quartile (1-4) of the metric within each session (as reactivation_AUC_coherence_PRE)
q = NaN(size(metric));
for s = unique(session_ids)'
    idx = find(session_ids == s);
    pcts = prctile(metric(idx), 0:25:100);
    for iBin = 1:4
        if iBin == 1
            in_q = metric(idx) >= pcts(iBin) & metric(idx) <= pcts(iBin+1);
        else
            in_q = metric(idx) > pcts(iBin) & metric(idx) <= pcts(iBin+1);
        end
        q(idx(in_q)) = iBin;
    end
end
end

function r = condition_auc(sel, meas, C, th, sess_events, nBoot)
% Hierarchical bootstrap AUC per condition (sessions, then events within session)
nC = numel(C); nS = numel(sess_events);
Cm = cat(1, C{:});
b_real = NaN(nBoot, nC); b_shuf = NaN(nBoot, nC);
parfor iBoot = 1:nBoot
    s = RandStream('philox4x32_10', 'Seed', iBoot);
    pick = randi(s, nS, nS, 1);
    idx = cell(1, nS); idx_perm = idx;
    for k = 1:nS
        ev = sess_events{pick(k)};
        e = ev(randi(s, numel(ev), numel(ev), 1));
        idx{k} = e;
        idx_perm{k} = e(randperm(s, numel(e)));   % within-session shuffle of the measure
    end
    idx = [idx{:}]; idx_perm = [idx_perm{:}];
    a = NaN(1, nC); as = a;
    for c = 1:nC
        a(c) = auc(sel(idx), meas(idx), Cm(c, idx), th{c});
        as(c) = auc(sel(idx), meas(idx_perm), Cm(c, idx), th{c});
    end
    b_real(iBoot, :) = a; b_shuf(iBoot, :) = as;
end
r.n = sum(Cm, 2)';
r.AUC_mean = mean(b_real, 1, 'omitnan'); r.AUC_CI = prctile(b_real, [2.5 97.5], 1);   % 2 x nC
r.AUC_shuffled_mean = mean(b_shuf, 1, 'omitnan'); r.AUC_shuffled_CI = prctile(b_shuf, [2.5 97.5], 1);
r.p_vs_shuffled = mean(b_real <= b_shuf);
if nC > 1
    d = b_real(:, end) - b_real(:, 1);             % last - first condition
    r.diff_boot = d;
    r.diff_mean = mean(d, 'omitnan'); r.diff_CI = prctile(d, [2.5 97.5]);
    r.diff_p = mean(d <= 0); r.diff_p_two = min(1, 2 * min(r.diff_p, 1 - r.diff_p));
else
    r.diff_boot = []; r.diff_mean = NaN; r.diff_CI = [NaN NaN]; r.diff_p = NaN; r.diff_p_two = NaN;
end
end

function a = auc(sel, meas, c, th)
d = NaN(1, numel(th));
for j = 1:numel(th)
    t1 = sel >= th(j) & c; t2 = sel <= -th(j) & c;
    if any(t1) && any(t2)
        d(j) = mean(meas(t1), 'omitnan') - mean(meas(t2), 'omitnan');
    end
end
a = trapz(th, d) / (max(th) - min(th));
end

function plot_static(res, a, v1_windows, save_path)
nC = numel(a.labels);
fig = figure('Name', a.title, 'Position', [100 100 260 + 160 * nC * numel(v1_windows) 380]);
tiledlayout(1, numel(v1_windows), 'TileSpacing', 'compact');
for w = 1:numel(v1_windows)
    ax = nexttile; hold on;
    r = res.(v1_windows(w).name);
    for c = 1:nC
        bar(c - 0.2, r.AUC_shuffled_mean(c), 0.35, 'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        errorbar(c - 0.2, r.AUC_shuffled_mean(c), r.AUC_shuffled_mean(c) - r.AUC_shuffled_CI(1, c), ...
            r.AUC_shuffled_CI(2, c) - r.AUC_shuffled_mean(c), 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'CapSize', 5);
        bar(c + 0.2, r.AUC_mean(c), 0.35, 'FaceColor', a.colors(c, :), 'FaceAlpha', 0.7, 'EdgeColor', 'none');
        errorbar(c + 0.2, r.AUC_mean(c), r.AUC_mean(c) - r.AUC_CI(1, c), r.AUC_CI(2, c) - r.AUC_mean(c), ...
            'k', 'LineWidth', 1.5, 'CapSize', 6);
    end
    yline(0, 'k-');
    set(ax, 'XTick', 1:nC, 'XTickLabel', a.labels, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
    xlim([0.4 nC + 0.6]); ylabel('Reactivation coherence AUC');
    if nC > 1
        sub = sprintf('%s - %s: %.3f [%.3f %.3f], p = %.3f', a.labels{end}, a.labels{1}, r.diff_mean, r.diff_CI(1), r.diff_CI(2), r.diff_p);
    else
        sub = sprintf('p vs shuffle = %.3f', r.p_vs_shuffled);
    end
    title({v1_windows(w).label, sub}, 'FontSize', 9, 'FontWeight', 'normal');
end
linkaxes(findobj(fig, 'Type', 'axes'), 'y');
sgtitle({a.title, 'original POST decoder, hierarchical (session -> event) bootstrap (grey: shuffled)'}, 'FontSize', 11);
savefig(fig, [save_path '.fig']);
exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
end

function plot_temporal(A, a, save_path)
nC = numel(a.labels);
nT = nC + (nC > 1);
fig = figure('Name', [a.title ' temporal'], 'Position', [100 100 420 190 * nT + 60]);
tiledlayout(nT, 1, 'TileSpacing', 'compact');
t = A.time_bins;
for c = 1:nC
    nexttile; hold on;
    ci = squeeze(A.ci(:, c, :)); cs = squeeze(A.shuffled_ci(:, c, :));
    fill([t fliplr(t)], [cs(:, 1)' fliplr(cs(:, 2)')], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    plot(t, A.shuffled_mean(:, c), 'k', 'LineWidth', 1);
    fill([t fliplr(t)], [ci(:, 1)' fliplr(ci(:, 2)')], a.colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(t, A.mean(:, c), 'Color', a.colors(c, :), 'LineWidth', 2);
    xline(0, '--k'); yline(0, '--r'); xlim([-0.5 0.5]); ylim([-0.1 0.3]);
    title(a.labels{c}, 'FontWeight', 'normal', 'FontSize', 9); ylabel('AUC');
    set(gca, 'TickDir', 'out', 'Box', 'off');
end
if nC > 1
    nexttile; hold on;
    fill([t fliplr(t)], [A.diff_ci(:, 1)' fliplr(A.diff_ci(:, 2)')], [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(t, A.diff_mean, 'k', 'LineWidth', 2);
    sig = A.diff_ci(:, 1)' > 0;
    if any(sig), plot(t(sig), 0.2 * ones(1, sum(sig)), 'k.', 'MarkerSize', 8); end
    xline(0, '--k'); yline(0, 'k-'); xlim([-0.5 0.5]); ylim([-0.15 0.25]);
    title(sprintf('%s - %s (dots: 95%% CI > 0)', a.labels{end}, a.labels{1}), 'FontWeight', 'normal', 'FontSize', 9);
    ylabel('AUC difference');
    set(gca, 'TickDir', 'out', 'Box', 'off');
end
if strcmp(a.direction, 'HPC->V1'), xlabel('V1 window centre (s rel. ripple onset)'); else, xlabel('V1 window start (s rel. ripple onset)'); end
sgtitle({[a.title ': temporal'], 'original POST decoder, hierarchical bootstrap (black: within-session shuffled)'}, 'FontSize', 10);
savefig(fig, [save_path '.fig']);
exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
end
