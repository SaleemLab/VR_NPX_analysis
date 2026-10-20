%% Split-half V1->HPC reactivation coherence AUC with hierarchical (session) bootstrap (POST)
% Same analyses, conditions and AUC as split_half_V1_AUC_POST, but resampling
% respects the nesting of events within sessions (R1.6):
%   each iteration draws the 22 sessions with replacement, then, within each drawn
%   session, its events with replacement (as many as the session has); AUCs of all
%   conditions are computed on this pooled resample (same resample for all
%   conditions, so condition differences are paired).
%   Shuffled control: within each drawn session, HPC bias is permuted across that
%   session's resampled events (selector bias and conditions kept), which breaks
%   the V1-HPC pairing without mixing sessions.
%   Thresholds on |selector| are fixed (0:10:90 prctile over all events).
% Directions (dominance half -> selector half): Full->Full, A->B, B->A.
% Conditions: spindle quartiles (spindle_percentile, 10-100), SO peak vs trough,
% SO trough bilateral vs unilateral, SO peak bilateral vs unilateral (dominant
% hemisphere from the sign of the dominance half's V1 bias).
% p: one-sided fraction of iterations with difference (last - first condition)
% <= 0; two-sided = 2 * min(p, 1 - p).
% Static AUC: V1 -0.2-0 s and 0-0.2 s. Temporal AUC: 0.1 s V1 window starting at
% -0.5:0.02:0.5 s.
% Input: KDE_reactivation_ripples_PSTH_POST_unique_V1_split_half.mat
% Output: nature_neuroscience_revision_analysis\split_half_V1_AUC_POST_hierarchical
%   which : analyses to run (default all: 'spindle', 'SO_phase',
%           'SO_trough_uni_bi', 'SO_peak_uni_bi')
function split_half_V1_AUC_POST_hierarchical(which)
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    out_dir = fullfile(rev_dir, 'split_half_V1_AUC_POST_hierarchical');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    nBoot = 1000;

    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique_V1_split_half.mat'), 'KDE_reactivation_ripples_PSTH');
    P = KDE_reactivation_ripples_PSTH; clear KDE_reactivation_ripples_PSTH
    mask = P.nan_mask';
    HC = P.HPC_z_logodds_ripples' + mask;
    V1.A = P.V1A_z_logodds_ripples' + mask;
    V1.B = P.V1B_z_logodds_ripples' + mask;
    V1.full = P.V1full_z_logodds_ripples' + mask;

    bc = (-1:0.01:1); bc = bc(1:end-1) + 0.005;
    HPC_b = mean(HC(bc > 0 & bc < 0.1, :), 'omitnan');
    sessions = unique(ripple_info.session_id(:))';
    sess_events = arrayfun(@(s) find(ripple_info.session_id(:)' == s), sessions, 'UniformOutput', false);
    fprintf('%d sessions, %d events\n', numel(sessions), numel(HPC_b));

    dirs = struct('name', {'Full->Full', 'A->B', 'B->A'}, 'dom', {'full', 'A', 'B'}, 'sel', {'full', 'B', 'A'});
    an(1).name = 'spindle'; an(1).title = 'Spindle-band power (dominant hemisphere, within-session percentile)';
    an(1).labels = {'Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)'}; an(1).colors = cool(4);
    an(2).name = 'SO_phase'; an(2).title = 'SO phase (dominant hemisphere)';
    an(2).labels = {'SO Peak', 'SO Trough'}; an(2).colors = [0.2 0.6 0.8; 0.8 0.2 0.2];
    an(3).name = 'SO_trough_uni_bi'; an(3).title = 'SO trough (dominant): bilateral vs unilateral';
    an(3).labels = {'Trough bilateral', 'Trough unilateral'}; an(3).colors = [0.95 0.6 0.6; 0.8 0.2 0.2];
    an(4).name = 'SO_peak_uni_bi'; an(4).title = 'SO peak (dominant): bilateral vs unilateral';
    an(4).labels = {'Peak bilateral', 'Peak unilateral'}; an(4).colors = [0.6 0.8 0.95; 0.2 0.6 0.8];
    if nargin >= 1, an = an(ismember({an.name}, which)); end
    F.spindle = ripple_info.spindle_percentile;   % bottom 10% NaN, as original
    F.SO_phase = ripple_info.SO_phase;
    F.SO_trough_uni_bi = ripple_info.SO_phase;
    F.SO_peak_uni_bi = ripple_info.SO_phase;

    %% Static AUC
    v1_windows = struct('name', {'pre', 'post'}, 'label', {'V1 -0.2-0 s (pre)', 'V1 0-0.2 s (post)'}, ...
        'bins', {bc > -0.2 & bc < 0, bc > 0 & bc < 0.2});
    T = table();
    for iA = 1:numel(an)
        res = struct();
        for d = 1:numel(dirs)
            for w = 1:numel(v1_windows)
                dom_b = mean(V1.(dirs(d).dom)(v1_windows(w).bins, :), 'omitnan');
                sel_b = mean(V1.(dirs(d).sel)(v1_windows(w).bins, :), 'omitnan');
                C = condition_masks(an(iA).name, F.(an(iA).name), dom_b);
                r = condition_auc(sel_b, HPC_b, C, sess_events, nBoot);
                res.(matlab.lang.makeValidName(dirs(d).name)).(v1_windows(w).name) = r;
                for c = 1:numel(C)
                    T = [T; table(string(an(iA).name), string(dirs(d).name), string(v1_windows(w).name), string(an(iA).labels{c}), ...
                        r.n(c), r.AUC_mean(c), r.AUC_CI(1, c), r.AUC_CI(2, c), r.AUC_shuffled_mean(c), r.AUC_shuffled_CI(1, c), ...
                        r.AUC_shuffled_CI(2, c), r.p_vs_shuffled(c), r.diff_mean, r.diff_CI(1), r.diff_CI(2), r.diff_p, r.diff_p_two, ...
                        'VariableNames', {'analysis', 'direction', 'V1_window', 'condition', 'n_events', 'AUC_mean', 'AUC_CI_lo', ...
                        'AUC_CI_hi', 'AUC_shuffled_mean', 'AUC_shuffled_CI_lo', 'AUC_shuffled_CI_hi', 'p_vs_shuffled', ...
                        'diff_last_minus_first', 'diff_CI_lo', 'diff_CI_hi', 'diff_p_one_sided', 'diff_p_two_sided'})]; %#ok<AGROW>
                end
                fprintf('%s %s %s: %s | %s - %s = %.3f [%.3f %.3f], p = %.3f (two-sided %.3f)\n', an(iA).name, dirs(d).name, ...
                    v1_windows(w).name, sprintf('%.3f ', r.AUC_mean), an(iA).labels{end}, an(iA).labels{1}, r.diff_mean, ...
                    r.diff_CI(1), r.diff_CI(2), r.diff_p, r.diff_p_two);
            end
        end
        save(fullfile(out_dir, ['hier_AUC_static_' an(iA).name '.mat']), 'res', 'v1_windows', 'dirs');
        plot_static(res, an(iA), dirs, v1_windows, fullfile(out_dir, ['hier_AUC_static_' an(iA).name]));
        writetable(T(T.analysis == an(iA).name, :), fullfile(out_dir, ['hier_AUC_static_' an(iA).name '.csv']));
    end
    plot_AUC_low_vs_high_controls('split_half_hierarchical');   % low vs high bar plot of all four analyses

    %% Temporal AUC
    time_bins = -0.5:0.02:0.5; win = 0.1;
    for iA = 1:numel(an)
        AUC = struct();
        for d = 1:numel(dirs)
            fprintf('Temporal %s %s\n', an(iA).name, dirs(d).name);
            nC = numel(an(iA).labels);
            A.mean = NaN(numel(time_bins), nC); A.ci = NaN(numel(time_bins), nC, 2);
            A.shuffled_mean = A.mean; A.shuffled_ci = A.ci;
            A.diff_mean = NaN(numel(time_bins), 1); A.diff_ci = NaN(numel(time_bins), 2); A.diff_p = A.diff_mean;
            for t = 1:numel(time_bins)
                in = bc >= time_bins(t) & bc < time_bins(t) + win;
                dom_b = mean(V1.(dirs(d).dom)(in, :), 'omitnan');
                sel_b = mean(V1.(dirs(d).sel)(in, :), 'omitnan');
                C = condition_masks(an(iA).name, F.(an(iA).name), dom_b);
                r = condition_auc(sel_b, HPC_b, C, sess_events, nBoot);
                A.mean(t, :) = r.AUC_mean; A.ci(t, :, :) = permute(r.AUC_CI, [3 2 1]);
                A.shuffled_mean(t, :) = r.AUC_shuffled_mean; A.shuffled_ci(t, :, :) = permute(r.AUC_shuffled_CI, [3 2 1]);
                A.diff_mean(t) = r.diff_mean; A.diff_ci(t, :) = r.diff_CI; A.diff_p(t) = r.diff_p;
            end
            A.labels = an(iA).labels; A.time_bins = time_bins; A.win_size = win; A.window_alignment = 'start';
            AUC.(matlab.lang.makeValidName(dirs(d).name)) = A;
        end
        save(fullfile(out_dir, ['hier_AUC_temporal_' an(iA).name '.mat']), 'AUC', 'dirs');
        plot_temporal(AUC, an(iA), dirs, fullfile(out_dir, ['hier_AUC_temporal_' an(iA).name]));
    end
    fprintf('Saved to %s\n', out_dir);
end

function C = condition_masks(name, Fv, dom_b)
% Per-event condition masks from the dominant (f) and non-dominant (g)
% hemisphere's feature
dom = NaN(size(dom_b)); dom(dom_b < 0) = 1; dom(dom_b > 0) = 2;
f = NaN(size(dom_b)); g = f;
f(dom == 1) = Fv(dom == 1, 1); g(dom == 1) = Fv(dom == 1, 2);
f(dom == 2) = Fv(dom == 2, 2); g(dom == 2) = Fv(dom == 2, 1);
peak = @(x) x >= -pi/2 & x <= pi/2;
trough = @(x) (x >= -pi & x < -pi/2) | (x > pi/2 & x <= pi);
switch name
    case 'spindle'
        C = {f <= 25, f > 25 & f <= 50, f > 50 & f <= 75, f > 75};
    case 'SO_phase'
        C = {peak(f), trough(f)};
    case 'SO_trough_uni_bi'
        C = {trough(f) & trough(g), trough(f) & peak(g)};
    case 'SO_peak_uni_bi'
        C = {peak(f) & peak(g), peak(f) & trough(g)};
end
end

function r = condition_auc(sel, meas, C, sess_events, nBoot)
% Hierarchical bootstrap AUC per condition (sessions, then events within session)
nC = numel(C); nS = numel(sess_events);
th = prctile(abs(sel), 0:10:100); th = th(1:end-1);
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
        idx_perm{k} = e(randperm(s, numel(e)));   % within-session shuffle of HPC
    end
    idx = [idx{:}]; idx_perm = [idx_perm{:}];
    a = NaN(1, nC); as = a;
    for c = 1:nC
        a(c) = auc(sel(idx), meas(idx), Cm(c, idx), th);
        as(c) = auc(sel(idx), meas(idx_perm), Cm(c, idx), th);
    end
    b_real(iBoot, :) = a; b_shuf(iBoot, :) = as;
end
r.n = sum(Cm, 2)';
r.AUC_mean = mean(b_real, 'omitnan'); r.AUC_CI = prctile(b_real, [2.5 97.5]);
r.AUC_shuffled_mean = mean(b_shuf, 'omitnan'); r.AUC_shuffled_CI = prctile(b_shuf, [2.5 97.5]);
r.p_vs_shuffled = mean(b_real <= b_shuf);
d = b_real(:, end) - b_real(:, 1);                 % last - first condition
r.diff_boot = d;
r.diff_mean = mean(d, 'omitnan'); r.diff_CI = prctile(d, [2.5 97.5]);
r.diff_p = mean(d <= 0);
r.diff_p_two = min(1, 2 * min(r.diff_p, 1 - r.diff_p));
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

function plot_static(res, a, dirs, v1_windows, save_path)
nC = numel(a.labels);
fig = figure('Name', a.title, 'Position', [100 100 330 * numel(v1_windows) * 1.2 300 * numel(dirs)]);
tiledlayout(numel(dirs), numel(v1_windows), 'TileSpacing', 'compact');
for d = 1:numel(dirs)
    for w = 1:numel(v1_windows)
        ax = nexttile; hold on;
        r = res.(matlab.lang.makeValidName(dirs(d).name)).(v1_windows(w).name);
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
        title({sprintf('%s, %s', dirs(d).name, v1_windows(w).label), ...
            sprintf('%s - %s: %.3f [%.3f %.3f], p = %.3f', a.labels{end}, a.labels{1}, r.diff_mean, r.diff_CI(1), r.diff_CI(2), r.diff_p)}, ...
            'FontSize', 9, 'FontWeight', 'normal');
    end
end
linkaxes(findobj(fig, 'Type', 'axes'), 'y');
sgtitle({a.title, 'hierarchical (session -> event) bootstrap; dominance half -> selector half (grey: shuffled)'}, 'FontSize', 11);
savefig(fig, [save_path '.fig']);
exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
end

function plot_temporal(AUC, a, dirs, save_path)
nC = numel(a.labels);
fig = figure('Name', [a.title ' temporal'], 'Position', [100 100 360 * numel(dirs) 200 * (nC + 1)]);
tiledlayout(nC + 1, numel(dirs), 'TileSpacing', 'compact', 'TileIndexing', 'columnmajor');
for d = 1:numel(dirs)
    A = AUC.(matlab.lang.makeValidName(dirs(d).name));
    t = A.time_bins;
    for c = 1:nC
        nexttile; hold on;
        ci = squeeze(A.ci(:, c, :)); cs = squeeze(A.shuffled_ci(:, c, :));
        fill([t fliplr(t)], [cs(:, 1)' fliplr(cs(:, 2)')], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
        plot(t, A.shuffled_mean(:, c), 'k', 'LineWidth', 1);
        fill([t fliplr(t)], [ci(:, 1)' fliplr(ci(:, 2)')], a.colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        plot(t, A.mean(:, c), 'Color', a.colors(c, :), 'LineWidth', 2);
        xline(0, '--k'); yline(0, '--r'); xlim([-0.5 0.5]); ylim([-0.1 0.25]);
        title(sprintf('%s: %s', dirs(d).name, a.labels{c}), 'FontWeight', 'normal', 'FontSize', 9);
        if d == 1, ylabel('AUC'); end
        set(gca, 'TickDir', 'out', 'Box', 'off');
    end
    nexttile; hold on;
    fill([t fliplr(t)], [A.diff_ci(:, 1)' fliplr(A.diff_ci(:, 2)')], [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(t, A.diff_mean, 'k', 'LineWidth', 2);
    sig = A.diff_ci(:, 1)' > 0;
    if any(sig), plot(t(sig), 0.2 * ones(1, sum(sig)), 'k.', 'MarkerSize', 8); end
    xline(0, '--k'); yline(0, 'k-'); xlim([-0.5 0.5]); ylim([-0.15 0.25]);
    title(sprintf('%s - %s (dots: 95%% CI > 0)', a.labels{end}, a.labels{1}), 'FontWeight', 'normal', 'FontSize', 9);
    xlabel('V1 window start (s rel. ripple onset)');
    if d == 1, ylabel('AUC difference'); end
    set(gca, 'TickDir', 'out', 'Box', 'off');
end
sgtitle({[a.title ': temporal V1->HPC AUC, hierarchical bootstrap'], ...
    'columns: dominance half -> selector half (black: within-session shuffled)'}, 'FontSize', 11);
savefig(fig, [save_path '.fig']);
exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
end
