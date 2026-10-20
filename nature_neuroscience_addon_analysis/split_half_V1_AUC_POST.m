%% Split-half V1->HPC reactivation coherence AUC: spindle-band power and SO phase (POST)
% R1.3 control: the V1 half that defines the dominant hemisphere (whose spindle /
% SO feature sets the condition) is not the half used to select Track 1 / Track 2
% events for the coherence measure.
%   dominance half D : sign of D's V1 bias -> dominant hemisphere (D > 0 ->
%                      hemisphere 2, D < 0 -> hemisphere 1); the event's condition
%                      is that hemisphere's spindle quartile or SO phase.
%   selector half S  : Track 1 events S >= th, Track 2 events S <= -th (th =
%                      0:10:90 prctile of |S| over all events), restricted to
%                      events whose condition holds.
%   measure          : HPC bias (mean z log odds 0-0.1 s); AUC = normalised trapz
%                      of mean(HPC | T1) - mean(HPC | T2) across thresholds.
% Directions: A->B (A dominance, B selector), B->A, and Full->Full (reference,
% same decoder code as the halves; as the original analysis). Static AUC also
% saves 'Avg': 2-fold cross-validated mean of A->B and B->A, averaged per
% bootstrap resample (resamples are seeded, so both folds use the same events);
% CI, p vs shuffled and the High - Low p come from the averaged distribution.
% Bootstrap over events (1000, seeded; same resamples for all conditions, so
% condition differences are paired); shuffled control keeps selector and
% conditions while resampling HPC (as reactivation_coherence_auc).
% Spindle: within-session percentile of spindle amplitude
% (ripple_info.spindle_percentile, 10-100: bottom 10% NaN, excluded as in the
% original paper). Quartiles <= 25, 25-50, 50-75, > 75. SO: peak [-pi/2, pi/2]
% vs trough.
% Static AUC: V1 window -0.2-0 s (pre) and 0-0.2 s (post). Temporal AUC: 0.1 s V1
% window starting at each time bin (-0.5:0.02:0.5 s), as the original V1->HPC
% temporal analysis.
% Unilateral vs bilateral SO (as corticohippocampal_reactivation_temporal_KDE_diff_
% V1_detection.m): trough unilateral = dominant trough & non-dominant peak, trough
% bilateral = both trough; peak unilateral = dominant peak & non-dominant trough,
% peak bilateral = both peak. Difference = unilateral - bilateral.
% Ripple power (HPC->V1, as reactivation_AUC_coherence_PRE): within-session
% ripple power quartiles; HPC bias selects Track 1 / 2 (thresholds from |HPC| of
% each quartile's events) and the selector half's V1 bias is measured (Full,
% B for A->B, A for B->A; no dominance half). Temporal: V1 window centred on
% each bin, as the original HPC->V1 temporal analysis.
% Input: KDE_reactivation_ripples_PSTH_POST_unique_V1_split_half.mat
% Output: nature_neuroscience_revision_analysis\split_half_V1_AUC_POST
%   which       : analyses to run (default / empty all: 'spindle', 'SO_phase',
%                 'SO_trough_uni_bi', 'SO_peak_uni_bi', 'ripple_power')
%   do_temporal : run the temporal AUC (default true)
function split_half_V1_AUC_POST(which, do_temporal)
    if nargin < 1, which = []; end
    if nargin < 2, do_temporal = true; end
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    out_dir = fullfile(rev_dir, 'split_half_V1_AUC_POST');
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

    dirs = struct('name', {'Full->Full', 'A->B', 'B->A'}, 'dom', {'full', 'A', 'B'}, 'sel', {'full', 'B', 'A'});
    an(1).name = 'spindle'; an(1).title = 'Spindle-band power (dominant hemisphere, within-session percentile)';
    an(1).labels = {'Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)'}; an(1).colors = cool(4);
    an(2).name = 'SO_phase'; an(2).title = 'SO phase (dominant hemisphere)';
    an(2).labels = {'SO Peak', 'SO Trough'}; an(2).colors = [0.2 0.6 0.8; 0.8 0.2 0.2];
    an(3).name = 'SO_trough_uni_bi'; an(3).title = 'SO trough (dominant): bilateral vs unilateral';
    an(3).labels = {'Trough bilateral', 'Trough unilateral'}; an(3).colors = [0.95 0.6 0.6; 0.8 0.2 0.2];
    an(4).name = 'SO_peak_uni_bi'; an(4).title = 'SO peak (dominant): bilateral vs unilateral';
    an(4).labels = {'Peak bilateral', 'Peak unilateral'}; an(4).colors = [0.6 0.8 0.95; 0.2 0.6 0.8];
    an(5).name = 'ripple_power'; an(5).title = 'Ripple power (within-session quartiles)';
    an(5).labels = {'Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)'}; an(5).colors = cool(4);
    [an.direction] = deal('V1->HPC'); an(5).direction = 'HPC->V1';
    if ~isempty(which), an = an(ismember({an.name}, which)); end
    F.spindle = ripple_info.spindle_percentile;
    F.SO_phase = ripple_info.SO_phase;
    F.SO_trough_uni_bi = ripple_info.SO_phase;
    F.SO_peak_uni_bi = ripple_info.SO_phase;
    rp_q = within_session_quartile(ripple_info.ripple_power(:), ripple_info.session_id(:))';

    %% Static AUC
    v1_windows = struct('name', {'pre', 'post'}, 'label', {'V1 -0.2-0 s (pre)', 'V1 0-0.2 s (post)'}, ...
        'bins', {bc > -0.2 & bc < 0, bc > 0 & bc < 0.2});
    T = table();
    dirs_static = [dirs, struct('name', 'Avg', 'dom', '', 'sel', '')];   % Avg = mean of A->B and B->A
    for iA = 1:numel(an)
        res = struct();
        for d = 1:numel(dirs_static)
            for w = 1:numel(v1_windows)
                if strcmp(dirs_static(d).name, 'Avg')
                    r = average_folds(res.(matlab.lang.makeValidName('A->B')).(v1_windows(w).name), ...
                        res.(matlab.lang.makeValidName('B->A')).(v1_windows(w).name));
                else
                    dom_b = mean(V1.(dirs(d).dom)(v1_windows(w).bins, :), 'omitnan');
                    sel_b = mean(V1.(dirs(d).sel)(v1_windows(w).bins, :), 'omitnan');
                    r = run_analysis(an(iA), F, rp_q, dom_b, sel_b, HPC_b, nBoot);
                end
                res.(matlab.lang.makeValidName(dirs_static(d).name)).(v1_windows(w).name) = r;
                for c = 1:numel(r.n)
                    T = [T; table(string(an(iA).name), string(dirs_static(d).name), string(v1_windows(w).name), string(an(iA).labels{c}), ...
                        r.n(c), r.AUC_mean(c), r.AUC_CI(1, c), r.AUC_CI(2, c), r.AUC_shuffled_mean(c), r.AUC_shuffled_CI(1, c), ...
                        r.AUC_shuffled_CI(2, c), r.p_vs_shuffled(c), ...
                        'VariableNames', {'analysis', 'direction', 'V1_window', 'condition', 'n_events', 'AUC_mean', 'AUC_CI_lo', ...
                        'AUC_CI_hi', 'AUC_shuffled_mean', 'AUC_shuffled_CI_lo', 'AUC_shuffled_CI_hi', 'p_vs_shuffled'})]; %#ok<AGROW>
                end
                fprintf('%s %s %s: %s | %s - %s = %.3f [%.3f %.3f], p = %.3f\n', an(iA).name, dirs_static(d).name, v1_windows(w).name, ...
                    sprintf('%.3f ', r.AUC_mean), an(iA).labels{end}, an(iA).labels{1}, r.diff_mean, r.diff_CI(1), r.diff_CI(2), r.diff_p);
            end
        end
        S = struct('res', res); S.v1_windows = v1_windows; S.dirs = dirs_static;   % saved as 'dirs' incl. Avg
        save(fullfile(out_dir, ['split_half_AUC_static_' an(iA).name '.mat']), '-struct', 'S');
        plot_static(res, an(iA), dirs_static, v1_windows, fullfile(out_dir, ['split_half_AUC_static_' an(iA).name]));
    end
    for iA = 1:numel(an)
        writetable(T(T.analysis == an(iA).name, :), fullfile(out_dir, ['split_half_AUC_static_' an(iA).name '.csv']));
    end
    if isempty(which) % combined table only when all analyses are run
        writetable(T, fullfile(out_dir, 'split_half_AUC_static.csv'));
    end
    plot_AUC_low_vs_high_controls('split_half');   % low vs high bar plot of all four analyses
    if ~do_temporal, fprintf('Saved to %s\n', out_dir); return; end

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
            if strcmp(an(iA).direction, 'HPC->V1')
                t_off = -win / 2;   % centred, as original HPC->V1
            else
                t_off = 0;          % starting, as original V1->HPC
            end
            for t = 1:numel(time_bins)
                in = bc >= time_bins(t) + t_off & bc < time_bins(t) + t_off + win;
                dom_b = mean(V1.(dirs(d).dom)(in, :), 'omitnan');
                sel_b = mean(V1.(dirs(d).sel)(in, :), 'omitnan');
                r = run_analysis(an(iA), F, rp_q, dom_b, sel_b, HPC_b, nBoot);
                A.mean(t, :) = r.AUC_mean; A.ci(t, :, :) = permute(r.AUC_CI, [3 2 1]);
                A.shuffled_mean(t, :) = r.AUC_shuffled_mean; A.shuffled_ci(t, :, :) = permute(r.AUC_shuffled_CI, [3 2 1]);
                A.diff_mean(t) = r.diff_mean; A.diff_ci(t, :) = r.diff_CI; A.diff_p(t) = r.diff_p;
            end
            A.labels = an(iA).labels; A.time_bins = time_bins; A.win_size = win;
            if t_off < 0, A.window_alignment = 'centred'; else, A.window_alignment = 'start'; end
            AUC.(matlab.lang.makeValidName(dirs(d).name)) = A;
        end
        save(fullfile(out_dir, ['split_half_AUC_temporal_' an(iA).name '.mat']), 'AUC', 'dirs');
        plot_temporal(AUC, an(iA), dirs, fullfile(out_dir, ['split_half_AUC_temporal_' an(iA).name]));
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

function r = run_analysis(a, F, rp_q, dom_b, sel_b, HPC_b, nBoot)
% V1->HPC: selector-half V1 selects, HPC measured, conditions from the dominance
% half. HPC->V1 (ripple power): HPC selects, selector-half V1 measured, per-
% quartile thresholds.
if strcmp(a.direction, 'HPC->V1')
    C = arrayfun(@(q) rp_q == q, 1:4, 'UniformOutput', false);
    r = condition_auc(HPC_b, sel_b, C, nBoot, true);
else
    C = condition_masks(a.name, F.(a.name), dom_b);
    r = condition_auc(sel_b, HPC_b, C, nBoot, false);
end
end

function q = within_session_quartile(metric, session_ids)
% Quartile (1-4) of the metric within each session (as reactivation_AUC_coherence_PRE)
q = NaN(size(metric));
for s = unique(session_ids)'
    idx = find(session_ids == s);
    p = prctile(metric(idx), 0:25:100);
    for b = 1:4
        if b == 1
            m = metric(idx) >= p(b) & metric(idx) <= p(b+1);
        else
            m = metric(idx) > p(b) & metric(idx) <= p(b+1);
        end
        q(idx(m)) = b;
    end
end
end

function r = condition_auc(sel, meas, C, nBoot, per_cond_th)
% Bootstrapped AUC per condition (same resamples for all conditions).
% Thresholds: 0:10:90 prctile of |sel| over all events, or over each
% condition's events if per_cond_th (as the original HPC->V1 quartiles)
nC = numel(C); n_ev = numel(sel);
Cm = cat(1, C{:});
th = NaN(nC, 10);
for c = 1:nC
    if per_cond_th, x = abs(sel(Cm(c, :))); else, x = abs(sel); end
    t = prctile(x, 0:10:100); th(c, :) = t(1:end-1);
end
b_real = NaN(nBoot, nC); b_shuf = NaN(nBoot, nC);
parfor iBoot = 1:nBoot
    s = RandStream('philox4x32_10', 'Seed', iBoot);
    idx = randi(s, n_ev, n_ev, 1)';
    a = NaN(1, nC); as = a;
    for c = 1:nC
        a(c) = auc(sel(idx), meas(idx), Cm(c, idx), th(c, :));
        as(c) = auc(sel, meas(idx), Cm(c, :), th(c, :));
    end
    b_real(iBoot, :) = a; b_shuf(iBoot, :) = as;
end
r = boot_summary(b_real, b_shuf);
r.n = sum(Cm, 2)';
end

function r = average_folds(r1, r2)
% 2-fold cross-validated AUC: mean of the two folds (A->B, B->A) per bootstrap
% resample (same seeded resamples), then summarised as one distribution
r = boot_summary((r1.b_real + r2.b_real) / 2, (r1.b_shuf + r2.b_shuf) / 2);
r.n = NaN(size(r1.n));   % event sets differ between folds
end

function r = boot_summary(b_real, b_shuf)
r.b_real = b_real; r.b_shuf = b_shuf;
r.AUC_mean = mean(b_real, 'omitnan'); r.AUC_CI = prctile(b_real, [2.5 97.5]);
r.AUC_shuffled_mean = mean(b_shuf, 'omitnan'); r.AUC_shuffled_CI = prctile(b_shuf, [2.5 97.5]);
r.p_vs_shuffled = mean(b_real <= b_shuf);
d = b_real(:, end) - b_real(:, 1);                 % High - Low / Trough - Peak
r.diff_mean = mean(d, 'omitnan'); r.diff_CI = prctile(d, [2.5 97.5]); r.diff_p = mean(d <= 0);
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
sgtitle({a.title, 'dominance half -> selector half (grey: shuffled)'}, 'FontSize', 11);
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
    xlabel(sprintf('V1 window %s (s rel. ripple onset)', A.window_alignment));
    if d == 1, ylabel('AUC difference'); end
    set(gca, 'TickDir', 'out', 'Box', 'off');
end
sgtitle({[a.title ': temporal ' a.direction ' AUC'],'columns: dominance half -> selector half (black: shuffled)'}, 'FontSize', 11);
savefig(fig, [save_path '.fig']);
exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
end
