%% Per-session or per-animal bootstrapped reactivation coherence AUC, original POST decoder
% Original POST PLS-KDE V1 and HPC z log odds (KDE_reactivation_ripples_PSTH_POST_unique.mat,
% non-duplicated events). For each unit (session or animal) the unit's events are
% bootstrapped (nBoot = 1000, event resampling within the unit, seeded as
% reactivation_coherence_auc) giving the mean AUC and its 95% CI, and the
% shuffled AUC (selector and its conditions un-resampled, measure resampled)
% with its 95% CI. Same AUC as reactivation_coherence_auc: thresholds 0:10:90
% prctile of |selector| of the unit's condition events.
% Conditions are defined within session (also when pooled per animal):
%   all_ripples  (HPC->V1): all events
%   ripple_power (HPC->V1): within-session ripple power Q1 (Low) vs Q4 (High)
%   spindle      (V1->HPC): hemisphere-matched spindle_percentile Q1 vs Q4
%                           (bottom 10% NaN, Q1 = 10-25th, as original paper)
%   SO_phase     (V1->HPC): hemisphere-matched SO peak vs trough
%   SO_trough_uni_bi (V1->HPC): dominant hemisphere in SO trough, non-dominant in
%                           peak (unilateral) vs also in trough (bilateral)
% V1->HPC conditions use all events of the unit with the same resamples (paired
% Trough - Peak / High - Low / Bilateral - Unilateral difference); ripple power
% Q1 and Q4 are resampled separately.
% Static: V1 -0.2-0 s (pre) and 0-0.2 s (post). Temporal: 0.1 s V1 window,
% -0.5:0.02:0.5 s, centred on each bin for HPC->V1 and starting at each bin for
% V1->HPC (as hierarchical_AUC_POST_original). HPC selector/measure 0-0.1 s.
% p_vs_shuffle: fraction of resamples with AUC <= shuffled AUC (one-sided).
% Outputs: nature_neuroscience_revision_analysis\POST_original_bootAUC_per_<unit>
%   bootAUC_<analysis>.mat         static and temporal results (unit x ...)
%   bootAUC_static.csv, bootAUC_temporal.csv
%   summary_temporal               mean bootstrapped AUC trace of every unit
%                                  (Low vs High / Peak vs Trough, no shuffle)
%   summary_temporal_diff          each unit's bootstrapped mean difference trace
%                                  (High - Low, Trough - Peak, Unilateral - Bilateral)
%   static\<unit>_<id>_static      per-unit bars (real vs shuffled, 95% CI), optional
%   temporal\<unit>_<id>_temporal  per-unit traces (real vs shuffled, 95% CI), optional
%   unit         : 'session' (default) or 'animal'
%   compute      : analyses to (re)compute (default or [] all; {} = none). Tables and
%                  figures use every analysis with a saved bootAUC_*.mat.
%   unit_figures : true = also draw the per-unit figures (default false)
%   decoder      : 'original' (default) or 'rate_lt10' (PLS-KDE decoders from neurons
%                  < 10 Hz, KDE_reactivation_ripples_PSTH_POST_unique_rate_lt10.mat;
%                  output POST_rate_lt10_bootAUC_per_<unit>)
function boot_AUC_per_unit_POST_original(unit, compute, unit_figures, decoder)
    all_analyses = {'all_ripples', 'ripple_power', 'spindle', 'SO_phase', 'SO_trough_uni_bi'};
    if nargin < 1 || isempty(unit), unit = 'session'; end
    if nargin < 2 || (isnumeric(compute) && isempty(compute)), compute = all_analyses; end   % {} = none
    if nargin < 3 || isempty(unit_figures), unit_figures = false; end
    if nargin < 4 || isempty(decoder), decoder = 'original'; end
    assert(ismember(decoder, {'original', 'rate_lt10'}), 'decoder must be ''original'' or ''rate_lt10''');
    rev_dir = addon_setup();
    out_dir = fullfile(rev_dir, ['POST_' decoder '_bootAUC_per_' unit]);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    if ~isempty(compute)
        compute_boot_AUC(rev_dir, out_dir, unit, compute, decoder);
    end
    saved = all_analyses(cellfun(@(a) exist(fullfile(out_dir, ['bootAUC_' a '.mat']), 'file') > 0, all_analyses));
    for iA = 1:numel(saved)
        RES(iA) = load(fullfile(out_dir, ['bootAUC_' saved{iA} '.mat']), 'static', 'temporal', 'info'); %#ok<AGROW>
    end

    %% Tables
    TS = table(); TT = table();
    for iA = 1:numel(RES)
        [ts, tt] = result_tables(RES(iA));
        TS = [TS; ts]; TT = [TT; tt]; %#ok<AGROW>
    end
    writetable(TS, fullfile(out_dir, 'bootAUC_static.csv'));
    writetable(TT, fullfile(out_dir, 'bootAUC_temporal.csv'));

    %% Figures: summary of temporal traces; per unit (static, temporal) if asked
    plot_summary_temporal(RES, unit, fullfile(out_dir, 'summary_temporal'));
    plot_summary_temporal_diff(RES, unit, fullfile(out_dir, 'summary_temporal_diff'));
    if unit_figures
        for d = {fullfile(out_dir, 'static'), fullfile(out_dir, 'temporal')}
            if ~exist(d{1}, 'dir'), mkdir(d{1}); end
        end
        in = RES(1).info;
        for k = 1:numel(in.unit_ids)
            name = sprintf('%s_%d', unit, in.unit_ids(k));
            if strcmp(unit, 'session'), name_t = sprintf('Session %d (animal %d)', in.unit_ids(k), in.animal_ids(k));
            else, name_t = sprintf('Animal %d (%d sessions)', in.unit_ids(k), in.n_sessions(k)); end
            plot_unit_static(RES, k, name_t, fullfile(out_dir, 'static', [name '_static']));
            plot_unit_temporal(RES, k, name_t, fullfile(out_dir, 'temporal', [name '_temporal']));
        end
    end
    fprintf('Saved to %s\n', out_dir);
end

function compute_boot_AUC(rev_dir, out_dir, unit, which, decoder)
% Bootstrapped AUC per unit, analysis and V1 window; saves bootAUC_<analysis>.mat
    nBoot = 1000;
    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    if strcmp(decoder, 'original')
        psth_file = 'KDE_reactivation_ripples_PSTH_POST_unique.mat';
    else
        psth_file = 'KDE_reactivation_ripples_PSTH_POST_unique_rate_lt10.mat';
    end
    load(fullfile(rev_dir, psth_file), 'KDE_reactivation_ripples_PSTH');
    P = KDE_reactivation_ripples_PSTH; clear KDE_reactivation_ripples_PSTH
    HC = (P.HPC_z_logodds_ripples + P.nan_mask)';
    V1 = (P.V1_z_logodds_ripples + P.nan_mask)';
    clear P

    switch unit
        case 'session', uid = ripple_info.session_id(:)';
        case 'animal',  uid = ripple_info.subject_id(:)';
        otherwise, error('unit must be ''session'' or ''animal''');
    end
    units = unique(uid); nU = numel(units);
    sess_animal = arrayfun(@(u) ripple_info.subject_id(find(uid == u, 1)), units);
    n_sessions = arrayfun(@(u) numel(unique(ripple_info.session_id(uid == u))), units);

    % V1 per window: static (pre, post) and temporal (centred / starting)
    bc = (-1:0.01:1); bc = bc(1:end-1) + 0.005;
    HPC_b = mean(HC(bc > 0 & bc < 0.1, :), 'omitnan');
    v1_windows = struct('name', {'pre', 'post'}, 'label', {'V1 -0.2-0 s (pre)', 'V1 0-0.2 s (post)'}, ...
        'range', {[-0.2 0], [0 0.2]});
    V_static = cell2mat(arrayfun(@(w) mean(V1(bc > w.range(1) & bc < w.range(2), :), 'omitnan'), v1_windows', 'UniformOutput', false));
    time_bins = -0.5:0.02:0.5; win = 0.1;
    V_centred = cell2mat(arrayfun(@(t) mean(V1(bc >= t - win/2 & bc < t + win/2, :), 'omitnan'), time_bins', 'UniformOutput', false));
    V_start = cell2mat(arrayfun(@(t) mean(V1(bc >= t & bc < t + win, :), 'omitnan'), time_bins', 'UniformOutput', false));
    clear V1 HC
    nW = numel(v1_windows); nT = numel(time_bins);

    an = analysis_list(ripple_info);
    an = an(ismember({an.name}, which));
    fprintf('%s: %d units, %d events, %d analyses\n', unit, nU, numel(uid), numel(an));

    %% Bootstrapped AUC per unit, analysis and V1 window (static rows, then temporal rows)
    for iA = 1:numel(an)
        a = an(iA); nC = numel(a.labels);
        is_hpc = strcmp(a.direction, 'HPC->V1');
        if is_hpc, V = [V_static; V_centred]; else, V = [V_static; V_start]; end
        nR = size(V, 1);
        R.mean = NaN(nU, nR, nC); R.ci = NaN(nU, nR, nC, 2);
        R.shuffled_mean = R.mean; R.shuffled_ci = R.ci; R.p_vs_shuffled = R.mean;
        R.diff_mean = NaN(nU, nR); R.diff_ci = NaN(nU, nR, 2); R.diff_p = R.diff_mean;
        R.n_events = NaN(nU, nW, nC);
        tic
        for k = 1:nU
            in_u = find(uid == units(k));
            ev = cell(1, nC); c1 = ev; c2 = ev; I = ev;
            for c = 1:nC
                ev{c} = find(a.mask{c}(in_u));
                c1{c} = a.c1{c}(in_u(ev{c})); c2{c} = a.c2{c}(in_u(ev{c}));
                I{c} = boot_indices(numel(ev{c}), nBoot);
            end
            [br, bs] = unit_boot(V(:, in_u), HPC_b(in_u), is_hpc, ev, c1, c2, I);
            for c = 1:nC
                for w = 1:nW   % events in the condition (V1->HPC: dominant hemisphere by V1 sign)
                    if is_hpc
                        R.n_events(k, w, c) = numel(ev{c});
                    else
                        v = V(w, in_u(ev{c}));
                        R.n_events(k, w, c) = sum((v > 0 & c1{c}) | (v < 0 & c2{c}));
                    end
                end
                x = br(:, :, c); xs = bs(:, :, c);
                R.mean(k, :, c) = mean(x, 2, 'omitnan');
                R.ci(k, :, c, :) = prctile(x, [2.5 97.5], 2);
                R.shuffled_mean(k, :, c) = mean(xs, 2, 'omitnan');
                R.shuffled_ci(k, :, c, :) = prctile(xs, [2.5 97.5], 2);
                p = mean(x <= xs, 2); p(all(isnan(x), 2)) = NaN;
                R.p_vs_shuffled(k, :, c) = p;
            end
            if nC > 1
                d = br(:, :, end) - br(:, :, 1);
                R.diff_mean(k, :) = mean(d, 2, 'omitnan');
                R.diff_ci(k, :, :) = prctile(d, [2.5 97.5], 2);
                p = mean(d <= 0, 2); p(all(isnan(d), 2)) = NaN;
                R.diff_p(k, :) = p;
            end
            fprintf('  %s %s %d (%d/%d) done\n', a.name, unit, units(k), k, nU);
        end
        fprintf('%s: %.0f s\n', a.name, toc);

        static = split_rows(R, 1:nW); static.n_events = R.n_events;
        static.windows = v1_windows;
        temporal = split_rows(R, nW + (1:nT));
        temporal.time_bins = time_bins; temporal.win_size = win;
        temporal.alignment = ternary(is_hpc, 'centred', 'starting');
        info = struct('unit', unit, 'unit_ids', units, 'animal_ids', sess_animal, 'n_sessions', n_sessions, ...
            'analysis', a.name, 'title', a.title, 'direction', a.direction, 'labels', {a.labels}, ...
            'colors', a.colors, 'nBoot', nBoot, 'decoder', decoder);
        save(fullfile(out_dir, ['bootAUC_' a.name '.mat']), 'static', 'temporal', 'info');
    end
end

function an = analysis_list(info)
% Condition masks (events entering the AUC) and Track 1 / Track 2 conditions
nEv = numel(info.session_id); all_true = true(1, nEv);
S = AUC_low_high_style();
q = within_session_quartile(info.ripple_power(:), info.session_id(:))';
Hs = hemisphere_conditions('spindle', info);
Hso = hemisphere_conditions('SO', info);
Hub = hemisphere_conditions('SO_trough_bilateral', info);   % 1 unilateral, 2 bilateral
row = @(x) x(:)';
an = struct('name', {}, 'title', {}, 'direction', {}, 'labels', {}, 'colors', {}, 'mask', {}, 'c1', {}, 'c2', {});
an(1).name = 'all_ripples'; an(1).title = 'All ripples (HPC->V1)'; an(1).direction = 'HPC->V1';
an(1).labels = {'All'}; an(1).colors = [231 41 138] / 256;
an(1).mask = {all_true}; an(1).c1 = {all_true}; an(1).c2 = {all_true};
an(2).name = 'ripple_power'; an(2).title = 'Ripple power Q1 vs Q4 (HPC->V1)'; an(2).direction = 'HPC->V1';
an(2).labels = {'Low (Q1)', 'High (Q4)'}; an(2).colors = S.ripple_power.colors;
an(2).mask = {q == 1, q == 4}; an(2).c1 = {all_true, all_true}; an(2).c2 = an(2).c1;
an(3).name = 'spindle'; an(3).title = 'Spindle Q1 vs Q4 (V1->HPC)'; an(3).direction = 'V1->HPC';
an(3).labels = {'Low (Q1)', 'High (Q4)'}; an(3).colors = S.spindle.colors;
an(3).mask = {all_true, all_true}; an(3).c1 = {row(Hs(1).c1), row(Hs(4).c1)}; an(3).c2 = {row(Hs(1).c2), row(Hs(4).c2)};
an(4).name = 'SO_phase'; an(4).title = 'SO peak vs trough (V1->HPC)'; an(4).direction = 'V1->HPC';
an(4).labels = {'Peak', 'Trough'}; an(4).colors = S.SO_phase.colors;
an(4).mask = {all_true, all_true}; an(4).c1 = {row(Hso(1).c1), row(Hso(2).c1)}; an(4).c2 = {row(Hso(1).c2), row(Hso(2).c2)};
an(5).name = 'SO_trough_uni_bi'; an(5).title = 'Dominant SO trough: unilateral vs bilateral (V1->HPC)'; an(5).direction = 'V1->HPC';
an(5).labels = {'Unilateral trough', 'Bilateral trough'}; an(5).colors = S.SO_trough_uni_bi.colors;
an(5).mask = {all_true, all_true}; an(5).c1 = {row(Hub(1).c1), row(Hub(2).c1)}; an(5).c2 = {row(Hub(1).c2), row(Hub(2).c2)};
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

function I = boot_indices(n, nBoot)
% nBoot x n event resamples, seeded per iteration as reactivation_coherence_auc
I = zeros(nBoot, n, 'int32');
if n == 0, return; end
for b = 1:nBoot
    s = RandStream('philox4x32_10', 'Seed', b);
    I(b, :) = randi(s, n, n, 1)';
end
end

function [br, bs] = unit_boot(V, hpc, is_hpc, ev, c1, c2, I)
% Bootstrapped real / shuffled AUC (nRows x nBoot x nC) for one unit, one row per V1 window
nR = size(V, 1); nC = numel(ev); nBoot = size(I{1}, 1);
br = NaN(nR, nBoot, nC); bs = br;
parfor r = 1:nR
    v = V(r, :);
    a = NaN(nBoot, nC); as = a;
    for c = 1:nC
        if is_hpc
            [a(:, c), as(:, c)] = boot_auc(hpc(ev{c}), v(ev{c}), c1{c}, c2{c}, I{c});
        else
            [a(:, c), as(:, c)] = boot_auc(v(ev{c}), hpc(ev{c}), c1{c}, c2{c}, I{c});
        end
    end
    br(r, :, :) = reshape(a, [1 nBoot nC]);
    bs(r, :, :) = reshape(as, [1 nBoot nC]);
end
end

function [a, as] = boot_auc(sel, meas, c1, c2, I)
% Vectorised reactivation_coherence_auc: real (sel, meas, conditions resampled
% together) and shuffled (only meas resampled) AUC per resample
nBoot = size(I, 1); n = numel(sel);
a = NaN(nBoot, 1); as = a;
if n < 10, return; end
th = prctile(abs(sel), 0:10:100); th = th(1:end-1);
S = sel(I); M = meas(I);
ok = ~isnan(M); M(~ok) = 0;
d = track_diff(S, M, ok & c1(I), ok & c2(I), th);
ds = track_diff(repmat(sel, nBoot, 1), M, ok & c1, ok & c2, th);
a = trapz(th, d, 2) / (max(th) - min(th));
as = trapz(th, ds, 2) / (max(th) - min(th));
end

function d = track_diff(S, M, w1, w2, th)
% mean(M | S >= th(j) & w1) - mean(M | S <= -th(j) & w2) for every threshold and
% resample (rows); events are binned by the number of thresholds they pass, then
% reverse cumulative sums over bins give the nested threshold sets
[nBoot, n] = size(M); nT = numel(th);
jp = zeros(nBoot, n); jn = jp;
for j = 1:nT
    jp = jp + (S >= th(j)); jn = jn + (S <= -th(j));
end
row = repmat((1:nBoot)', 1, n);
pass = @(jj, x) flip(cumsum(flip(reshape(accumarray(row(:) + nBoot * jj(:), x(:), [nBoot * (nT + 1) 1]), ...
    nBoot, nT + 1), 2), 2), 2);
s1 = pass(jp, M .* w1); n1 = pass(jp, double(w1 & true(nBoot, 1)));
s2 = pass(jn, M .* w2); n2 = pass(jn, double(w2 & true(nBoot, 1)));
d = s1(:, 2:end) ./ n1(:, 2:end) - s2(:, 2:end) ./ n2(:, 2:end);
end

function X = split_rows(R, rows)
X.mean = R.mean(:, rows, :); X.ci = R.ci(:, rows, :, :);
X.shuffled_mean = R.shuffled_mean(:, rows, :); X.shuffled_ci = R.shuffled_ci(:, rows, :, :);
X.p_vs_shuffled = R.p_vs_shuffled(:, rows, :);
X.diff_mean = R.diff_mean(:, rows); X.diff_ci = R.diff_ci(:, rows, :); X.diff_p = R.diff_p(:, rows);
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function s = decoder_label(in)
% Figure label of the decoder (results saved before the decoder option: original)
if isfield(in, 'decoder') && strcmp(in.decoder, 'rate_lt10')
    s = 'POST PLS-KDE decoder, neurons < 10 Hz';
else
    s = 'original POST decoder';
end
end

function [ts, tt] = result_tables(X)
% Long tables: unit x V1 window x condition (static), unit x time x condition (temporal)
in = X.info; nU = numel(in.unit_ids); nC = numel(in.labels);
if strcmp(in.unit, 'session'), second = in.animal_ids; second_name = 'animal_id';
else, second = in.n_sessions; second_name = 'n_sessions'; end
ts = table(); tt = table();
for k = 1:nU
    for c = 1:nC
        for w = 1:numel(X.static.windows)
            ts = [ts; table(in.unit_ids(k), second(k), string(in.analysis), string(in.direction), ...
                string(X.static.windows(w).name), string(in.labels{c}), X.static.n_events(k, w, c), ...
                X.static.mean(k, w, c), X.static.ci(k, w, c, 1), X.static.ci(k, w, c, 2), ...
                X.static.shuffled_mean(k, w, c), X.static.shuffled_ci(k, w, c, 1), X.static.shuffled_ci(k, w, c, 2), ...
                X.static.p_vs_shuffled(k, w, c), X.static.ci(k, w, c, 1) > X.static.shuffled_ci(k, w, c, 2), ...
                X.static.diff_mean(k, w), X.static.diff_ci(k, w, 1), X.static.diff_ci(k, w, 2), X.static.diff_p(k, w), ...
                'VariableNames', {[in.unit '_id'], second_name, 'analysis', 'direction', 'V1_window', 'condition', 'n_events', ...
                'AUC_boot_mean', 'AUC_CI_lo', 'AUC_CI_hi', 'shuffle_mean', 'shuffle_CI_lo', 'shuffle_CI_hi', ...
                'p_vs_shuffle', 'CI_lo_above_shuffle_CI_hi', 'diff_last_minus_first', 'diff_CI_lo', 'diff_CI_hi', ...
                'diff_p_one_sided'})]; %#ok<AGROW>
        end
        t = X.temporal.time_bins(:); nT = numel(t);
        tt = [tt; table(repmat(in.unit_ids(k), nT, 1), repmat(second(k), nT, 1), repmat(string(in.analysis), nT, 1), ...
            repmat(string(in.labels{c}), nT, 1), t, X.temporal.mean(k, :, c)', X.temporal.ci(k, :, c, 1)', ...
            X.temporal.ci(k, :, c, 2)', X.temporal.shuffled_mean(k, :, c)', X.temporal.shuffled_ci(k, :, c, 1)', ...
            X.temporal.shuffled_ci(k, :, c, 2)', X.temporal.p_vs_shuffled(k, :, c)', ...
            'VariableNames', {[in.unit '_id'], second_name, 'analysis', 'condition', 'time_bin', 'AUC_boot_mean', ...
            'AUC_CI_lo', 'AUC_CI_hi', 'shuffle_mean', 'shuffle_CI_lo', 'shuffle_CI_hi', 'p_vs_shuffle'})]; %#ok<AGROW>
    end
end
end

function save_fig(fig, save_path)
for f = {'.fig', '.png', '.pdf'}
    try   % a file open in a viewer is locked: warn instead of stopping
        switch f{1}
            case '.fig', savefig(fig, [save_path f{1}]);
            case '.png', exportgraphics(fig, [save_path f{1}], 'Resolution', 300);
            case '.pdf', exportgraphics(fig, [save_path f{1}], 'ContentType', 'vector');
        end
    catch err
        warning('Could not write %s%s (open elsewhere?): %s', save_path, f{1}, err.message);
    end
end
close(fig);
end

function plot_unit_static(RES, k, name_t, save_path)
% Rows = analyses, columns = V1 windows; shuffled (grey) and real (colour) bars with 95% CI
nA = numel(RES); W = RES(1).static.windows; nW = numel(W);
fig = figure('Visible', 'off', 'Position', [100 100 260 * nW + 120 230 * nA]);
tiledlayout(nA, nW, 'TileSpacing', 'compact');
for iA = 1:nA
    X = RES(iA).static; in = RES(iA).info; nC = numel(in.labels);
    for w = 1:nW
        ax = nexttile; hold on;
        for c = 1:nC
            m = X.shuffled_mean(k, w, c); ci = squeeze(X.shuffled_ci(k, w, c, :));
            bar(c - 0.2, m, 0.35, 'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
            errorbar(c - 0.2, m, m - ci(1), ci(2) - m, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'CapSize', 5);
            m = X.mean(k, w, c); ci = squeeze(X.ci(k, w, c, :));
            bar(c + 0.2, m, 0.35, 'FaceColor', in.colors(c, :), 'FaceAlpha', 0.8, 'EdgeColor', 'none');
            errorbar(c + 0.2, m, m - ci(1), ci(2) - m, 'k', 'LineWidth', 1.5, 'CapSize', 6);
        end
        yline(0, 'k-');
        lab = cellfun(@(l, n) sprintf('%s\\newline n = %d', l, n), in.labels, num2cell(squeeze(X.n_events(k, w, :))'), 'UniformOutput', false);
        set(ax, 'XTick', 1:nC, 'XTickLabel', lab, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 8);
        xlim([0.4 nC + 0.6]); ylabel('AUC');
        if nC > 1
            sub = sprintf('%s - %s: %.3f [%.3f %.3f]', in.labels{end}, in.labels{1}, X.diff_mean(k, w), X.diff_ci(k, w, 1), X.diff_ci(k, w, 2));
        else
            sub = sprintf('p vs shuffle = %.3f', X.p_vs_shuffled(k, w, 1));
        end
        title({sprintf('%s, %s', in.title, W(w).label), sub}, 'FontSize', 8, 'FontWeight', 'normal');
    end
end
linkaxes(findobj(fig, 'Type', 'axes'), 'y');
sgtitle({name_t, [decoder_label(RES(1).info) ', bootstrapped AUC (95% CI; grey: shuffled)']}, 'FontSize', 10);
save_fig(fig, save_path);
end

function plot_unit_temporal(RES, k, name_t, save_path)
% Rows = analyses, columns = conditions; shuffled (grey) and real (colour) traces with 95% CI
nA = numel(RES); nCol = max(arrayfun(@(X) numel(X.info.labels), RES));
fig = figure('Visible', 'off', 'Position', [100 100 330 * nCol + 60 200 * nA + 60]);
tiledlayout(nA, nCol, 'TileSpacing', 'compact');
for iA = 1:nA
    X = RES(iA).temporal; in = RES(iA).info; t = X.time_bins;
    for c = 1:numel(in.labels)
        nexttile((iA - 1) * nCol + c); hold on;
        ok = isfinite(X.mean(k, :, c));
        tk = t(ok);
        cs = squeeze(X.shuffled_ci(k, ok, c, :)); ci = squeeze(X.ci(k, ok, c, :));
        if any(ok)
            fill([tk fliplr(tk)], [cs(:, 1)' fliplr(cs(:, 2)')], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
            plot(tk, X.shuffled_mean(k, ok, c), 'k', 'LineWidth', 1);
            fill([tk fliplr(tk)], [ci(:, 1)' fliplr(ci(:, 2)')], in.colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.35);
            plot(tk, X.mean(k, ok, c), 'Color', in.colors(c, :) * 0.75, 'LineWidth', 2);   % darker: Low colours are pale
        end
        xline(0, '--k'); yline(0, '--r'); xlim([t(1) t(end)]);
        title(sprintf('%s: %s', in.title, in.labels{c}), 'FontWeight', 'normal', 'FontSize', 8);
        xlabel(sprintf('V1 window %s (s)', ternary(strcmp(X.alignment, 'centred'), 'centre', 'start')));
        ylabel('AUC'); set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 8);
    end
end
linkaxes(findobj(fig, 'Type', 'axes'), 'xy');
sgtitle({name_t, [decoder_label(RES(1).info) ', bootstrapped temporal AUC (95% CI; black: shuffled)']}, 'FontSize', 10);
save_fig(fig, save_path);
end

function plot_summary_temporal(RES, unit, save_path)
% One tile per analysis: mean bootstrapped AUC trace of every unit (thin) and
% the across-unit mean +/- SEM (thick), per condition; no shuffle
nA = numel(RES); nU = numel(RES(1).info.unit_ids);
fig = figure('Visible', 'off', 'Position', [100 100 420 * min(nA, 2) 330 * ceil(nA / 2)]);
tiledlayout(ceil(nA / 2), min(nA, 2), 'TileSpacing', 'compact');
for iA = 1:nA
    X = RES(iA).temporal; in = RES(iA).info; t = X.time_bins;
    nexttile; hold on;
    h = gobjects(1, numel(in.labels));
    for c = 1:numel(in.labels)
        Y = X.mean(:, :, c);
        plot(t, Y', '-', 'Color', in.colors(c, :) * 0.85, 'LineWidth', 0.6);   % opaque (as rendered in Affinity)
        m = mean(Y, 1, 'omitnan'); se = std(Y, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(Y), 1));
        ok = isfinite(m) & isfinite(se);
        fill([t(ok) fliplr(t(ok))], [m(ok) - se(ok) fliplr(m(ok) + se(ok))], in.colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        h(c) = plot(t, m, 'Color', in.colors(c, :) * 0.7, 'LineWidth', 2.5);
    end
    xline(0, '--k'); yline(0, 'k-'); xlim([t(1) t(end)]);
    legend(h, in.labels, 'Location', 'northwest', 'Box', 'off', 'FontSize', 8);
    title(in.title, 'FontWeight', 'normal', 'FontSize', 9);
    xlabel(sprintf('V1 window %s (s rel. ripple onset)', ternary(strcmp(X.alignment, 'centred'), 'centre', 'start')));
    ylabel('Bootstrapped mean AUC');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 9);
end
sgtitle({sprintf('Per-%s temporal AUC (n = %d %ss), %s', unit, nU, unit, decoder_label(RES(1).info)), ...
    sprintf('thin: each %s''s bootstrapped mean; thick: mean \\pm SEM across %ss', unit, unit)}, 'FontSize', 10);
save_fig(fig, save_path);
end

function plot_summary_temporal_diff(RES, unit, save_path)
% One tile per two-condition analysis: bootstrapped mean AUC difference of every
% unit (thin) and the across-unit mean +/- SEM (thick): High - Low, Trough - Peak,
% Unilateral - Bilateral (sign flipped from the saved Bilateral - Unilateral)
RES = RES(arrayfun(@(X) numel(X.info.labels), RES) > 1);
nA = numel(RES); nU = numel(RES(1).info.unit_ids);
fig = figure('Visible', 'off', 'Position', [100 100 420 * min(nA, 2) 330 * ceil(nA / 2)]);
tiledlayout(ceil(nA / 2), min(nA, 2), 'TileSpacing', 'compact');
for iA = 1:nA
    X = RES(iA).temporal; in = RES(iA).info; t = X.time_bins;
    if strcmp(in.analysis, 'SO_trough_uni_bi'), sgn = -1; lab = {in.labels{1}, in.labels{2}};
    else, sgn = 1; lab = {in.labels{end}, in.labels{1}}; end
    col = in.colors(end, :);
    Y = sgn * X.diff_mean;
    nexttile; hold on;
    plot(t, Y', '-', 'Color', col * 0.85, 'LineWidth', 0.6);
    m = mean(Y, 1, 'omitnan'); se = std(Y, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(Y), 1));
    ok = isfinite(m) & isfinite(se);
    fill([t(ok) fliplr(t(ok))], [m(ok) - se(ok) fliplr(m(ok) + se(ok))], col, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(t, m, 'Color', col * 0.7, 'LineWidth', 2.5);
    xline(0, '--k'); yline(0, 'k-'); xlim([t(1) t(end)]);
    title({in.title, sprintf('%s - %s', lab{:})}, 'FontWeight', 'normal', 'FontSize', 9);
    xlabel(sprintf('V1 window %s (s rel. ripple onset)', ternary(strcmp(X.alignment, 'centred'), 'centre', 'start')));
    ylabel('\Delta bootstrapped mean AUC');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 9);
end
sgtitle({sprintf('Per-%s temporal AUC difference (n = %d %ss), %s', unit, nU, unit, decoder_label(RES(1).info)), ...
    sprintf('thin: each %s''s bootstrapped mean difference; thick: mean \\pm SEM across %ss', unit, unit)}, 'FontSize', 10);
save_fig(fig, save_path);
end
