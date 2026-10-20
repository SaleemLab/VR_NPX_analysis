function run_Lseq_RUN_decoding_CV(varargin)
% 10-fold cross-validated sequence-likelihood (L_seq; Huh et al. 2026, Nat Commun
% 17:8304) decoding of RUN1 laps with Track 1 and Track 2 templates, per region (all
% spatial cells; RUN1_decoding_input.mat).
% Templates (Lseq_template.m) are built from running spikes (speed > threshold) of the
% training laps of each track (laps = trials). Held-out laps are scored with both
% templates (Lseq_score.m, cell-identity shuffles) as
%   'lap'    : all running spikes of the lap (whole-track sequence)
%   'pseudo' : replay-like events of 'pseudo_spikes' spikes sampled at random from the
%              lap and kept in temporal order (as the paper's simulations), >= min_cells
% Default-sequence check: each template is split into a shared component
% W_shared = (W_T1 + W_T2)/2 and a track-specific component W_specific = (W_T1 - W_T2)/2
% (L_seq(T1) = L_shared + L_specific, L_seq(T2) = L_shared - L_specific). A sequence
% common to both tracks (e.g. intrinsic co-firing order) gives positive L_seq with the
% other track's template and with W_shared. Whole-session diagnostics also compare
% templates: split-half reliability within a track vs similarity across tracks, a
% stationary-period template (speed < still_speed, epochs >= still_min_dur s), and the
% order predicted by rate-map peaks.
% Track assignment: an event is assigned to a track whose template gives a two-sided
% shuffle p < 0.05; if both do, the larger |z| wins. Track bias = |z_T1| - |z_T2|.
% Options: 'sessions' (1:22), 'regions' ({'HPC','V1','V1_L','V1_R'}),
% 'speed_threshold' (5), 'still_speed' (2), 'still_min_dur' (1), 'nShuffle' (500),
% 'n_pseudo' (20 per lap), 'pseudo_spikes' (10), 'min_cells' (5), 'min_pairs' (20 spike
% pairs for a cell pair to enter template correlations), 'out_dir'
% (nature_neuroscience_revision_analysis\Lseq_RUN_decoding), 'overwrite', 'summary_only',
% 'replot' (redraw session figures of sessions already done).
% Output: out_dir\<animal>\session_<date>.mat (+ figures), per-animal and all-animal
% summaries, out_dir\metrics.csv.

p = inputParser;
addParameter(p, 'sessions', 1:22, @isnumeric);
addParameter(p, 'regions', {'HPC', 'V1', 'V1_L', 'V1_R'}, @iscell);
addParameter(p, 'speed_threshold', 5, @isnumeric);
addParameter(p, 'still_speed', 2, @isnumeric);
addParameter(p, 'still_min_dur', 1, @isnumeric);
addParameter(p, 'nShuffle', 500, @isnumeric);
addParameter(p, 'n_pseudo', 20, @isnumeric);
addParameter(p, 'pseudo_spikes', 10, @isnumeric);
addParameter(p, 'min_cells', 5, @isnumeric);
addParameter(p, 'min_pairs', 20, @isnumeric);
addParameter(p, 'out_dir', '', @ischar);
addParameter(p, 'overwrite', false, @islogical);
addParameter(p, 'summary_only', false, @islogical);
addParameter(p, 'replot', false, @islogical);
parse(p, varargin{:});
P = p.Results;
P.pos_edges = 0:10:140;
P.nFolds = 10;
P.event_types = {'lap', 'pseudo'};
P.template_names = {'T1', 'T2', 'shared', 'specific'};

revision_folder = get_revision_folder();
if isempty(P.out_dir), P.out_dir = fullfile(revision_folder, 'Lseq_RUN_decoding'); end

if ~P.summary_only
    S = load(fullfile(revision_folder, 'RUN1_decoding_input.mat'), 'D'); D = S.D;
    warning('off', 'stats:glmfit:IterationLimit'); warning('off', 'stats:glmfit:BadScaling');
    warning('off', 'stats:glmfit:PerfectSeparation');
    for nsession = P.sessions
        score_session(D(nsession), nsession, P);
    end
end
summarise_all(P);
end

function score_session(Ds, nsession, P)
out_file = fullfile(P.out_dir, Ds.subject, sprintf('session_%s.mat', Ds.session));
if exist(out_file, 'file') && ~P.overwrite
    fprintf('Session %d exists, skipping\n', nsession);
    if P.replot
        S = load(out_file, 'out'); plot_session(S.out, out_file, P);
    end
    return
end
fprintf('\n=== Session %d: %s %s ===\n', nsession, Ds.subject, Ds.session);
[beh, spk, folds] = RUN1_session_behaviour(Ds, nsession, P);
nLaps = size(Ds.lap_times, 1);
nPos = numel(P.pos_edges) - 1;
still_epochs = stationary_epochs(beh, P);

out = struct();
out.subject = Ds.subject; out.session = Ds.session; out.session_index = nsession; out.params = P;
out.lap_times = Ds.lap_times; out.lap_track = Ds.lap_track; out.lap_fold = folds;
out.n_still_epochs = size(still_epochs, 1);
out.still_duration = sum(diff(still_epochs, 1, 2));
for nr = 1:numel(P.regions)
    tic
    cells = region_cells(Ds.region, P.regions{nr});
    nC = sum(cells);
    new_idx = zeros(numel(cells), 1); new_idx(cells) = 1:nC;
    rs = RandStream('mrg32k3a', 'Seed', 1000 * nsession + nr);
    R = struct('name', P.regions{nr}, 'cell_id', Ds.cell_id(cells), 'n_cells', nC);

    % Ordered pair counts per lap (running spikes)
    in_region = cells(spk.cell);
    lap_cells = cell(nLaps, 1); lap_t = lap_cells;
    Mlap = zeros(nC, nC, nLaps);
    for k = 1:nLaps
        s = spk.running & in_region & spk.lap == k;
        [lap_t{k}, order] = sort(spk.t(s));
        c = new_idx(spk.cell(s)); lap_cells{k} = c(order);
        Mlap(:, :, k) = Lseq_pair_counts(spk.t(s), new_idx(spk.cell(s)), nC);
    end

    % Whole-session template diagnostics
    R.diagnostics = template_diagnostics(Mlap, Ds.lap_track, spk, in_region, new_idx, nC, still_epochs, beh, P, nPos);

    % Cross-validated scoring of held-out laps
    nW = 4;
    ev = struct('lap', zeros(0, 1), 'track', zeros(0, 1), 'fold', zeros(0, 1), 'type', zeros(0, 1), ...
        'n_spikes', zeros(0, 1), 'n_cells', zeros(0, 1), 'L', zeros(0, nW), 'z', zeros(0, nW), ...
        'p_fwd', zeros(0, nW), 'p_rev', zeros(0, nW));
    for f = 1:P.nFolds
        train = folds ~= f;
        T1 = Lseq_template(sum(Mlap(:, :, train & Ds.lap_track == 1), 3));
        T2 = Lseq_template(sum(Mlap(:, :, train & Ds.lap_track == 2), 3));
        Ws = {T1.W, T2.W, (T1.W + T2.W) / 2, (T1.W - T2.W) / 2};
        for k = find(folds == f)'
            c = lap_cells{k}; tk = lap_t{k};
            if numel(unique(c)) >= P.min_cells
                ev = add_event(ev, Lseq_score(tk, c, Ws, P.nShuffle, rs), k, Ds.lap_track(k), f, 1);
            end
            if numel(c) < P.pseudo_spikes, continue; end
            for e = 1:P.n_pseudo
                for attempt = 1:50
                    pick = sort(randperm(rs, numel(c), P.pseudo_spikes));
                    if numel(unique(c(pick))) >= P.min_cells, break; end
                end
                if numel(unique(c(pick))) < P.min_cells, continue; end
                ev = add_event(ev, Lseq_score(tk(pick), c(pick), Ws, P.nShuffle, rs), k, Ds.lap_track(k), f, 2);
            end
        end
    end
    R.events = ev;
    R.confusion = confusion_summary(ev, P);
    out.region(nr) = R;
    for t = 1:numel(P.event_types)
        C = R.confusion(t);
        fprintf('  %-4s (%3d cells) %-6s n=%4d | mean z [T1|T1 %.1f, T1|T2 %.1f; T2|T1 %.1f, T2|T2 %.1f] | sig fwd [%.2f %.2f; %.2f %.2f] | assigned %.2f, acc %.2f, bias %.2f\n', ...
            P.regions{nr}, nC, P.event_types{t}, C.n_events(1) + C.n_events(2), C.mean_z(1, 1), C.mean_z(1, 2), ...
            C.mean_z(2, 1), C.mean_z(2, 2), C.frac_fwd(1, 1), C.frac_fwd(1, 2), C.frac_fwd(2, 1), C.frac_fwd(2, 2), ...
            mean(C.frac_assigned), C.track_acc, mean(C.mean_bias));
    end
    d = R.diagnostics;
    fprintf('       template r: T1 split %.2f, T2 split %.2f, T1 vs T2 (cross-split) %.2f, still vs T1 %.2f, still vs T2 %.2f (%.0f s)\n', ...
        d.r_T1_split, d.r_T2_split, d.r_T1_T2_cross_split, d.r_still_T1, d.r_still_T2, toc);
end
folder = fileparts(out_file);
if ~exist(folder, 'dir'), mkdir(folder); end
save(out_file, 'out', '-v7.3');
plot_session(out, out_file, P);
end

function plot_session(out, out_file, P)
try
    plot_Lseq_summary({out.region.confusion}, {out.region.diagnostics}, {out.region.name}, [out.region.n_cells], ...
        sprintf('%s %s (session %d): L_{seq}, 10-fold CV', out.subject, out.session, out.session_index), strrep(out_file, '.mat', ''), P);
    plot_template_scatter(out, P, strrep(out_file, '.mat', '_template_scatter'));
catch err
    warning('Session %d plotting failed: %s', out.session_index, err.message);
end
close all
end

function ev = add_event(ev, S, lap, track, fold, type)
ev.lap(end+1, 1) = lap; ev.track(end+1, 1) = track; ev.fold(end+1, 1) = fold; ev.type(end+1, 1) = type;
ev.n_spikes(end+1, 1) = S.n_spikes; ev.n_cells(end+1, 1) = S.n_cells;
ev.L(end+1, :) = S.L; ev.z(end+1, :) = S.z; ev.p_fwd(end+1, :) = S.p_fwd; ev.p_rev(end+1, :) = S.p_rev;
end

function cells = region_cells(region, name)
switch name
    case 'HPC', cells = contains(region, 'HPC');
    case 'V1', cells = contains(region, 'V1');
    otherwise, cells = region == string(name);
end
end

function E = stationary_epochs(beh, P)
% [start end] of contiguous periods with speed < still_speed lasting >= still_min_dur
still = beh.speed < P.still_speed;
d = diff([false; still; false]);
on = find(d == 1); off = find(d == -1) - 1;
E = [beh.tvec(on) beh.tvec(off)];
E = E(E(:, 2) - E(:, 1) >= P.still_min_dur, :);
end

function C = confusion_summary(ev, P)
% Per event type: rows = true track, columns = template (T1, T2)
for t = 1:numel(P.event_types)
    C(t).type = P.event_types{t}; %#ok<AGROW>
    for a = 1:2
        sel = ev.type == t & ev.track == a & all(isfinite(ev.z(:, 1:2)), 2);
        C(t).n_events(a) = sum(sel);
        for b = 1:2
            C(t).mean_z(a, b) = mean(ev.z(sel, b));
            C(t).median_z(a, b) = median(ev.z(sel, b));
            C(t).mean_L(a, b) = mean(ev.L(sel, b));
            C(t).frac_fwd(a, b) = mean(ev.p_fwd(sel, b) < 0.05);
            C(t).frac_rev(a, b) = mean(ev.p_rev(sel, b) < 0.05);
        end
        % Track assignment: two-sided shuffle p < 0.05 on a track template; if both
        % tracks are significant, the larger |z| wins (others unassigned).
        % Track bias = |z_T1| - |z_T2| (positive = Track 1), here signed to the true track.
        p2 = min(1, 2 * min(ev.p_fwd(sel, 1:2), ev.p_rev(sel, 1:2)));
        az = abs(ev.z(sel, 1:2));
        sig = p2 < 0.05;
        assigned = zeros(sum(sel), 1);
        assigned(sig(:, 1) & (~sig(:, 2) | az(:, 1) >= az(:, 2))) = 1;
        assigned(sig(:, 2) & (~sig(:, 1) | az(:, 2) > az(:, 1))) = 2;
        C(t).classification(a, :) = [mean(assigned == 1) mean(assigned == 2)];   % fractions of all events
        C(t).frac_assigned(a) = mean(assigned > 0);
        C(t).acc_assigned(a) = mean(assigned(assigned > 0) == a);
        C(t).mean_bias(a) = mean((3 - 2 * a) * (az(:, 1) - az(:, 2)));
        % Shared and specific components; specific signed towards the true track
        sgn = 3 - 2 * a;
        C(t).mean_z_shared(a) = mean(ev.z(sel, 3));
        C(t).mean_z_specific(a) = mean(sgn * ev.z(sel, 4));
        C(t).frac_fwd_shared(a) = mean(ev.p_fwd(sel, 3) < 0.05);
        if a == 1, C(t).frac_sig_specific(a) = mean(ev.p_fwd(sel, 4) < 0.05);
        else, C(t).frac_sig_specific(a) = mean(ev.p_rev(sel, 4) < 0.05); end
    end
    C(t).track_acc = mean(C(t).acc_assigned, 'omitnan');   % among assigned events
end
end

function d = template_diagnostics(Mlap, lap_track, spk, in_region, new_idx, nC, still_epochs, beh, P, nPos)
% Similarity of whole-session templates (W, upper triangle of cell pairs with at
% least min_pairs ordered spike pairs in both templates)
laps1 = find(lap_track == 1); laps2 = find(lap_track == 2);
T1 = Lseq_template(sum(Mlap(:, :, laps1), 3));
T2 = Lseq_template(sum(Mlap(:, :, laps2), 3));
T1o = Lseq_template(sum(Mlap(:, :, laps1(1:2:end)), 3)); T1e = Lseq_template(sum(Mlap(:, :, laps1(2:2:end)), 3));
T2o = Lseq_template(sum(Mlap(:, :, laps2(1:2:end)), 3)); T2e = Lseq_template(sum(Mlap(:, :, laps2(2:2:end)), 3));
% Stationary-period template: each stationary epoch is a trial
Ms = zeros(nC);
s_all = in_region;
for k = 1:size(still_epochs, 1)
    s = s_all & spk.t >= still_epochs(k, 1) & spk.t <= still_epochs(k, 2);
    if sum(s) > 1, Ms = Ms + Lseq_pair_counts(spk.t(s), new_idx(spk.cell(s)), nC); end
end
Tst = Lseq_template(Ms);
d.r_T1_T2 = wcorr(T1, T2, P);
d.r_T1_split = wcorr(T1o, T1e, P);
d.r_T2_split = wcorr(T2o, T2e, P);
d.r_T1_T2_cross_split = mean([wcorr(T1o, T2e, P) wcorr(T1e, T2o, P)]);
d.r_still_T1 = wcorr(Tst, T1, P);
d.r_still_T2 = wcorr(Tst, T2, P);
d.n_still_pairs = sum(Ms(:));
% Order predicted by rate-map peaks (running samples, all laps)
r = beh.running;
occ = accumarray([beh.track(r) beh.pos_bin(r)], beh.dt(r), [2 nPos]);
s = spk.running & in_region;
cnt = accumarray([new_idx(spk.cell(s)) spk.pos_bin(s) spk.track(s)], 1, [nC nPos 2]);
peak = NaN(nC, 2);
for track = 1:2
    rate = cnt(:, :, track) ./ occ(track, :);
    [mx, peak(:, track)] = max(rate, [], 2);
    peak(mx == 0, track) = NaN;
end
d.peak_bin = peak;
U = triu(true(nC), 1);
order1 = sign(peak(:, 1)' - peak(:, 1)); order2 = sign(peak(:, 2)' - peak(:, 2));   % (a,b): +1 if b peaks after a
ok = U & order1 ~= 0 & order2 ~= 0 & ~isnan(order1) & ~isnan(order2);
d.peak_order_agreement = mean(order1(ok) == order2(ok));   % 0.5 = unrelated orders across tracks
d.r_W_T1_peak_T1 = corr_ok(T1.W, order1, T1.N >= P.min_pairs & U & order1 ~= 0);
d.r_W_T2_peak_T2 = corr_ok(T2.W, order2, T2.N >= P.min_pairs & U & order2 ~= 0);
d.r_W_T1_peak_T2 = corr_ok(T1.W, order2, T1.N >= P.min_pairs & U & order2 ~= 0);
d.r_W_T2_peak_T1 = corr_ok(T2.W, order1, T2.N >= P.min_pairs & U & order1 ~= 0);
d.W_T1 = single(T1.W); d.W_T2 = single(T2.W); d.W_still = single(Tst.W);
d.N_T1 = single(T1.N); d.N_T2 = single(T2.N); d.N_still = single(Tst.N);
end

function r = wcorr(A, B, P)
U = triu(true(size(A.W)), 1);
ok = U & A.N >= P.min_pairs & B.N >= P.min_pairs;
r = corr_ok(A.W, B.W, ok);
end

function r = corr_ok(x, y, ok)
if sum(ok(:)) < 10, r = NaN; return; end
r = corr(x(ok), y(ok));
end

function plot_Lseq_summary(Cs, Ds, names, n_cells, title_str, save_path, P)
% Rows = regions. Columns: 2 x 2 (true track x template) mean shuffle z and fraction
% significant forward for whole laps and pseudo-events, pseudo-event track
% classification, shared vs specific components, and template correlations.
nR = numel(Cs);
fig = figure('Name', title_str, 'Position', [30 30 2000 300 * nR], 'Color', 'w');
tl = tiledlayout(nR, 7, 'TileSpacing', 'compact', 'Padding', 'compact');
for nr = 1:nR
    C = Cs{nr};
    for t = 1:2
        zmax = max(abs([C(1).mean_z(:); C(2).mean_z(:)]), [], 'omitnan');
        nexttile; heat2x2(C(t).mean_z, [-zmax zmax], diverging_map(), '%.1f');
        title(sprintf('%s (n=%d) %s: mean z', names{nr}, n_cells(nr), C(t).type), 'FontWeight', 'normal', 'Interpreter', 'none');
        nexttile; heat2x2(C(t).frac_fwd, [0 1], flipud(gray), '%.2f');
        title(sprintf('%s: frac. sig. forward', C(t).type), 'FontWeight', 'normal');
    end
    nexttile; heat2x2(C(2).classification, [0 1], flipud(gray), '%.2f');
    xlabel({'Assigned track (p < 0.05, larger |z|)', 'fraction of all events'});
    title(sprintf('pseudo: assigned %.2f, acc %.2f, bias %.2f', mean(C(2).frac_assigned), C(2).track_acc, mean(C(2).mean_bias)), 'FontWeight', 'normal');
    nexttile; hold on;
    y = [C(1).mean_z_shared(:) C(1).mean_z_specific(:)];
    b = bar(y, 'grouped'); b(1).FaceColor = [0.65 0.65 0.65]; b(2).FaceColor = [0.15 0.15 0.15];
    set(gca, 'XTick', 1:2, 'XTickLabel', {'T1 laps', 'T2 laps'}, 'TickDir', 'out');
    ylabel('Mean z (whole laps)');
    title({'Shared (grey) vs specific (black)', 'component, specific signed to true track'}, 'FontWeight', 'normal');
    nexttile; hold on;
    d = Ds{nr};
    v = [d.r_T1_split d.r_T2_split d.r_T1_T2_cross_split d.r_still_T1 d.r_still_T2];
    bar(v, 'FaceColor', [0.5 0.5 0.5]);
    set(gca, 'XTick', 1:5, 'XTickLabel', {'T1 half', 'T2 half', 'T1-T2', 'still-T1', 'still-T2'}, ...
        'XTickLabelRotation', 45, 'TickDir', 'out');
    ylabel('Template correlation (W)'); ylim([min([0 v]) - 0.05, 1]);
    title(sprintf('peak-order agreement %.2f', d.peak_order_agreement), 'FontWeight', 'normal');
end
title(tl, title_str);
save_figure(fig, save_path);
end

function heat2x2(M, cl, cmap, fmt)
imagesc(M, cl); colormap(gca, cmap); axis square; colorbar;
for a = 1:2
    for b = 1:2
        text(b, a, sprintf(fmt, M(a, b)), 'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', [0.1 0.6 0.1], 'FontWeight', 'bold');
    end
end
set(gca, 'XTick', 1:2, 'XTickLabel', {'T1 template', 'T2 template'}, 'YTick', 1:2, 'YTickLabel', {'T1 data', 'T2 data'}, 'TickDir', 'out');
end

function cmap = diverging_map()
n = 128;
up = [linspace(0.13, 1, n)' linspace(0.4, 1, n)' linspace(0.67, 1, n)'];
dn = [linspace(1, 0.7, n)' linspace(1, 0.09, n)' linspace(1, 0.17, n)'];
cmap = [up; dn];
end

function plot_template_scatter(out, P, save_path)
% Pairwise W of Track 1 vs Track 2 templates (whole session) and of the stationary
% template vs Track 1, per region
nR = numel(out.region);
fig = figure('Name', 'template scatter', 'Position', [50 50 300 * nR 600], 'Color', 'w');
tl = tiledlayout(2, nR, 'TileSpacing', 'compact');
for row = 1:2
    for nr = 1:nR
        d = out.region(nr).diagnostics;
        U = triu(true(size(d.W_T1)), 1);
        if row == 1
            x = d.W_T1; y = d.W_T2; ok = U & d.N_T1 >= P.min_pairs & d.N_T2 >= P.min_pairs;
            xl = 'W Track 1'; yl = 'W Track 2'; r = d.r_T1_T2;
        else
            x = d.W_T1; y = d.W_still; ok = U & d.N_T1 >= P.min_pairs & d.N_still >= P.min_pairs;
            xl = 'W Track 1'; yl = 'W stationary'; r = d.r_still_T1;
        end
        idx = find(ok);
        if numel(idx) > 5000, idx = idx(round(linspace(1, numel(idx), 5000))); end
        nexttile; hold on;
        scatter(x(idx), y(idx), 4, 'k', 'filled', 'MarkerFaceAlpha', 0.2);
        plot([-5 5], [-5 5], 'r--'); axis square; xlim([-5 5]); ylim([-5 5]);
        xlabel(xl); ylabel(yl);
        title(sprintf('%s: r = %.2f (%d pairs)', out.region(nr).name, r, sum(ok(:))), 'FontWeight', 'normal', 'Interpreter', 'none');
        set(gca, 'TickDir', 'out');
    end
end
title(tl, sprintf('%s %s: pairwise order log odds W (whole session templates)', out.subject, out.session));
save_figure(fig, save_path);
print(fig, [save_path '.svg'], '-dsvg', '-vector');
end

function summarise_all(P)
files = dir(fullfile(P.out_dir, 'M*', 'session_*.mat'));
if isempty(files), return; end
sessions = struct('subject', {}, 'C', {}, 'D', {}, 'names', {}, 'n_cells', {});
rows = {};
for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name), 'out'); out = S.out;
    sessions(k).subject = out.subject;
    sessions(k).C = {out.region.confusion};
    sessions(k).D = cellfun(@(d) rmfield(d, {'W_T1', 'W_T2', 'W_still', 'N_T1', 'N_T2', 'N_still', 'peak_bin'}), ...
        {out.region.diagnostics}, 'UniformOutput', false);
    sessions(k).names = {out.region.name};
    sessions(k).n_cells = [out.region.n_cells];
    for nr = 1:numel(out.region)
        d = sessions(k).D{nr};
        for t = 1:numel(P.event_types)
            C = out.region(nr).confusion(t);
            rows(end+1, :) = {out.subject, out.session, out.region(nr).name, out.region(nr).n_cells, C.type, ...
                sum(C.n_events), C.mean_z(1, 1), C.mean_z(1, 2), C.mean_z(2, 1), C.mean_z(2, 2), ...
                C.frac_fwd(1, 1), C.frac_fwd(1, 2), C.frac_fwd(2, 1), C.frac_fwd(2, 2), C.track_acc, ...
                mean(C.frac_assigned), mean(C.mean_bias), ...
                mean(C.mean_z_shared), mean(C.mean_z_specific), d.r_T1_split, d.r_T2_split, d.r_T1_T2_cross_split, ...
                d.r_still_T1, d.r_still_T2, d.peak_order_agreement}; %#ok<AGROW>
        end
    end
end
names = sessions(1).names;
animals = unique({sessions.subject});
animal_C = cell(1, numel(animals)); animal_D = animal_C;
for na = 1:numel(animals)
    in = strcmp({sessions.subject}, animals{na});
    animal_C{na} = average_structs({sessions(in).C});
    animal_D{na} = average_structs({sessions(in).D});
    plot_Lseq_summary(animal_C{na}, animal_D{na}, names, round(mean(vertcat(sessions(in).n_cells), 1)), ...
        sprintf('%s: L_{seq}, %d sessions (mean n cells)', animals{na}, sum(in)), ...
        fullfile(P.out_dir, animals{na}, sprintf('%s_summary', animals{na})), P);
    close all
end
plot_Lseq_summary(average_structs(animal_C), average_structs(animal_D), names, round(mean(vertcat(sessions.n_cells), 1)), ...
    sprintf('All animals (%d animals, %d sessions): L_{seq}, mean of animal means', numel(animals), numel(sessions)), ...
    fullfile(P.out_dir, 'all_animals_summary'), P);
close all
T = cell2table(rows, 'VariableNames', {'subject', 'session', 'region', 'n_cells', 'event_type', 'n_events', ...
    'z_T1data_T1tmpl', 'z_T1data_T2tmpl', 'z_T2data_T1tmpl', 'z_T2data_T2tmpl', 'fwd_T1data_T1tmpl', ...
    'fwd_T1data_T2tmpl', 'fwd_T2data_T1tmpl', 'fwd_T2data_T2tmpl', 'track_acc', 'frac_assigned', 'track_bias', ...
    'z_shared', 'z_specific', ...
    'r_T1_split', 'r_T2_split', 'r_T1_T2_cross_split', 'r_still_T1', 'r_still_T2', 'peak_order_agreement'});
writetable(T, fullfile(P.out_dir, 'metrics.csv'));
plot_session_metrics(T, P);
end

function A = average_structs(list)
% Element-wise mean over a cell array of (cell arrays of) structs with numeric fields
nR = numel(list{1});
A = cell(1, nR);
for nr = 1:nR
    items = cellfun(@(x) x{nr}, list, 'UniformOutput', false);
    A{nr} = items{1};
    for e = 1:numel(A{nr})
        for fn = fieldnames(A{nr})'
            v = cellfun(@(x) x(e).(fn{1}), items, 'UniformOutput', false);
            if isnumeric(v{1}) || islogical(v{1})
                A{nr}(e).(fn{1}) = mean(cat(3, v{:}), 3, 'omitnan');
            end
        end
    end
end
end

function plot_session_metrics(T, P)
% Per session (dots coloured by animal): whole-lap z with own vs other-track
% template, shared and specific components, track accuracy, and template correlations
animals = unique(T.subject);
cmap = lines(numel(animals));
Tp = T(strcmp(T.event_type, 'lap'), :);
measures = {{'z_T1data_T1tmpl', 'z_T2data_T2tmpl'}, {'z_T1data_T2tmpl', 'z_T2data_T1tmpl'}, {'z_shared'}, {'z_specific'}, ...
    {'track_acc'}, {'r_T1_split', 'r_T2_split'}, {'r_T1_T2_cross_split'}, {'r_still_T1', 'r_still_T2'}};
labels = {'z, own template', 'z, other-track template', 'z, shared component', 'z, specific component', ...
    'Track accuracy', 'r within track (split half)', 'r T1 vs T2 (cross split)', 'r stationary vs track'};
fig = figure('Name', 'Lseq sessions', 'Position', [50 50 1800 650], 'Color', 'w');
tl = tiledlayout(2, 4, 'TileSpacing', 'compact');
for k = 1:numel(measures)
    nexttile; hold on;
    for nr = 1:numel(P.regions)
        in = strcmp(Tp.region, P.regions{nr});
        y = mean(Tp{in, measures{k}}, 2);
        [~, ia] = ismember(Tp.subject(in), animals);
        scatter(nr + 0.15 * randn(sum(in), 1), y, 20, cmap(ia, :), 'filled', 'MarkerFaceAlpha', 0.7);
        plot(nr + [-0.3 0.3], median(y, 'omitnan') * [1 1], 'k-', 'LineWidth', 2);
    end
    plot([0.5 numel(P.regions) + 0.5], (k == 5) * 0.5 * [1 1], 'k--');
    set(gca, 'XTick', 1:numel(P.regions), 'XTickLabel', P.regions, 'TickDir', 'out', 'TickLabelInterpreter', 'none');
    ylabel(labels{k});
end
title(tl, 'L_{seq} RUN decoding (10-fold CV, held-out whole laps; template r from whole session): sessions (colour = animal), bar = median');
save_figure(fig, fullfile(P.out_dir, 'session_metrics'));
close(fig);
end

function save_figure(fig, save_path)
savefig(fig, [save_path '.fig']);
exportgraphics(fig, [save_path '.png'], 'Resolution', 150);
exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
end

function folder = get_revision_folder()
if exist('D:\corticohippocampal_replay', 'dir') > 0
    folder = 'D:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    folder = 'P:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
end
end
