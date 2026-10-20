function run_Lseq_replay_decoding_POST(varargin)
% Sequence-likelihood (L_seq; Huh et al. 2026, Nat Commun 17:8304) replay scoring of
% POST sleep ripples, on the same events and windows as run_Bayesian_replay_decoding_POST
% (read from Bayesian_replay_decoding\<animal>\session_<date>.mat):
%   'post' = onset to onset + n_post_bins x 20 ms (max(0.1 s, ripple duration))
%   'pre'  = the same duration immediately before onset
% Templates per region (HPC, whole V1; all spatial cells) from running spikes of all
% RUN1 laps of each track (Lseq_template.m). Track-specific templates remove the part
% of each track's template shared with the other track:
%   W1_specific = W1 - b12 W2, b12 = <W1,W2>/<W2,W2> (over all cell pairs), and vice versa
% so a Track 2 sequence (forward or reverse) scores ~0 with W1_specific, and the sign of
% L_seq gives forward (+) or reverse (-). Each event (spikes in the window, temporal
% order) is scored with W1, W2, W1_specific and W2_specific (Lseq_score definition:
% k = 2 sequence likelihood, forward - reverse, per spike-pair term) against nShuffle
% cell-identity shuffles (permutation of the participating cells; same permutations for
% all templates). Two-sided p = min(1, 2 min(p_forward, p_reverse)); significant if
% p < 0.05. Seeded per session, region and window, RandStream substream = ripple_info row.
% Track assignment (summary) uses the raw T1/T2 templates: significant on one track, or
% the larger |z| if both; track bias = |z_T1| - |z_T2|. The overlap-removed templates are
% kept for reference only.
% Events need >= min_cells participating cells.
% Options: 'sessions' (1:22), 'nShuffle' (1000), 'min_cells' (5), 'speed_threshold' (5),
% 'seed' (2026), 'out_dir' (...\Lseq_replay_decoding), 'overwrite',
% 'summary_only'.
% Output: out_dir\<animal>\session_<date>.mat (+ examples figure); out_dir\
% Lseq_replay_events_all.mat/.csv, Lseq_replay_session_summary.csv, Lseq_replay_summary.

p = inputParser;
addParameter(p, 'sessions', 1:22, @isnumeric);
addParameter(p, 'nShuffle', 1000, @isnumeric);
addParameter(p, 'min_cells', 5, @isnumeric);
addParameter(p, 'speed_threshold', 5, @isnumeric);
addParameter(p, 'seed', 2026, @isnumeric);
addParameter(p, 'out_dir', '', @ischar);
addParameter(p, 'overwrite', false, @islogical);
addParameter(p, 'summary_only', false, @islogical);
parse(p, varargin{:});
P = p.Results;
P.pos_edges = 0:10:140; P.nPos = 14; P.nFolds = 10;
P.regions = {'HPC', 'V1'};
P.windows = {'post', 'pre'};
P.templates = {'T1', 'T2', 'T1_specific', 'T2_specific'};

revision_folder = get_revision_folder();
P.bayes_dir = fullfile(revision_folder, 'Bayesian_replay_decoding');
if isempty(P.out_dir), P.out_dir = fullfile(revision_folder, 'Lseq_replay_decoding'); end
if ~exist(P.out_dir, 'dir'), mkdir(P.out_dir); end

if ~P.summary_only
    addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
    S = load(fullfile(revision_folder, 'RUN1_decoding_input.mat'), 'D'); D = S.D;
    SUBJECTS = {'M24016', 'M24017', 'M24018', 'M24062', 'M24064', 'M24065'};
    experiment_info = subject_session_stimuli_mapping(SUBJECTS, 'bilateral');
    experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
    for nsession = P.sessions
        score_session(D(nsession), experiment_info(nsession), nsession, P);
    end
end
summarise_all(P);
end

function score_session(Ds, info, nsession, P)
out_file = fullfile(P.out_dir, Ds.subject, sprintf('session_%s.mat', Ds.session));
if exist(out_file, 'file') && ~P.overwrite
    fprintf('Session %d exists, skipping\n', nsession);
    return
end
tic
B = load(fullfile(P.bayes_dir, Ds.subject, sprintf('session_%s.mat', Ds.session)), 'out'); B = B.out;
Es = B.events; nEv = numel(Es.onset);
tb = B.params.time_bin;
win = {[Es.onset, Es.onset + Es.n_post_bins * tb], [Es.onset - Es.n_post_bins * tb, Es.onset]};
fprintf('\n=== Session %d: %s %s, %d ripples ===\n', nsession, Ds.subject, Ds.session, nEv);

[~, spk] = RUN1_session_behaviour(Ds, nsession, P);
is_post = strcmp(info.StimulusName, 'SleepChronic');
options = info.session(is_post).probe(1);
S = load(fullfile(options.ANALYSIS_DATAPATH, '..', 'session_clusters_Sleep.mat'), 'session_clusters');
sc = S.session_clusters;
keep = ismember(sc.cluster_id, Ds.cell_id);
t = vertcat(sc.spike_times{keep}); id = vertcat(sc.spike_id{keep});
[~, cell_idx] = ismember(id, Ds.cell_id);
[t, o] = sort(t); cell_idx = cell_idx(o);

out = struct();
out.subject = Ds.subject; out.session = Ds.session; out.session_index = nsession; out.params = P;
out.events = Es; out.time_bin = tb; out.cell_id = Ds.cell_id; out.region_labels = Ds.region; out.rate_maps = B.rate_maps;
nLaps = size(Ds.lap_times, 1);
for nr = 1:numel(P.regions)
    cells = region_cells(Ds.region, P.regions{nr});
    nC = sum(cells);
    new_idx = zeros(numel(cells), 1); new_idx(cells) = 1:nC;
    % RUN templates (all laps, running spikes)
    in_region = cells(spk.cell);
    M = zeros(nC, nC, 2);
    for k = 1:nLaps
        s = spk.running & in_region & spk.lap == k;
        tr = Ds.lap_track(k);
        M(:, :, tr) = M(:, :, tr) + Lseq_pair_counts(spk.t(s), new_idx(spk.cell(s)), nC);
    end
    T1 = Lseq_template(M(:, :, 1)); T2 = Lseq_template(M(:, :, 2));
    W1 = T1.W; W2 = T2.W;
    b12 = sum(W1 .* W2, 'all') / sum(W2 .^ 2, 'all');
    b21 = sum(W1 .* W2, 'all') / sum(W1 .^ 2, 'all');
    Ws = {W1, W2, W1 - b12 * W2, W2 - b21 * W1};
    U = triu(true(nC), 1);
    R = struct('name', P.regions{nr}, 'cell_id', Ds.cell_id(cells), 'n_cells', nC);
    R.template = struct('W', {Ws}, 'b12', b12, 'b21', b21, 'r_T1_T2', corr(W1(U), W2(U)), ...
        'r_specific', corr(Ws{3}(U), Ws{4}(U)), 'N_T1', T1.N, 'N_T2', T2.N);
    % sleep spikes of the region
    in = cells(cell_idx);
    ts = t(in); cs = new_idx(cell_idx(in));
    for w = 1:numel(P.windows)
        X = struct('n_spikes', zeros(nEv, 1), 'n_cells', zeros(nEv, 1), 'scored', false(nEv, 1), ...
            'L', NaN(nEv, 4), 'z', NaN(nEv, 4), 'p_fwd', NaN(nEv, 4), 'p_rev', NaN(nEv, 4));
        lo = first_ge(ts, win{w}(:, 1)); hi = first_ge(ts, win{w}(:, 2)) - 1;
        for e = 1:nEv
            c = cs(lo(e):hi(e));
            X.n_spikes(e) = numel(c); X.n_cells(e) = numel(unique(c));
            if X.n_cells(e) < P.min_cells, continue; end
            rs = RandStream('mrg32k3a', 'Seed', P.seed + 1000 * nsession + 10 * nr + w);
            rs.Substream = Es.ripple_info_row(e);
            Sc = score_event(ts(lo(e):hi(e)), c, Ws, P.nShuffle, rs);
            X.scored(e) = true;
            X.L(e, :) = Sc.L; X.z(e, :) = Sc.z; X.p_fwd(e, :) = Sc.p_fwd; X.p_rev(e, :) = Sc.p_rev;
        end
        X.p = min(1, 2 * min(X.p_fwd, X.p_rev));          % two-sided
        R.(P.windows{w}) = X;
        ok = X.scored;
        fprintf('  %-3s (%3d cells, r(W1,W2) %.2f) %-4s: %4d/%4d scored | sig (2-sided) T1 %.3f T2 %.3f | specific T1 %.3f T2 %.3f either %.3f\n', ...
            P.regions{nr}, nC, R.template.r_T1_T2, P.windows{w}, sum(ok), nEv, mean(X.p(ok, 1) < 0.05), mean(X.p(ok, 2) < 0.05), ...
            mean(X.p(ok, 3) < 0.05), mean(X.p(ok, 4) < 0.05), mean(min(X.p(ok, 3:4), [], 2) < 0.05));
    end
    out.region(nr) = R;
end
folder = fileparts(out_file);
if ~exist(folder, 'dir'), mkdir(folder); end
save(out_file, 'out', '-v7.3');
fprintf('  saved (%.0f s)\n', toc);
try
    plot_Lseq_replay_examples('sessions', nsession);   % 4-row rasters (HPC/V1 x T1/T2 order)
catch err
    warning('Session %d example plot failed: %s', nsession, err.message);
end
close all
end

function S = score_event(t, cells, Ws, nShuffle, rs)
% As Lseq_score (k = 2 likelihood, forward - reverse, per spike-pair term, spike pairs
% at exactly the same time excluded), with the shuffle sums vectorised over the
% non-zero spike-pair counts of the event. t is sorted (spikes of the window in order).
t = t(:); cells = cells(:);
[u, ~, ci] = unique(cells);
m = numel(u);
[~, ~, g] = unique(t);                          % time groups (ties share a group)
G = full(sparse(g, ci, 1, max(g), m));
before = cumsum(G, 1) - G;                      % spikes of each cell strictly earlier
step = diff(t) > 0;                             % consecutive spikes not tied
A = before' * G + accumarray([ci([step; false]) ci([false; step])], 1, [m m]);
A(1:m + 1:end) = 0;
A = A / 2;
nz = find(A); [ia, ib] = ind2sub([m m], nz); a = A(nz);
n_terms = sum(a);
[~, perm] = sort(rand(rs, nShuffle, m), 2);
nW = numel(Ws);
S.L = NaN(1, nW); S.z = NaN(1, nW); S.p_fwd = NaN(1, nW); S.p_rev = NaN(1, nW);
if n_terms == 0, return; end
idx = (perm(:, ib) - 1) * m + perm(:, ia);           % nShuffle x nnz
for w = 1:nW
    Wu = Ws{w}(u, u);
    S.L(w) = (Wu(nz)' * a) / n_terms;
    Ls = reshape(Wu(idx), nShuffle, []) * a / n_terms;
    S.z(w) = (S.L(w) - mean(Ls)) / std(Ls);
    S.p_fwd(w) = (1 + sum(Ls >= S.L(w))) / (1 + nShuffle);
    S.p_rev(w) = (1 + sum(Ls <= S.L(w))) / (1 + nShuffle);
end
end

function idx = first_ge(t, x)
idx = zeros(size(x));
n = numel(t);
for k = 1:numel(x)
    lo = 1; hi = n + 1;
    while lo < hi
        mid = floor((lo + hi) / 2);
        if t(mid) < x(k), lo = mid + 1; else, hi = mid; end
    end
    idx(k) = lo;
end
end

function cells = region_cells(region, name)
switch name
    case 'HPC', cells = contains(region, 'HPC');
    case 'V1', cells = contains(region, 'V1');
end
end

function summarise_all(P)
files = dir(fullfile(P.out_dir, 'M*', 'session_*.mat'));
if isempty(files), return; end
T = table();
for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name), 'out'); out = S.out;
    B = load(fullfile(P.bayes_dir, out.subject, sprintf('session_%s.mat', out.session)), 'out'); B = B.out;
    nEv = numel(out.events.onset);
    Tk = table(repmat(string(out.subject), nEv, 1), repmat(out.session_index, nEv, 1), out.events.ripple_info_row, ...
        'VariableNames', {'subject', 'session_index', 'ripple_info_row'});
    for nr = 1:2
        for w = 1:2
            X = out.region(nr).(P.windows{w});
            pre = sprintf('%s_%s_', P.regions{nr}, P.windows{w});
            Tk.([pre 'scored']) = X.scored; Tk.([pre 'n_cells']) = X.n_cells;
            for q = 1:4
                Tk.([pre 'z_' P.templates{q}]) = X.z(:, q); Tk.([pre 'p_' P.templates{q}]) = X.p(:, q);
            end
            [tr, bi] = assign_track(X.z(:, 1:2), X.p(:, 1:2), X.scored);
            Tk.([pre 'track']) = tr; Tk.([pre 'bias']) = bi;
        end
        Tk.([P.regions{nr} '_bayes_logodds_post']) = B.region(nr).post.logodds;
        Tk.([P.regions{nr} '_bayes_p_post']) = min(B.region(nr).post.p, [], 2);
    end
    T = [T; Tk]; %#ok<AGROW>
end
save(fullfile(P.out_dir, 'Lseq_replay_events_all.mat'), 'T', '-v7.3');
writetable(T, fullfile(P.out_dir, 'Lseq_replay_events_all.csv'));

% Per session measures
sessions = unique(T.session_index); nS = numel(sessions);
names = {'HPC_rate_pre', 'HPC_rate_post', 'V1_rate_pre', 'V1_rate_post', 'HPC_specrate_pre', 'HPC_specrate_post', ...
    'HPC_absz_pre', 'HPC_absz_post', 'V1_absz_pre', 'V1_absz_post', ...
    'V1z_same_pre', 'V1z_same_post', 'V1z_other_pre', 'V1z_other_post', ...
    'V1sig_same_pre', 'V1sig_same_post', 'V1sig_other_pre', 'V1sig_other_post', ...
    'n_HPC_sig', 'bayes_agree', 'bayes_sig_overlap'};
for f = names, M.(f{1}) = NaN(nS, 1); end
subj = strings(nS, 1);
for k = 1:nS
    Ts = T(T.session_index == sessions(k), :); subj(k) = Ts.subject(1);
    for reg = ["HPC", "V1"]
        for w = ["pre", "post"]
            ok = Ts.(reg + "_" + w + "_scored");
            zs = [Ts.(reg + "_" + w + "_z_T1") Ts.(reg + "_" + w + "_z_T2")];
            M.(reg + "_rate_" + w)(k) = mean(Ts.(reg + "_" + w + "_track")(ok) > 0);
            M.(reg + "_absz_" + w)(k) = mean(max(abs(zs(ok, :)), [], 2));
        end
    end
    for w = ["pre", "post"]   % for reference: overlap-removed track-specific templates
        ok = Ts.("HPC_" + w + "_scored");
        ps = [Ts.("HPC_" + w + "_p_T1_specific") Ts.("HPC_" + w + "_p_T2_specific")];
        M.("HPC_specrate_" + w)(k) = mean(min(ps(ok, :), [], 2) < 0.05);
    end
    % HPC replay content: assigned track (post window)
    htr = Ts.HPC_post_track;
    sig = Ts.HPC_post_scored & htr > 0;
    htr(htr == 0) = 1;   % placeholder for unassigned events (excluded by sig below)
    M.n_HPC_sig(k) = sum(sig);
    n = height(Ts);
    for w = ["pre", "post"]
        z = abs([Ts.("V1_" + w + "_z_T1") Ts.("V1_" + w + "_z_T2")]);
        pv = [Ts.("V1_" + w + "_p_T1") Ts.("V1_" + w + "_p_T2")];
        ok = sig & Ts.("V1_" + w + "_scored");
        same = sub2ind([n 2], (1:n)', htr); other = sub2ind([n 2], (1:n)', 3 - htr);
        M.("V1z_same_" + w)(k) = mean(z(same(ok))); M.("V1z_other_" + w)(k) = mean(z(other(ok)));
        M.("V1sig_same_" + w)(k) = mean(pv(same(ok)) < 0.05); M.("V1sig_other_" + w)(k) = mean(pv(other(ok)) < 0.05);
    end
    lo = Ts.HPC_bayes_logodds_post;
    ok = sig & isfinite(lo);
    M.bayes_agree(k) = mean((lo(ok) > 0) == (htr(ok) == 1));             % Bayesian log odds favours the Lseq track
    M.bayes_sig_overlap(k) = mean(Ts.HPC_bayes_p_post(sig) < 0.05);       % also Bayesian-significant
end
Ms = struct2table(M); Ms.subject = subj; Ms.session_index = sessions;
writetable(Ms, fullfile(P.out_dir, 'Lseq_replay_session_summary.csv'));
plot_summary(Ms, P);
end

function [track, bias] = assign_track(z, p, scored)
% Track assignment from the raw T1/T2 templates: two-sided p < 0.05 on a track; if both
% tracks are significant, the larger |z| wins; 0 = unassigned. Bias = |z_T1| - |z_T2|
% (positive = Track 1).
az = abs(z);
sig = p < 0.05 & scored;
track = zeros(size(z, 1), 1);
track(sig(:, 1) & (~sig(:, 2) | az(:, 1) >= az(:, 2))) = 1;
track(sig(:, 2) & (~sig(:, 1) | az(:, 2) > az(:, 1))) = 2;
bias = az(:, 1) - az(:, 2);
bias(~scored) = NaN;
end

function plot_summary(Ms, P)
animals = unique(Ms.subject);
cmap = lines(numel(animals));
[~, ia] = ismember(Ms.subject, animals);
panels = {{'HPC_rate_pre', 'HPC_rate_post'}, {'HPC_specrate_pre', 'HPC_specrate_post'}, {'V1_rate_pre', 'V1_rate_post'}, ...
    {'HPC_absz_pre', 'HPC_absz_post'}, {'V1_absz_pre', 'V1_absz_post'}, ...
    {'V1z_same_pre', 'V1z_same_post', 'V1z_other_pre', 'V1z_other_post'}, ...
    {'V1sig_same_pre', 'V1sig_same_post', 'V1sig_other_pre', 'V1sig_other_post'}, {'bayes_agree', 'bayes_sig_overlap'}};
labels = {{'pre', 'post'}, {'pre', 'post'}, {'pre', 'post'}, {'pre', 'post'}, {'pre', 'post'}, ...
    {'same pre', 'same post', 'other pre', 'other post'}, {'same pre', 'same post', 'other pre', 'other post'}, ...
    {'track agrees', 'Bayes sig'}};
ylabs = {'HPC sig (either track)', 'HPC sig, overlap-removed templates (ref.)', 'V1 sig (either track)', ...
    'HPC mean max |z|', 'V1 mean max |z|', 'V1 |z|, HPC-sig events', ...
    'V1 sig, HPC-sig events', 'HPC-sig events: Bayesian agreement'};
chance = [1 - 0.95^2, 1 - 0.95^2, 1 - 0.95^2, NaN, NaN, NaN, 0.05, 0.5];
fig = figure('Name', 'Lseq replay summary', 'Position', [30 30 1900 430], 'Color', 'w');
tl = tiledlayout(1, numel(panels), 'TileSpacing', 'compact');
for k = 1:numel(panels)
    nexttile; hold on;
    Y = Ms{:, panels{k}};
    if size(Y, 2) > 1 && k < 8
        for c = 1:2:size(Y, 2)
            for s = 1:size(Y, 1), plot([c c + 1], Y(s, [c c + 1]), '-', 'Color', [cmap(ia(s), :) 0.5]); end
        end
    end
    for c = 1:size(Y, 2)
        scatter(c * ones(size(Y, 1), 1) + 0.05 * randn(size(Y, 1), 1), Y(:, c), 16, cmap(ia, :), 'filled');
        plot(c + [-0.25 0.25], mean(Y(:, c), 'omitnan') * [1 1], 'k-', 'LineWidth', 2);
    end
    if isfinite(chance(k)), yline(chance(k), 'k--'); end
    if k < 8
        ttl = sprintf('post vs pre p = %.2g', signrank(Y(:, 2), Y(:, 1)));
        if size(Y, 2) == 4
            ttl = {sprintf('same post vs pre p = %.2g', signrank(Y(:, 2), Y(:, 1))), ...
                sprintf('(same-other) post-pre p = %.2g', signrank((Y(:, 2) - Y(:, 1)) - (Y(:, 4) - Y(:, 3))))};
        end
    else
        ttl = sprintf('agree vs 0.5 p = %.2g', signrank(Y(:, 1) - 0.5));
    end
    title(ttl, 'FontWeight', 'normal', 'FontSize', 8);
    set(gca, 'XTick', 1:size(Y, 2), 'XTickLabel', labels{k}, 'TickDir', 'out', 'XTickLabelRotation', 30);
    xlim([0.5 size(Y, 2) + 0.5]); ylabel(ylabs{k}, 'FontSize', 8);
end
title(tl, sprintf(['POST ripples, L_{seq} (cell-ID shuffle %d, two-sided; track = raw template with p < 0.05, ' ...
    'larger |z| if both). Same window as Bayesian. Dots = sessions (colour = animal); dashed = chance'], P.nShuffle), 'FontSize', 9);
save_figure(fig, fullfile(P.out_dir, 'Lseq_replay_summary'));
close(fig);
end

function save_figure(fig, save_path)
for attempt = 1:10
    try
        savefig(fig, [save_path '.fig']);
        exportgraphics(fig, [save_path '.png'], 'Resolution', 150);
        exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
        return
    catch err
        if attempt == 10, rethrow(err); end
        pause(3);
    end
end
end

function folder = get_revision_folder()
if exist('D:\corticohippocampal_replay', 'dir') > 0
    folder = 'D:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    folder = 'P:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
end
end
