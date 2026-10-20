function run_Bayesian_replay_decoding_POST(varargin)
% Bayesian replay decoding of POST sleep ripples (non-duplicated SWS ripples of
% ripple_info.mat) with a multinomial naive Bayes decoder (bayes_decode_two_track.m):
% 20 ms bins from -0.5 to 0.5 s around ripple onset, 10 cm position bins of both tracks
% (28 joint states, posterior of each bin sums to 1 over both tracks), RUN1 rate maps
% from all laps (running samples > 5 cm/s, unsmoothed, 0.1 Hz floor), spike counts
% unsmoothed. Regions: HPC (HPC_L + HPC_R) and V1 (V1_L + V1_R), all spatial cells
% (RUN1_decoding_input.mat). Bins with no spikes have a flat posterior and are set to
% NaN / left out of the sequence scores.
% Per event and region:
%   log odds per bin and per window = log(sum P(Track 1) / sum P(Track 2))
%   scoring windows: 'post' = onset to max(onset + 0.1 s, offset) (whole 20 ms bins,
%   at most 0.5 s), 'pre' = the same number of bins immediately before onset
%   per window and track (posterior restricted to the track, renormalised per bin,
%   empty bins left out): weighted correlation (position x time) and line fit
%   (Davidson et al. 2009: mean posterior within +/- 1 position bin of the best line)
%   significance of the weighted correlation: cell-identity shuffle (rate maps
%   randomly reassigned among all decoder cells, re-decoded), nShuffle per event,
%   p = (1 + #(|r_shuffle| >= |r|)) / (1 + nShuffle); significant if p < 0.05;
%   z = (|r| - mean |r_shuffle|) / std |r_shuffle|. Scores need >= min_cells active
%   cells and >= min_bins non-empty bins in the window. Shuffles are seeded per
%   session, region, window and event (RandStream substream = ripple_info row), so
%   any event is reproducible on its own.
% Options: 'sessions' (1:22), 'time_bin' (0.02), 'window' ([-0.5 0.5]), 'min_post'
% (0.1 s), 'nShuffle' (1000), 'min_cells' (5), 'min_bins' (3), 'line_band' (1 bin),
% 'rate_floor' (0.1), 'speed_threshold' (5), 'seed' (2026), 'n_examples' (8),
% 'cell_filter' (decoder cells within the spatial-cell summary: 'all' (default),
%   'rate_lt10' = mean RUN1 firing rate < 10 Hz, 'not_flagged' = interneurons_all false),
% 'out_dir' (nature_neuroscience_revision_analysis\Bayesian_replay_decoding, with
%   _<cell_filter> appended unless 'all'), 'overwrite', 'summary_only',
% 'replot' (redraw example figures of sessions already decoded), 'refit_lines' (recompute
% the line fit of sessions already decoded from their stored posteriors; other scores kept).
% 'shuffle_types' (default {}): additional null distributions for the sessions already
%   decoded, any of {'cell_id', 'spike', 'field', 'time'}; when given, only this post-hoc
%   mode runs (no decoding, no summarise_all). Same score (|weighted correlation| per track),
%   nShuffle per event, scored events only, p = (1 + #(|r_sh| >= |r|)) / (1 + nShuffle):
%     cell_id: the decoding shuffle above, recomputed with its own seeds (checked against
%              the stored p) to get its empirical chance
%     spike:   spike train circular shift, each active cell's counts independently shifted
%              in time within the scoring window, re-decoded (empty bins left out)
%     field:   place field circular shift, nShuffle decoders per session and region with
%              every decoder cell's rate map shifted independently on T1 and T2 (0 to
%              nPos - 1 bins), multinomial log likelihoods recomputed, shared by all events
%     time:    time bin permutation of the real posterior, all bins of the window permuted
%              (empty bins move and stay out of the score); exhaustive (exact p) when there
%              are <= nShuffle arrangements of the non-empty bins
%   plus event log odds z-scored against the per-bin log odds of all decoded bins (-0.5 to
%   0.5 s of every ripple, session and region; finite bins only, +/-Inf capped at the
%   0.5/99.5 prctile). Empirical chance per event = fraction of null replicates whose own p
%   (within the same null) is < 0.05, per track and on either track.
%   Saved to out_dir\<animal>\shuffles_session_<date>.mat (struct sh: sh.<region>.<window>.
%   p_<type>, z_<type>, mean_<type>, sd_<type>, chance_<type> (nEv x 2), chance_either_<type>,
%   logodds_z (nEv x 1); rows = out.events; only the types run are replaced), followed by
%   the shuffle summary.
% 'shuffle_summary' (false): only redraw the shuffle summary, out_dir\shuffle_comparison:
%   event table, one figure per significance criterion ('significance_criteria', default
%   {{'cell_id'}, {'spike', 'field'}, {'field', 'time'}}; see replay_significance.m) with
%   replay rates pre vs post and proportion of significant events vs mean log odds
%   difference (T1 - T2 events) across 'alphas' (95% bootstrap CI, 'nBoot' resamples of
%   events), and one shuffle agreement figure per region.
% Output: out_dir\<animal>\session_<date>.mat (+ example event figure), and in out_dir
% replay_events_all.mat (event table), replay_summary figures and CSV, and across-session
% example figures replay_best_examples1-2 / replay_nonsignificant_examples1-2 (10 events each).

p = inputParser;
addParameter(p, 'sessions', 1:22, @isnumeric);
addParameter(p, 'time_bin', 0.02, @isnumeric);
addParameter(p, 'window', [-0.5 0.5], @isnumeric);
addParameter(p, 'min_post', 0.1, @isnumeric);
addParameter(p, 'nShuffle', 1000, @isnumeric);
addParameter(p, 'min_cells', 5, @isnumeric);
addParameter(p, 'min_bins', 3, @isnumeric);
addParameter(p, 'line_band', 1, @isnumeric);
addParameter(p, 'rate_floor', 0.1, @isnumeric);
addParameter(p, 'speed_threshold', 5, @isnumeric);
addParameter(p, 'seed', 2026, @isnumeric);
addParameter(p, 'n_examples', 8, @isnumeric);
addParameter(p, 'out_dir', '', @ischar);
addParameter(p, 'overwrite', false, @islogical);
addParameter(p, 'summary_only', false, @islogical);
addParameter(p, 'replot', false, @islogical);
addParameter(p, 'refit_lines', false, @islogical);
addParameter(p, 'cell_filter', 'all', @(x) ismember(x, {'all', 'rate_lt10', 'not_flagged'}));
addParameter(p, 'shuffle_types', {}, @(x) iscell(x) && all(ismember(x, {'cell_id', 'spike', 'field', 'time'})));
addParameter(p, 'shuffle_summary', false, @islogical);
addParameter(p, 'significance_criteria', {{'cell_id'}, {'spike', 'field'}, {'field', 'time'}}, @iscell);
addParameter(p, 'alphas', [0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001], @isnumeric);
addParameter(p, 'nBoot', 1000, @isnumeric);
parse(p, varargin{:});
P = p.Results;
P.pos_edges = 0:10:140;
P.nPos = numel(P.pos_edges) - 1;
P.nFolds = 10;                                     % unused here (RUN1_session_behaviour)
P.edges = P.window(1):P.time_bin:P.window(2) + 1e-9;
P.nBins = numel(P.edges) - 1;
P.bin_centres = P.edges(1:end-1) + P.time_bin / 2;
P.onset_bin = find(P.edges(1:end-1) >= -1e-9, 1);  % first bin starting at onset
P.regions = {'HPC', 'V1'};
P.windows = {'post', 'pre'};
P.ripple_colours = [0.9 0.5 0; 0.5 0.2 0.8];          % example-plot shading: probe 1, probe 2

revision_folder = get_revision_folder();
if isempty(P.out_dir)
    P.out_dir = fullfile(revision_folder, 'Bayesian_replay_decoding');
    if ~strcmp(P.cell_filter, 'all'), P.out_dir = [P.out_dir '_' P.cell_filter]; end
end
if ~exist(P.out_dir, 'dir'), mkdir(P.out_dir); end

if P.refit_lines
    refit_lines_all(P);
end
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
if ~isempty(P.shuffle_types) || P.shuffle_summary
    if ~isempty(P.shuffle_types), extra_shuffles_all(P); end
    summarise_shuffles(P);
    return
end
if ~P.summary_only
    S = load(fullfile(revision_folder, 'RUN1_decoding_input.mat'), 'D'); D = S.D;
    E = ripple_events(revision_folder);
    for nsession = P.sessions
        in = E.session == nsession;
        decode_session(D(nsession), structfun(@(x) x(in), E, 'UniformOutput', false), session_experiment_info(nsession), nsession, P);
    end
end
summarise_all(P);
end

function E = ripple_events(revision_folder)
% Times of the non-duplicated POST SWS ripples (ripple_info.event_id = row of the full
% SWS PSTH: probe 1 SWS ripples of sessions 1-22, then probe 2)
S = load(fullfile(revision_folder, 'ripple_info.mat'), 'ripple_info'); info = S.ripple_info;
S = load(fullfile(revision_folder, 'ripples_all_POST.mat'), 'ripples_all'); RA = S.ripples_all;
sws = {find(RA(1).SWS_index == 1), find(RA(2).SWS_index == 1)};
eid = info.event_id(:);
E.ripple_info_row = (1:numel(eid))';
E.probe = 1 + (eid > numel(sws{1}));
E.onset = NaN(size(eid)); E.offset = E.onset; E.peak = E.onset; E.session = E.onset;
for pr = 1:2
    in = E.probe == pr;
    idx = sws{pr}(eid(in) - (pr - 1) * numel(sws{1}));
    E.onset(in) = RA(pr).onset(idx); E.offset(in) = RA(pr).offset(idx);
    E.peak(in) = RA(pr).peaktimes(idx); E.session(in) = RA(pr).session_count(idx);
end
assert(isequal(E.session, info.session_id(:)), 'ripple_info session_id does not match ripples_all');
end

function decode_session(Ds, Es, info, nsession, P)
out_file = fullfile(P.out_dir, Ds.subject, sprintf('session_%s.mat', Ds.session));
if exist(out_file, 'file') && ~P.overwrite
    fprintf('Session %d exists, skipping\n', nsession);
    if P.replot
        S = load(out_file, 'out');
        [t, cell_idx] = load_sleep_spikes(info, Ds.cell_id);
        try
            plot_session_examples(S.out, t, cell_idx, P, strrep(out_file, '.mat', '_examples'));
        catch err
            warning('Session %d example plot failed: %s', nsession, err.message);
        end
        close all
    end
    return
end
fprintf('\n=== Session %d: %s %s, %d ripples ===\n', nsession, Ds.subject, Ds.session, numel(Es.onset));
tic
[~, order] = sort(Es.onset);
Es = structfun(@(x) x(order), Es, 'UniformOutput', false);
nEv = numel(Es.onset);
Es.duration = Es.offset - Es.onset;
Es.n_post_bins = min(P.nBins - P.onset_bin + 1, ceil(round(max(P.min_post, Es.duration) / P.time_bin, 6)));

% RUN1 rate maps, all laps (running samples)
[beh, spk] = RUN1_session_behaviour(Ds, nsession, P);
nC = numel(Ds.cell_id);
r = beh.running;
occ = accumarray([beh.track(r) beh.pos_bin(r)], beh.dt(r), [2 P.nPos]);
r = spk.running;
cnt = accumarray([spk.cell(r) spk.pos_bin(r) spk.track(r)], 1, [nC P.nPos 2]);
maps = [cnt(:, :, 1) ./ occ(1, :) cnt(:, :, 2) ./ occ(2, :)];
maps(~isfinite(maps)) = 0;

% Decoder cells: spatial cells of the summary, optionally filtered
run_rate = accumarray(spk.cell, 1, [nC 1]) / (Ds.tvec(end) - Ds.tvec(1));   % mean RUN1 rate (Hz)
switch P.cell_filter
    case 'all', use = true(nC, 1);
    case 'rate_lt10', use = run_rate < 10;
    case 'not_flagged', use = ~Ds.interneuron(:);
end

% POST sleep spikes of the same cells
[t, cell_idx] = load_sleep_spikes(info, Ds.cell_id);

out = struct();
out.subject = Ds.subject; out.session = Ds.session; out.session_index = nsession; out.params = P;
out.events = Es;
out.rate_maps = maps; out.cell_id = Ds.cell_id; out.region_labels = Ds.region;
out.run_rate = run_rate; out.decoder_cell = use;
for nr = 1:numel(P.regions)
    cells = region_cells(Ds.region, P.regions{nr}) & use;
    nCr = sum(cells);
    new_idx = zeros(nC, 1); new_idx(cells) = 1:nCr;
    in = cells(cell_idx);
    tr = t(in); cr = new_idx(cell_idx(in));
    mr = max(maps(cells, :), P.rate_floor);
    logp = log(mr ./ sum(mr, 1));                     % nCells x 28 (multinomial)
    R = struct('name', P.regions{nr}, 'cell_id', Ds.cell_id(cells), 'n_cells', nCr);
    R.posterior = NaN(2 * P.nPos, P.nBins, nEv, 'single');
    R.n_spikes = zeros(nEv, P.nBins, 'single');
    R.n_active_bin = zeros(nEv, P.nBins, 'single');
    R.logodds_bin = NaN(nEv, P.nBins, 'single');
    for w = 1:numel(P.windows)
        R.(P.windows{w}) = empty_scores(nEv);
    end
    lo = first_ge(tr, Es.onset + P.window(1));
    hi = first_ge(tr, Es.onset + P.window(2)) - 1;
    for e = 1:nEv
        s = (lo(e):hi(e))';
        b = floor((tr(s) - (Es.onset(e) + P.window(1))) / P.time_bin) + 1;
        ok = b >= 1 & b <= P.nBins;
        counts = accumarray([cr(s(ok)) b(ok)], 1, [nCr P.nBins]);
        nonempty = sum(counts, 1) > 0;
        post = bayes_decode_two_track(counts, mr, P.time_bin, 'multinomial', P.rate_floor);
        post(:, ~nonempty) = NaN;
        R.posterior(:, :, e) = post;
        R.n_spikes(e, :) = sum(counts, 1);
        R.n_active_bin(e, :) = sum(counts > 0, 1);
        R.logodds_bin(e, :) = log(sum(post(1:P.nPos, :), 1) ./ sum(post(P.nPos + 1:end, :), 1));
        nb = Es.n_post_bins(e);
        win_bins = {P.onset_bin + (0:nb - 1), P.onset_bin - nb + (0:nb - 1)};
        for w = 1:numel(P.windows)
            rs = RandStream('mrg32k3a', 'Seed', P.seed + 1000 * nsession + 10 * nr + w);
            rs.Substream = Es.ripple_info_row(e);
            sc_w = score_window(counts(:, win_bins{w}), logp, P, rs);
            R.(P.windows{w}) = fill_scores(R.(P.windows{w}), e, sc_w);
        end
    end
    out.region(nr) = R;
    for w = 1:numel(P.windows)
        X = R.(P.windows{w});
        ok = X.scored;
        fprintf('  %-3s (%3d cells) %-4s: %4d/%4d scored | sig T1 %.3f, T2 %.3f, either %.3f | mean |wc| %.2f\n', ...
            P.regions{nr}, nCr, P.windows{w}, sum(ok), nEv, mean(X.p(ok, 1) < 0.05), mean(X.p(ok, 2) < 0.05), ...
            mean(min(X.p(ok, :), [], 2) < 0.05), mean(max(abs(X.wc(ok, :)), [], 2)));
    end
end
folder = fileparts(out_file);
if ~exist(folder, 'dir'), mkdir(folder); end
save(out_file, 'out', '-v7.3');
fprintf('  saved (%.0f s)\n', toc);
try
    plot_session_examples(out, t, cell_idx, P, strrep(out_file, '.mat', '_examples'));
catch err
    warning('Session %d example plot failed: %s', nsession, err.message);
end
close all
end

function X = empty_scores(nEv)
X.n_active = zeros(nEv, 1); X.n_nonempty = zeros(nEv, 1); X.n_spikes = zeros(nEv, 1);
X.scored = false(nEv, 1); X.logodds = NaN(nEv, 1);
X.wc = NaN(nEv, 2); X.p = NaN(nEv, 2); X.z = NaN(nEv, 2);
X.shuffle_mean = NaN(nEv, 2); X.shuffle_sd = NaN(nEv, 2);
X.line_score = NaN(nEv, 2); X.line_slope = NaN(nEv, 2); X.line_y = NaN(nEv, 2, 2);
end

function X = fill_scores(X, e, s)
for f = fieldnames(s)'
    if strcmp(f{1}, 'line_y'), X.line_y(e, :, :) = reshape(s.line_y, 1, 2, 2); else, X.(f{1})(e, :) = s.(f{1}); end
end
end

function s = score_window(cw, logp, P, rs)
% Log odds, weighted correlation (+ cell-ID shuffle) and line fit of one window
nPos = P.nPos;
nbw = size(cw, 2);
nonempty = sum(cw, 1) > 0;
tb = find(nonempty);                 % bin positions within the window (gaps = empty bins)
cw = cw(:, nonempty);
A = find(sum(cw, 2) > 0);
s.n_active = numel(A); s.n_nonempty = numel(tb); s.n_spikes = sum(cw(:));
s.scored = false; s.logodds = NaN;
s.wc = NaN(1, 2); s.p = NaN(1, 2); s.z = NaN(1, 2); s.shuffle_mean = NaN(1, 2); s.shuffle_sd = NaN(1, 2);
s.line_score = NaN(1, 2); s.line_slope = NaN(1, 2); s.line_y = NaN(2, 2);
if isempty(A), return; end
L = logp(A, :)' * cw(A, :);
post = softmax1(L);
s.logodds = log(sum(post(1:nPos, :), 'all') / sum(post(nPos + 1:end, :), 'all'));
if s.n_active < P.min_cells || s.n_nonempty < P.min_bins, return; end
s.scored = true;
% cell-identity shuffles: random rate maps (out of all decoder cells) for the active cells
nS = P.nShuffle; nA = numel(A);
[~, perm] = sort(rand(rs, nS, size(logp, 1)), 2);
idx = perm(:, 1:nA)';                                      % nA x nS
G = reshape(logp(idx(:), :), nA, nS, 2 * nPos);
G = reshape(permute(G, [2 3 1]), nS * 2 * nPos, nA);
Lsh = reshape(G * cw(A, :), nS, 2 * nPos, numel(tb));
Lsh = exp(Lsh - max(Lsh, [], 2));
Lsh = Lsh ./ sum(Lsh, 2);
for k = 1:2
    rows = (k - 1) * nPos + (1:nPos);
    Pk = post(rows, :) ./ sum(post(rows, :), 1);
    s.wc(k) = wcorr_batch(reshape(Pk, [1 nPos numel(tb)]), tb);
    if ~isfinite(s.wc(k)), continue; end
    Wk = Lsh(:, rows, :); Wk = Wk ./ sum(Wk, 2);
    r_sh = abs(wcorr_batch(Wk, tb));
    r_sh(~isfinite(r_sh)) = 0;
    s.p(k) = (1 + sum(r_sh >= abs(s.wc(k)))) / (1 + nS);
    s.shuffle_mean(k) = mean(r_sh); s.shuffle_sd(k) = std(r_sh);
    s.z(k) = (abs(s.wc(k)) - s.shuffle_mean(k)) / s.shuffle_sd(k);
    [s.line_score(k), s.line_slope(k), s.line_y(k, :)] = line_fit(Pk, tb, nbw, P);
end
end

function post = softmax1(L)
post = exp(L - max(L, [], 1));
post = post ./ sum(post, 1);
end

function r = wcorr_batch(W, t)
% Weighted correlation between position and time; W: nS x nPos x nT (weights), t: bin times
nX = size(W, 2);
x = reshape(1:nX, 1, nX, 1);
tt = reshape(t, 1, 1, []);
sw = sum(W, [2 3]);
mx = sum(W .* x, [2 3]) ./ sw;
mt = sum(W .* tt, [2 3]) ./ sw;
dx = x - mx; dt = tt - mt;
cxy = sum(W .* dx .* dt, [2 3]) ./ sw;
vx = sum(W .* dx.^2, [2 3]) ./ sw;
vt = sum(W .* dt.^2, [2 3]) ./ sw;
r = cxy ./ sqrt(vx .* vt);
end

function [score, slope, y] = line_fit(Pk, tb, nbw, P)
% Best line through the posterior (Davidson et al. 2009): lines from position y0 at the
% first window bin to y1 at the last; score = mean over non-empty bins of the posterior
% within +/- line_band bins of the line. y0 and y1 range over one track length beyond
% each end (bins 1-nPos+1 to 2nPos-1, step 0.5), so a fast sequence that runs off the
% track before the window ends (e.g. last bins empty) is not forced flat; parts of a
% line outside the track capture nothing.
if nbw < 2, score = NaN; slope = NaN; y = [NaN NaN]; return; end
nPos = size(Pk, 1);
yr = (2 - nPos):0.5:(2 * nPos - 1);
[Y0, Y1] = ndgrid(yr, yr);
Y0 = Y0(:); Y1 = Y1(:);
Yt = Y0 + (Y1 - Y0) .* ((tb - 1) / (nbw - 1));             % nLines x nT
x = reshape(1:nPos, 1, 1, nPos);
band = abs(x - Yt) <= P.line_band;                           % nLines x nT x nPos
S = sum(band .* reshape(Pk', 1, numel(tb), nPos), 3);
[score, k] = max(mean(S, 2));
y = [Y0(k) Y1(k)];
slope = (Y1(k) - Y0(k)) * (P.pos_edges(2) - P.pos_edges(1)) / ((nbw - 1) * P.time_bin);   % cm/s
end

function refit_lines_all(P)
% Recompute line_score, line_slope and line_y of every scored event (both regions,
% windows and tracks) from the stored posteriors with the current line_fit
files = dir(fullfile(P.out_dir, 'M*', 'session_*.mat'));
for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name);
    S = load(f, 'out'); out = S.out;
    if ~ismember(out.session_index, P.sessions), continue; end
    tic
    ob = out.params.onset_bin; nPos = out.params.nPos;
    changed = 0; n = 0;
    for r = 1:numel(out.region)
        for w = 1:numel(P.windows)
            X = out.region(r).(P.windows{w});
            for e = find(X.scored)'
                nb = out.events.n_post_bins(e);
                bins = ob + (w == 1) * (0:nb - 1) + (w == 2) * (-nb:-1);
                post = double(out.region(r).posterior(:, bins, e));
                ne = all(isfinite(post), 1); tb = find(ne);
                for tr = 1:2
                    rows = (tr - 1) * nPos + (1:nPos);
                    Pk = post(rows, ne) ./ sum(post(rows, ne), 1);
                    old = X.line_y(e, tr, :);
                    [X.line_score(e, tr), X.line_slope(e, tr), y] = line_fit(Pk, tb, nb, P);
                    X.line_y(e, tr, :) = reshape(y, 1, 1, 2);
                    changed = changed + any(old(:) ~= y(:)); n = n + 1;
                end
            end
            out.region(r).(P.windows{w}) = X;
        end
    end
    save(f, 'out', '-v7.3');
    fprintf('Session %d: line fit recomputed for %d event-tracks (%d changed, %.0f s)\n', out.session_index, n, changed, toc);
end
end

function idx = first_ge(t, x)
% Index of the first element of sorted t that is >= x (numel(t)+1 if none)
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
    otherwise, cells = region == string(name);
end
end

function [track, p, z] = replay_track(X, e)
% Replayed track of event e in window X: the track with the lower p (ties: larger z)
[p, track] = min(X.p(e, :));
if X.p(e, 1) == X.p(e, 2), [~, track] = max(X.z(e, :)); end
z = X.z(e, track);
end

function [t, cell_idx] = load_sleep_spikes(info, cell_id)
% POST sleep (SleepChronic) spikes of the given clusters: times and index into cell_id
is_post = strcmp(info.StimulusName, 'SleepChronic');
assert(sum(is_post) == 1, 'Expected one SleepChronic stimulus');
options = info.session(is_post).probe(1);
S = load(fullfile(options.ANALYSIS_DATAPATH, '..', 'session_clusters_Sleep.mat'), 'session_clusters');
sc = S.session_clusters;
keep = ismember(sc.cluster_id, cell_id);
t = vertcat(sc.spike_times{keep}); id = vertcat(sc.spike_id{keep});
[~, cell_idx] = ismember(id, cell_id);
[t, o] = sort(t); cell_idx = cell_idx(o);
end

function plot_session_examples(out, t, cell_idx, P, save_path)
% Most significant HPC replays (post window) of one session, one column per event
% (plot_example_column)
X = out.region(1).post;
cand = find(X.scored & min(X.p, [], 2) < 0.05);
if isempty(cand), return; end
zbest = max(X.z(cand, :), [], 2);
[~, o] = sort(zbest, 'descend');
ex = cand(o(1:min(P.n_examples, numel(o))));
nE = numel(ex);
[fig, tl] = example_figure(nE);
for k = 1:nE
    plot_example_column(tl, k, nE, out, ex(k), t, cell_idx, P, '');
end
title(tl, {sprintf('%s %s (session %d): most significant HPC replays', out.subject, out.session, ...
    out.session_index), example_legend(P)}, 'Interpreter', 'none', 'FontSize', 9);
save_figure(fig, save_path);
end

function [fig, tl] = example_figure(nE)
fig = figure('Name', 'examples', 'Position', [10 10 250 * nE 1700], 'Color', 'w');
tl = tiledlayout(8, nE, 'TileSpacing', 'compact', 'Padding', 'compact');
end

function s = example_legend(P)
s = sprintf(['Multinomial, %d ms, cell-ID shuffle. Joint = posterior over both tracks; renorm = one track ' ...
    'renormalised per bin (scored); red = best line; bold = p < 0.05; shaded = detected ripples (orange probe 1, ' ...
    'purple probe 2)'], round(P.time_bin * 1000));
end

function plot_example_column(tl, k, nE, out, e, t, cell_idx, P, prefix)
% One event in column k. Rows: HPC raster (decoder cells active in the window, ordered
% by place-field peak on the replayed track); HPC joint posterior (both tracks, as
% decoded); HPC Track 1 and Track 2 posteriors renormalised within the track (what the
% weighted correlation and line fit score) with that track's best line; the same three
% posterior rows for V1; per-bin log odds (T1/T2) of HPC and V1. Dashed: scoring window.
% Shaded: all detected ripples of both probes in the window. Titles: wc, p, line score
% and slope per track. prefix: text before the first title (e.g. animal and session).
H = out.region(1); V = out.region(2);
[track, p] = replay_track(H.post, e);
onset = out.events.onset(e); nb = out.events.n_post_bins(e);
hpc_cells = find(ismember(out.cell_id, H.cell_id));
% raster
nexttile(tl, k); hold on;
win = [-nb nb] * P.time_bin + [-0.05 0.05];
s = cell_idx > 0 & t >= onset + win(1) & t < onset + win(2) & ismember(cell_idx, hpc_cells);
active = unique(cell_idx(s));
[~, pk] = max(out.rate_maps(active, (track - 1) * P.nPos + (1:P.nPos)), [], 2);
[~, ord] = sort(pk); yrow = zeros(numel(out.cell_id), 1); yrow(active(ord)) = 1:numel(active);
scatter(t(s) - onset, yrow(cell_idx(s)), 12, 'k', '|', 'LineWidth', 1.2);
xline([0 nb * P.time_bin], 'r--'); xlim(win); ylim([0 numel(active) + 1]);
shade_ripples(out.session_index, onset, P);
ylabel(sprintf('HPC (T%d peak order)', track), 'FontSize', 7); set(gca, 'TickDir', 'out', 'FontSize', 7);
title(sprintf('%sev %d (row %d)\nHPC T%d p %.3f', prefix, e,out.events.ripple_info_row(e), track, p), ...
    'FontWeight', 'normal', 'FontSize', 7, 'Interpreter', 'none');
regs = {H, V};
for nr = 1:2
    Y = regs{nr}.post;
    nexttile(tl, (1 + 3 * (nr - 1)) * nE + k);
    plot_posterior(regs{nr}, e, nb, P);
    shade_ripples(out.session_index, onset, P);
    title(sprintf('%s joint | LO %.1f', regs{nr}.name, Y.logodds(e)), 'FontWeight', 'normal', 'FontSize', 7);
    for tr = 1:2
        nexttile(tl, (1 + 3 * (nr - 1) + tr) * nE + k);
        plot_track_posterior(regs{nr}, e, nb, P, tr);
        shade_ripples(out.session_index, onset, P);
        if Y.p(e, tr) < 0.05, fw = 'bold'; else, fw = 'normal'; end
        title(sprintf('%s T%d renorm: wc %.2f p %.3f\nline %.2f, %.0f cm/s', regs{nr}.name, tr, Y.wc(e, tr), Y.p(e, tr), ...
            Y.line_score(e, tr), Y.line_slope(e, tr)), 'FontWeight', fw, 'FontSize', 7);
    end
end
% per-bin log odds, HPC vs V1
nexttile(tl, 7 * nE + k); hold on;
plot(P.bin_centres, H.logodds_bin(e, :), 'k.-');
plot(P.bin_centres, V.logodds_bin(e, :), '.-', 'Color', [0.2 0.6 0.2]);
xline([0 nb * P.time_bin], 'r--'); yline(0, 'k:');
xlim(P.window); ylabel('log odds T1/T2', 'FontSize', 7);
shade_ripples(out.session_index, onset, P);
title(sprintf('window LO: HPC %.1f, V1 %.1f', H.post.logodds(e), V.post.logodds(e)), 'FontWeight', 'normal', 'FontSize', 7);
if k == 1, legend({'HPC', 'V1'}, 'Box', 'off', 'Location', 'best', 'FontSize', 7); end
set(gca, 'TickDir', 'out', 'FontSize', 7);
xlabel('Time from ripple onset (s)');
end

function shade_ripples(session_index, onset, P)
% Shade all detected POST ripples of the session (both probes, any state) that overlap
% the -0.5 to 0.5 s window of the event in the current axes (time relative to onset)
persistent RIP
if isempty(RIP)
    S = load(fullfile(get_revision_folder(), 'ripples_all_POST.mat'), 'ripples_all');
    for pr = 1:2
        RIP(pr).session = S.ripples_all(pr).session_count(:);
        RIP(pr).onset = S.ripples_all(pr).onset(:);
        RIP(pr).offset = S.ripples_all(pr).offset(:);
    end
end
yl = ylim;
for pr = 1:2
    in = RIP(pr).session == session_index & RIP(pr).offset > onset + P.window(1) & RIP(pr).onset < onset + P.window(2);
    on = RIP(pr).onset(in) - onset; off = RIP(pr).offset(in) - onset;
    for j = 1:numel(on)
        patch([on(j) off(j) off(j) on(j)], yl([1 1 2 2]), P.ripple_colours(pr, :), 'FaceAlpha', 0.2, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end
ylim(yl);
end

function plot_posterior(R, e, nb, P)
% Joint posterior over both tracks (T1 bottom, T2 top), -0.5 to 0.5 s
post = double(R.posterior(:, :, e));
imagesc(P.bin_centres, 1:2 * P.nPos, post); axis xy; hold on;
colormap(gca, flipud(gray)); clim([0 max([post(:); 0.1], [], 'omitnan')]);
yline(P.nPos + 0.5, 'b-', 'LineWidth', 1);
xline([0 nb * P.time_bin], 'r--');
set(gca, 'YTick', [P.nPos / 2 1.5 * P.nPos] + 0.5, 'YTickLabel', {'T1', 'T2'}, 'TickDir', 'out', 'FontSize', 7);
xlim(P.window);
end

function plot_track_posterior(R, e, nb, P, track)
% One track's posterior renormalised per bin (as scored), -0.5 to 0.5 s, with the best
% line of the post window
post = double(R.posterior((track - 1) * P.nPos + (1:P.nPos), :, e));
post = post ./ sum(post, 1);
imagesc(P.bin_centres, P.pos_edges(1:end-1) + 5, post); axis xy; hold on;
colormap(gca, flipud(gray)); clim([0 1]);
xline([0 nb * P.time_bin], 'r--');
y = squeeze(R.post.line_y(e, track, :));
if all(isfinite(y))
    plot([0.5 nb - 0.5] * P.time_bin, (y' - 1) * 10 + 5, 'r-', 'LineWidth', 1.5);
end
set(gca, 'YTick', [5 70 135], 'YTickLabel', {'0', '70', '140'}, 'TickDir', 'out', 'FontSize', 7);
ylabel(sprintf('T%d (cm)', track), 'FontSize', 7);
xlim(P.window);
end

function summarise_all(P)
% Event table across sessions, replay rates, V1 sequence pre vs post by HPC replay
% content, HPC-V1 log odds coherence, and the best example events
files = dir(fullfile(P.out_dir, 'M*', 'session_*.mat'));
if isempty(files), return; end
T = table(); ex = struct('score', {}, 'pmin', {}, 'file', {}, 'e', {});
for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name), 'out'); out = S.out;
    nEv = numel(out.events.onset);
    Tk = table(repmat(string(out.subject), nEv, 1), repmat(string(out.session), nEv, 1), repmat(out.session_index, nEv, 1), ...
        out.events.ripple_info_row, out.events.probe, out.events.onset, out.events.duration, out.events.n_post_bins, ...
        'VariableNames', {'subject', 'session', 'session_index', 'ripple_info_row', 'probe', 'onset', 'duration', 'n_post_bins'});
    for nr = 1:numel(out.region)
        for w = 1:numel(P.windows)
            X = out.region(nr).(P.windows{w});
            pre = sprintf('%s_%s_', out.region(nr).name, P.windows{w});
            Tk.([pre 'n_active']) = X.n_active; Tk.([pre 'scored']) = X.scored; Tk.([pre 'logodds']) = X.logodds;
            for tr = 1:2
                Tk.(sprintf('%swc_T%d', pre, tr)) = X.wc(:, tr); Tk.(sprintf('%sp_T%d', pre, tr)) = X.p(:, tr);
                Tk.(sprintf('%sz_T%d', pre, tr)) = X.z(:, tr); Tk.(sprintf('%sline_T%d', pre, tr)) = X.line_score(:, tr);
                Tk.(sprintf('%sslope_T%d', pre, tr)) = X.line_slope(:, tr);
            end
        end
    end
    T = [T; Tk]; %#ok<AGROW>
    X = out.region(1).post;
    cand = find(X.scored & X.n_active >= 8);
    for c = cand(:)'
        ex(end+1) = struct('score', max(X.z(c, :)), 'pmin', min(X.p(c, :)), ...
            'file', fullfile(files(k).folder, files(k).name), 'e', c); %#ok<AGROW>
    end
end
save(fullfile(P.out_dir, 'replay_events_all.mat'), 'T', '-v7.3');
writetable(T, fullfile(P.out_dir, 'replay_events_all.csv'));

% HPC replay content (post window): the more significant track of HPC-significant events
hp = [T.HPC_post_p_T1 T.HPC_post_p_T2]; hz = [T.HPC_post_z_T1 T.HPC_post_z_T2];
[hp_min, htrack] = min(hp, [], 2);
tie = hp(:, 1) == hp(:, 2); [~, zt] = max(hz, [], 2); htrack(tie) = zt(tie);
T.HPC_sig = T.HPC_post_scored & hp_min < 0.05;
T.HPC_track = htrack;
sessions = unique(T.session_index);
nS = numel(sessions);
M = struct();
names = {'HPC_rate_post', 'HPC_rate_pre', 'V1_rate_post', 'V1_rate_pre', ...
    'V1z_same_post', 'V1z_same_pre', 'V1z_other_post', 'V1z_other_pre', ...
    'V1sig_same_post', 'V1sig_same_pre', 'V1sig_other_post', 'V1sig_other_pre', 'r_logodds_post', 'n_HPC_sig', 'n_scored'};
for f = names, M.(f{1}) = NaN(nS, 1); end
subj = strings(nS, 1);
for k = 1:nS
    Ts = T(T.session_index == sessions(k), :);
    subj(k) = Ts.subject(1);
    for reg = ["HPC", "V1"]
        for w = ["post", "pre"]
            ok = Ts.(reg + "_" + w + "_scored");
            pm = min([Ts.(reg + "_" + w + "_p_T1") Ts.(reg + "_" + w + "_p_T2")], [], 2);
            M.(reg + "_rate_" + w)(k) = mean(pm(ok) < 0.05);
        end
    end
    M.n_scored(k) = sum(Ts.HPC_post_scored);
    sig = Ts.HPC_sig; M.n_HPC_sig(k) = sum(sig);
    for w = ["post", "pre"]
        z = [Ts.("V1_" + w + "_z_T1") Ts.("V1_" + w + "_z_T2")];
        pv = [Ts.("V1_" + w + "_p_T1") Ts.("V1_" + w + "_p_T2")];
        ok = sig & Ts.("V1_" + w + "_scored");
        same = sub2ind(size(z), (1:height(Ts))', Ts.HPC_track); other = sub2ind(size(z), (1:height(Ts))', 3 - Ts.HPC_track);
        M.("V1z_same_" + w)(k) = mean(z(same(ok))); M.("V1z_other_" + w)(k) = mean(z(other(ok)));
        M.("V1sig_same_" + w)(k) = mean(pv(same(ok)) < 0.05); M.("V1sig_other_" + w)(k) = mean(pv(other(ok)) < 0.05);
    end
    ok = isfinite(Ts.HPC_post_logodds) & isfinite(Ts.V1_post_logodds);
    M.r_logodds_post(k) = corr(Ts.HPC_post_logodds(ok), Ts.V1_post_logodds(ok), 'Type', 'Spearman');
end
Ms = struct2table(M); Ms.subject = subj; Ms.session_index = sessions;
writetable(Ms, fullfile(P.out_dir, 'replay_session_summary.csv'));
try
    plot_summary(Ms, P);
catch err
    warning('Summary figure failed: %s', err.message);
end
try
    plot_best_examples(ex, P);
catch err
    warning('Best-examples figure failed: %s', err.message);
end
end

function plot_summary(Ms, P)
animals = unique(Ms.subject);
cmap = lines(numel(animals));
[~, ia] = ismember(Ms.subject, animals);
fig = figure('Name', 'replay summary', 'Position', [40 40 1700 420], 'Color', 'w');
tl = tiledlayout(1, 5, 'TileSpacing', 'compact');
panels = {{'HPC_rate_pre', 'HPC_rate_post'}, {'V1_rate_pre', 'V1_rate_post'}, ...
    {'V1z_same_pre', 'V1z_same_post', 'V1z_other_pre', 'V1z_other_post'}, ...
    {'V1sig_same_pre', 'V1sig_same_post', 'V1sig_other_pre', 'V1sig_other_post'}};
labels = {{'pre', 'post'}, {'pre', 'post'}, {'same pre', 'same post', 'other pre', 'other post'}, {'same pre', 'same post', 'other pre', 'other post'}};
ylabs = {'HPC significant replay, either track (fraction)', 'V1 significant sequence, either track (fraction)', ...
    'V1 shuffle z, HPC replay events', 'V1 significant, HPC replay events'};
for k = 1:numel(panels)
    nexttile; hold on;
    Y = Ms{:, panels{k}};
    for c = 1:2:size(Y, 2)
        for s = 1:size(Y, 1)
            plot([c c + 1], Y(s, [c c + 1]), '-', 'Color', [cmap(ia(s), :) 0.5]);
        end
    end
    for c = 1:size(Y, 2)
        scatter(c * ones(size(Y, 1), 1), Y(:, c), 18, cmap(ia, :), 'filled');
        plot(c + [-0.25 0.25], mean(Y(:, c), 'omitnan') * [1 1], 'k-', 'LineWidth', 2);
    end
    % chance: p < 0.05 on either of two tracks (panels 1-2) or on one track (panel 4)
    switch k
        case {1, 2}, yline(1 - 0.95^2, 'k--', 'chance (either track)', 'FontSize', 7);
        case 3, yline(0, 'k--');
        case 4, yline(0.05, 'k--', 'chance', 'FontSize', 7);
    end
    pv = signrank(Y(:, 2), Y(:, 1));
    ttl = sprintf('post vs pre p = %.3g', pv);
    if size(Y, 2) == 4
        ttl = {sprintf('same: post vs pre p = %.2g', pv), sprintf('(same-other) post-pre p = %.2g', ...
            signrank((Y(:, 2) - Y(:, 1)) - (Y(:, 4) - Y(:, 3))))};
    end
    title(ttl, 'FontWeight', 'normal', 'FontSize', 8);
    set(gca, 'XTick', 1:size(Y, 2), 'XTickLabel', labels{k}, 'TickDir', 'out'); xlim([0.5 size(Y, 2) + 0.5]);
    ylabel(ylabs{k});
end
nexttile; hold on;
scatter(ia + 0.1 * randn(size(ia)), Ms.r_logodds_post, 25, cmap(ia, :), 'filled');
yline(0, 'k--');
set(gca, 'XTick', 1:numel(animals), 'XTickLabel', animals, 'TickDir', 'out', 'XTickLabelRotation', 45);
ylabel('Spearman r, HPC vs V1 log odds (post window)');
title(sprintf('mean r = %.2f, signrank p = %.2g', mean(Ms.r_logodds_post), signrank(Ms.r_logodds_post)), 'FontWeight', 'normal', 'FontSize', 8);
title(tl, sprintf(['POST ripples, multinomial %d ms decoding (sessions = dots, colour = animal). Same/other = V1 on the ' ...
    'track HPC replays / the other track (HPC-significant events). Scoring window = onset to max(0.1 s, offset); pre = same n bins before onset'], ...
    round(P.time_bin * 1000)), 'FontSize', 9);
save_figure(fig, fullfile(P.out_dir, 'replay_summary'));
close(fig);
end

function plot_best_examples(ex, P)
% Across sessions (scored HPC post window, >= 8 active cells): the n_best HPC replays
% with the largest shuffle z -> replay_best_examples1, 2, ...; n_best non-significant
% events (p >= 0.05 on both tracks) drawn at random -> replay_nonsignificant_examples1, 2, ...
% Same columns as the per-session example figures, n_per_page events per figure.
n_best = 20;
if isempty(ex), return; end
best = ex([ex.pmin] < 0.05);
[~, o] = sort([best.score], 'descend');
best = best(o(1:min(n_best, numel(o))));
plot_example_pages(best, P, 'replay_best_examples', 'Most significant HPC replays across sessions (largest shuffle z)');
nonsig = ex([ex.pmin] >= 0.05);
rs = RandStream('mrg32k3a', 'Seed', P.seed);
nonsig = nonsig(sort(randperm(rs, numel(nonsig), min(n_best, numel(nonsig)))));
plot_example_pages(nonsig, P, 'replay_nonsignificant_examples', ...
    'Non-significant HPC events across sessions (p >= 0.05 on both tracks, random draw)');
end

function plot_example_pages(ex, P, name, description)
n_per_page = 10;
if isempty(ex), return; end
% Load each session once; keep only the spikes around each event
files = unique({ex.file}, 'stable');
outs = cell(size(files)); spk = cell(numel(ex), 2);
for f = 1:numel(files)
    S = load(files{f}, 'out'); outs{f} = S.out;
    [t, cell_idx] = load_sleep_spikes(session_experiment_info(S.out.session_index), S.out.cell_id);
    for k = find(strcmp({ex.file}, files{f}))
        onset = S.out.events.onset(ex(k).e);
        s = t >= onset - 1 & t < onset + 1;
        spk(k, :) = {t(s), cell_idx(s)};
    end
end
[~, file_of] = ismember({ex.file}, files);
nPages = ceil(numel(ex) / n_per_page);
for pg = 1:nPages
    idx = (pg - 1) * n_per_page + 1:min(pg * n_per_page, numel(ex));
    nE = numel(idx);
    [fig, tl] = example_figure(nE);
    for k = 1:nE
        out = outs{file_of(idx(k))};
        plot_example_column(tl, k, nE, out, ex(idx(k)).e, spk{idx(k), 1}, spk{idx(k), 2}, P, ...
            sprintf('%s %s (s%d)\n', out.subject, out.session, out.session_index));
    end
    title(tl, {sprintf('%s, %d-%d of %d', description, idx(1), idx(end), numel(ex)), example_legend(P)}, ...
        'Interpreter', 'none', 'FontSize', 9);
    save_figure(fig, fullfile(P.out_dir, sprintf('%s%d', name, pg)));
    close(fig);
end
end

function extra_shuffles_all(P)
% Additional shuffles (P.shuffle_types) of every scored event of the decoded sessions ->
% <animal>\shuffles_session_<date>.mat next to each session file
files = dir(fullfile(P.out_dir, 'M*', 'session_*.mat'));
for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name);
    S = load(f, 'out'); out = S.out;
    if ~ismember(out.session_index, P.sessions), continue; end
    sh_file = fullfile(files(k).folder, ['shuffles_' files(k).name]);
    sh = struct();
    if exist(sh_file, 'file'), S = load(sh_file, 'sh'); sh = S.sh; end
    tic
    sh = session_shuffles(out, sh, P);
    save(sh_file, 'sh', '-v7.3');
    fprintf('Session %d (%s %s): %s done (%.0f s)\n', out.session_index, out.subject, out.session, ...
        strjoin(P.shuffle_types, ', '), toc);
end
end

function sh = session_shuffles(out, sh, P)
% Null distributions of the stored scored events of one session (see shuffle_types)
Q = out.params; nPos = Q.nPos; nS = Q.nShuffle;
types = P.shuffle_types;
type_seed = struct('cell_id', 0, 'spike', 100000, 'field', 200000, 'time', 300000);
[t, cell_idx] = load_sleep_spikes(session_experiment_info(out.session_index), out.cell_id);
nEv = numel(out.events.onset);
sh.subject = out.subject; sh.session = out.session; sh.session_index = out.session_index;
sh.params = get_field(sh, 'params', struct());
sh.params.nShuffle = nS;
sh.params.seed = Q.seed;
sh.params.types_run = union(get_field(sh.params, 'types_run', {}), types);
sh.params.note = ['p = (1 + #(|r_sh| >= |r|)) / (1 + nShuffle), time exhaustive: fraction of arrangements; ' ...
    'chance = fraction of null replicates with p < 0.05 within the same null; field: shared decoders per ' ...
    'session and region, T1/T2 shifted independently; time: all window bins permuted (empty bins move)'];
for nr = 1:numel(out.region)
    R = out.region(nr);
    cells = region_cells(out.region_labels, R.name);
    cells = cells(:) & get_field(out, 'decoder_cell', true(numel(out.cell_id), 1));   % absent = all cells (decoded before cell_filter)
    assert(isequal(out.cell_id(cells), R.cell_id), 'Decoder cells of %s do not match the stored decoder', R.name);
    nCr = sum(cells);
    new_idx = zeros(numel(out.cell_id), 1); new_idx(cells) = 1:nCr;
    in = cells(cell_idx);
    tr = t(in); cr = new_idx(cell_idx(in));
    mr = max(out.rate_maps(cells, :), Q.rate_floor);
    logp = log(mr ./ sum(mr, 1));
    logp_field = [];
    if ismember('field', types)
        rs = RandStream('mrg32k3a', 'Seed', Q.seed + 1000 * out.session_index + 10 * nr + type_seed.field);
        logp_field = field_shift_decoders(mr, nPos, nS, rs);
    end
    % event log odds z-scored against the per-bin log odds of all decoded bins
    lb = double(R.logodds_bin(:)); lb = lb(isfinite(lb));
    mu = mean(lb); sd = std(lb);
    Rs = get_field(sh, R.name, struct());
    Rs.logodds_bin_mean = mu; Rs.logodds_bin_sd = sd;
    W = cell(1, numel(Q.windows));
    for w = 1:numel(Q.windows)
        Y = get_field(Rs, Q.windows{w}, struct());
        z = (double(R.(Q.windows{w}).logodds) - mu) / sd;
        fin = isfinite(z);
        if any(fin)
            z(z == Inf) = prctile(z(fin), 99.5); z(z == -Inf) = prctile(z(fin), 0.5);
        end
        Y.logodds_z = z;
        for ty = types
            Y.(['p_' ty{1}]) = NaN(nEv, 2); Y.(['z_' ty{1}]) = NaN(nEv, 2);
            Y.(['mean_' ty{1}]) = NaN(nEv, 2); Y.(['sd_' ty{1}]) = NaN(nEv, 2);
            Y.(['chance_' ty{1}]) = NaN(nEv, 2); Y.(['chance_either_' ty{1}]) = NaN(nEv, 1);
        end
        if ismember('time', types), Y.n_perm_time = NaN(nEv, 1); end
        W{w} = Y;
    end
    lo = first_ge(tr, out.events.onset + Q.window(1));
    hi = first_ge(tr, out.events.onset + Q.window(2)) - 1;
    max_dwc = 0; max_dp = 0;
    for e = 1:nEv
        scored = cellfun(@(wn) R.(wn).scored(e), Q.windows);
        if ~any(scored), continue; end
        s = (lo(e):hi(e))';
        b = floor((tr(s) - (out.events.onset(e) + Q.window(1))) / Q.time_bin) + 1;
        ok = b >= 1 & b <= Q.nBins;
        counts = accumarray([cr(s(ok)) b(ok)], 1, [nCr Q.nBins]);
        nb = out.events.n_post_bins(e);
        win_bins = {Q.onset_bin + (0:nb - 1), Q.onset_bin - nb + (0:nb - 1)};
        for w = find(scored)
            X = R.(Q.windows{w});
            cw = counts(:, win_bins{w});
            r = real_wc(cw, logp, nPos);
            max_dwc = max([max_dwc abs(r - X.wc(e, :))]);
            for ty = types
                rs = RandStream('mrg32k3a', 'Seed', Q.seed + 1000 * out.session_index + 10 * nr + w + type_seed.(ty{1}));
                rs.Substream = out.events.ripple_info_row(e);
                [r_null, exact] = shuffle_null(ty{1}, cw, logp, logp_field, nPos, nS, rs);
                [pv, zv, mv, sv, ch, che] = null_stats(r, r_null, exact, strcmp(ty{1}, 'cell_id'));
                Y = W{w};
                Y.(['p_' ty{1}])(e, :) = pv; Y.(['z_' ty{1}])(e, :) = zv;
                Y.(['mean_' ty{1}])(e, :) = mv; Y.(['sd_' ty{1}])(e, :) = sv;
                Y.(['chance_' ty{1}])(e, :) = ch; Y.(['chance_either_' ty{1}])(e) = che;
                if strcmp(ty{1}, 'time'), Y.n_perm_time(e) = size(r_null, 1); end
                if strcmp(ty{1}, 'cell_id'), max_dp = max([max_dp abs(pv - X.p(e, :))]); end
                W{w} = Y;
            end
        end
    end
    assert(max_dwc < 1e-6, '%s: re-decoded |r| differs from the stored wc (max %.2g)', R.name, max_dwc);
    if ismember('cell_id', types) && max_dp > 0
        warning('%s session %d: recomputed cell_id p differs from the stored p (max %.3g)', R.name, out.session_index, max_dp);
    end
    for w = 1:numel(Q.windows), Rs.(Q.windows{w}) = W{w}; end
    sh.(R.name) = Rs;
end
end

function v = get_field(s, name, default)
if isfield(s, name), v = s.(name); else, v = default; end
end

function r = real_wc(cw, logp, nPos)
% |weighted correlation| per track of one window (empty bins left out), as score_window
r = NaN(1, 2);
ne = sum(cw, 1) > 0;
A = find(sum(cw, 2) > 0);
if isempty(A), return; end
post = softmax1(logp(A, :)' * cw(A, ne));
for k = 1:2
    Pk = post((k - 1) * nPos + (1:nPos), :);
    r(k) = wcorr_batch(reshape(Pk ./ sum(Pk, 1), [1 nPos sum(ne)]), find(ne));
end
end

function logp_field = field_shift_decoders(mr, nPos, nS, rs)
% nS place field circular shift decoders: every cell's map shifted independently on each
% track by 0 to nPos - 1 bins, multinomial log likelihoods (nCells x 2 nPos x nS)
nC = size(mr, 1);
shift = randi(rs, nPos, nC, 2, nS) - 1;
M = zeros(nC, 2 * nPos, nS);
rowi = repmat((1:nC)', 1, nPos, nS);
for k = 1:2
    cols = (k - 1) * nPos + (1:nPos);
    Mk = mr(:, cols);
    col = mod((0:nPos - 1) - reshape(shift(:, k, :), nC, 1, nS), nPos) + 1;
    M(:, cols, :) = Mk(sub2ind([nC nPos], rowi, col));
end
logp_field = log(M ./ sum(M, 1));
end

function [r_null, exact] = shuffle_null(type, cw, logp, logp_field, nPos, nS, rs)
% |weighted correlation| per track (columns) of each null replicate (rows); exact = true
% when the null is the complete set of arrangements (time, few bins)
exact = false;
nbw = size(cw, 2);
ne = sum(cw, 1) > 0; tb = find(ne);
A = find(sum(cw, 2) > 0); nA = numel(A);
switch type
    case {'cell_id', 'field'}
        if strcmp(type, 'cell_id')
            % identical to score_window (same random draws)
            [~, perm] = sort(rand(rs, nS, size(logp, 1)), 2);
            idx = perm(:, 1:nA)';
            G = reshape(logp(idx(:), :), nA, nS, 2 * nPos);
            G = reshape(permute(G, [2 3 1]), nS * 2 * nPos, nA);
        else
            G = reshape(permute(logp_field(A, :, :), [3 2 1]), nS * 2 * nPos, nA);
        end
        L = reshape(G * cw(A, ne), nS, 2 * nPos, numel(tb));
        L = exp(L - max(L, [], 2)); L = L ./ sum(L, 2);
        r_null = zeros(nS, 2);
        for k = 1:2
            Wk = L(:, (k - 1) * nPos + (1:nPos), :); Wk = Wk ./ sum(Wk, 2);
            r_null(:, k) = abs(wcorr_batch(Wk, tb));
        end
    case 'spike'
        shift = randi(rs, nbw, nA, 1, nS) - 1;
        col = mod((0:nbw - 1) - shift, nbw) + 1;                 % nA x nbw x nS
        cwA = cw(A, :);
        C = cwA(sub2ind([nA nbw], repmat((1:nA)', 1, nbw, nS), col));
        L = reshape(logp(A, :)' * reshape(C, nA, nbw * nS), 2 * nPos, nbw, nS);
        post = exp(L - max(L, [], 1)); post = post ./ sum(post, 1);
        has = sum(C, 1) > 0;                                     % non-empty bins of each shuffle
        r_null = zeros(nS, 2);
        for k = 1:2
            Pk = post((k - 1) * nPos + (1:nPos), :, :);
            Pk = Pk ./ sum(Pk, 1) .* has;                        % empty bins: zero weight
            r_null(:, k) = abs(wcorr_batch(permute(Pk, [3 1 2]), 1:nbw));
        end
    case 'time'
        nt = numel(tb);
        post = softmax1(logp(A, :)' * cw(A, ne));
        n_arr = exp(gammaln(nbw + 1) - gammaln(nbw - nt + 1));   % ordered placements of the non-empty bins
        if n_arr <= nS
            c = nchoosek(1:nbw, nt); q = perms(1:nt);
            T = zeros(size(c, 1) * size(q, 1), nt);
            for i = 1:size(q, 1)
                T((i - 1) * size(c, 1) + (1:size(c, 1)), :) = c(:, q(i, :));
            end
            exact = true;
        else
            T = zeros(nS, nt);
            for j = 1:nS, T(j, :) = randperm(rs, nbw, nt); end
        end
        nP = size(T, 1);
        r_null = zeros(nP, 2);
        for k = 1:2
            Pk = post((k - 1) * nPos + (1:nPos), :); Pk = Pk ./ sum(Pk, 1);
            W2 = zeros(nP * nbw, nPos);
            for i = 1:nt
                W2(sub2ind([nP nbw], (1:nP)', T(:, i)), :) = repmat(Pk(:, i)', nP, 1);
            end
            r_null(:, k) = abs(wcorr_batch(permute(reshape(W2, nP, nbw, nPos), [1 3 2]), 1:nbw));
        end
end
r_null(~isfinite(r_null)) = 0;
end

function [pv, zv, mv, sv, ch, che] = null_stats(r, r_null, exact, strict)
% p, z, null mean/sd per track, and the empirical chance of p < 0.05 per track / either
% track (each null replicate's p within the same null). strict = no tie tolerance (cell_id,
% to reproduce the stored p exactly)
tol = 1e-12 * ~strict;
nN = size(r_null, 1);
pv = NaN(1, 2); zv = pv; ch = pv;
mv = mean(r_null, 1); sv = std(r_null, 0, 1);
p_rep = zeros(nN, 2);
for k = 1:2
    if isfinite(r(k))
        n_ge = sum(r_null(:, k) >= abs(r(k)) - tol);
        if exact, pv(k) = n_ge / nN; else, pv(k) = (1 + n_ge) / (1 + nN); end
        zv(k) = (abs(r(k)) - mv(k)) / sv(k);
    end
    [~, ~, rk] = unique(-round(r_null(:, k), 12));
    p_grp = cumsum(accumarray(rk, 1)) / nN;
    p_rep(:, k) = p_grp(rk);
    ch(k) = mean(p_rep(:, k) < 0.05);
end
che = mean(min(p_rep, [], 2) < 0.05);
end

function summarise_shuffles(P)
% Event table (cell_id from the session files, other shuffles from the shuffles_ files)
% and the shuffle comparison figures -> out_dir\shuffle_comparison
types_all = {'cell_id', 'spike', 'field', 'time'};
cmp_dir = fullfile(P.out_dir, 'shuffle_comparison');
if ~exist(cmp_dir, 'dir'), mkdir(cmp_dir); end
files = dir(fullfile(P.out_dir, 'M*', 'session_*.mat'));
T = table();
for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name), 'out'); out = S.out;
    sh = struct();
    sh_file = fullfile(files(k).folder, ['shuffles_' files(k).name]);
    if exist(sh_file, 'file'), S = load(sh_file, 'sh'); sh = S.sh; end
    nEv = numel(out.events.onset);
    Tk = table(repmat(string(out.subject), nEv, 1), repmat(out.session_index, nEv, 1), (1:nEv)', ...
        out.events.ripple_info_row, 'VariableNames', {'subject', 'session_index', 'event', 'ripple_info_row'});
    for nr = 1:numel(out.region)
        R = out.region(nr);
        for w = 1:numel(P.windows)
            X = R.(P.windows{w});
            Y = get_field(get_field(sh, R.name, struct()), P.windows{w}, struct());
            pre = sprintf('%s_%s_', R.name, P.windows{w});
            Tk.([pre 'scored']) = X.scored; Tk.([pre 'n_active']) = X.n_active;
            Tk.([pre 'logodds']) = X.logodds; Tk.([pre 'logodds_z']) = get_field(Y, 'logodds_z', NaN(nEv, 1));
            for ty = types_all
                if strcmp(ty{1}, 'cell_id')
                    pv = X.p; zv = X.z;
                else
                    pv = get_field(Y, ['p_' ty{1}], NaN(nEv, 2)); zv = get_field(Y, ['z_' ty{1}], NaN(nEv, 2));
                end
                for tr = 1:2
                    Tk.(sprintf('%sp_%s_T%d', pre, ty{1}, tr)) = pv(:, tr);
                    Tk.(sprintf('%sz_%s_T%d', pre, ty{1}, tr)) = zv(:, tr);
                end
                Tk.([pre 'chance_either_' ty{1}]) = get_field(Y, ['chance_either_' ty{1}], NaN(nEv, 1));
            end
        end
    end
    T = [T; Tk]; %#ok<AGROW>
end
save(fullfile(cmp_dir, 'replay_events_shuffles.mat'), 'T', '-v7.3');
writetable(T, fullfile(cmp_dir, 'replay_events_shuffles.csv'));
for c = 1:numel(P.significance_criteria)
    try
        plot_criterion(T, P.significance_criteria{c}, P, cmp_dir);
    catch err
        warning('Criterion figure %s failed: %s', strjoin(P.significance_criteria{c}, '+'), err.message);
    end
end
for reg = P.regions
    try
        plot_agreement(T, reg{1}, P, cmp_dir);
    catch err
        warning('Agreement figure %s failed: %s', reg{1}, err.message);
    end
end
end

function E = shuffle_p(T, reg, win, types, field)
% p (or z) values of the given shuffle types from the event table, E.(type) = n x 2
if nargin < 5, field = 'p'; end
E = struct();
for ty = types
    E.(ty{1}) = [T.(sprintf('%s_%s_%s_%s_T1', reg, win, field, ty{1})) T.(sprintf('%s_%s_%s_%s_T2', reg, win, field, ty{1}))];
end
end

function plot_criterion(T, crit, P, cmp_dir)
% Replay rates (pre vs post) and proportion of significant events vs mean log odds
% difference across alpha levels, for one significance criterion, HPC and V1
fig = figure('Color', 'w', 'Position', [40 40 1500 820]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for reg = P.regions
    nexttile; plot_rate_panel(T, reg{1}, crit);
    nexttile; plot_alpha_curve(T, reg{1}, crit, 'logodds', P);
    nexttile; plot_alpha_curve(T, reg{1}, crit, 'logodds_z', P);
end
title(tl, {sprintf('Significance criterion: %s (p < alpha on every listed shuffle, same track)', strjoin(crit, ' + ')), ...
    ['Scored events (>= 5 active cells, >= 3 non-empty bins). Left: alpha = 0.05, sessions = dots (colour = animal), ' ...
    'dashed = mean empirical chance.'], ['Middle/right: events pooled, alpha 0.2 to 0.001 (larger marker = stricter), ' ...
    'shaded = 95% bootstrap CI (events resampled); log odds z = relative to all decoded bins of the session']}, ...
    'Interpreter', 'none', 'FontSize', 9);
save_figure(fig, fullfile(cmp_dir, ['criterion_' strjoin(crit, '_')]));
close(fig);
end

function plot_rate_panel(T, reg, crit)
hold on
sess = unique(T.session_index); nS = numel(sess);
subj = strings(nS, 1);
for k = 1:nS, subj(k) = T.subject(find(T.session_index == sess(k), 1)); end
animals = unique(subj); cmap = lines(numel(animals)); [~, ia] = ismember(subj, animals);
wins = {'pre', 'post'};
Y = NaN(nS, 2); CH = NaN(nS, 2);
for c = 1:2
    ok = T.(sprintf('%s_%s_scored', reg, wins{c}));
    sig = replay_significance(shuffle_p(T, reg, wins{c}, crit), crit, 0.05);
    either = any(sig, 2);
    ch = NaN(height(T), 1);
    if isscalar(crit), ch = T.(sprintf('%s_%s_chance_either_%s', reg, wins{c}, crit{1})); end
    for k = 1:nS
        in = ok & T.session_index == sess(k);
        Y(k, c) = mean(either(in)); CH(k, c) = mean(ch(in), 'omitnan');
    end
end
for s = 1:nS, plot([1 2], Y(s, :), '-', 'Color', [cmap(ia(s), :) 0.5]); end
for c = 1:2
    scatter(c * ones(nS, 1), Y(:, c), 18, cmap(ia, :), 'filled');
    plot(c + [-0.25 0.25], mean(Y(:, c), 'omitnan') * [1 1], 'k-', 'LineWidth', 2);
    if any(isfinite(CH(:, c))), plot(c + [-0.35 0.35], mean(CH(:, c), 'omitnan') * [1 1], 'k--'); end
end
if all(isnan(CH(:)))
    if isscalar(crit), note = 'empirical chance not computed for this shuffle'; else, note = 'no empirical chance for combined criteria'; end
    text(0.55, 0.97, note, 'Units', 'normalized', 'FontSize', 7, 'VerticalAlignment', 'top', 'Interpreter', 'none');
end
pv = NaN;
if sum(all(isfinite(Y), 2)) > 1, pv = signrank(Y(:, 2), Y(:, 1)); end
title(sprintf('%s: post %.3f vs pre %.3f, p = %.3g', reg, mean(Y(:, 2), 'omitnan'), mean(Y(:, 1), 'omitnan'), pv), ...
    'FontWeight', 'normal', 'FontSize', 8);
set(gca, 'XTick', 1:2, 'XTickLabel', wins, 'TickDir', 'out'); xlim([0.5 2.5]);
ylabel(sprintf('%s significant, either track (fraction)', reg));
end

function plot_alpha_curve(T, reg, crit, lo_name, P)
% Proportion of significant events (multi-track counted once) vs mean log odds of T1
% events - mean log odds of T2 events (multi-track events in both), across alpha levels
hold on
cols = {[0.4 0.4 0.4], [0.85 0.33 0.1]}; wins = {'pre', 'post'};
nA = numel(P.alphas); msz = linspace(12, 80, nA);
rs = RandStream('mrg32k3a', 'Seed', P.seed);
h = gobjects(1, 2);
for c = 1:2
    ok = T.(sprintf('%s_%s_scored', reg, wins{c}));
    sig = replay_significance(shuffle_p(T(ok, :), reg, wins{c}, crit), crit, P.alphas);
    lo = T.(sprintf('%s_%s_%s', reg, wins{c}, lo_name))(ok);
    [prop, dif] = alpha_stats(sig, lo);
    n = sum(ok);
    propB = NaN(P.nBoot, nA); difB = propB;
    for b = 1:P.nBoot
        idx = randi(rs, n, n, 1);
        [propB(b, :), difB(b, :)] = alpha_stats(sig(idx, :, :), lo(idx));
    end
    cp = prctile(propB, [2.5 97.5]); cd = prctile(difB, [2.5 97.5]);
    for a = 1:nA
        if all(isfinite([cp(:, a); cd(:, a)]))
            patch(cd([1 2 2 1], a), cp([1 1 2 2], a), cols{c}, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        end
    end
    h(c) = plot(dif, prop, '-', 'Color', cols{c});
    scatter(dif, prop, msz, cols{c}, 'filled', 'MarkerEdgeColor', 'k');
    if c == 2
        for a = 1:nA
            if isfinite(dif(a)), text(dif(a), prop(a), sprintf('  %g', P.alphas(a)), 'FontSize', 6); end
        end
    end
end
xline(0, 'k:');
xlabel(sprintf('mean %s difference (T1 - T2 events)', strrep(lo_name, '_', ' ')));
ylabel('proportion of significant events');
legend(h, wins, 'Location', 'northeast', 'Box', 'off');
title(sprintf('%s, %s', reg, strrep(lo_name, '_', ' ')), 'FontWeight', 'normal', 'FontSize', 8);
set(gca, 'TickDir', 'out');
end

function [prop, dif] = alpha_stats(sig, lo)
% sig: n x 2 x nAlpha, lo: n x 1 log odds
fin = isfinite(lo);
prop = reshape(mean(any(sig, 2), 1), 1, []);
dif = NaN(1, size(sig, 3));
for a = 1:size(sig, 3)
    t1 = sig(:, 1, a) & fin; t2 = sig(:, 2, a) & fin;
    if any(t1) && any(t2), dif(a) = mean(lo(t1)) - mean(lo(t2)); end
end
end

function plot_agreement(T, reg, P, cmp_dir)
% Agreement between single shuffles at alpha = 0.05 (either track), post and pre windows
types_all = {'cell_id', 'spike', 'field', 'time'};
avail = types_all(cellfun(@(ty) any(isfinite(T.(sprintf('%s_post_p_%s_T1', reg, ty)))), types_all));
nT = numel(avail);
if nT < 2, return; end
fig = figure('Color', 'w', 'Position', [40 40 1900 780]);
tl = tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
for win = P.windows
    ok = T.(sprintf('%s_%s_scored', reg, win{1}));
    n = sum(ok);
    Ep = shuffle_p(T(ok, :), reg, win{1}, avail); Ez = shuffle_p(T(ok, :), reg, win{1}, avail, 'z');
    sig = false(n, nT); trk = zeros(n, nT); pp = cell(1, nT);
    for i = 1:nT
        pv = Ep.(avail{i}); zv = Ez.(avail{i});
        sig(:, i) = any(pv < 0.05, 2);
        [~, tk] = min(pv, [], 2);
        tie = pv(:, 1) == pv(:, 2); [~, tz] = max(zv, [], 2); tk(tie) = tz(tie);
        trk(:, i) = tk; pp{i} = pv(:);
    end
    cond = NaN(nT); jac = NaN(nT); rho = NaN(nT); same = NaN(nT);
    for i = 1:nT
        for j = 1:nT
            cond(i, j) = mean(sig(sig(:, i), j));
            jac(i, j) = sum(sig(:, i) & sig(:, j)) / sum(sig(:, i) | sig(:, j));
            rho(i, j) = corr(pp{i}, pp{j}, 'Type', 'Spearman', 'Rows', 'complete');
            both = sig(:, i) & sig(:, j);
            same(i, j) = mean(trk(both, i) == trk(both, j));
        end
    end
    rates = strjoin(cellfun(@(ty, r) sprintf('%s %.3f', ty, r), avail, num2cell(mean(sig, 1)), 'UniformOutput', false), ', ');
    nexttile; plot_heat(cond, avail, {sprintf('%s %s (n = %d): P(sig column | sig row)', reg, win{1}, n), rates}, [0 1]);
    nexttile; plot_heat(jac, avail, 'Jaccard overlap of significant events', [0 1]);
    nexttile; plot_heat(rho, avail, 'Spearman r of p (per track, pooled)', [0 1]);
    nexttile; plot_heat(same, avail, 'same replayed track | both significant', [0 1]);
    nexttile;
    code = sig * (2 .^ (0:nT - 1))';
    cnt = accumarray(code + 1, 1, [2 ^ nT 1]); cnt = cnt(2:end);
    names = arrayfun(@(m) strjoin(avail(bitget(m, 1:nT) == 1), '+'), (1:2 ^ nT - 1)', 'UniformOutput', false);
    [cs, o] = sort(cnt, 'descend');
    bar(cs / n, 'FaceColor', [0.5 0.5 0.5]);
    set(gca, 'XTick', 1:numel(o), 'XTickLabel', names(o), 'XTickLabelRotation', 60, 'TickLabelInterpreter', 'none', ...
        'TickDir', 'out', 'FontSize', 7);
    ylabel('fraction of scored events');
    title(sprintf('significant on exactly these shuffles (any: %.3f)', mean(any(sig, 2))), 'FontWeight', 'normal', 'FontSize', 8);
end
title(tl, sprintf('%s: agreement between shuffles (p < 0.05 on either track; replayed track = lower p, ties: larger z)', reg), ...
    'Interpreter', 'none', 'FontSize', 9);
save_figure(fig, fullfile(cmp_dir, sprintf('shuffle_agreement_%s', reg)));
close(fig);
end

function plot_heat(M, lab, ttl, lims)
imagesc(M, lims); axis square; colorbar;
for i = 1:size(M, 1)
    for j = 1:size(M, 2)
        if isfinite(M(i, j))
            text(j, i, sprintf('%.2f', M(i, j)), 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [1 1 1] * (M(i, j) < 0.6));
        end
    end
end
set(gca, 'XTick', 1:numel(lab), 'XTickLabel', lab, 'YTick', 1:numel(lab), 'YTickLabel', lab, ...
    'TickLabelInterpreter', 'none', 'FontSize', 7);
title(ttl, 'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 8);
end

function info = session_experiment_info(nsession)
% experiment_info of the 22 bilateral sessions (index = session number)
persistent INFO
if isempty(INFO)
    SUBJECTS = {'M24016', 'M24017', 'M24018', 'M24062', 'M24064', 'M24065'};
    INFO = subject_session_stimuli_mapping(SUBJECTS, 'bilateral');
    INFO = INFO([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
end
info = INFO(nsession);
end

function save_figure(fig, save_path)
% Retries: files on D: are occasionally locked for a moment after being written
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
