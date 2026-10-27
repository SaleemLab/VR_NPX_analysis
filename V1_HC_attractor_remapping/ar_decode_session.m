function out = ar_decode_session(Ds, T, s, P, old_setups)
% Track decoding of one RUN1 session in the 100 ms bins of ar_bin_session (a
% whole-session grid plus bins locked to lap ends and onsets), for each region
% in P.regions:
%   Bayesian  multinomial (main) and Poisson, 28 joint states (2 tracks x 14
%             10 cm bins), flat prior; track log odds = log of the posterior mass
%             on T1 over T2 (position marginalised, so no position ground truth
%             is needed in the grey screen)
%   PLS-KDE   location-agnostic track decoder (paper decoder: PLS on good
%             positions, 3-D KDE per track); log odds = log KDE(T1) - log KDE(T2)
% Cross-validation by held-out blocks: windows in the segment of lap k (the
% grey screen before lap k and lap k itself) are decoded by a decoder trained on
% laps outside blocks b-h..b+h (b = block of lap k), with the number of training
% laps equalised between tracks (majority track randomly subsampled):
%   LO  : h = 1 ('lobo3'): the nearest same-track training blocks (b+-2) are closer
%         in time than the nearest other-track blocks (b+-3), so slow drift would
%         favour the current block's track (liberal for a block-track bias)
%   LO5 : h = 2 ('lobo5'): nearest other-track blocks (b+-3) are closer than
%         same-track blocks (b+-4), so drift would favour the other track
%         (conservative); blocks without training laps of both tracks are NaN
% A decoder trained on all laps (LO_all) is kept for grey-screen sensitivity
% checks only (lap windows are then in its training data).
% old_setups (optional): struct with a PLS setup per region name (from an
% earlier run on the same training bins) to skip the shuffle-based setup.
% Output out: subject, session, session_index, G (bins without counts, plus
% n_spikes per region), region(r): name, n_cells, cell_id, and per decoder
% (bayes_multinomial, bayes_poisson, pls) LO, LO5, LO_all, and for Bayesian the
% posteriors (28 x nBins) of lobo3 and lobo5, for PLS proj (nBins x 3, lobo3) and setup.
if nargin < 5, old_setups = struct(); end

G = ar_bin_session(Ds, T, P);
nWin = numel(G.t_start);
nLaps = height(T);
nPos = numel(P.pos_edges) - 1;
nBlocks = max(T.block);
win_block = T.block(G.segment);
[occ, cnt] = ar_bayes_lap_stats(Ds, T, P);
% training laps per block for each fold set (h = 1: lobo3, h = 2: lobo5)
train_laps = cell(nBlocks, 2);
for h = 1:2
    for b = 1:nBlocks
        train_laps{b, h} = ar_balanced_train_laps(T, b, h, 1000 * s + 10 * b + h);
    end
end

out.subject = Ds.subject; out.session = Ds.session; out.session_index = s;
out.params = P;
Gs = rmfield(G, 'counts');
for r = 1:numel(P.regions)
    name = P.regions{r};
    cells = region_cells(Ds.region, name);
    R = struct('name', name, 'n_cells', sum(cells), 'cell_id', Ds.cell_id(cells));
    Gs.n_spikes.(name) = sum(G.counts(:, cells), 2);
    if sum(cells) < P.min_cells
        fprintf('  %-5s %3d cells: skipped\n', name, sum(cells));
        out.region(r) = fill_empty(R);
        continue
    end
    C = double(G.counts(:, cells))';    % cells x windows

    % ---- Bayesian ----
    for method = {'multinomial', 'poisson'}
        LO = NaN(nWin, 2); post = {NaN(2 * nPos, nWin, 'single'), NaN(2 * nPos, nWin, 'single')};
        for h = 1:2
            for b = 1:nBlocks
                w = win_block == b;
                if ~any(w) || isempty(train_laps{b, h}), continue; end
                train = ismember((1:nLaps)', train_laps{b, h});
                maps = ar_rate_maps(occ, cnt(cells, :, :), train, T.track, nPos);
                [LO(w, h), post{h}(:, w)] = ar_bayes_track(C(:, w), maps, P.bin_width, method{1}, P.rate_floor);
            end
        end
        maps = ar_rate_maps(occ, cnt(cells, :, :), true(nLaps, 1), T.track, nPos);
        LO_all = ar_bayes_track(C, maps, P.bin_width, method{1}, P.rate_floor);
        R.(['bayes_' method{1}]) = struct('LO', single(LO(:, 1)), 'LO5', single(LO(:, 2)), 'LO_all', single(LO_all), ...
            'posterior', post{1}, 'posterior5', post{2});
    end

    % ---- PLS-KDE ----
    X = C';
    base = G.nonoverlap & G.speed >= P.pls_speed;
    mu1 = mean(X(base, :), 1); sd1 = std(X(base, :), 0, 1);
    Z = (X - mu1) ./ sd1; Z(:, sd1 == 0) = 0;
    Z = Z - mean(Z, 2);
    trn = G.nonoverlap & G.lap > 0 & G.speed >= P.pls_speed;
    pos5 = discretize(G.pos, P.pls_pos_edges);
    trn = trn & ~isnan(pos5);
    Y = T.track(G.lap(trn));
    if isfield(old_setups, name) && ~isempty(old_setups.(name))
        setup = old_setups.(name);
    else
        setup = ar_pls_setup(Z(trn, :), Y, pos5(trn), P.pls_nShuffle);
    end
    good = trn & ismember(pos5, setup.good_position);
    LO = NaN(nWin, 2); proj = NaN(nWin, 3, 'single'); LO_all = NaN(nWin, 1);
    for h = 1:2 * ~isnan(setup.bandwidth)
        for b = 1:nBlocks
            w = win_block == b;
            if ~any(w) || isempty(train_laps{b, h}), continue; end
            [lo, pj] = pls_fold(Z, good & ismember(G.lap, train_laps{b, h}), T.track, G.lap, w, setup.bandwidth);
            LO(w, h) = lo;
            if h == 1, proj(w, :) = pj; end
        end
    end
    if ~isnan(setup.bandwidth)
        LO_all = pls_fold(Z, good, T.track, G.lap, true(nWin, 1), setup.bandwidth);
    end
    R.pls = struct('LO', single(LO(:, 1)), 'LO5', single(LO(:, 2)), 'LO_all', single(LO_all), 'proj', proj, 'setup', setup);

    % In-lap cross-validated track accuracy (running windows) for the log
    run = G.kind == 0 & G.lap > 0 & G.speed > P.bayes_speed;
    true1 = T.track(G.lap(run)) == 1;
    acc = @(lo) track_accuracy(lo(run), true1);
    fprintf('  %-5s %3d cells: in-lap track acc lobo3/lobo5: multinomial %.2f/%.2f, poisson %.2f/%.2f, PLS %.2f/%.2f (bw %.1f, %d good pos)\n', ...
        name, R.n_cells, acc(R.bayes_multinomial.LO), acc(R.bayes_multinomial.LO5), acc(R.bayes_poisson.LO), ...
        acc(R.bayes_poisson.LO5), acc(R.pls.LO), acc(R.pls.LO5), setup.bandwidth, numel(setup.good_position));
    out.region(r) = R;
end
out.G = Gs;
end

function R = fill_empty(R)
e = struct('LO', [], 'LO5', [], 'LO_all', [], 'posterior', [], 'posterior5', []);
R.bayes_multinomial = e; R.bayes_poisson = e;
R.pls = struct('LO', [], 'LO5', [], 'LO_all', [], 'proj', [], 'setup', []);
end

function a = track_accuracy(lo, true1)
% fraction of decoded windows (finite LO) on the true track; LO = 0 counts as wrong
ok = isfinite(lo);
a = mean((lo(ok) > 0) == true1(ok) & lo(ok) ~= 0);
end

function cells = region_cells(region, name)
switch name
    case 'HPC', cells = contains(region, 'HPC');
    case 'V1', cells = contains(region, 'V1');
    otherwise, cells = region == string(name);
end
end

function [LO, proj] = pls_fold(Z, train, lap_track, win_lap, test, bw)
% PLS (3 components) on the training windows (z-scored with their own mean/std),
% KDE per track on the training projections, log odds of the test windows
Xtr = Z(train, :);
Ytr = lap_track(win_lap(train));
mu2 = mean(Xtr, 1); sd2 = std(Xtr, 0, 1); sd2(sd2 == 0) = Inf;
Xtr = (Xtr - mu2) ./ sd2;
XL = plsregress(Xtr, Ytr, 3);
[W, ~] = qr(XL, 0);
P1 = Xtr(Ytr == 1, :) * W; P2 = Xtr(Ytr == 2, :) * W;
proj = ((Z(test, :) - mu2) ./ sd2) * W;
LO = ar_log_kde(P1, proj, bw) - ar_log_kde(P2, proj, bw);
proj = single(proj);
end
