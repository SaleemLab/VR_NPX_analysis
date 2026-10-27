function R = ar_grey_lag(varargin)
% Lead / lag between HPC and V1 track bias using only bins that lie inside the
% grey screen and away from its visual events: the first P_offset s after the
% lap end (V1 latency and offset response to the screen turning grey) are
% excluded, and bins end at track onset.
% Bins: the whole-session 100 ms grid (ar_bin_session kind 0), kept if they lie
% fully within [t_off(k-1) + offset, t_on(k)] of a valid transition.
% Windows (per trial): 'whole' (offset s after lap end to onset), 'post_offset'
% (offset to 1.5 s after lap end), 'late' (last 1.5 s before onset).
% Per session, region pair and decoder:
%   residual = zLO minus the mean time course (from lap end) of trials with the
%   same previous x upcoming track pair, minus the trial mean; optionally
%   ('ctrl') with log spike count and speed regressed out per region first
%   lagged correlation r(L) = corr(HPC(t), V1(t + L)) over bin pairs within the
%   same trial and window (L > 0: V1 follows HPC), L = -5..5 bins (+-500 ms)
%   null: HPC paired with V1 from another trial of the same track pair
%   asymmetry = mean r(L = 1..3) - mean r(L = -3..-1) (> 0: HPC leads)
% Output: results\grey_lag\ (mat, csv, figure)
p = inputParser;
addParameter(p, 'decoders', {'bayes_multinomial', 'bayes_multinomial_all', 'pls'});
addParameter(p, 'offset', 0.3);
addParameter(p, 'nShuffle', 200);
addParameter(p, 'plot_only', false);   % redraw the figure from the saved results
parse(p, varargin{:});
O = p.Results;
P = ar_setup();
Lt = load(fullfile(P.results, 'data', 'lap_tables.mat'), 'L'); L = Lt.L;
out_dir = fullfile(P.results, 'grey_lag');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
if O.plot_only
    S = load(fullfile(out_dir, 'grey_lag.mat'), 'R', 'lags', 'O'); R = S.R;
    summarise(R, S.lags, P, S.O, out_dir);
    return
end
bad = find(ismember(P.subjects, P.V1_exclude_subjects));
lags = -5:5;
wins = {'whole', 'post_offset', 'late'};
rows = {};
rng(13);
for s = 1:numel(L)
    if ismember(L(s).animal_index, bad), continue; end
    f = fullfile(P.results, 'decoding', sprintf('session_%02d.mat', s));
    S = load(f, 'out'); out = S.out; G = out.G; T = L(s).laps;
    trans = find(T.valid_transition);
    grid = G.kind == 0;
    % grey bins per transition (grid bins fully inside the grey, after the offset)
    bins = cell(numel(trans), 1); tfe = bins; tbo = bins;
    for i = 1:numel(trans)
        k = trans(i);
        b = find(grid & G.t_start >= T.t_off(k - 1) + O.offset & G.t_start + P.bin_width <= T.t_on(k));
        bins{i} = b;
        tfe{i} = G.t_start(b) - T.t_off(k - 1);               % from lap end
        tbo{i} = T.t_on(k) - (G.t_start(b) + P.bin_width);    % to onset (>= 0)
    end
    pair = T.prev_track(trans) * 10 + T.track(trans);
    rH = find(strcmp({out.region.name}, 'HPC')); rV = find(strcmp({out.region.name}, 'V1'));
    for d = 1:numel(O.decoders)
        [zH, okH] = zlo(out.region(rH), O.decoders{d}, grid);
        [zV, okV] = zlo(out.region(rV), O.decoders{d}, grid);
        if ~okH || ~okV, continue; end
        for ctrl = [false true]
            h = zH; v = zV;
            if ctrl
                h = regress_out(h, G, 'HPC', grid); v = regress_out(v, G, 'V1', grid);
            end
            for w = 1:numel(wins)
                sel = cellfun(@(a, e) select_window(a, e, wins{w}), tfe, tbo, 'UniformOutput', false);
                [Hs, Vs] = deal(cell(numel(trans), 1));
                for i = 1:numel(trans)
                    Hs{i} = h(bins{i}(sel{i})); Vs{i} = v(bins{i}(sel{i}));
                end
                [Hs, Vs] = residualise(Hs, Vs, pair);
                rr = lagged(Hs, Vs, lags);
                nul = NaN(O.nShuffle, numel(lags));
                for sh = 1:O.nShuffle
                    perm = (1:numel(trans))';
                    for g = unique(pair)'
                        ix = find(pair == g); perm(ix) = ix(randperm(numel(ix)));
                    end
                    nul(sh, :) = lagged(Hs, Vs(perm), lags);
                end
                asym = mean(rr(lags >= 1 & lags <= 3)) - mean(rr(lags <= -1 & lags >= -3));
                asym_nul = mean(nul(:, lags >= 1 & lags <= 3), 2) - mean(nul(:, lags <= -1 & lags >= -3), 2);
                rows(end + 1, :) = {string(O.decoders{d}), ctrl, string(wins{w}), s, L(s).animal_index, ...
                    median(cellfun(@numel, Hs)), rr, mean(nul, 1), std(nul, 0, 1), asym, mean(asym_nul), std(asym_nul)}; %#ok<AGROW>
            end
        end
    end
    fprintf('Session %d done\n', s);
end
R = cell2table(rows, 'VariableNames', {'decoder', 'ctrl', 'window', 'session', 'animal', 'median_bins', 'r', ...
    'null_mean', 'null_sd', 'asym', 'asym_null_mean', 'asym_null_sd'});
save(fullfile(out_dir, 'grey_lag.mat'), 'R', 'lags', 'O');
summarise(R, lags, P, O, out_dir);
end

%% ------------------------------------------------------------------------
function [z, ok] = zlo(Rg, dec, grid)
parts = regexp(dec, '^(bayes_multinomial|bayes_poisson|pls)_?(lobo5|all)?$', 'tokens', 'once');
if isempty(Rg.(parts{1}).LO), z = []; ok = false; return; end
switch parts{2}
    case '', LO = double(Rg.(parts{1}).LO);
    case 'lobo5', LO = double(Rg.(parts{1}).LO5);
    case 'all', LO = double(Rg.(parts{1}).LO_all);
end
r = isfinite(LO) & grid;
z = (LO - mean(LO(r))) / std(LO(r));
lim = prctile(z(r), [0.5 99.5]);
z = min(max(z, lim(1)), lim(2));
ok = true;
end

function z = regress_out(z, G, region, grid)
% remove the part of the bin log odds explained by log spike count and speed
X = [ones(numel(z), 1) log1p(double(G.n_spikes.(region))) G.speed];
r = grid & all(isfinite(X), 2) & isfinite(z);
b = X(r, :) \ z(r);
z = z - X(:, 2:3) * b(2:3);
end

function s = select_window(tfe, tbo, w)
switch w
    case 'whole', s = true(size(tfe));
    case 'post_offset', s = tfe < 1.5;
    case 'late', s = tbo < 1.5;
end
end

function [Hs, Vs] = residualise(Hs, Vs, pair)
% remove the mean time course (aligned to the first bin, i.e. from lap end for
% 'whole' / 'post_offset' and from the window start for 'late') per track pair,
% then each trial's mean
for g = unique(pair)'
    ix = find(pair == g);
    n = max(cellfun(@numel, Hs(ix)));
    if n == 0, continue; end
    mh = zeros(n, 1); mv = zeros(n, 1); c = zeros(n, 1);
    for i = ix'
        m = numel(Hs{i});
        mh(1:m) = mh(1:m) + Hs{i}; mv(1:m) = mv(1:m) + Vs{i}; c(1:m) = c(1:m) + 1;
    end
    mh = mh ./ c; mv = mv ./ c;
    for i = ix'
        m = numel(Hs{i});
        Hs{i} = Hs{i} - mh(1:m); Vs{i} = Vs{i} - mv(1:m);
    end
end
for i = 1:numel(Hs)
    Hs{i} = Hs{i} - mean(Hs{i}); Vs{i} = Vs{i} - mean(Vs{i});
end
end

function r = lagged(Hs, Vs, lags)
% r(L) = corr(H(t), V(t + L)) over bin pairs within the same trial; trials with
% unequal lengths (shuffled pairs) use their common length
r = NaN(1, numel(lags));
for k = 1:numel(lags)
    Lg = lags(k); a = []; b = [];
    for i = 1:numel(Hs)
        n = min(numel(Hs{i}), numel(Vs{i}));
        if n <= abs(Lg) + 2, continue; end
        h = Hs{i}(1:n); v = Vs{i}(1:n);
        if Lg >= 0, a = [a; h(1:n - Lg)]; b = [b; v(1 + Lg:n)]; %#ok<AGROW>
        else, a = [a; h(1 - Lg:n)]; b = [b; v(1:n + Lg)]; end %#ok<AGROW>
    end
    ok = isfinite(a) & isfinite(b);
    r(k) = corr(a(ok), b(ok));
end
end

%% ------------------------------------------------------------------------
function summarise(R, lags, P, O, out_dir)
acol = lines(6);
wins = {'whole', 'post_offset', 'late'};
wl = {sprintf('whole grey (%.1f s after end to onset)', O.offset), sprintf('%.1f-1.5 s after lap end', O.offset), 'last 1.5 s before onset'};
stats = {};
fig = figure('Position', [50 50 1700 760], 'Color', 'w');
tl = tiledlayout(2, numel(wins) + 1, 'TileSpacing', 'compact');
for c = [false true]
    for w = 1:numel(wins)
        nexttile(tl); hold on;
        for d = 1:numel(O.decoders)
            x = R(R.decoder == O.decoders{d} & R.window == wins{w} & R.ctrl == c, :);
            if isempty(x), continue; end
            dr = x.r - x.null_mean;
            m = mean(dr, 1, 'omitnan'); se = std(dr, 0, 1, 'omitnan') / sqrt(height(x));
            col = [0 0 0; 0.55 0.55 0.55; 0.2 0.45 0.75]; col = col(d, :);
            errorbar(lags * 100 + (d - 2) * 8, m, se, 'o-', 'Color', col, 'MarkerFaceColor', col, 'CapSize', 0, 'LineWidth', 1.2);
            for k = 1:numel(lags)
                stats(end + 1, :) = {string(O.decoders{d}), c, string(wins{w}), lags(k) * 100, m(k), se(k), ...
                    signrank(dr(:, k)), sum(dr(:, k) > 0), height(x)}; %#ok<AGROW>
            end
        end
        xline(0, 'k:'); yline(0, 'k:');
        xlabel('Lag (ms; > 0: V1 follows HPC)'); ylabel('HPC-V1 r minus shuffled-trial null');
        ttl = {wl{w}, 'as decoded'}; if c, ttl{2} = 'spikes & speed removed'; end
        title(ttl, 'FontWeight', 'normal');
        if w == 1 && ~c, legend(strrep(O.decoders, '_', ' '), 'Box', 'off', 'Location', 'southwest'); end
    end
    % asymmetry per session (main decoder)
    nexttile(tl); hold on;
    for w = 1:numel(wins)
        x = R(R.decoder == O.decoders{1} & R.window == wins{w} & R.ctrl == c, :);
        a = x.asym - x.asym_null_mean;
        scatter(w + 0.1 * randn(height(x), 1), a, 20, acol(x.animal, :), 'filled', 'MarkerFaceAlpha', 0.7);
        errorbar(w + 0.3, mean(a), std(a) / sqrt(numel(a)), 'ko', 'MarkerFaceColor', 'k', 'CapSize', 0, 'LineWidth', 1.5);
        [ua, ~, ia] = unique(x.animal); am = accumarray(ia, a, [], @mean);
        text(w, max(a) + 0.01, sprintf('p %.2g\n%d/%d anim', signrank(a), sum(am > 0), numel(am)), 'HorizontalAlignment', 'center', 'FontSize', 7);
        stats(end + 1, :) = {string(O.decoders{1}), c, string(wins{w}), NaN, mean(a), std(a) / sqrt(numel(a)), signrank(a), sum(a > 0), numel(a)}; %#ok<AGROW>
    end
    yline(0, 'k:');
    set(gca, 'XTick', 1:numel(wins), 'XTickLabel', {'whole', 'post-offset', 'late'}); xlim([0.5 numel(wins) + 0.5]);
    ylabel('Asymmetry (HPC leads > 0), minus null');
    title({sprintf('%s asymmetry', strrep(O.decoders{1}, '_', ' ')), '(+100..300 vs -100..-300 ms)'}, 'FontWeight', 'normal');
end
title(tl, 'HPC-V1 lead/lag inside the grey screen only (mean +- SEM over sessions; V1 animals)');
ar_save_figure(fig, fullfile(out_dir, 'grey_lag'));
close(fig);
T = cell2table(stats, 'VariableNames', {'decoder', 'ctrl', 'window', 'lag_ms', 'mean_minus_null', 'sem', 'p_signrank', 'n_pos', 'n'});
writetable(T, fullfile(out_dir, 'grey_lag_stats.csv'));
end
