function diagnose_replay_chance(varargin)
% Read-only diagnostics of the replay significance rates of run_Bayesian_replay_decoding_POST
% (why are HPC replay rates below the nominal 'either track' chance of 1 - 0.95^2?).
% Uses only the stored session outputs (posteriors and scores); nothing is re-decoded or
% overwritten. Per scored event (pre and post window, one region):
%   cell-ID shuffle p (stored), |r| vs shuffle mean |r| and z (stored), track dependence
%   time-bin permutation null on the stored posterior (non-empty bins permuted among the
%   non-empty bin positions; exhaustive for <= 6 bins, else nPerm random): permutation p
%   per track, and the empirical chance of p < 0.05 per track / on either track (fraction of
%   null replicates that would pass), which accounts for the discreteness of the null and
%   for the dependence between the two track tests (same spikes)
%   positional spread of the per-track posterior (weighted position SD over the window;
%   low = stationary / coherent location)
% Options: 'out_dir' (Bayesian_replay_decoding_rate_lt10), 'region' ('HPC'), 'nPerm' (1000),
% 'seed' (2026). Output: out_dir\chance_diagnostics (figures, event table, summary CSV).

p = inputParser;
addParameter(p, 'out_dir', 'D:\corticohippocampal_replay\nature_neuroscience_revision_analysis\Bayesian_replay_decoding_rate_lt10', @ischar);
addParameter(p, 'region', 'HPC', @ischar);
addParameter(p, 'nPerm', 1000, @isnumeric);
addParameter(p, 'seed', 2026, @isnumeric);
parse(p, varargin{:});
P = p.Results;
diag_dir = fullfile(P.out_dir, 'chance_diagnostics');
if ~exist(diag_dir, 'dir'), mkdir(diag_dir); end

files = dir(fullfile(P.out_dir, 'M*', 'session_*.mat'));
windows = {'post', 'pre'};
T = table();
for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name), 'out'); out = S.out;
    R = out.region(strcmp({out.region.name}, P.region));
    prm = out.params; nPos = prm.nPos;
    fprintf('%s %s: %d events\n', out.subject, out.session, numel(out.events.onset));
    for w = 1:2
        X = R.(windows{w});
        ev = find(X.scored);
        rs = RandStream('mrg32k3a', 'Seed', P.seed + 100 * out.session_index + w);
        n = numel(ev);
        C = struct('wc_check', NaN(n, 2), 'p_time', NaN(n, 2), 'chance_track', NaN(n, 2), ...
            'chance_either', NaN(n, 1), 'spread', NaN(n, 2), 'n_perm', NaN(n, 1));
        for i = 1:n
            e = ev(i);
            nb = out.events.n_post_bins(e);
            if w == 1, bins = prm.onset_bin + (0:nb - 1); else, bins = prm.onset_bin - nb + (0:nb - 1); end
            post = double(R.posterior(:, bins, e));
            tb = find(all(isfinite(post), 1));
            post = post(:, tb);
            nt = numel(tb);
            if nt <= 6
                perm_idx = perms(1:nt);
            else
                [~, perm_idx] = sort(rand(rs, P.nPerm, nt), 2);
                perm_idx = [1:nt; perm_idx]; %#ok<AGROW>   % real order = replicate 1
            end
            nP = size(perm_idx, 1);
            r_null = NaN(nP, 2); p_rep = NaN(nP, 2);
            for tr = 1:2
                Pk = post((tr - 1) * nPos + (1:nPos), :);
                Pk = Pk ./ sum(Pk, 1);
                C.wc_check(i, tr) = wcorr_batch(reshape(Pk, [1 nPos nt]), tb);
                W = permute(reshape(Pk(:, perm_idx(:)), nPos, nP, nt), [2 1 3]);
                r_null(:, tr) = abs(wcorr_batch(W, tb));
                r_null(~isfinite(r_null(:, tr)), tr) = 0;
                C.p_time(i, tr) = mean(r_null(:, tr) >= abs(C.wc_check(i, tr)) - 1e-12);
                % p of every null replicate within the same null distribution
                [~, ~, rk] = unique(-round(r_null(:, tr), 12));
                p_grp = cumsum(accumarray(rk, 1)) / nP;
                p_rep(:, tr) = p_grp(rk);
                C.chance_track(i, tr) = mean(p_rep(:, tr) < 0.05);
                x = (1:nPos)';
                wx = sum(Pk, 2) / sum(Pk(:));
                C.spread(i, tr) = sqrt(sum(wx .* (x - sum(wx .* x)).^2));
            end
            C.chance_either(i) = mean(min(p_rep, [], 2) < 0.05);
            C.n_perm(i) = nP;
        end
        Tk = table(repmat(string(out.subject), n, 1), repmat(out.session_index, n, 1), repmat(string(windows{w}), n, 1), ...
            ev, X.n_active(ev), X.n_nonempty(ev), X.n_spikes(ev), X.wc(ev, :), C.wc_check, X.p(ev, :), X.z(ev, :), ...
            X.shuffle_mean(ev, :), X.shuffle_sd(ev, :), C.p_time, C.chance_track, C.chance_either, C.spread, C.n_perm, ...
            'VariableNames', {'subject', 'session_index', 'window', 'event', 'n_active', 'n_nonempty', 'n_spikes', 'wc', ...
            'wc_check', 'p', 'z', 'shuffle_mean', 'shuffle_sd', 'p_time', 'chance_track', 'chance_either', 'spread', 'n_perm'});
        T = [T; Tk]; %#ok<AGROW>
    end
end
fprintf('max |wc stored - wc from stored posterior| = %.2g\n', max(abs(T.wc(:) - T.wc_check(:)), [], 'omitnan'));
save(fullfile(diag_dir, sprintf('%s_chance_diagnostics_events.mat', P.region)), 'T', 'P', '-v7.3');

% ---- pooled summary ----
lines = {};
for w = windows
    Tw = T(T.window == w{1}, :);
    pm = min(Tw.p, [], 2);
    a1 = mean(Tw.p(:, 1) < 0.05); a2 = mean(Tw.p(:, 2) < 0.05);
    lines{end+1} = sprintf(['%s %s (n = %d scored): cell-ID shuffle sig T1 %.3f, T2 %.3f, either %.3f ' ...
        '(independence prediction from marginals %.3f) | Spearman(p1, p2) = %.2f, P(p2<.05 | p1<.05) = %.2f'], ...
        P.region, w{1}, height(Tw), a1, a2, mean(pm < 0.05), 1 - (1 - a1) * (1 - a2), ...
        corr(Tw.p(:, 1), Tw.p(:, 2), 'Type', 'Spearman'), mean(Tw.p(Tw.p(:, 1) < 0.05, 2) < 0.05)); %#ok<AGROW>
    lines{end+1} = sprintf(['   |r| < shuffle mean |r|: %.3f of tests | median z %.2f | mean p %.3f (uniform 0.5) | ' ...
        'p > 0.95: %.3f, p < 0.05: %.3f'], mean(abs(Tw.wc(:)) < Tw.shuffle_mean(:)), median(Tw.z(:), 'omitnan'), ...
        mean(Tw.p(:)), mean(Tw.p(:) > 0.95), mean(Tw.p(:) < 0.05)); %#ok<AGROW>
    lines{end+1} = sprintf(['   time-bin perm null: empirical chance per track %.3f, either %.3f | observed time-perm ' ...
        'sig either %.3f | events with <= 4 non-empty bins (cannot reach p<.05 by time perm): %.3f'], ...
        mean(Tw.chance_track(:)), mean(Tw.chance_either), mean(min(Tw.p_time, [], 2) < 0.05), mean(Tw.n_nonempty <= 4)); %#ok<AGROW>
end
fprintf('\n%s\n', strjoin(lines, newline));
fid = fopen(fullfile(diag_dir, sprintf('%s_chance_diagnostics_summary.txt', P.region)), 'w');
fprintf(fid, '%s\n', strjoin(lines, newline)); fclose(fid);

% per-session table
sess = unique(T.session_index);
M = table(sess, 'VariableNames', {'session_index'});
for w = windows
    for k = 1:numel(sess)
        Tw = T(T.window == w{1} & T.session_index == sess(k), :);
        M.("rate_either_" + w{1})(k) = mean(min(Tw.p, [], 2) < 0.05);
        M.("rate_T1_" + w{1})(k) = mean(Tw.p(:, 1) < 0.05);
        M.("rate_T2_" + w{1})(k) = mean(Tw.p(:, 2) < 0.05);
        M.("rate_time_either_" + w{1})(k) = mean(min(Tw.p_time, [], 2) < 0.05);
        M.("chance_either_" + w{1})(k) = mean(Tw.chance_either);
        M.("median_z_" + w{1})(k) = median(Tw.z(:), 'omitnan');
        M.("frac_below_shuffle_mean_" + w{1})(k) = mean(abs(Tw.wc(:)) < Tw.shuffle_mean(:));
        M.("median_nonempty_" + w{1})(k) = median(Tw.n_nonempty);
    end
end
writetable(M, fullfile(diag_dir, sprintf('%s_chance_diagnostics_sessions.csv', P.region)));

plot_diagnostics(T, M, P, diag_dir);
end

function plot_diagnostics(T, M, P, diag_dir)
col = struct('post', [0.85 0.33 0.1], 'pre', [0.3 0.3 0.3]);
windows = {'pre', 'post'};
fig = figure('Color', 'w', 'Position', [50 50 1500 820]);
tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s replay significance diagnostics (cell-ID shuffle, stored scores) - %s', P.region, ...
    P.out_dir), 'FontSize', 10, 'Interpreter', 'none');
edges = 0:0.05:1;

% a: p per track (pooled)
nexttile; hold on
for w = windows
    Tw = T(T.window == w{1}, :);
    histogram(Tw.p(:), edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', col.(w{1}), 'LineWidth', 1.5);
end
yline(1, 'k--'); xlabel('cell-ID shuffle p (per track)'); ylabel('density');
legend(windows, 'Location', 'north', 'Box', 'off'); title('p per track (uniform = 1)', 'FontWeight', 'normal');

% b: min p over tracks
nexttile; hold on
for w = windows
    Tw = T(T.window == w{1}, :);
    histogram(min(Tw.p, [], 2), edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', col.(w{1}), 'LineWidth', 1.5);
end
x = linspace(0, 1, 100); plot(x, 2 * (1 - x), 'k--'); plot(x, ones(size(x)), 'k:');
xlabel('min(p_{T1}, p_{T2})'); ylabel('density');
title({'min p over tracks', '-- independent nulls, : fully dependent'}, 'FontWeight', 'normal');

% c: p1 vs p2
nexttile; hold on
Tw = T(T.window == "post", :);
scatter(Tw.p(:, 1), Tw.p(:, 2), 4, col.post, 'filled', 'MarkerFaceAlpha', 0.2);
xlabel('p T1'); ylabel('p T2'); axis square; xlim([0 1]); ylim([0 1]);
title(sprintf('post: Spearman r = %.2f', corr(Tw.p(:, 1), Tw.p(:, 2), 'Type', 'Spearman')), 'FontWeight', 'normal');

% d: |r| vs shuffle mean |r|
nexttile; hold on
for w = windows
    Tw = T(T.window == w{1}, :);
    scatter(Tw.shuffle_mean(:), abs(Tw.wc(:)), 4, col.(w{1}), 'filled', 'MarkerFaceAlpha', 0.15);
end
plot([0 1], [0 1], 'k--'); axis square; xlim([0 1]); ylim([0 1]);
xlabel('shuffle mean |r|'); ylabel('real |r|');
Tw = T(T.window == "post", :);
title(sprintf('post: %.0f%% of tests below shuffle mean', 100 * mean(abs(Tw.wc(:)) < Tw.shuffle_mean(:))), 'FontWeight', 'normal');

% e: z distribution
nexttile; hold on
for w = windows
    Tw = T(T.window == w{1}, :);
    histogram(Tw.z(:), -4:0.2:6, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', col.(w{1}), 'LineWidth', 1.5);
end
xline(0, 'k--'); xlabel('shuffle z'); ylabel('density');
title(sprintf('median z: pre %.2f, post %.2f', median(T.z(T.window == "pre", :), 'all', 'omitnan'), ...
    median(T.z(T.window == "post", :), 'all', 'omitnan')), 'FontWeight', 'normal');

% f: rates by number of non-empty bins
nexttile; hold on
grp = [3 4 5 6 7 9 12 inf]; lab = {'3', '4', '5', '6', '7-8', '9-11', '12+'};
for w = windows
    Tw = T(T.window == w{1}, :);
    g = discretize(Tw.n_nonempty, grp);
    rate = splitapply(@(p) mean(min(p, [], 2) < 0.05), Tw.p, g);
    chance = splitapply(@mean, Tw.chance_either, g);
    gi = unique(g);
    plot(gi, rate, '-o', 'Color', col.(w{1}), 'MarkerFaceColor', col.(w{1}));
    plot(gi, chance, ':s', 'Color', col.(w{1}));
end
yline(1 - 0.95^2, 'k--', 'nominal 0.0975', 'FontSize', 7);
set(gca, 'XTick', 1:numel(lab), 'XTickLabel', lab); xlabel('non-empty bins in window');
ylabel('fraction sig, either track'); title({'solid = cell-ID rate', 'dotted = time-perm empirical chance'}, 'FontWeight', 'normal');

% g: z vs positional spread
nexttile; hold on
Tw = T(T.window == "post", :);
scatter(Tw.spread(:), Tw.z(:), 4, col.post, 'filled', 'MarkerFaceAlpha', 0.15);
yline(0, 'k--'); xlabel('posterior position SD (bins)'); ylabel('shuffle z');
ok = isfinite(Tw.z(:));
title(sprintf('post: Spearman r = %.2f (low SD = stationary)', corr(Tw.spread(ok), Tw.z(ok), 'Type', 'Spearman')), ...
    'FontWeight', 'normal');

% h: per-session rates vs empirical chance
nexttile; hold on
for w = windows
    xo = 1 + strcmp(w{1}, 'post');
    jit = xo + 0.15 * (rand(height(M), 1) - 0.5);
    plot(jit - 0.2, M.("rate_either_" + w{1}), 'o', 'Color', col.(w{1}));
    plot(jit + 0.2, M.("chance_either_" + w{1}), 's', 'Color', col.(w{1}));
end
yline(1 - 0.95^2, 'k--');
set(gca, 'XTick', [0.8 1.2 1.8 2.2], 'XTickLabel', {'pre rate', 'pre chance', 'post rate', 'post chance'});
ylabel('fraction (per session)'); title('cell-ID rate vs time-perm chance', 'FontWeight', 'normal');

save_figure(fig, fullfile(diag_dir, sprintf('%s_chance_diagnostics', P.region)));
end

function r = wcorr_batch(W, t)
% Weighted correlation between position and time (as in run_Bayesian_replay_decoding_POST)
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

function save_figure(fig, save_path)
savefig(fig, [save_path '.fig']);
exportgraphics(fig, [save_path '.png'], 'Resolution', 150);
exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
print(fig, [save_path '.svg'], '-dsvg', '-vector');
end
