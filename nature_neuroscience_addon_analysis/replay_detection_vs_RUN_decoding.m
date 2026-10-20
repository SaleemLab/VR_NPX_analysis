function replay_detection_vs_RUN_decoding(varargin)
% Read-only: per-session replay detection (stored cell-ID shuffle scores of
% run_Bayesian_replay_decoding_POST) vs RUN1 10-fold cross-validated decoding
% (run_Bayesian_RUN_decoding_CV, multinomial, 100 ms), and the number of time bins of the
% ripple scoring windows.
% Replay measures per session and region (scored events, alpha = 0.05, either track):
%   rate post, rate pre, rate post - pre, mean log odds of T1 - T2 significant events (post,
%   multi-track events in both), median shuffle z (post, better track), n scored.
% RUN measures: track accuracy (posterior > 0.5 on the true track), mean soft position error
%   with the own-track template (unweighted mean over true positions, from template_soft),
%   median error (own template), number of decoder cells.
% Bins: distributions of window length (n_post_bins, all ripples) and non-empty bins /
%   active cells of scored events; cell-ID replay rate by number of non-empty bins.
% Options: 'replay_dir' (Bayesian_replay_decoding_rate_lt10), 'run_dir'
% (Bayesian_RUN_decoding_rate_lt10), 'method' ('multinomial'). Output: replay_dir\run_cv_vs_replay.

p = inputParser;
root = 'D:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
addParameter(p, 'replay_dir', fullfile(root, 'Bayesian_replay_decoding_rate_lt10'), @ischar);
addParameter(p, 'run_dir', fullfile(root, 'Bayesian_RUN_decoding_rate_lt10'), @ischar);
addParameter(p, 'method', 'multinomial', @ischar);
parse(p, varargin{:});
P = p.Results;
out_dir = fullfile(P.replay_dir, 'run_cv_vs_replay');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
regions = {'HPC', 'V1'}; windows = {'post', 'pre'};

% ---- replay sessions: per-session measures and per-event bin counts ----
files = dir(fullfile(P.replay_dir, 'M*', 'session_*.mat'));
S = table(); EV = table();
for k = 1:numel(files)
    L = load(fullfile(files(k).folder, files(k).name), 'out'); out = L.out;
    fprintf('%s %s\n', out.subject, out.session);
    row = table(string(out.subject), string(out.session), out.session_index, 'VariableNames', {'subject', 'session', 'session_index'});
    nEv = numel(out.events.onset);
    ev = table(repmat(out.session_index, nEv, 1), out.events.n_post_bins, out.events.duration, ...
        'VariableNames', {'session_index', 'n_post_bins', 'duration'});
    for nr = 1:numel(out.region)
        R = out.region(nr); reg = R.name;
        row.([reg '_n_cells']) = R.n_cells;
        for w = 1:2
            X = R.(windows{w}); ok = X.scored;
            pmin = min(X.p, [], 2);
            pre = [reg '_' windows{w} '_'];
            row.([pre 'rate']) = mean(pmin(ok) < 0.05);
            row.([pre 'n_scored']) = sum(ok);
            ev.([pre 'scored']) = ok; ev.([pre 'n_nonempty']) = X.n_nonempty; ev.([pre 'n_active']) = X.n_active;
            ev.([pre 'sig']) = ok & pmin < 0.05;
            if w == 1
                s1 = ok & X.p(:, 1) < 0.05 & isfinite(X.logodds); s2 = ok & X.p(:, 2) < 0.05 & isfinite(X.logodds);
                row.([pre 'logodds_diff']) = mean(X.logodds(s1)) - mean(X.logodds(s2));
                row.([pre 'median_z']) = median(max(X.z(ok, :), [], 2), 'omitnan');
            end
        end
        row.([reg '_rate_diff']) = row.([reg '_post_rate']) - row.([reg '_pre_rate']);
    end
    S = [S; row]; %#ok<AGROW>
    EV = [EV; ev]; %#ok<AGROW>
end

% ---- RUN cross-validated decoding ----
M = readtable(fullfile(P.run_dir, 'metrics_100ms.csv'), 'TextType', 'string');
M = M(M.method == P.method, :);
for reg = regions
    acc = NaN(height(S), 1); med = acc; soft = acc;
    for k = 1:height(S)
        r = M(M.subject == S.subject(k) & string(M.session) == S.session(k) & M.region == reg{1}, :);
        if height(r) ~= 1, continue; end
        acc(k) = r.track_acc; med(k) = mean([r.err_T1data_T1tmpl r.err_T2data_T2tmpl]);
        f = fullfile(P.run_dir, P.method, S.subject(k), sprintf('session_%s_100ms.mat', S.session(k)));
        if exist(f, 'file')
            L = load(f, 'out');
            Rr = L.out.region(strcmp({L.out.region.name}, reg{1}));
            if ~isempty(Rr), soft(k) = soft_error(Rr.metrics, 10); end
        end
    end
    S.([reg{1} '_run_track_acc']) = acc; S.([reg{1} '_run_median_err']) = med; S.([reg{1} '_run_soft_err']) = soft;
end
writetable(S, fullfile(out_dir, 'session_run_cv_vs_replay.csv'));
save(fullfile(out_dir, 'run_cv_vs_replay.mat'), 'S', 'EV', 'P');

plot_cv_vs_replay(S, regions, out_dir);
plot_bins(EV, regions, out_dir);
end

function e = soft_error(Mt, bin_cm)
% Mean expected |decoded - true| position (cm), own-track templates, unweighted over positions
e = NaN(1, 2);
for a = 1:2
    C = Mt.template_soft{a, a};                 % rows true pos, columns decoded (soft)
    n = size(C, 1);
    D = abs((1:n)' - (1:n)) * bin_cm;
    e(a) = mean(sum(C .* D, 2), 'omitnan');
end
e = mean(e);
end

function plot_cv_vs_replay(S, regions, out_dir)
animals = unique(S.subject); cmap = lines(numel(animals)); [~, ia] = ismember(S.subject, animals);
xs = {'run_track_acc', 'run_soft_err', 'n_cells'};
xl = {'RUN CV track accuracy', 'RUN CV mean position error, own track (cm)', 'decoder cells'};
ys = {'post_rate', 'pre_rate', 'rate_diff', 'post_logodds_diff', 'post_median_z'};
yl = {'replay rate post (either track)', 'replay rate pre', 'rate post - pre', ...
    'log odds T1 - T2 sig. events (post)', 'median shuffle z (post)'};
for nr = 1:numel(regions)
    reg = regions{nr};
    fig = figure('Color', 'w', 'Position', [30 30 1700 950]);
    tl = tiledlayout(numel(xs), numel(ys), 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:numel(xs)
        x = S.([reg '_' xs{i}]);
        for j = 1:numel(ys)
            nexttile; hold on
            y = S.([reg '_' ys{j}]);
            scatter(x, y, 30, cmap(ia, :), 'filled');
            for k = 1:height(S), text(x(k), y(k), sprintf(' %d', S.session_index(k)), 'FontSize', 6); end
            if any(strcmp(ys{j}, {'post_rate', 'pre_rate'})), yline(1 - 0.95^2, 'k--'); end
            if any(strcmp(ys{j}, {'rate_diff', 'post_logodds_diff', 'post_median_z'})), yline(0, 'k:'); end
            ok = isfinite(x) & isfinite(y);
            [r, pv] = corr(x(ok), y(ok), 'Type', 'Spearman');
            title(sprintf('Spearman r = %.2f, p = %.3g', r, pv), 'FontWeight', 'normal', 'FontSize', 8);
            if i == numel(xs), xlabel(xl{i}); else, xlabel(xl{i}, 'FontSize', 7); end
            if j == 1 || true, ylabel(yl{j}, 'FontSize', 7); end
            set(gca, 'TickDir', 'out', 'FontSize', 7);
        end
    end
    lg = arrayfun(@(a) sprintf('\\color[rgb]{%.2f %.2f %.2f}%s', cmap(a, :), animals(a)), 1:numel(animals), 'UniformOutput', false);
    title(tl, {sprintf('%s: RUN1 cross-validated decoding (multinomial, 100 ms, 10-fold) vs POST ripple replay (cell-ID shuffle, alpha 0.05)', reg), ...
        ['sessions = dots (number = session index), colour = animal: ' strjoin(lg, '  ')]}, 'FontSize', 9);
    save_fig(fig, fullfile(out_dir, sprintf('run_cv_vs_replay_%s', reg)));
    close(fig);
end
end

function plot_bins(EV, regions, out_dir)
fig = figure('Color', 'w', 'Position', [30 30 1700 800]);
tl = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
cols = struct('post', [0.85 0.33 0.1], 'pre', [0.4 0.4 0.4]);
for nr = 1:numel(regions)
    reg = regions{nr};
    % window length of all ripples (same for both regions)
    nexttile; hold on
    histogram(EV.n_post_bins, 0.5:1:25.5, 'FaceColor', [0.6 0.6 0.6]);
    sc = EV.([reg '_post_scored']);
    histogram(EV.n_post_bins(sc), 0.5:1:25.5, 'FaceColor', cols.post);
    xlabel('window length (20 ms bins) = max(0.1 s, ripple duration)'); ylabel('ripples');
    title(sprintf('all ripples (n = %d), %s scored post (n = %d, %.0f%%)', height(EV), reg, sum(sc), 100 * mean(sc)), ...
        'FontWeight', 'normal', 'FontSize', 8);
    legend({'all', 'scored'}, 'Box', 'off');
    % non-empty bins of scored events
    nexttile; hold on
    for w = {'pre', 'post'}
        s = EV.([reg '_' w{1} '_scored']);
        histogram(EV.([reg '_' w{1} '_n_nonempty'])(s), 0.5:1:25.5, 'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', 'EdgeColor', cols.(w{1}), 'LineWidth', 1.5);
    end
    s = EV.([reg '_post_scored']); nn = EV.([reg '_post_n_nonempty'])(s);
    xlabel('non-empty bins (scored)'); ylabel('fraction of events');
    title(sprintf('%s post: median %d, <= 4 bins %.0f%%, >= 6 bins %.0f%%', reg, median(nn), 100 * mean(nn <= 4), 100 * mean(nn >= 6)), ...
        'FontWeight', 'normal', 'FontSize', 8);
    legend({'pre', 'post'}, 'Box', 'off');
    % active cells
    nexttile; hold on
    for w = {'pre', 'post'}
        s = EV.([reg '_' w{1} '_scored']);
        histogram(EV.([reg '_' w{1} '_n_active'])(s), 0.5:1:40.5, 'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', 'EdgeColor', cols.(w{1}), 'LineWidth', 1.5);
    end
    na = EV.([reg '_post_n_active'])(EV.([reg '_post_scored']));
    xlabel('active decoder cells (scored)'); ylabel('fraction of events');
    title(sprintf('%s post: median %d active cells', reg, median(na)), 'FontWeight', 'normal', 'FontSize', 8);
    % replay rate by non-empty bins
    nexttile; hold on
    edges = [3 4 5 6 7 9 12 inf]; lab = {'3', '4', '5', '6', '7-8', '9-11', '12+'};
    for w = {'pre', 'post'}
        s = EV.([reg '_' w{1} '_scored']);
        g = discretize(EV.([reg '_' w{1} '_n_nonempty'])(s), edges);
        sig = double(EV.([reg '_' w{1} '_sig'])(s));
        rate = accumarray(g, sig, [numel(lab) 1], @mean, NaN); n = accumarray(g, 1, [numel(lab) 1]);
        plot(1:numel(lab), rate, '-o', 'Color', cols.(w{1}), 'MarkerFaceColor', cols.(w{1}));
        if strcmp(w{1}, 'post')
            for i = 1:numel(lab), text(i, rate(i), sprintf('  n=%d', n(i)), 'FontSize', 6); end
        end
    end
    yline(1 - 0.95^2, 'k--', 'nominal 0.0975', 'FontSize', 7);
    set(gca, 'XTick', 1:numel(lab), 'XTickLabel', lab); xlabel('non-empty bins'); ylabel('cell-ID replay rate (either track)');
    title(sprintf('%s: replay rate by non-empty bins', reg), 'FontWeight', 'normal', 'FontSize', 8);
    legend({'pre', 'post'}, 'Box', 'off', 'Location', 'northwest');
end
title(tl, 'Ripple scoring windows: length, non-empty bins and active cells (scored = >= 5 active cells and >= 3 non-empty bins)', 'FontSize', 9);
save_fig(fig, fullfile(out_dir, 'ripple_bins'));
close(fig);
end

function save_fig(fig, f)
savefig(fig, [f '.fig']);
exportgraphics(fig, [f '.png'], 'Resolution', 150);
exportgraphics(fig, [f '.pdf'], 'ContentType', 'vector');
end
