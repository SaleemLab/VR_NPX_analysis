function replay_logodds_KDE_check(varargin)
% Read-only: proportion of significant replay events vs mean log odds difference (T1 - T2
% events) across alpha levels, as in the shuffle_comparison criterion figures of
% run_Bayesian_replay_decoding_POST, but with the track log odds taken from the PLS-KDE
% decoder (neurons < 10 Hz, KDE_reactivation_ripples_PSTH_POST_unique_rate_lt10.mat)
% instead of the Bayesian decoder used for detection.
% KDE event log odds = mean over the 10 ms PSTH bins (-1 to 1 s from onset, nan_mask
% applied) whose centres fall in the Bayesian scoring window: post = onset to onset +
% n_post_bins * 20 ms, pre = the same length before onset. Raw and z log odds (z per
% session and probe, as in the PSTH file). Columns per detecting region: Bayesian z log
% odds (reference), KDE raw, KDE z (same region), KDE z of the other region.
% Options: 'replay_dir' (Bayesian_replay_decoding_rate_lt10), 'criteria'
% ({{'cell_id'}, {'spike','field'}, {'field','time'}}), 'alphas', 'nBoot' (1000).
% Output: replay_dir\shuffle_comparison\KDE_logodds\.

p = inputParser;
rev = 'D:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
addParameter(p, 'replay_dir', fullfile(rev, 'Bayesian_replay_decoding_rate_lt10'), @ischar);
addParameter(p, 'kde_file', fullfile(rev, 'KDE_reactivation_ripples_PSTH_POST_unique_rate_lt10.mat'), @ischar);
addParameter(p, 'criteria', {{'cell_id'}, {'spike', 'field'}, {'field', 'time'}}, @iscell);
addParameter(p, 'alphas', [0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001], @isnumeric);
addParameter(p, 'nBoot', 1000, @isnumeric);
addParameter(p, 'seed', 2026, @isnumeric);
parse(p, varargin{:});
P = p.Results;
out_dir = fullfile(P.replay_dir, 'shuffle_comparison', 'KDE_logodds');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
regions = {'HPC', 'V1'}; windows = {'post', 'pre'};

L = load(fullfile(P.replay_dir, 'shuffle_comparison', 'replay_events_shuffles.mat'), 'T'); T = L.T;
L = load(fullfile(P.replay_dir, 'replay_events_all.mat'), 'T'); A = L.T;
assert(isequal(T.ripple_info_row, A.ripple_info_row), 'Event tables are not aligned');
T.n_post_bins = A.n_post_bins;
K = load(P.kde_file, 'KDE_reactivation_ripples_PSTH'); K = K.KDE_reactivation_ripples_PSTH;
bc = -1:0.01:1; bc = bc(1:end - 1) + 0.005;          % PSTH bin centres (s from onset)
mask = K.nan_mask;
for reg = regions
    for w = windows
        for m = {'logodds', 'z_logodds'}
            X = K.(sprintf('%s_%s_ripples', reg{1}, m{1})) + mask;
            X(~isfinite(X)) = NaN;
            v = NaN(height(T), 1);
            for nb = unique(T.n_post_bins)'
                if strcmp(w{1}, 'post'), cols = bc > 0 & bc < nb * 0.02; else, cols = bc > -nb * 0.02 & bc < 0; end
                rows = T.n_post_bins == nb;
                v(rows) = mean(X(T.ripple_info_row(rows), cols), 2, 'omitnan');
            end
            T.(sprintf('%s_%s_kde_%s', reg{1}, w{1}, m{1})) = v;
        end
    end
end
save(fullfile(out_dir, 'replay_events_KDE_logodds.mat'), 'T', 'P', '-v7.3');

% agreement between the Bayesian and KDE event log odds (scored post events)
agree = struct();
for reg = regions
    ok = T.([reg{1} '_post_scored']) & isfinite(T.([reg{1} '_post_logodds_z'])) & isfinite(T.([reg{1} '_post_kde_z_logodds']));
    agree.(reg{1}) = [corr(T.([reg{1} '_post_logodds_z'])(ok), T.([reg{1} '_post_kde_z_logodds'])(ok), 'Type', 'Spearman') ...
        mean(sign(T.([reg{1} '_post_logodds_z'])(ok)) == sign(T.([reg{1} '_post_kde_z_logodds'])(ok)))];
    fprintf('%s post: Bayesian vs KDE event z log odds Spearman r = %.2f, sign agreement %.1f%%\n', reg{1}, ...
        agree.(reg{1})(1), 100 * agree.(reg{1})(2));
end

for c = 1:numel(P.criteria)
    crit = P.criteria{c};
    fig = figure('Color', 'w', 'Position', [30 30 1900 820]);
    tl = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for nr = 1:2
        reg = regions{nr}; other = regions{3 - nr};
        lo_cols = {[reg '_%s_logodds_z'], [reg '_%s_kde_logodds'], [reg '_%s_kde_z_logodds'], [other '_%s_kde_z_logodds']};
        lo_lab = {sprintf('%s Bayesian z log odds (reference)', reg), sprintf('%s KDE log odds', reg), ...
            sprintf('%s KDE z log odds', reg), sprintf('%s KDE z log odds (other region)', other)};
        for j = 1:numel(lo_cols)
            nexttile;
            plot_alpha_curve(T, reg, crit, lo_cols{j}, lo_lab{j}, P);
        end
    end
    [~, kde_name] = fileparts(P.kde_file); [~, det_name] = fileparts(P.replay_dir);
    title(tl, {sprintf(['Significance criterion: %s (detection = Bayesian weighted correlation, %s, row = detecting region). ' ...
        'Track log odds from the PLS-KDE decoder (%s) over the same scoring window'], strjoin(crit, ' + '), det_name, kde_name), ...
        sprintf(['Bayesian vs KDE event z log odds (scored post): HPC Spearman r = %.2f (sign agreement %.0f%%), V1 r = %.2f (%.0f%%). ' ...
        'Events pooled; alpha 0.2 to 0.001 (larger marker = stricter); shaded = 95%% bootstrap CI'], agree.HPC(1), 100 * agree.HPC(2), ...
        agree.V1(1), 100 * agree.V1(2))}, 'Interpreter', 'none', 'FontSize', 9);
    save_fig(fig, fullfile(out_dir, ['criterion_' strjoin(crit, '_') '_KDE_logodds']));
    close(fig);
end
end

function plot_alpha_curve(T, reg, crit, lo_pattern, lo_label, P)
% Proportion of significant events (multi-track counted once) vs mean log odds of T1
% events - T2 events (multi-track events in both), pre and post windows
hold on
cols = {[0.4 0.4 0.4], [0.85 0.33 0.1]}; wins = {'pre', 'post'};
nA = numel(P.alphas); msz = linspace(12, 80, nA);
rs = RandStream('mrg32k3a', 'Seed', P.seed);
h = gobjects(1, 2);
for c = 1:2
    ok = T.(sprintf('%s_%s_scored', reg, wins{c}));
    E = struct();
    for ty = crit
        E.(ty{1}) = [T.(sprintf('%s_%s_p_%s_T1', reg, wins{c}, ty{1}))(ok) T.(sprintf('%s_%s_p_%s_T2', reg, wins{c}, ty{1}))(ok)];
    end
    sig = replay_significance(E, crit, P.alphas);
    lo = T.(sprintf(lo_pattern, wins{c}))(ok);
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
xlabel('mean log odds difference (T1 - T2 events)');
ylabel(sprintf('proportion of significant %s events', reg));
legend(h, wins, 'Location', 'northeast', 'Box', 'off');
title(lo_label, 'FontWeight', 'normal', 'FontSize', 8, 'Interpreter', 'none');
set(gca, 'TickDir', 'out');
end

function [prop, dif] = alpha_stats(sig, lo)
fin = isfinite(lo);
prop = reshape(mean(any(sig, 2), 1), 1, []);
dif = NaN(1, size(sig, 3));
for a = 1:size(sig, 3)
    t1 = sig(:, 1, a) & fin; t2 = sig(:, 2, a) & fin;
    if any(t1) && any(t2), dif(a) = mean(lo(t1)) - mean(lo(t2)); end
end
end

function save_fig(fig, f)
savefig(fig, [f '.fig']);
exportgraphics(fig, [f '.png'], 'Resolution', 150);
exportgraphics(fig, [f '.pdf'], 'ContentType', 'vector');
end
