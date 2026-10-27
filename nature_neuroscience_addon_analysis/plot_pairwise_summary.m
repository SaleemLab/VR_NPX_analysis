function plot_pairwise_summary(gain_control, cell_set, sleep, measure)
% Summary figures of the pairwise V1-HC coupling for the revision.
% measure: 'partial' (default; partial correlation, population MUA removed) or 'raw'
% (Pearson correlation of the counts without removing MUA; equal to the correlation of
% z-scored counts). Raw figures are saved with the suffix _raw.
% Reads session_summary.csv and group_stats.csv of plot_pairwise_V1_HC_coupling and the
% per-session bootstrapped AUC of boot_AUC_per_unit_POST_original('session')
% (POST_original_bootAUC_per_session\bootAUC_all_ripples.mat, static mean AUC,
% all ripples, HPC -> V1).
%   summary_pairwise_partial         : sessions in one column per V1 window (pre, post),
%                                      mean +- SEM across sessions, signrank p, number of
%                                      sessions / animals (animal means) > 0
%   summary_pairwise_vs_KDE_logodds_r: per session, pairwise effect vs Pearson r between
%                                      PLS-KDE V1 and HC log odds across ripples
%   summary_pairwise_vs_KDE_bootAUC  : per session, pairwise effect vs mean bootstrapped
%                                      PLS-KDE reactivation coherence AUC
% Pairwise effect = Spearman rho(partial coupling, z_V1 * z_HC) - permutation null mean.
% cell_set: 'lt10' (default, cells < 10 Hz) or 'all_cells'; sleep: 'POST' (default) or
% 'PRE' (pairwise_cell_set_folders). For PRE the per-session bootstrapped AUC (all
% ripples, HPC -> V1) is taken from compare_session_bootAUC_vs_shuffleCI
% (session_AUC_vs_shuffle\session_bootAUC_vs_shuffleCI.mat, AUC_boot_mean / shuffle_mean).
% Output: nature_neuroscience_revision_analysis\pairwise_V1_HC_coupling\gain_<gain_control>[_all_cells][_PRE]

if nargin < 1, gain_control = 'total_V1'; end
if nargin < 2, cell_set = 'lt10'; end
if nargin < 3 || isempty(sleep), sleep = 'POST'; end
if nargin < 4 || isempty(measure), measure = 'partial'; end
assert(ismember(measure, {'partial', 'raw'}), 'measure must be ''partial'' or ''raw''');
sfx = ''; if strcmp(measure, 'raw'), sfx = '_raw'; end
rev_dir = addon_setup();
[~, out_name, cells_txt] = pairwise_cell_set_folders(gain_control, cell_set, sleep);
out_dir = fullfile(rev_dir, 'pairwise_V1_HC_coupling', out_name);
SS = readtable(fullfile(out_dir, 'session_summary.csv'));
G = readtable(fullfile(out_dir, 'group_stats.csv'), 'TextType', 'string');
% Per-session bootstrapped AUC (sessions x V1 window [pre post]) and its shuffled mean
has_auc = true;
AUC = NaN(height(SS), 2); AUC_shuf = NaN(height(SS), 2);
if strcmp(sleep, 'POST')
    B = load(fullfile(rev_dir, 'POST_original_bootAUC_per_session', 'bootAUC_all_ripples.mat'), 'static', 'info');
    [~, auc_row] = ismember(SS.session_id, B.info.unit_ids);
    AUC = B.static.mean(auc_row, :, 1); AUC_shuf = B.static.shuffled_mean(auc_row, :, 1);
else
    L = load(fullfile(rev_dir, 'session_AUC_vs_shuffle', 'session_bootAUC_vs_shuffleCI.mat'), 'T');
    A = L.T(string(L.T.sleep) == sleep & string(L.T.analysis) == "all_ripples", :);
    wn = {'pre', 'post'};
    for w = 1:2
        Aw = A(string(A.V1_window) == wn{w}, :);
        [ok, r] = ismember(SS.session_id, Aw.session_id);
        AUC(ok, w) = Aw.AUC_boot_mean(r(ok)); AUC_shuf(ok, w) = Aw.shuffle_mean(r(ok));
    end
end

wins = {'V1_pre', 'V1_post'};
win_labels = {'V1 -0.2-0 s', 'V1 0-0.2 s'};
ylab = {'Pairwise \rho(coupling, z_{V1} z_{HC}) - null'};
animals = unique(SS.animal_id);
cols = lines(numel(animals));

%% 1. Pairwise effect, all sessions in one column per window
col = @(w) ['rho_minus_null_' measure '_' wins{w}];
f1 = figure('Name', ['summary_pairwise_' measure], 'Position', [100 100 380 420], 'Color', 'w');
hold on
all_val = SS{:, {col(1), col(2)}};
y_lo = min(all_val, [], 'all'); y_hi = max(all_val, [], 'all');
y_text = y_hi + 0.08 * (y_hi - y_lo);
for w = 1:2
    val = SS.(col(w));
    g = G(G.measure == measure & G.V1_window == wins{w} & G.statistic == "rho", :);
    rs = RandStream('twister', 'Seed', w);
    x = w - 0.15 + 0.12 * (rand(rs, numel(val), 1) - 0.5);
    for a = 1:numel(animals)
        idx = SS.animal_id == animals(a);
        scatter(x(idx), val(idx), 28, cols(a, :), 'filled', 'MarkerFaceAlpha', 0.7);
    end
    ok = isfinite(val);
    errorbar(w + 0.15, mean(val(ok)), std(val(ok)) / sqrt(sum(ok)), 'ko', ...
        'MarkerFaceColor', 'k', 'MarkerSize', 7, 'LineWidth', 1.5, 'CapSize', 0);
    text(w, y_text, ...
        {sprintf('signrank p=%s', p_str(g.p_signrank)), sprintf('%d/%d sessions > 0', g.n_sessions_pos, g.n_sessions), ...
        sprintf('%d/%d animals > 0', g.n_animals_pos, g.n_animals)}, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
end
yline(0, '--');
xlim([0.4 2.6]);
ylim([y_lo - 0.05 * (y_hi - y_lo), y_hi + 0.35 * (y_hi - y_lo)]);
set(gca, 'XTick', 1:2, 'XTickLabel', win_labels, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
ylabel(ylab);
if strcmp(gain_control, 'total_V1'), gain_txt = 'total V1 MUA'; else, gain_txt = 'own-hemisphere V1 MUA'; end
if strcmp(measure, 'raw')
    meas_txt = 'raw correlation (MUA not removed), HC 0-0.1 s';
else
    meas_txt = sprintf('partial out %s + HC MUA, HC 0-0.1 s', gain_txt);
end
title({'Pairwise V1-HC coupling vs shared track preference', meas_txt, cells_txt}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% 2. Pairwise effect vs PLS-KDE log-odds correlation; 3. vs PLS-KDE bootstrapped AUC
xvars = {@(w) SS.(['KDE_logodds_r_' wins{w}])};
xlabs = {{'PLS-KDE log-odds correlation', 'r(V1, HC) across ripples'}};
names = {['summary_pairwise_vs_KDE_logodds_r' sfx]};
if has_auc
    xvars{2} = @(w) AUC(:, w);
    xlabs{2} = {'PLS-KDE reactivation coherence AUC', '(bootstrap mean, all ripples)'};
    names{2} = ['summary_pairwise_vs_KDE_bootAUC' sfx];
end
figs = gobjects(1, numel(names));
for k = 1:numel(names)
    figs(k) = figure('Name', names{k}, 'Position', [150 + 40 * k, 150, 900, 420], 'Color', 'w');
    for w = 1:2
        subplot(1, 2, w); hold on
        xv = xvars{k}(w);
        yv = SS.(col(w));
        ok = isfinite(xv) & isfinite(yv);
        c = polyfit(xv(ok), yv(ok), 1);
        xx = [min(xv(ok)) max(xv(ok))];
        plot(xx, polyval(c, xx), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        for a = 1:numel(animals)
            idx = SS.animal_id == animals(a);
            scatter(xv(idx), yv(idx), 30, cols(a, :), 'filled', 'MarkerFaceAlpha', 0.8);
        end
        [rr, pp] = corr(xv(ok), yv(ok), 'type', 'Spearman');
        yline(0, '--');
        if k == 1
            xline(0, '--');
        else
            xline(mean(AUC_shuf(:, w), 'omitnan'), '--', 'shuffle', ...
                'LabelVerticalAlignment', 'bottom', 'FontSize', 8);
        end
        xlabel(xlabs{k}); ylabel(ylab);
        title({sprintf('%s, HC 0-0.1 s', win_labels{w}), sprintf('Spearman r = %.2f, p = %s (n = %d sessions)', ...
            rr, p_str(pp), sum(ok))}, 'FontWeight', 'normal', 'FontSize', 9);
        set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
    end
end

for f = [f1 figs]
    name = get(f, 'Name');
    savefig(f, fullfile(out_dir, [name '.fig']));
    exportgraphics(f, fullfile(out_dir, [name '.pdf']), 'ContentType', 'vector');
end
fprintf('Saved summary figures in %s\n', out_dir);
end

function s = p_str(p)
if p < 0.001
    s = sprintf('%.1e', p);
else
    s = sprintf('%.3f', p);
end
end
