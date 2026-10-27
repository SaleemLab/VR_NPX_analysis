function plot_pairwise_PRE_vs_POST(cell_set)
% PRE vs POST pairwise V1-HC coupling (partial correlation, gain 'total_V1') in the sessions
% that have PRE sleep. Reads session_summary.csv of plot_pairwise_V1_HC_coupling for both
% sleeps (pairwise_cell_set_folders). Session value = Spearman rho(partial coupling,
% z_V1 * z_HC) - permutation null mean.
% Statistics (two-sided Wilcoxon signed-rank across the matched sessions): PRE vs 0,
% POST vs 0, POST - PRE; number of sessions / animals (animal means) > 0.
% cell_set: 'all_cells' (default) or 'lt10'
% Output: pairwise_V1_HC_coupling\PRE_vs_POST_<cell_set>\ (stats.csv, PRE_vs_POST.fig/.pdf)

if nargin < 1, cell_set = 'all_cells'; end
rev_dir = addon_setup();
base = fullfile(rev_dir, 'pairwise_V1_HC_coupling');
[~, post_name, cells_txt] = pairwise_cell_set_folders('total_V1', cell_set, 'POST');
[~, pre_name] = pairwise_cell_set_folders('total_V1', cell_set, 'PRE');
out_dir = fullfile(base, ['PRE_vs_POST_' cell_set]);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
POST = readtable(fullfile(base, post_name, 'session_summary.csv'));
PRE = readtable(fullfile(base, pre_name, 'session_summary.csv'));
[sessions, ip, iq] = intersect(PRE.session_id, POST.session_id);
animal = PRE.animal_id(ip);
animals = unique(animal);
cols = lines(6);
wins = {'V1_pre', 'V1_post'};
win_labels = {'V1 -0.2-0 s, HC 0-0.1 s', 'V1 0-0.2 s, HC 0-0.1 s'};

T = table();
f = figure('Name', 'PRE_vs_POST', 'Position', [100 100 700 430], 'Color', 'w');
for w = 1:2
    col = ['rho_minus_null_partial_' wins{w}];
    a = PRE.(col)(ip); b = POST.(col)(iq); d = b - a;
    am = @(x) arrayfun(@(k) mean(x(animal == k)), animals);
    T = [T; table(string(wins{w}), numel(sessions), mean(a), mean(b), signrank(a), signrank(b), signrank(d), ...
        sum(a > 0), sum(b > 0), sum(d > 0), sum(am(a) > 0), sum(am(b) > 0), sum(am(d) > 0), numel(animals), ...
        'VariableNames', {'V1_window', 'n_sessions', 'mean_PRE', 'mean_POST', 'p_PRE_vs0', 'p_POST_vs0', ...
        'p_POST_minus_PRE', 'n_sess_PRE_pos', 'n_sess_POST_pos', 'n_sess_POST_gt_PRE', ...
        'n_anim_PRE_pos', 'n_anim_POST_pos', 'n_anim_POST_gt_PRE', 'n_animals'})]; %#ok<AGROW>

    subplot(1, 2, w); hold on
    plot([1 2], [a b]', '-', 'Color', [0.8 0.8 0.8]);
    for k = 1:numel(animals)
        idx = animal == animals(k);
        scatter([ones(sum(idx), 1); 2 * ones(sum(idx), 1)], [a(idx); b(idx)], 28, cols(animals(k), :), ...
            'filled', 'MarkerFaceAlpha', 0.75);
    end
    errorbar([0.75 2.25], [mean(a) mean(b)], [std(a) std(b)] / sqrt(numel(a)), 'ko-', ...
        'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'CapSize', 0);
    yline(0, '--');
    g = T(end, :);
    title({win_labels{w}, sprintf('PRE vs 0 p=%s | POST vs 0 p=%s', p_str(g.p_PRE_vs0), p_str(g.p_POST_vs0)), ...
        sprintf('POST-PRE p=%s (%d/%d sessions, %d/%d animals)', p_str(g.p_POST_minus_PRE), ...
        g.n_sess_POST_gt_PRE, g.n_sessions, g.n_anim_POST_gt_PRE, g.n_animals)}, 'FontWeight', 'normal', 'FontSize', 9);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'PRE', 'POST'}, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
    xlim([0.5 2.5]);
    if w == 1, ylabel('Pairwise \rho(coupling, z_{V1} z_{HC}) - null'); end
end
sgtitle({'Pairwise V1-HC coupling, PRE vs POST sleep (matched sessions)', ...
    sprintf('partial out total V1 MUA + HC MUA; %s', strrep(cells_txt, ', PRE sleep', ''))}, 'FontSize', 11);
writetable(T, fullfile(out_dir, 'stats.csv'));
disp(T);
savefig(f, fullfile(out_dir, 'PRE_vs_POST.fig'));
exportgraphics(f, fullfile(out_dir, 'PRE_vs_POST.pdf'), 'ContentType', 'vector');
fprintf('Saved %s\n', out_dir);
end

function s = p_str(p)
if p < 0.001
    s = sprintf('%.1e', p);
else
    s = sprintf('%.3f', p);
end
end
