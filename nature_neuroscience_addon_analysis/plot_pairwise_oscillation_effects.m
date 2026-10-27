function plot_pairwise_oscillation_effects(out_name, gain_folder)
% Group statistics and figures of pairwise_coupling_oscillation_POST
% (pairwise_V1_HC_coupling\<gain_folder>\<out_name>; defaults gain_total_V1, oscillation_effects).
% Session value = rho(partial coupling, z_V1 * z_HC) - null within each event subset.
%   oscillation_effects_paired          : per condition and V1 window, A -> B per session
%   oscillation_effects_paired_thinned  : the same after binomial spike thinning
%   oscillation_effects_delta           : B - A per session, original vs thinned
% Statistics (signrank across sessions): A vs 0, B vs 0, B - A vs 0; number of sessions
% and animals (animal mean) with B > A. Table: group_stats.csv

if nargin < 1 || isempty(out_name), out_name = 'oscillation_effects'; end
if nargin < 2 || isempty(gain_folder), gain_folder = 'gain_total_V1'; end
rev_dir = addon_setup();
out_dir = fullfile(rev_dir, 'pairwise_V1_HC_coupling', gain_folder, out_name);
load(fullfile(out_dir, 'oscillation_coupling_POST.mat'), 'T', 'P');
if ~isfield(P, 'hemisphere_mode'), P.hemisphere_mode = 'own'; end
if ~isfield(P, 'counts_name'), P.counts_name = 'counts'; end
hemi_txt = struct('own', 'own hemisphere', 'dominant', 'dominant hemisphere', ...
    'preference', 'preferred-track hemisphere');
ht = hemi_txt.(P.hemisphere_mode);
cells_txt = 'cells < 10 Hz';
if strcmp(P.counts_name, 'counts_all_rates'), cells_txt = 'all spatial cells'; end
if ~isfield(P, 'measure'), P.measure = 'correlation'; end
if strcmp(P.measure, 'covariance')
    y_lab = 'slope(partial cov, z_{V1} z_{HC}) - null';
    d_lab = '\Delta slope (B - A)';
    cells_txt = [cells_txt ', partial covariance'];
else
    y_lab = '\rho(coupling, z_{V1} z_{HC}) - null';
    d_lab = '\Delta\rho (B - A)';
end
conds = {'ripple_power', 'spindle', 'SO_phase', 'SO_trough_uni_bi'};
cond_titles = {'Ripple power', sprintf('Spindle-band power (%s)', ht), sprintf('SO phase (%s)', ht), ...
    sprintf('SO trough (%s)', ht)};
wins = {'V1_pre', 'V1_post'};
win_labels = {'V1 -0.2-0 s', 'V1 0-0.2 s'};
versions = {'', '_thinned'};
version_labels = {'original', 'spike-thinned'};
animals = unique(T.animal_id);
cols = lines(numel(animals));

%% Group statistics
G = table();
for v = 1:2
    for w = 1:2
        for c = 1:numel(conds)
            S = T(T.V1_window == wins{w} & T.condition == conds{c}, :);
            a = S.(['rho_A' versions{v}]); b = S.(['rho_B' versions{v}]); d = b - a;
            ok = isfinite(d);
            dm = arrayfun(@(x) mean(d(ok & S.animal_id == x)), animals);
            G = [G; table(string(version_labels{v}), string(wins{w}), string(conds{c}), S.label_A(1), S.label_B(1), ...
                sum(ok), mean(a(ok)), mean(b(ok)), signrank(a(ok)), signrank(b(ok)), mean(d(ok)), signrank(d(ok)), ...
                sum(d(ok) > 0), sum(dm > 0), sum(isfinite(dm)), ...
                mean(S.V1_count_A(ok)), mean(S.V1_count_B(ok)), mean(S.HC_count_A(ok)), mean(S.HC_count_B(ok)), ...
                'VariableNames', {'version', 'V1_window', 'condition', 'label_A', 'label_B', 'n_sessions', ...
                'mean_A', 'mean_B', 'p_A_vs0', 'p_B_vs0', 'mean_delta', 'p_delta', 'n_sessions_B_gt_A', ...
                'n_animals_B_gt_A', 'n_animals', 'V1_count_A', 'V1_count_B', 'HC_count_A', 'HC_count_B'})]; %#ok<AGROW>
        end
    end
end
writetable(G, fullfile(out_dir, 'group_stats.csv'));
disp(G(:, {'version', 'V1_window', 'condition', 'mean_A', 'mean_B', 'p_delta', 'n_sessions_B_gt_A', 'n_animals_B_gt_A'}));

%% Paired A -> B per session
figs = gobjects(1, 3);
for v = 1:2
    figs(v) = figure('Name', ['oscillation_effects_paired' versions{v}], 'Position', [40 40 1300 650], 'Color', 'w');
    for w = 1:2
        for c = 1:numel(conds)
            subplot(2, numel(conds), (w - 1) * numel(conds) + c); hold on
            S = T(T.V1_window == wins{w} & T.condition == conds{c}, :);
            a = S.(['rho_A' versions{v}]); b = S.(['rho_B' versions{v}]);
            ok = isfinite(a) & isfinite(b);
            plot([1 2], [a(ok) b(ok)]', '-', 'Color', [0.8 0.8 0.8]);
            for k = 1:numel(animals)
                idx = ok & S.animal_id == animals(k);
                scatter(ones(sum(idx), 1), a(idx), 20, cols(k, :), 'filled', 'MarkerFaceAlpha', 0.7);
                scatter(2 * ones(sum(idx), 1), b(idx), 20, cols(k, :), 'filled', 'MarkerFaceAlpha', 0.7);
            end
            m = [mean(a(ok)) mean(b(ok))];
            se = [std(a(ok)) std(b(ok))] / sqrt(sum(ok));
            errorbar([0.75 2.25], m, se, 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'CapSize', 0);
            yline(0, '--');
            g = G(G.version == version_labels{v} & G.V1_window == wins{w} & G.condition == conds{c}, :);
            title({cond_titles{c}, sprintf('%s | B-A signrank p=%s', win_labels{w}, p_str(g.p_delta)), ...
                sprintf('B>A: %d/%d sessions, %d/%d animals', g.n_sessions_B_gt_A, g.n_sessions, ...
                g.n_animals_B_gt_A, g.n_animals)}, 'FontWeight', 'normal', 'FontSize', 8);
            set(gca, 'XTick', [1 2], 'XTickLabel', {S.label_A(1), S.label_B(1)}, 'TickDir', 'out', ...
                'Box', 'off', 'FontSize', 9);
            xlim([0.5 2.5]);
            if c == 1, ylabel(y_lab); end
        end
    end
    sgtitle({sprintf('Pairwise V1-HC coupling by oscillation condition (%s, %s)', version_labels{v}, cells_txt), ...
        sprintf('partial out total V1 MUA + HC MUA, HC 0-0.1 s; spindle / SO from the %s', ht)}, ...
        'FontSize', 11);
end

%% Delta B - A, original vs thinned
figs(3) = figure('Name', 'oscillation_effects_delta', 'Position', [80 80 1100 420], 'Color', 'w');
short = {'Ripple power', 'Spindle', 'SO trough-peak', 'Uni-bilateral'};
for w = 1:2
    subplot(1, 2, w); hold on
    ptxt = cell(numel(conds), 2);
    for c = 1:numel(conds)
        S = T(T.V1_window == wins{w} & T.condition == conds{c}, :);
        for v = 1:2
            d =S.(['rho_B' versions{v}]) - S.(['rho_A' versions{v}]);
            ok = isfinite(d);
            x0 = c + (v - 1.5) * 0.35;
            rs = RandStream('twister', 'Seed', c * 10 + v);
            xj = x0 + 0.1 * (rand(rs, sum(ok), 1) - 0.5);
            scatter(xj, d(ok), 14, [0.5 0.5 0.5] + (v == 2) * [-0.3 0.1 0.4], 'filled', 'MarkerFaceAlpha', 0.5);
            errorbar(x0, mean(d(ok)), std(d(ok)) / sqrt(sum(ok)), 'ko', 'MarkerFaceColor', 'k', ...
                'LineWidth', 1.5, 'CapSize', 0);
            g = G(G.version == version_labels{v} & G.V1_window == wins{w} & G.condition == conds{c}, :);
            ptxt{c, v} = p_str(g.p_delta);
        end
    end
    yline(0, '--');
    yl = ylim;
    for c = 1:numel(conds)
        for v = 1:2
            text(c + (v - 1.5) * 0.35, yl(2), ptxt{c, v}, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', 7);
        end
    end
    ylim([yl(1) yl(2) + 0.1 * diff(yl)]);
    set(gca, 'XTick', 1:numel(conds), 'XTickLabel', short, 'TickDir', 'out', 'Box', 'off', 'FontSize', 9);
    ylabel(d_lab);
    title({win_labels{w}, 'left: original, right: spike-thinned; numbers: signrank p'}, 'FontWeight', 'normal', 'FontSize', 9);
end
sgtitle(sprintf('Oscillation effect on pairwise V1-HC coupling (B - A per session; %s, %s)', ht, cells_txt), ...
    'FontSize', 11);

for f = figs
    name = get(f, 'Name');
    savefig(f, fullfile(out_dir, [name '.fig']));
    exportgraphics(f, fullfile(out_dir, [name '.pdf']), 'ContentType', 'vector');
end
fprintf('Saved figures and group_stats.csv in %s\n', out_dir);
end

function s = p_str(p)
if p < 0.001
    s = sprintf('%.1e', p);
else
    s = sprintf('%.3f', p);
end
end
