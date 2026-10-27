function plot_cluster_MUA_correlation(cell_set, window, example_session)
% Simple illustration of the cluster-MUA correlation and its effect on raw vs partial V1-HC
% pair coupling (POST ripples, gain control 'total_V1').
% Cluster-MUA correlation m of a cluster = Pearson correlation, across ripples, of its spike
% count with the MUA of its own region without its own spikes (V1: V1_L + V1_R in the V1
% window; HC: HPC MUA 0-0.1 s).
%   a : two example V1 clusters (high / low m): firing relative to the cluster's mean in
%       quintiles of V1 MUA
%   b : distribution of m (V1, HC clusters)
%   c : pairs split into within-session quintiles of the V1 cluster's m_V1, of the HC
%       cluster's m_HC, and of their product m_V1 * m_HC; mean pair correlation per
%       quintile (raw, partial), mean +- SEM across sessions
%   d : per session Spearman rho(pair correlation, m_V1), rho(pair correlation, m_HC) and
%       rho(pair correlation, m_V1 * m_HC), raw vs partial (two-sided signrank vs 0)
% cell_set: 'all_cells' (default), 'lt10', 'rate_0p5_5'; window: 'V1_post' (default) or
% 'V1_pre'; example_session: session for panel a (default 17)
% Output: pairwise_V1_HC_coupling\gain_total_V1[_<cell_set>]\gain_coupling_check\
%   cluster_MUA_correlation_<window>.fig/.pdf

if nargin < 1 || isempty(cell_set), cell_set = 'all_cells'; end
if nargin < 2 || isempty(window), window = 'V1_post'; end
if nargin < 3 || isempty(example_session), example_session = 17; end
rev_dir = addon_setup();
base = fullfile(rev_dir, 'pairwise_V1_HC_coupling');
[counts_name, out_name, cells_txt] = pairwise_cell_set_folders('total_V1', cell_set);
out_dir = fullfile(base, out_name, 'gain_coupling_check');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
load(fullfile(base, out_name, 'coupling_POST.mat'), 'S');
win_label = struct('V1_pre', 'V1 -0.2-0 s', 'V1_post', 'V1 0-0.2 s').(window);
nq = 5;

m_V1_all = []; m_HC_all = [];
bin_raw = NaN(numel(S), nq, 3); bin_par = NaN(numel(S), nq, 3);
rho_raw = NaN(numel(S), 3); rho_par = NaN(numel(S), 3);
for s = 1:numel(S)
    L = load(fullfile(base, counts_name, sprintf('session_%02d.mat', S(s).session_count)), 'C');
    C = L.C;
    W = S(s).(window);
    [~, iv] = ismember(W.V1_cluster_id, C.V1.cluster_id);
    [~, ih] = ismember(W.HPC_cluster_id, C.HPC.cluster_id);
    v = C.V1.counts.(window)(:, iv); h = C.HPC.counts.HC(:, ih);
    V1_tot = C.MUA.V1_L.(window) + C.MUA.V1_R.(window);
    m_v = col_corr(v, V1_tot - v); m_h = col_corr(h, C.MUA.HPC.HC - h);
    m_V1_all = [m_V1_all; m_v(:)]; m_HC_all = [m_HC_all; m_h(:)]; %#ok<AGROW>
    g = m_v(:) * m_h(:)';
    ok = isfinite(g) & isfinite(W.r_raw) & isfinite(W.r_partial);
    rr = W.r_raw(ok); rp = W.r_partial(ok);
    [mv_pair, mh_pair] = ndgrid(m_v(:), m_h(:));
    x_group = {mv_pair(ok), mh_pair(ok), g(ok)};      % V1 cluster, HC cluster, product
    for k = 1:3
        rho_raw(s, k) = corr(rr, x_group{k}, 'type', 'Spearman');
        rho_par(s, k) = corr(rp, x_group{k}, 'type', 'Spearman');
        q = discretize(x_group{k}, unique(prctile(x_group{k}, linspace(0, 100, nq + 1))));
        if max(q) < nq, continue; end                  % too few distinct values for quintiles
        bin_raw(s, :, k) = accumarray(q, rr, [nq 1], @mean, NaN)';
        bin_par(s, :, k) = accumarray(q, rp, [nq 1], @mean, NaN)';
    end
    if S(s).session_count == example_session
        ex = struct('v', v, 'mua', V1_tot - v, 'm', m_v, 'session', example_session);
    end
end
animal = [S.subject_id]';
cols = lines(6);
col_raw = [0.3 0.3 0.3]; col_par = [0.85 0.33 0.1];

f = figure('Name', ['cluster_MUA_correlation_' window], 'Position', [60 40 1150 1050], 'Color', 'w');

% a: example clusters
subplot(3, 3, 1); hold on
active = mean(ex.v, 1) >= 0.3;                       % enough spikes for a clear curve
cand = find(active);
[~, o] = sort(ex.m(cand));
pick = [cand(o(end - round(0.1 * numel(o)))), cand(o(max(1, round(0.1 * numel(o)))))];
ec = [0.15 0.15 0.15; 0.6 0.6 0.6];
qm = discretize(ex.mua(:, 1), prctile(ex.mua(:, 1), linspace(0, 100, nq + 1)));
for k = 1:2
    i = pick(k);
    qi = discretize(ex.mua(:, i), prctile(ex.mua(:, i), linspace(0, 100, nq + 1)));
    if any(isnan(qi)), qi = qm; end
    y = accumarray(qi, ex.v(:, i), [nq 1], @mean) / mean(ex.v(:, i));
    plot(1:nq, y, 'o-', 'Color', ec(k, :), 'MarkerFaceColor', ec(k, :), 'LineWidth', 1.5);
end
yline(1, ':');
legend({sprintf('high (m = %.2f)', ex.m(pick(1))), sprintf('low (m = %.2f)', ex.m(pick(2)))}, ...
    'Box', 'off', 'Location', 'northwest', 'FontSize', 8);
set(gca, 'XTick', 1:nq, 'XTickLabel', {'low', '', '', '', 'high'}, 'TickDir', 'out', 'Box', 'off', 'FontSize', 9);
xlabel('V1 MUA without the cluster (quintile)'); ylabel('firing / cluster''s mean');
title({'a  Cluster-MUA correlation', sprintf('two V1 clusters, session %d', ex.session)}, 'FontWeight', 'normal', 'FontSize', 10);

% b: distribution
subplot(3, 3, 2); hold on
edges = -0.2:0.05:0.9;
histogram(m_V1_all, edges, 'Normalization', 'probability', 'FaceColor', col_raw, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
histogram(m_HC_all, edges, 'Normalization', 'probability', 'FaceColor', col_par, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
xline(0, '--');
legend({sprintf('V1 (median %.2f)', median(m_V1_all, 'omitnan')), sprintf('HC (median %.2f)', median(m_HC_all, 'omitnan'))}, ...
    'Box', 'off', 'FontSize', 8);
xlabel('cluster-MUA correlation m'); ylabel('fraction of clusters');
title({'b  Distribution', sprintf('%d V1, %d HC clusters x sessions', numel(m_V1_all), numel(m_HC_all))}, ...
    'FontWeight', 'normal', 'FontSize', 10);
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 9);

% definitions
subplot(3, 3, 3); axis off
text(0, 0.95, {'\bfCluster-MUA correlation m\rm', ...
    'Pearson r, across ripples, between a cluster''s', 'spike count and the MUA of its own region', ...
    'without that cluster', '', ...
    'm_{V1}: V1 cluster vs V1 MUA (L + R, V1 window)', 'm_{HC}: HC cluster vs HC MUA (0-0.1 s)', '', ...
    'pair r: V1-HC cluster pair correlation,', 'raw or partial (V1 and HC MUA removed)'}, ...
    'VerticalAlignment', 'top', 'FontSize', 9);

% c: pair correlation by quintile of the V1 cluster's m, the HC cluster's m and their
% product (same y axis); d: per-session rho with each (same y axis)
x_labs = {'m_{V1} of the V1 cluster (quintile)', 'm_{HC} of the HC cluster (quintile)', ...
    'm_{V1} \times m_{HC} (quintile)'};
d_ylabs = {'\rho(pair r, m_{V1})', '\rho(pair r, m_{HC})', '\rho(pair r, m_{V1} \times m_{HC})'};
col_titles = {'Grouped by V1 cluster', 'Grouped by HC cluster', 'Grouped by both (product)'};
c_tags = {'c  ', '', ''}; d_tags = {'d  ', '', ''};
ax_c = gobjects(1, 3); ax_d = gobjects(1, 3);
for g_i = 1:3
    ax_c(g_i) = subplot(3, 3, 3 + g_i); hold on
    for k = 1:2
        B = bin_raw(:, :, g_i); colr = col_raw;
        if k == 2, B = bin_par(:, :, g_i); colr = col_par; end
        n_ok = sum(all(isfinite(B), 2));
        errorbar(1:nq, mean(B, 1, 'omitnan'), std(B, 0, 1, 'omitnan') / sqrt(n_ok), 'o-', ...
            'Color', colr, 'MarkerFaceColor', colr, 'LineWidth', 1.5, 'CapSize', 0);
    end
    yline(0, '--');
    if g_i == 1
        legend({'raw', 'partial (MUA removed)'}, 'Box', 'off', 'Location', 'northwest', 'FontSize', 8);
        ylabel('V1-HC pair correlation r');
    end
    set(gca, 'XTick', 1:nq, 'XTickLabel', {'low', '', '', '', 'high'}, 'TickDir', 'out', 'Box', 'off', 'FontSize', 9);
    xlabel(x_labs{g_i});
    title({[c_tags{g_i} col_titles{g_i}], sprintf('mean \\pm SEM, %d sessions', n_ok)}, ...
        'FontWeight', 'normal', 'FontSize', 10);

    ax_d(g_i) = subplot(3, 3, 6 + g_i); hold on
    plot([1 2], [rho_raw(:, g_i) rho_par(:, g_i)]', '-', 'Color', [0.85 0.85 0.85]);
    for a = unique(animal)'
        idx = animal == a;
        scatter([ones(sum(idx), 1); 2 * ones(sum(idx), 1)], [rho_raw(idx, g_i); rho_par(idx, g_i)], 22, cols(a, :), ...
            'filled', 'MarkerFaceAlpha', 0.75);
    end
    errorbar([0.75 2.25], [mean(rho_raw(:, g_i)) mean(rho_par(:, g_i))], ...
        [std(rho_raw(:, g_i)) std(rho_par(:, g_i))] / sqrt(numel(S)), 'ko', ...
        'MarkerFaceColor', 'k', 'LineWidth', 1.4, 'CapSize', 0);
    yline(0, '--');
    xlim([0.5 2.5]);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'raw', 'partial'}, 'TickDir', 'out', 'Box', 'off', 'FontSize', 9);
    ylabel(d_ylabs{g_i});
    title({[d_tags{g_i} 'Per session'], sprintf('raw p = %s, partial p = %s (signrank)', ...
        p_str(signrank(rho_raw(:, g_i))), p_str(signrank(rho_par(:, g_i))))}, 'FontWeight', 'normal', 'FontSize', 10);
end
linkaxes(ax_c, 'y'); linkaxes(ax_d, 'y');

sgtitle(sprintf('Cluster-MUA correlation and V1-HC pair correlation (POST ripples, %s, HC 0-0.1 s; %s)', ...
    win_label, cells_txt), 'FontSize', 11);
savefig(f, fullfile(out_dir, ['cluster_MUA_correlation_' window '.fig']));
exportgraphics(f, fullfile(out_dir, ['cluster_MUA_correlation_' window '.pdf']), 'ContentType', 'vector');
fprintf('%s: rho(raw r) mean %s, rho(partial r) mean %s for m_V1, m_HC, product\n', window, ...
    mat2str(mean(rho_raw), 3), mat2str(mean(rho_par), 3));
fprintf('signrank p raw %s, partial %s\n', mat2str(arrayfun(@(k) signrank(rho_raw(:, k)), 1:3), 3), ...
    mat2str(arrayfun(@(k) signrank(rho_par(:, k)), 1:3), 3));
fprintf('Saved %s\n', fullfile(out_dir, ['cluster_MUA_correlation_' window '.pdf']));
end

function c = col_corr(X, Y)
X = X - mean(X, 1); Y = Y - mean(Y, 1);
c = sum(X .* Y, 1) ./ sqrt(sum(X .^ 2, 1) .* sum(Y .^ 2, 1));
end

function s = p_str(p)
if p < 0.001
    s = sprintf('%.1e', p);
else
    s = sprintf('%.3f', p);
end
end
