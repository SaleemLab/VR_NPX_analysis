function pairwise_gain_coupling_check(cell_set)
% Why the raw (MUA not removed) pairwise V1-HC coupling does not show the track-preference
% effect: tests whether each pair's raw correlation is set by how strongly its two clusters
% follow the MUA of their own region (cluster-MUA correlation m).
% Per session and V1 window (pairs and r of pairwise_V1_HC_coupling_POST('total_V1')):
%   m_V1(i) : corr across ripples of V1 cluster i's count with the V1 MUA (V1_L + V1_R)
%             without its own spikes, same window
%   m_HC(j) : corr of HC cluster j's count (0-0.1 s) with the HC MUA without its own spikes
%   gain term g_ij = m_V1(i) * m_HC(j)
% Statistics per session (Spearman across pairs):
%   gain_vs_raw     : rho(r_raw, g)          -> raw coupling is set by cluster-MUA correlation
%   gain_vs_partial : rho(r_partial, g)      -> removed by partialling the MUA
%   gain_vs_pref    : rho(g, z_V1 * z_HC)    -> cluster-MUA correlation unrelated to preference
%   R2_gain         : variance of r_raw explained by a linear fit on g
%   pref_raw, pref_raw_minus_gain, pref_partial : rho(coupling, z_V1 * z_HC) - permutation
%                     null (nPerm, preferences permuted within region) for r_raw, for r_raw
%                     with the linear gain term regressed out across pairs, and for r_partial
% Two-sided signrank across sessions. cell_set: 'all_cells' (default), 'lt10', 'rate_0p5_5'.
% Output: pairwise_V1_HC_coupling\gain_total_V1[_<cell_set>]\gain_coupling_check\

if nargin < 1, cell_set = 'all_cells'; end
nPerm = 1000;
rev_dir = addon_setup();
base = fullfile(rev_dir, 'pairwise_V1_HC_coupling');
[counts_name, out_name, cells_txt] = pairwise_cell_set_folders('total_V1', cell_set);
out_dir = fullfile(base, out_name, 'gain_coupling_check');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
load(fullfile(base, out_name, 'coupling_POST.mat'), 'S');
wins = {'V1_pre', 'V1_post'};
win_labels = {'V1 -0.2-0 s, HC 0-0.1 s', 'V1 0-0.2 s, HC 0-0.1 s'};

T = table(); ex = struct();
for s = 1:numel(S)
    L = load(fullfile(base, counts_name, sprintf('session_%02d.mat', S(s).session_count)), 'C');
    C = L.C;
    h_all = C.HPC.counts.HC; H = C.MUA.HPC.HC;
    for w = 1:2
        W = S(s).(wins{w});
        [~, iv] = ismember(W.V1_cluster_id, C.V1.cluster_id);
        [~, ih] = ismember(W.HPC_cluster_id, C.HPC.cluster_id);
        v = C.V1.counts.(wins{w})(:, iv); h = h_all(:, ih);
        V1_tot = C.MUA.V1_L.(wins{w}) + C.MUA.V1_R.(wins{w});
        m_v = col_corr(v, V1_tot - v);
        m_h = col_corr(h, H - h);
        g = m_v(:) * m_h(:)';
        pz = W.z_V1(:) * W.z_HC(:)';
        ok = isfinite(W.r_raw) & isfinite(W.r_partial) & isfinite(g);
        rr = W.r_raw(ok); rp = W.r_partial(ok); gg = g(ok); pzz = pz(ok);
        X = [ones(numel(gg), 1) gg];
        r_res = rr - X * (X \ rr);
        R2 = 1 - sum(r_res .^ 2) / sum((rr - mean(rr)) .^ 2);
        [pv, ph] = ndgrid(1:numel(iv), 1:numel(ih)); pv = pv(ok); ph = ph(ok);
        rng(S(s).session_count * 10 + w, 'twister');
        T = [T; table(S(s).session_count, S(s).subject_id, string(wins{w}), nnz(ok), mean(m_v), mean(m_h), ...
            corr(rr, gg, 'type', 'Spearman'), corr(rp, gg, 'type', 'Spearman'), corr(gg, pzz, 'type', 'Spearman'), R2, ...
            pref_minus_null(rr, pv, ph, W.z_V1, W.z_HC, nPerm), ...
            pref_minus_null(r_res, pv, ph, W.z_V1, W.z_HC, nPerm), ...
            pref_minus_null(rp, pv, ph, W.z_V1, W.z_HC, nPerm), ...
            'VariableNames', {'session_id', 'animal_id', 'V1_window', 'n_pairs', 'mean_m_V1', 'mean_m_HC', ...
            'gain_vs_raw', 'gain_vs_partial', 'gain_vs_pref', 'R2_gain', ...
            'pref_raw', 'pref_raw_minus_gain', 'pref_partial'})]; %#ok<AGROW>
        if s == find([S.session_count] == 17, 1) || (s == 1 && ~any([S.session_count] == 17))
            ex.(wins{w}) = struct('g', gg, 'raw', rr, 'partial', rp, 'session', S(s).session_count);
        end
    end
    fprintf('Session %d done\n', S(s).session_count);
end
writetable(T, fullfile(out_dir, 'session_stats.csv'));

%% Group statistics
metrics = {'gain_vs_raw', 'gain_vs_partial', 'gain_vs_pref', 'R2_gain', 'pref_raw', 'pref_raw_minus_gain', 'pref_partial'};
animals = unique(T.animal_id);
G = table();
for w = 1:2
    Tw = T(T.V1_window == wins{w}, :);
    for i = 1:numel(metrics)
        x = Tw.(metrics{i});
        am = arrayfun(@(a) mean(x(Tw.animal_id == a)), animals);
        G = [G; table(string(wins{w}), string(metrics{i}), height(Tw), mean(x), std(x) / sqrt(numel(x)), ...
            signrank(x), sum(x > 0), sum(am > 0), numel(animals), ...
            'VariableNames', {'V1_window', 'metric', 'n_sessions', 'mean', 'SEM', 'p_signrank', ...
            'n_sessions_pos', 'n_animals_pos', 'n_animals'})]; %#ok<AGROW>
    end
    d = Tw.pref_raw_minus_gain - Tw.pref_raw;
    G = [G; table(string(wins{w}), "pref_raw_minus_gain - pref_raw", height(Tw), mean(d), std(d) / sqrt(numel(d)), ...
        signrank(d), sum(d > 0), sum(arrayfun(@(a) mean(d(Tw.animal_id == a)), animals) > 0), numel(animals), ...
        'VariableNames', G.Properties.VariableNames)]; %#ok<AGROW>
end
writetable(G, fullfile(out_dir, 'group_stats.csv'));
disp(G);

%% Figure
cols = lines(6);
f = figure('Name', 'gain_coupling_check', 'Position', [60 60 1250 700], 'Color', 'w');
for w = 1:2
    Tw = T(T.V1_window == wins{w}, :);
    % (a) example session: raw and partial coupling vs gain term
    subplot(2, 3, (w - 1) * 3 + 1); hold on
    E = ex.(wins{w});
    scatter(E.g, E.raw, 4, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.25);
    scatter(E.g, E.partial, 4, [0.85 0.33 0.1], 'filled', 'MarkerFaceAlpha', 0.25);
    xlabel('m_{V1} \times m_{HC} (cluster-MUA correlation)'); ylabel('pair coupling r');
    legend({sprintf('raw (\\rho=%.2f)', corr(E.raw, E.g, 'type', 'Spearman')), ...
        sprintf('partial (\\rho=%.2f)', corr(E.partial, E.g, 'type', 'Spearman'))}, 'Box', 'off', 'Location', 'northwest');
    title({sprintf('Example session %d', E.session), win_labels{w}}, 'FontWeight', 'normal', 'FontSize', 9);
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 9);
    % (b) per session: rho(coupling, gain term) raw vs partial, and gain vs preference
    subplot(2, 3, (w - 1) * 3 + 2); hold on
    vals = [Tw.gain_vs_raw Tw.gain_vs_partial Tw.gain_vs_pref];
    paired_panel(vals, Tw.animal_id, cols, {'raw', 'partial', 'gain vs pref'});
    ylabel('\rho with m_{V1} \times m_{HC}');
    g = G(G.V1_window == wins{w}, :);
    title({sprintf('raw p=%s, partial p=%s, pref p=%s', p_str(g.p_signrank(g.metric == "gain_vs_raw")), ...
        p_str(g.p_signrank(g.metric == "gain_vs_partial")), p_str(g.p_signrank(g.metric == "gain_vs_pref"))), ...
        sprintf('mean R^2 of raw r on gain term = %.2f', mean(Tw.R2_gain))}, 'FontWeight', 'normal', 'FontSize', 9);
    % (c) preference effect: raw, raw with the gain term removed, partial
    subplot(2, 3, (w - 1) * 3 + 3); hold on
    vals = [Tw.pref_raw Tw.pref_raw_minus_gain Tw.pref_partial];
    paired_panel(vals, Tw.animal_id, cols, {'raw', 'raw - gain', 'partial'});
    ylabel('\rho(coupling, z_{V1} z_{HC}) - null');
    title({sprintf('vs 0: raw p=%s, raw-gain p=%s, partial p=%s', p_str(g.p_signrank(g.metric == "pref_raw")), ...
        p_str(g.p_signrank(g.metric == "pref_raw_minus_gain")), p_str(g.p_signrank(g.metric == "pref_partial"))), ...
        sprintf('raw-gain vs raw: p=%s', p_str(g.p_signrank(g.metric == "pref_raw_minus_gain - pref_raw")))}, ...
        'FontWeight', 'normal', 'FontSize', 9);
end
sgtitle({'Raw pair coupling depends on the clusters'' cluster-MUA correlation; partialling the MUA removes it', cells_txt}, 'FontSize', 11);
savefig(f, fullfile(out_dir, 'gain_coupling_check.fig'));
exportgraphics(f, fullfile(out_dir, 'gain_coupling_check.pdf'), 'ContentType', 'vector');
fprintf('Saved %s\n', out_dir);
end

function c = col_corr(X, Y)
% Column-wise Pearson correlation of X(:, k) with Y(:, k)
X = X - mean(X, 1); Y = Y - mean(Y, 1);
c = sum(X .* Y, 1) ./ sqrt(sum(X .^ 2, 1) .* sum(Y .^ 2, 1));
end

function val = pref_minus_null(r, pv, ph, zV1, zHC, nPerm)
% Spearman rho(r, z_V1 * z_HC) across pairs minus the mean over permutations of the
% preferences among V1 cells and among HPC cells (as pairwise_V1_HC_coupling_POST)
rr = tiedrank(r);
obs = corr(rr, tiedrank(zV1(pv) .* zHC(ph)));
null = zeros(nPerm, 1);
for k = 1:nPerm
    zv = zV1(randperm(numel(zV1))); zh = zHC(randperm(numel(zHC)));
    null(k) = corr(rr, tiedrank(zv(pv) .* zh(ph)));
end
val = obs - mean(null);
end

function paired_panel(vals, animal, cols, labels)
n = size(vals, 2);
plot(1:n, vals', '-', 'Color', [0.85 0.85 0.85]);
for a = unique(animal)'
    idx = animal == a;
    for k = 1:n
        scatter(k * ones(sum(idx), 1), vals(idx, k), 18, cols(a, :), 'filled', 'MarkerFaceAlpha', 0.7);
    end
end
errorbar((1:n) + 0.25, mean(vals), std(vals) / sqrt(size(vals, 1)), 'ko', 'MarkerFaceColor', 'k', ...
    'LineWidth', 1.3, 'CapSize', 0);
yline(0, '--');
xlim([0.5 n + 0.6]);
set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'TickDir', 'out', 'Box', 'off', 'FontSize', 9);
end

function s = p_str(p)
if p < 0.001
    s = sprintf('%.1e', p);
else
    s = sprintf('%.3f', p);
end
end
