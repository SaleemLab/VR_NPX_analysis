%% Interhemispheric SO synchrony at ripple times: PRE vs POST sleep
% Hypothesis: after the (hemifield-lateralised) task, SO become less bilaterally
% synchronised, i.e. fewer ripples occur with both hemispheres in SO trough.
% SO phase of each hemisphere at each non-duplicated ripple event (SO_phase n x 2),
% peak [-pi/2, pi/2] vs trough as in hemisphere_conditions.
% Per session and sleep:
%   PLV            : |mean(exp(i*(phase1 - phase2)))|, interhemispheric phase locking
%   phi            : correlation of the two hemispheres' trough indicators
%                    (synchrony controlling for each hemisphere's trough fraction)
%   p_both_trough  : fraction of events with both hemispheres in trough
%   p_mixed        : fraction with one hemisphere in trough and the other in peak
%   p_bilateral_given_trough : both trough / at least one trough
%   p_trough_h1/h2 : fraction in trough per hemisphere (marginal control, two-sided only)
% PRE vs POST: paired signed-rank and exact sign-flip permutation (12 sessions with both),
% unpaired rank-sum (all sessions). One-sided direction: POST less synchronised
% (PLV, phi, p_both_trough, p_bilateral_given_trough lower; p_mixed higher).
% Repeated for ripples in the first 30 min of cumulative SWS and for equal SWS in PRE
% and POST (sws_epoch_masks): *_by_SWS_epoch.csv, *_metrics_matched_first30 /
% _equalSWS and SO_bilateral_synchrony_by_SWS_epoch figures.
% Output: nature_neuroscience_revision_analysis\SO_bilateral_synchrony_PRE_vs_POST
function analyse_SO_bilateral_synchrony_PRE_vs_POST()
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    out_dir = fullfile(rev_dir, 'SO_bilateral_synchrony_PRE_vs_POST');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    load(fullfile(rev_dir, 'ripple_info_PRE_unique.mat'), 'ripple_info_PRE');
    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    D(1).name = "PRE";  D(1).info = ripple_info_PRE;
    D(2).name = "POST"; D(2).info = ripple_info;

    % --- Event selection by cumulative SWS epoch (all / first30 / matched, see sws_epoch_masks) ---
    sws_cap_min = 30;
    epochs = {'all', 'first30', 'matched'};
    epoch_labels = {'All SWS', sprintf('First %d min SWS', sws_cap_min), 'Equal SWS (PRE = POST)'};
    epoch_files = {'', '_first30', '_equalSWS'};
    E = sws_epoch_masks(D(1).info, D(2).info, analysis_folder, sws_cap_min);

    % --- Per-session metrics and PRE vs POST tests, per SWS epoch ---
    edges_d = linspace(-pi, pi, 25);                 % phase-difference histogram
    edges_2d = linspace(-pi, pi, 13);                % joint phase histogram
    metrics = {'PLV', 'phi', 'p_both_trough', 'p_mixed', 'p_bilateral_given_trough', 'p_trough_h1', 'p_trough_h2'};
    metric_labels = {'Interhemispheric SO PLV', 'Trough correlation (phi)', 'P(both trough)', ...
        'P(one trough, one peak)', 'P(both trough | any trough)', 'P(trough), hemisphere 1', 'P(trough), hemisphere 2'};
    tails = {'left', 'left', 'left', 'right', 'left', 'both', 'both'};
    TE = table(); SE = table(); R = struct();
    for e = 1:numel(epochs)
        [Te, He] = session_metrics(D, {E(1).sel.(epochs{e}), E(2).sel.(epochs{e})}, edges_d, edges_2d);
        Se = table();
        for m = 1:numel(metrics)
            [x, y] = pre_post_matrix(Te, metrics{m});
            Se = [Se; prepost_tests(metrics{m}, x, y, tails{m})]; %#ok<AGROW>
        end
        R.(epochs{e}) = struct('T', Te, 'S', Se, 'H', He);
        TE = [TE; [table(repmat(string(epochs{e}), height(Te), 1), 'VariableNames', {'SWS_epoch'}) Te]]; %#ok<AGROW>
        SE = [SE; [table(repmat(string(epochs{e}), height(Se), 1), 'VariableNames', {'SWS_epoch'}) Se]]; %#ok<AGROW>
    end
    disp(SE);
    writetable(TE, fullfile(out_dir, 'SO_bilateral_synchrony_per_session_by_SWS_epoch.csv'));
    writetable(SE, fullfile(out_dir, 'SO_bilateral_synchrony_stats_by_SWS_epoch.csv'));

    % All SWS: outputs as before
    T = R.all.T; S = R.all.S; H = R.all.H;
    writetable(T, fullfile(out_dir, 'SO_bilateral_synchrony_per_session.csv'));
    writetable(S, fullfile(out_dir, 'SO_bilateral_synchrony_stats.csv'));
    save(fullfile(out_dir, 'SO_bilateral_synchrony.mat'), 'T', 'S', 'H', 'R', 'sws_cap_min');

    % --- Plot 1: per-session metrics, PRE vs POST ---
    sleep_colors = {[0.5 0.5 0.5], [231, 41, 138] / 256};
    fig = figure('Position', [100 100 1500 700]);
    tiledlayout(2, 4, 'TileSpacing', 'compact');
    for m = 1:numel(metrics)
        ax = nexttile; hold on;
        [x, y] = pre_post_matrix(T, metrics{m});
        ok = isfinite(x) & isfinite(y);
        plot([1 2], [x(ok) y(ok)]', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
        v = {x(isfinite(x)), y(isfinite(y))};
        for iS = 1:2
            scatter(iS * ones(size(v{iS})), v{iS}, 25, sleep_colors{iS}, 'filled', 'MarkerFaceAlpha', 0.7);
            errorbar(iS + 0.3, mean(v{iS}), std(v{iS}) / sqrt(numel(v{iS})), 'o', 'Color', 'k', ...
                'MarkerFaceColor', sleep_colors{iS}, 'LineWidth', 1.5, 'CapSize', 6);
        end
        st = S(S.metric == metrics{m}, :);
        set(ax, 'XTick', [1 2], 'XTickLabel', {'PRE', 'POST'}, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
        xlim([0.4 2.6]);
        ylabel(metric_labels{m});
        if strcmp(tails{m}, 'both')
            p_line = sprintf('perm p = %.3g, signrank p = %.3g (two-sided)', st.p_perm_two_sided, st.p_signrank_two_sided);
        else
            p_line = sprintf('%s: perm p = %.3g, signrank p = %.3g', st.direction, st.p_perm_one_sided, st.p_signrank_one_sided);
        end
        title({p_line, sprintf('n = %d paired; rank-sum p = %.3g', st.n_paired, st.p_ranksum_one_sided)}, ...
            'FontSize', 8, 'FontWeight', 'normal');
    end
    sgtitle('Interhemispheric SO synchrony at ripple times: PRE vs POST', 'FontSize', 12);
    savefig(fig, fullfile(out_dir, 'SO_bilateral_synchrony_metrics.fig'));
    exportgraphics(fig, fullfile(out_dir, 'SO_bilateral_synchrony_metrics.png'), 'Resolution', 300);

    % --- Plot 1b: matched sessions only (the points used by the paired tests), per SWS epoch ---
    for e = 1:numel(epochs)
        Te = R.(epochs{e}).T; Se = R.(epochs{e}).S;
        fig = figure('Position', [100 100 1500 700]);
        tiledlayout(2, 4, 'TileSpacing', 'compact');
        for m = 1:numel(metrics)
            ax = nexttile;
            [x, y] = pre_post_matrix(Te, metrics{m});
            n = plot_matched_pre_post(ax, x, y, [1 2]);
            set(ax, 'XTick', [1 2], 'XTickLabel', {'PRE', 'POST'}, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
            xlim(ax, [0.4 2.6]);
            ylabel(ax, metric_labels{m});
            title(ax, {p_line_text(Se(Se.metric == metrics{m}, :), tails{m}), sprintf('n = %d matched sessions', n)}, ...
                'FontSize', 8, 'FontWeight', 'normal');
        end
        sgtitle(sprintf('Interhemispheric SO synchrony at ripple times: PRE vs POST, matched sessions (%s)', epoch_labels{e}), ...
            'FontSize', 12);
        savefig(fig, fullfile(out_dir, ['SO_bilateral_synchrony_metrics_matched' epoch_files{e} '.fig']));
        exportgraphics(fig, fullfile(out_dir, ['SO_bilateral_synchrony_metrics_matched' epoch_files{e} '.png']), 'Resolution', 300);
    end

    % --- Plot 1c: key synchrony metrics by SWS epoch (rows), matched sessions ---
    key = {'PLV', 'phi', 'p_both_trough', 'p_mixed', 'p_bilateral_given_trough'};
    fig = figure('Position', [50 50 1750 1000]);
    tiledlayout(numel(epochs), numel(key), 'TileSpacing', 'compact');
    for e = 1:numel(epochs)
        Te = R.(epochs{e}).T; Se = R.(epochs{e}).S;
        for kk = 1:numel(key)
            m = find(strcmp(metrics, key{kk}));
            ax = nexttile;
            [x, y] = pre_post_matrix(Te, metrics{m});
            n = plot_matched_pre_post(ax, x, y, [1 2]);
            set(ax, 'XTick', [1 2], 'XTickLabel', {'PRE', 'POST'}, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
            xlim(ax, [0.4 2.6]);
            ylabel(ax, metric_labels{m}, 'FontSize', 9);
            title(ax, {epoch_labels{e}, p_line_text(Se(Se.metric == metrics{m}, :), tails{m}), ...
                sprintf('n = %d matched sessions', n)}, 'FontSize', 8, 'FontWeight', 'normal');
        end
    end
    sgtitle('Interhemispheric SO synchrony at ripple times by SWS epoch: PRE vs POST, matched sessions', 'FontSize', 12);
    savefig(fig, fullfile(out_dir, 'SO_bilateral_synchrony_by_SWS_epoch.fig'));
    exportgraphics(fig, fullfile(out_dir, 'SO_bilateral_synchrony_by_SWS_epoch.png'), 'Resolution', 300);

    % --- Plot 2: phase-difference distribution and joint phase distribution ---
    pre_ids = [H([H.sleep] == "PRE").session_id];
    fig = figure('Position', [100 100 1500 420]);
    tiledlayout(1, 4, 'TileSpacing', 'compact');
    nexttile; hold on;
    centers = edges_d(1:end-1) + diff(edges_d) / 2;
    groups = {[H.sleep] == "PRE", [H.sleep] == "POST" & ismember([H.session_id], pre_ids), [H.sleep] == "POST"};
    group_names = {'PRE', 'POST (matched)', 'POST (all)'};
    group_colors = {sleep_colors{1}, sleep_colors{2}, sleep_colors{2} * 0.6 + 0.4};
    for g = 1:3
        M = vertcat(H(groups{g}).dphase_hist);
        mu = mean(M, 1); se = std(M, [], 1) / sqrt(size(M, 1));
        fill([centers fliplr(centers)], [mu - se fliplr(mu + se)], group_colors{g}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        plot(centers, mu, '-', 'Color', group_colors{g}, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%s (n = %d)', group_names{g}, size(M, 1)));
    end
    xlim([-pi pi]); set(gca, 'XTick', [-pi 0 pi], 'XTickLabel', {'-\pi', '0', '\pi'}, 'TickDir', 'out', 'Box', 'off');
    xlabel('SO phase difference (hemi 1 - hemi 2)'); ylabel('Fraction of ripples');
    legend('Location', 'best', 'Box', 'off'); title('Phase difference (mean \pm SEM over sessions)', 'FontWeight', 'normal');

    J_pre = mean(cat(3, H(groups{1}).joint_hist), 3);
    J_post = mean(cat(3, H(groups{2}).joint_hist), 3);
    J = {J_pre, J_post, J_post - J_pre};
    J_titles = {'PRE', 'POST (matched sessions)', 'POST - PRE'};
    clim_abs = max([J_pre(:); J_post(:)]);
    for k = 1:3
        nexttile;
        imagesc(edges_2d([1 end]), edges_2d([1 end]), J{k}');
        axis xy square;
        if k < 3
            colormap(gca, parula); caxis([0 clim_abs]);
        else
            lim = max(abs(J{3}(:))); colormap(gca, redblue_map()); caxis([-lim lim]);
        end
        colorbar;
        hold on;
        xline([-pi/2 pi/2], 'w:'); yline([-pi/2 pi/2], 'w:');
        set(gca, 'XTick', [-pi 0 pi], 'XTickLabel', {'-\pi', '0', '\pi'}, 'YTick', [-pi 0 pi], ...
            'YTickLabel', {'-\pi', '0', '\pi'}, 'TickDir', 'out');
        xlabel('SO phase hemi 1'); ylabel('SO phase hemi 2');
        title(J_titles{k}, 'FontWeight', 'normal');
    end
    sgtitle('SO phase at ripple times (trough = |phase| > \pi/2)', 'FontSize', 12);
    savefig(fig, fullfile(out_dir, 'SO_bilateral_synchrony_distributions.fig'));
    exportgraphics(fig, fullfile(out_dir, 'SO_bilateral_synchrony_distributions.png'), 'Resolution', 300);
    fprintf('Saved to %s\n', out_dir);
end

function [T, H] = session_metrics(D, sel, edges_d, edges_2d)
% Per-session synchrony metrics and phase histograms over the events in sel{iS}
T = table(); H = struct('sleep', {}, 'session_id', {}, 'dphase_hist', {}, 'joint_hist', {});
for iS = 1:2
    info = D(iS).info;
    sess = info.session_id(:);
    for s = unique(sess(sel{iS}))'
        in_s = sess == s & sel{iS};
        ph = info.SO_phase(in_s, :);
        ok = all(isfinite(ph), 2);
        ph = ph(ok, :);
        trough = (ph >= -pi & ph < -pi/2) | (ph > pi/2 & ph <= pi);
        dph = angle(exp(1i * (ph(:, 1) - ph(:, 2))));
        n_any = sum(any(trough, 2));
        subj = info.subject_id(find(in_s, 1));
        T = [T; table(s, subj, D(iS).name, size(ph, 1), ...
            abs(mean(exp(1i * dph))), angle(mean(exp(1i * dph))), ...
            corr(double(trough(:, 1)), double(trough(:, 2))), ...
            mean(all(trough, 2)), mean(xor(trough(:, 1), trough(:, 2))), ...
            sum(all(trough, 2)) / n_any, mean(trough(:, 1)), mean(trough(:, 2)), ...
            'VariableNames', {'session_id', 'subject_id', 'sleep', 'n_events', 'PLV', 'mean_phase_diff', ...
            'phi', 'p_both_trough', 'p_mixed', 'p_bilateral_given_trough', 'p_trough_h1', 'p_trough_h2'})]; %#ok<AGROW>
        H(end+1) = struct('sleep', D(iS).name, 'session_id', s, ...
            'dphase_hist', histcounts(dph, edges_d, 'Normalization', 'probability'), ...
            'joint_hist', histcounts2(ph(:, 1), ph(:, 2), edges_2d, edges_2d, 'Normalization', 'probability')); %#ok<AGROW>
    end
end
end

function s = p_line_text(st, tail)
if strcmp(tail, 'both')
    s = sprintf('perm p = %.3g, signrank p = %.3g (two-sided)', st.p_perm_two_sided, st.p_signrank_two_sided);
else
    s = sprintf('%s: perm p = %.3g, signrank p = %.3g', st.direction, st.p_perm_one_sided, st.p_signrank_one_sided);
end
end

function [x, y] = pre_post_matrix(T, metric)
% PRE (x) and POST (y) values indexed by session_id; NaN where missing
n = max(T.session_id);
x = NaN(n, 1); y = NaN(n, 1);
Tp = T(T.sleep == "PRE", :);  x(Tp.session_id) = Tp.(metric);
Tq = T(T.sleep == "POST", :); y(Tq.session_id) = Tq.(metric);
end

function row = prepost_tests(metric, x, y, tail)
% Paired signed-rank and permutation (sessions with both), unpaired rank-sum (all)
ok = isfinite(x) & isfinite(y);
d = y(ok) - x(ok);
ps2 = signrank(y(ok), x(ok));
[pp2, pp1] = paired_perm_test(d, tail_or_right(tail));
xa = x(isfinite(x)); ya = y(isfinite(y));
pr2 = ranksum(ya, xa);
if strcmp(tail, 'both')
    ps1 = NaN; pp1 = NaN; pr1 = NaN; direction = "two-sided";
else
    ps1 = signrank(y(ok), x(ok), 'tail', tail);
    pr1 = ranksum(ya, xa, 'tail', tail);
    if strcmp(tail, 'left'), direction = "POST < PRE"; else, direction = "POST > PRE"; end
end
row = table(string(metric), sum(ok), mean(x(ok)), mean(y(ok)), sum(d > 0), ...
    ps2, ps1, pp2, pp1, numel(xa), numel(ya), mean(xa), mean(ya), pr2, pr1, direction, ...
    'VariableNames', {'metric', 'n_paired', 'mean_PRE', 'mean_POST', 'n_POST_gt_PRE', ...
    'p_signrank_two_sided', 'p_signrank_one_sided', 'p_perm_two_sided', 'p_perm_one_sided', ...
    'n_pre_all', 'n_post_all', 'mean_PRE_all', 'mean_POST_all', 'p_ranksum_two_sided', 'p_ranksum_one_sided', 'direction'});
end

function t = tail_or_right(tail)
if strcmp(tail, 'both'), t = 'right'; else, t = tail; end
end

function cmap = redblue_map()
% Blue (negative) - white - red (positive) diverging colormap
n = 128;
r = [linspace(0.15, 1, n) ones(1, n)]';
g = [linspace(0.35, 1, n) linspace(1, 0.2, n)]';
b = [ones(1, n) linspace(1, 0.2, n)]';
cmap = [r g b];
end
