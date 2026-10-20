function AUC_low_high_panel(ax, R, labels, colors, window_labels, y_lim, diff_p, bracket_y)
% One panel of the low vs high AUC bar plot (plot_AUC_low_vs_high style).
% R(w): the two conditions' AUC in V1 window w, fields m (1x2 mean), ci (2x2,
% rows 2.5 / 97.5 prctile, columns conditions), ms and cis (shuffled, same shape).
% Per window: Shuffle (mean of the two shuffled means; CI pooled as min lower,
% max upper), condition A (pale) and B bars with 95% CI. '*': condition 95% CI
% does not overlap the shuffled 95% CI; 'n.s.' otherwise.
% diff_p (optional, 1 x windows): B vs A p shown above a bracket at height
% bracket_y (default: 0.12 above the higher condition CI of that window).
    hold(ax, 'on');
    group_x = [1 2 3];     % Shuffle, A, B within a window
    group_gap = 4;
    xt = []; xtl = {};
    for w = 1:numel(R)
        x = group_x + (w - 1) * group_gap;
        ci_s = [min(R(w).cis(1, :)) max(R(w).cis(2, :))];
        m_s = mean(R(w).ms);
        errorbar(ax, x(1), m_s, m_s - ci_s(1), ci_s(2) - m_s, 'k', 'LineWidth', 1.5, 'CapSize', 8);
        for k = 1:2
            bar(ax, x(k + 1), R(w).m(k), 0.7, 'FaceColor', colors(k, :), 'EdgeColor', 'none');
            errorbar(ax, x(k + 1), R(w).m(k), R(w).m(k) - R(w).ci(1, k), R(w).ci(2, k) - R(w).m(k), ...
                'k', 'LineWidth', 1.5, 'CapSize', 8);
            if R(w).ci(1, k) > ci_s(2) || R(w).ci(2, k) < ci_s(1)
                text(ax, x(k + 1), R(w).ci(2, k) + 0.008, '*', 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', 'FontSize', 16);
            else
                text(ax, x(k + 1), R(w).ci(2, k) + 0.008, 'n.s.', 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', 'FontSize', 11);
            end
        end
        if nargin >= 7 && ~isempty(diff_p)
            if nargin >= 8 && ~isempty(bracket_y)
                yb = bracket_y;
            else
                yb = max(R(w).ci(2, :)) + 0.12;   % clears the '*' / 'n.s.' markers
            end
            plot(ax, x([2 2 3 3]), yb + [-0.01 0 0 -0.01], 'k', 'LineWidth', 1);
            if diff_p(w) < 0.001
                p_str = 'p < 0.001';
            else
                p_str = sprintf('p = %.3f', diff_p(w));
            end
            text(ax, mean(x(2:3)), yb + 0.005, p_str, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', 10);
        end
        text(ax, mean(x), y_lim(2), window_labels{w}, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', 'FontSize', 11);
        xt = [xt x]; %#ok<AGROW>
        xtl = [xtl {'Shuffle'} labels]; %#ok<AGROW>
    end
    yline(ax, 0, 'k--');
    set(ax, 'XTick', xt, 'XTickLabel', xtl, 'XTickLabelRotation', 45);
    xlim(ax, [0.3, xt(end) + 0.7]);
    ylim(ax, y_lim);
    yticks(ax, -0.2:0.1:0.3);
    ylabel(ax, 'Reactivation coherence AUC');
    set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12);
end
