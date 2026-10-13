function n = plot_matched_pre_post(ax, pre, post, xpos)
% PRE vs POST for sessions with both values only: paired lines, points and
% mean +/- SEM of the matched sessions. pre/post indexed by session; xpos = [x_PRE x_POST].
sleep_colors = {[0.5 0.5 0.5], [231, 41, 138] / 256};
ok = isfinite(pre) & isfinite(post);
v = {pre(ok), post(ok)};
n = sum(ok);
hold(ax, 'on');
plot(ax, xpos, [v{1} v{2}]', '-', 'Color', [0 0 0 0.3], 'LineWidth', 1);
dx = [-0.3 0.3];
for iS = 1:2
    scatter(ax, xpos(iS) * ones(n, 1), v{iS}, 25, sleep_colors{iS}, 'filled', 'MarkerFaceAlpha', 0.7);
    errorbar(ax, xpos(iS) + dx(iS), mean(v{iS}), std(v{iS}) / sqrt(n), 'o', 'Color', 'k', ...
        'MarkerFaceColor', sleep_colors{iS}, 'LineWidth', 1.5, 'CapSize', 6);
end
end
