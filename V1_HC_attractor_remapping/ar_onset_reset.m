function T = ar_onset_reset(varargin)
% Does track onset reset the ongoing track-bias state, or add to it?
% For a reference time t0, predictor = mean residual zLO over [t0 - 1, t0) and
% response = mean over [t0, t0 + 0.2) (onset-locked 100 ms bins; residual =
% block x track-pair mean removed per bin). Spearman across transitions per
% session. t0 inside the grey (-0.5 .. -0.2 s: both windows in the grey) gives
% the state's own continuity; t0 = 0 gives continuity across track onset.
% If onset resets the state, rho(t0 = 0) falls below the within-grey rho.
% Output: results\trial_prediction\onset_reset.(csv, fig/png/pdf)
p = inputParser;
addParameter(p, 'decoders', {'bayes_multinomial', 'pls'});
parse(p, varargin{:});
O = p.Results;
P = ar_setup();
S = load(fullfile(P.results, 'data', 'epoch_table.mat'), 'TC'); TC = S.TC;
bad = find(ismember(P.subjects, P.V1_exclude_subjects));
t = TC.t_on;
t0s = -0.5:0.1:0;
rows = {};
for d = 1:numel(O.decoders)
    for reg = {'HPC', 'V1'}
        sel = find(strcmp({TC.on.region}, reg{1}) & strcmp({TC.on.decoder}, O.decoders{d}));
        for i = sel
            s = TC.on(i).session; Ds = TC.trans(TC.trans.session == s, :);
            if strcmp(reg{1}, 'V1') && ismember(Ds.animal(1), bad), continue; end
            z = double(TC.on(i).z);
            g = findgroups(Ds.block, Ds.prev_track, Ds.track);
            for k = unique(g)'
                z(g == k, :) = z(g == k, :) - mean(z(g == k, :), 1, 'omitnan');
            end
            rho = NaN(1, numel(t0s));
            for j = 1:numel(t0s)
                x = mean(z(:, t >= t0s(j) - 1 & t < t0s(j)), 2, 'omitnan');
                y = mean(z(:, t >= t0s(j) & t < t0s(j) + 0.2), 2, 'omitnan');
                rho(j) = corr(x, y, 'Type', 'Spearman', 'Rows', 'complete');
            end
            rows(end + 1, :) = {string(O.decoders{d}), string(reg{1}), s, Ds.animal(1), rho}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'decoder', 'region', 'session', 'animal', 'rho'});
% figure and stats: across-onset (t0 = 0) minus within-grey mean (t0 = -0.5..-0.2)
acol = lines(6);
fig = figure('Position', [50 50 1400 420], 'Color', 'w');
tl = tiledlayout(1, 3, 'TileSpacing', 'compact');
st = {};
regs = {'HPC', 'V1'};
for r = 1:2
    nexttile(tl); hold on;
    for d = 1:numel(O.decoders)
        x = T(T.decoder == O.decoders{d} & T.region == regs{r}, :);
        R = vertcat(x.rho);
        m = mean(R, 1, 'omitnan'); se = std(R, 0, 1, 'omitnan') / sqrt(height(x));
        c = [0 0 0; 0.2 0.45 0.75]; c = c(d, :);
        errorbar(t0s * 1000, m, se, 'o-', 'Color', c, 'MarkerFaceColor', c, 'CapSize', 0, 'LineWidth', 1.5);
        within = mean(R(:, t0s <= -0.2 + 1e-9), 2); across = R(:, end);
        dd = across - within;
        [ua, ~, ia] = unique(x.animal); am = accumarray(ia, dd, [], @mean);
        st(end + 1, :) = {string(O.decoders{d}), string(regs{r}), mean(within), mean(across), mean(dd), signrank(dd), ...
            sum(dd > 0), numel(dd), sum(am > 0), numel(am)}; %#ok<AGROW>
    end
    xline(0, 'k-'); yline(0, 'k:'); ylim([0 0.6]);
    xlabel('t0 (ms from onset): 1 s before t0 predicts the 200 ms after t0'); ylabel('Spearman rho (residual bias)');
    legend(strrep(O.decoders, '_', ' '), 'Box', 'off', 'Location', 'southwest');
    title(sprintf('%s: continuity of the track-bias state', regs{r}), 'FontWeight', 'normal');
end
nexttile(tl); hold on;
for r = 1:2
    x = T(T.decoder == O.decoders{1} & T.region == regs{r}, :);
    R = vertcat(x.rho); dd = R(:, end) - mean(R(:, t0s <= -0.2 + 1e-9), 2);
    scatter(r + 0.1 * randn(numel(dd), 1), dd, 18, acol(x.animal, :), 'filled', 'MarkerFaceAlpha', 0.7);
    errorbar(r + 0.3, mean(dd), std(dd) / sqrt(numel(dd)), 'ko', 'MarkerFaceColor', 'k', 'CapSize', 0, 'LineWidth', 1.5);
end
yline(0, 'k:'); set(gca, 'XTick', 1:2, 'XTickLabel', regs); xlim([0.5 2.5]);
ylabel('rho across onset - rho within grey');
title({'Does onset reset the state?', '(< 0: partial reset; multinomial Bayes)'}, 'FontWeight', 'normal');
out_dir = fullfile(P.results, 'trial_prediction');
ar_save_figure(fig, fullfile(out_dir, 'onset_reset'));
close(fig);
S2 = cell2table(st, 'VariableNames', {'decoder', 'region', 'rho_within_grey', 'rho_across_onset', 'difference', 'p_signrank', ...
    'n_pos', 'n', 'n_anim_pos', 'n_anim'});
writetable(S2, fullfile(out_dir, 'onset_reset.csv'));
disp(S2)
end
