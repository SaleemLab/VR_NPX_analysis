function ar_decoder_comparison()
% Key population results side by side for every decoder variant
% (results\population\<decoder>\population_stats.csv), one marker per animal:
%   grey-screen bias toward the previous track (early / late grey, all transitions)
%   onset: switch - laps 6+ bias toward the upcoming track
% lobo3 decoders: drift would favour the block's track; lobo5: drift would favour
% the other track; _all: trained on all laps (grey epochs only are out of sample).
% Output: results\population\decoder_comparison.(fig/png/pdf) and .csv
P = ar_setup();
decs = {'bayes_multinomial', 'bayes_multinomial_lobo5', 'bayes_multinomial_all', 'bayes_poisson', ...
    'bayes_poisson_lobo5', 'pls', 'pls_lobo5', 'pls_all'};
tests = {'toward_prev_all', 'grey_early'; 'toward_prev_all', 'grey_late'; ...
    'toward_prev_all_nolick', 'grey_early'; 'toward_next_switch_minus_late', 'onset'; ...
    'toward_next_switch_minus_late', 'grey_late'};
regions = {'HPC', 'V1'};
acol = lines(6);
rows = {};
fig = figure('Position', [50 50 1700 330 * numel(regions)], 'Color', 'w');
tl = tiledlayout(numel(regions), size(tests, 1), 'TileSpacing', 'compact');
for r = 1:numel(regions)
    for k = 1:size(tests, 1)
        nexttile(tl); hold on;
        for d = 1:numel(decs)
            f = fullfile(P.results, 'population', decs{d}, 'population_stats.csv');
            if ~exist(f, 'file'), continue; end
            T = readtable(f, 'TextType', 'string');
            x = T(T.region == regions{r} & T.analysis == tests{k, 1} & T.measure == tests{k, 2}, :);
            if isempty(x), continue; end
            am = str2double(split(x.animal_means(1)));
            scatter(d + 0.08 * (1:numel(am))' - 0.25, am, 15, acol(1:numel(am), :), 'filled');
            errorbar(d, x.mean(1), x.sem(1), 'ko', 'MarkerFaceColor', 'k', 'CapSize', 0);
            if x.lme_p(1) < 0.05
                text(d, max([am; x.mean(1) + x.sem(1)]) * 1.1 + 0.01, '*', 'HorizontalAlignment', 'center');
            end
            rows(end + 1, :) = {string(regions{r}), string(tests{k, 1}), string(tests{k, 2}), string(decs{d}), ...
                x.mean(1), x.sem(1), x.p_signrank(1), x.n_sess_pos(1), x.n_sess(1), x.n_anim_pos(1), x.n_anim(1), ...
                x.lme_est(1), x.lme_p(1)}; %#ok<AGROW>
        end
        yline(0, 'k:');
        set(gca, 'XTick', 1:numel(decs), 'XTickLabel', strrep(decs, 'bayes_', ''), 'TickLabelInterpreter', 'none', ...
            'TickDir', 'out', 'Box', 'off'); xtickangle(45);
        title(sprintf('%s: %s, %s', regions{r}, tests{k, 1}, tests{k, 2}), 'Interpreter', 'none', 'FontWeight', 'normal');
        ylabel('z');
    end
end
title(tl, 'Decoder variants: mean +- SEM over sessions (black), animal means (colour), * LME p < 0.05');
out = fullfile(P.results, 'population', 'decoder_comparison');
ar_save_figure(fig, out);
T = cell2table(rows, 'VariableNames', {'region', 'analysis', 'measure', 'decoder', 'mean', 'sem', 'p_signrank', ...
    'n_sess_pos', 'n_sess', 'n_anim_pos', 'n_anim', 'lme_est', 'lme_p'});
writetable(T, [out '.csv']);
end
