function plot_pairwise_V1_HC_coupling(gain_control, cell_set, sleep)
% Group statistics and figures of pairwise_V1_HC_coupling_POST(gain_control) for cell_set
% 'lt10' (default, cells < 10 Hz) or 'all_cells', and sleep 'POST' (default) or 'PRE'
% (pairwise_cell_set_folders).
% Per V1 window and measure (raw, partial), session statistic = rho (Spearman of pair
% coupling vs z_V1 * z_HC) minus the mean of its permutation null; same for the
% same - opposite preference difference.
%   group permutation : mean over sessions of the observed statistic vs the mean over
%                       sessions of each permutation (one-sided, > null)
%   signrank          : session values vs 0
%   hierarchical boot : animals, then sessions within animal, with replacement (two-sided,
%                       p floored at 1/nBoot; the group permutation p is floored at 1/(nPerm+1))
% Figures (.fig + vector PDF) and group_stats.csv in the gain_<gain_control>[_all_cells] folder.

if nargin < 2, cell_set = 'lt10'; end
if nargin < 3, sleep = 'POST'; end
rev_dir = addon_setup();
[counts_name, out_name, cells_txt] = pairwise_cell_set_folders(gain_control, cell_set, sleep);
out_dir = fullfile(rev_dir, 'pairwise_V1_HC_coupling', out_name);
load(fullfile(out_dir, 'coupling_POST.mat'), 'S');
nBoot = 10000;
wins = {'V1_pre', 'V1_post'};
win_labels = {'V1 -0.2-0 s, HC 0-0.1 s', 'V1 0-0.2 s, HC 0-0.1 s'};
measures = {'raw', 'partial'};
meas_labels = {'Raw correlation', sprintf('Partial correlation (%s)', strrep(gain_control, '_', ' '))};
if strcmp(gain_control, 'total_V1')
    gain_label = sprintf('gain: total V1 MUA + HC MUA (leave-one-out); %s', cells_txt);
else
    gain_label = sprintf('gain: own-hemisphere V1 MUA + HC MUA (leave-one-out); %s', cells_txt);
end

animal = [S.subject_id]';
animals = unique(animal);
cols = lines(numel(animals));
pop = population_coherence(rev_dir, S, counts_name, sleep);

% Per-session values used by the figures (also read by plot_pairwise_summary)
SS = table();
for w = 1:2
    for m = 1:2
        st = arrayfun(@(s) s.(wins{w}).(measures{m}), S);
        SS.(sprintf('rho_minus_null_%s_%s', measures{m}, wins{w})) = ...
            [st.rho]' - arrayfun(@(x) mean(x.rho_null, 'omitnan'), st)';
    end
    SS.(sprintf('KDE_logodds_r_%s', wins{w})) = pop.(wins{w});
end
SS = [table([S.session_count]', animal, 'VariableNames', {'session_id', 'animal_id'}) SS];
writetable(SS, fullfile(out_dir, 'session_summary.csv'));

G = table();
f1 = figure('Name', 'pairwise_preference_rho', 'Position', [50 50 900 750], 'Color', 'w');
f2 = figure('Name', 'pairwise_same_vs_opposite', 'Position', [80 80 900 750], 'Color', 'w');
for m = 1:2
    for w = 1:2
        st = arrayfun(@(s) s.(wins{w}).(measures{m}), S);
        stat_names = {'rho', 'diff'};
        for k = 1:2
            obs = [st.(stat_names{k})]';
            null = [st.([stat_names{k} '_null'])];             % nPerm x sessions
            val = obs - mean(null, 1, 'omitnan')';
            ok = isfinite(val);
            null_group = mean(null(:, ok), 2);
            p_perm = (1 + sum(null_group >= mean(obs(ok)))) / (1 + numel(null_group));
            p_sr = signrank(val(ok));
            [bm, bci, p_boot] = hier_boot(val(ok), animal(ok), nBoot);
            animal_mean = arrayfun(@(a) mean(val(ok & animal == a)), animals);
            G = [G; table(string(gain_control), string(measures{m}), string(wins{w}), string(stat_names{k}), ...
                sum(ok), mean(val(ok)), bm, bci(1), bci(2), p_boot, p_perm, p_sr, sum(val(ok) > 0), ...
                sum(animal_mean > 0), numel(animals), ...
                'VariableNames', {'gain_control', 'measure', 'V1_window', 'statistic', 'n_sessions', ...
                'mean_minus_null', 'boot_mean', 'boot_CI_lo', 'boot_CI_hi', 'p_hier_boot', ...
                'p_group_perm', 'p_signrank', 'n_sessions_pos', 'n_animals_pos', 'n_animals'})]; %#ok<AGROW>
            if k == 1, figure(f1); else, figure(f2); end
            subplot(2, 2, (m - 1) * 2 + w); hold on
            if k == 1
                ylab = '\rho(coupling, z_{V1} z_{HC}) - null';
            else
                ylab = 'Same - opposite preference r (- null)';
            end
            plot_sessions(val, animal, animals, cols, bm, bci);
            ylabel(ylab);
            title({meas_labels{m}, win_labels{w}, sprintf('boot p=%.3g, perm p=%.3g, signrank p=%.3g', ...
                p_boot, p_perm, p_sr)}, 'FontWeight', 'normal', 'FontSize', 9);
        end
    end
end
figure(f1); sgtitle({'Pair coupling vs shared track preference, POST ripples', gain_label}, 'FontSize', 11);
figure(f2); sgtitle({'Same- vs opposite-preference V1-HC pairs, POST ripples', gain_label}, 'FontSize', 11);
writetable(G, fullfile(out_dir, 'group_stats.csv'));
disp(G);

% Binned coupling vs preference product (pooled pairs, session-centred)
f3 = figure('Name', 'pairwise_binned_preference', 'Position', [110 110 900 400], 'Color', 'w');
nbin = 10;
for w = 1:2
    subplot(1, 2, w); hold on
    for m = 1:2
        r_all = []; pz_all = [];
        for s = 1:numel(S)
            W = S(s).(wins{w});
            r = W.(['r_' measures{m}]);
            pz = W.z_V1(:) * W.z_HC(:)';
            ok = isfinite(r);
            r_all = [r_all; r(ok) - mean(r(ok))]; %#ok<AGROW>
            pz_all = [pz_all; pz(ok)]; %#ok<AGROW>
        end
        edges = prctile(pz_all, linspace(0, 100, nbin + 1));
        b = discretize(pz_all, edges);
        mu = accumarray(b, r_all, [nbin 1], @mean);
        se = accumarray(b, r_all, [nbin 1], @(x) std(x) / sqrt(numel(x)));
        xc = accumarray(b, pz_all, [nbin 1], @median);
        errorbar(xc, mu, se, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4, 'CapSize', 0);
    end
    xline(0, '--'); yline(0, '--');
    xlabel('z_{V1} \times z_{HC} (decile median)'); ylabel('Pair coupling r (session-centred)');
    legend(measures, 'Box', 'off', 'Location', 'northwest');
    title(win_labels{w}, 'FontWeight', 'normal');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
end
sgtitle({'Pooled pairs, 10 quantile bins', gain_label}, 'FontSize', 11);

% Session pairwise effect vs population log-odds coherence
f4 = figure('Name', 'pairwise_vs_population_coherence', 'Position', [140 140 900 400], 'Color', 'w');
for w = 1:2
    subplot(1, 2, w); hold on
    st = arrayfun(@(s) s.(wins{w}).partial, S);
    val = [st.rho]' - arrayfun(@(x) mean(x.rho_null, 'omitnan'), st)';
    pc = pop.(wins{w});
    for a = 1:numel(animals)
        idx = animal == animals(a);
        scatter(pc(idx), val(idx), 30, cols(a, :), 'filled');
    end
    ok = isfinite(val) & isfinite(pc);
    [rr, pp] = corr(pc(ok), val(ok), 'type', 'Spearman');
    xline(0, '--'); yline(0, '--');
    xlabel('Population coherence r(KDE V1, KDE HC)'); ylabel('Pairwise \rho - null (partial)');
    title(sprintf('%s: Spearman r=%.2f, p=%.3g', win_labels{w}, rr, pp), 'FontWeight', 'normal', 'FontSize', 9);
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
end
sgtitle({'Per session: pairwise effect vs population reactivation coherence', gain_label}, 'FontSize', 11);

figs = [f1 f2 f3 f4];
for f = figs
    name = get(f, 'Name');
    savefig(f, fullfile(out_dir, [name '.fig']));
    exportgraphics(f, fullfile(out_dir, [name '.pdf']), 'ContentType', 'vector');
end
fprintf('Saved figures and group_stats.csv in %s\n', out_dir);
end

function plot_sessions(val, animal, animals, cols, bm, bci)
for a = 1:numel(animals)
    idx = animal == animals(a) & isfinite(val);
    x = a + 0.15 * (rand(sum(idx), 1) - 0.5);
    scatter(x, val(idx), 25, cols(a, :), 'filled', 'MarkerFaceAlpha', 0.6);
    plot(a + [-0.25 0.25], mean(val(idx)) * [1 1], '-', 'Color', cols(a, :), 'LineWidth', 2);
end
xp = numel(animals) + 1.2;
errorbar(xp, bm, bm - bci(1), bci(2) - bm, 'ko', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'CapSize', 0);
yline(0, '--');
xlim([0.4 xp + 0.6]);
set(gca, 'XTick', [1:numel(animals) xp], 'XTickLabel', [compose('A%d', 1:numel(animals)) {'All'}], ...
    'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
end

function [bm, ci, p] = hier_boot(val, animal, nBoot)
% Mean over resampled animals -> resampled sessions within animal
animals = unique(animal);
nA = numel(animals);
by_animal = arrayfun(@(a) val(animal == a), animals, 'UniformOutput', false);
rs = RandStream('twister', 'Seed', 1);
b = NaN(nBoot, 1);
for i = 1:nBoot
    pick = randi(rs, nA, nA, 1);
    v = [];
    for a = pick'
        x = by_animal{a};
        v = [v; x(randi(rs, numel(x), numel(x), 1))]; %#ok<AGROW>
    end
    b(i) = mean(v);
end
bm = mean(val);
ci = prctile(b, [2.5 97.5]);
p = max(1 / nBoot, min(1, 2 * min(mean(b <= 0), mean(b >= 0))));  % floored at 1/nBoot
end

function pop = population_coherence(rev_dir, S, counts_name, sleep)
% Per session Pearson r between KDE V1 log odds (pre / post window) and KDE HC log odds
% (0-0.1 s) over the same ripples (rows of ripple_info, or ripple_info_PRE_unique for PRE)
if strcmp(sleep, 'PRE')
    L = load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_PRE_unique.mat'), 'KDE_reactivation_ripples_PSTH_PRE');
    K = L.KDE_reactivation_ripples_PSTH_PRE;
else
    L = load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
    K = L.KDE_reactivation_ripples_PSTH;
end
HC = (K.HPC_z_logodds_ripples + K.nan_mask)';
V1 = (K.V1_z_logodds_ripples + K.nan_mask)';
bc = -1 + 0.005:0.01:1;
HPC_b = mean(HC(bc > 0 & bc < 0.1, :), 'omitnan')';
V1_b.V1_pre = mean(V1(bc > -0.2 & bc < 0, :), 'omitnan')';
V1_b.V1_post = mean(V1(bc > 0 & bc < 0.2, :), 'omitnan')';
counts_dir = fullfile(rev_dir, 'pairwise_V1_HC_coupling', counts_name);
for wn = {'V1_pre', 'V1_post'}
    pop.(wn{1}) = NaN(numel(S), 1);
    for s = 1:numel(S)
        C = load(fullfile(counts_dir, sprintf('session_%02d.mat', S(s).session_count)), 'C');
        rows = C.C.ripple_info_rows;
        pop.(wn{1})(s) = corr(V1_b.(wn{1})(rows), HPC_b(rows), 'rows', 'complete');
    end
end
end
