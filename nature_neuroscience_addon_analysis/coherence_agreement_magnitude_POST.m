function coherence_agreement_magnitude_POST()
% Splits the POST reactivation coherence of each oscillation condition into an agreement
% part and a magnitude part, to ask whether a condition makes V1 and HC agree more on the
% track or only express their (agreeing or not) track bias more strongly.
% Log odds: original PLS-KDE (KDE_reactivation_ripples_PSTH_POST_unique.mat), HPC = mean z
% log odds 0-0.1 s, V1 = mean z log odds pre (-0.2-0 s) or post (0-0.2 s).
% Conditions (A -> B), as compare_session_geomean_coherence / boot_AUC_per_unit_POST_original:
%   ripple_power     : within-session ripple power Q1 -> Q4
%   spindle          : dominant-hemisphere spindle percentile <= 25 -> > 75
%   SO_phase         : dominant-hemisphere SO peak -> trough
%   SO_trough_uni_bi : dominant trough with the other hemisphere in trough (bilateral) ->
%                      in peak (unilateral)
%   dominant hemisphere from the V1 sign in that window (hemisphere_conditions: V1 > 0 ->
%   hemisphere 2, V1 < 0 -> hemisphere 1)
% Per session and subset:
%   geomean   : mean of sign(HC * V1) * sqrt(|HC * V1|) (the paper's event coherence)
%   agreement : fraction of events with sign(V1) = sign(HC)
%   spearman  : Spearman r between V1 and HC log odds across events
%   abs_V1, abs_HC : mean |log odds| (magnitude of the expressed track bias)
% Paired two-sided signrank of B - A across sessions; sessions and animals with B > A.
% Output: nature_neuroscience_revision_analysis\coherence_agreement_vs_magnitude\
%   session_metrics.csv, group_stats.csv, agreement_vs_magnitude_<window>.fig/.pdf

rev_dir = addon_setup();
out_dir = fullfile(rev_dir, 'coherence_agreement_vs_magnitude');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
L = load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
HC = (L.KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples + L.KDE_reactivation_ripples_PSTH.nan_mask)';
V1 = (L.KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples + L.KDE_reactivation_ripples_PSTH.nan_mask)';
clear L
bc = -1 + 0.005:0.01:1;
HPC_b = mean(HC(bc > 0 & bc < 0.1, :), 'omitnan')';
wins = struct('name', {'V1_pre', 'V1_post'}, 'label', {'V1 -0.2-0 s', 'V1 0-0.2 s'}, ...
    'bins', {bc > -0.2 & bc < 0, bc > 0 & bc < 0.2});
sess = ripple_info.session_id(:);
animal = ripple_info.subject_id(:);
sessions = unique(sess)';
q = NaN(size(sess));
for s = sessions
    q(sess == s) = within_session_quartile(ripple_info.ripple_power(sess == s));
end
Hs = hemisphere_conditions('spindle', ripple_info);
Hso = hemisphere_conditions('SO', ripple_info);
Hub = hemisphere_conditions('SO_trough_bilateral', ripple_info);   % 1 unilateral, 2 bilateral
conds = {'ripple_power', 'spindle', 'SO_phase', 'SO_trough_uni_bi'};
labels = {{'Low (Q1)', 'High (Q4)'}, {'Low (Q1)', 'High (Q4)'}, {'Peak', 'Trough'}, ...
    {'Bilateral trough', 'Unilateral trough'}};
metrics = {'geomean', 'agreement', 'spearman', 'abs_V1', 'abs_HC'};
metric_labels = {'Geomean coherence', 'Sign agreement', 'Spearman r(V1, HC)', 'Mean |V1 log odds|', 'Mean |HC log odds|'};

T = table();
for w = 1:numel(wins)
    V1_b = mean(V1(wins(w).bins, :), 'omitnan')';
    matched = @(H) (V1_b > 0 & H.c1(:)) | (V1_b < 0 & H.c2(:));
    masks = {{q == 1, q == 4}, {matched(Hs(1)), matched(Hs(4))}, {matched(Hso(1)), matched(Hso(2))}, ...
        {matched(Hub(2)), matched(Hub(1))}};
    ok_ev = isfinite(V1_b) & isfinite(HPC_b);
    for c = 1:numel(conds)
        for s = sessions
            row = table(s, animal(find(sess == s, 1)), string(wins(w).name), string(conds{c}), ...
                'VariableNames', {'session_id', 'animal_id', 'V1_window', 'condition'});
            for k = 1:2
                e = masks{c}{k} & sess == s & ok_ev;
                m = event_metrics(V1_b(e), HPC_b(e));
                ab = 'AB';
                row.(['n_' ab(k)]) = sum(e);
                for i = 1:numel(metrics)
                    row.([metrics{i} '_' ab(k)]) = m.(metrics{i});
                end
            end
            T = [T; row]; %#ok<AGROW>
        end
    end
end
writetable(T, fullfile(out_dir, 'session_metrics.csv'));

%% Group statistics
animals = unique(T.animal_id);
G = table();
for w = 1:numel(wins)
    for c = 1:numel(conds)
        S = T(T.V1_window == wins(w).name & T.condition == conds{c}, :);
        for i = 1:numel(metrics)
            a = S.([metrics{i} '_A']); b = S.([metrics{i} '_B']); d = b - a; ok = isfinite(d);
            dm = arrayfun(@(x) mean(d(ok & S.animal_id == x)), animals);
            G = [G; table(string(wins(w).name), string(conds{c}), string(labels{c}{1}), string(labels{c}{2}), ...
                string(metrics{i}), sum(ok), mean(a(ok)), mean(b(ok)), mean(d(ok)), signrank(d(ok)), ...
                sum(d(ok) > 0), sum(dm > 0), numel(animals), ...
                'VariableNames', {'V1_window', 'condition', 'label_A', 'label_B', 'metric', 'n_sessions', ...
                'mean_A', 'mean_B', 'mean_delta', 'p_delta', 'n_sessions_B_gt_A', 'n_animals_B_gt_A', 'n_animals'})]; %#ok<AGROW>
        end
    end
end
writetable(G, fullfile(out_dir, 'group_stats.csv'));
disp(G(:, {'V1_window', 'condition', 'metric', 'mean_A', 'mean_B', 'p_delta', 'n_sessions_B_gt_A', 'n_animals_B_gt_A'}));

%% Figures: metrics (rows) x conditions (columns), one figure per V1 window
cols = lines(numel(animals));
for w = 1:numel(wins)
    f = figure('Name', ['agreement_vs_magnitude_' wins(w).name], 'Position', [40 40 1200 1100], 'Color', 'w');
    for i = 1:numel(metrics)
        for c = 1:numel(conds)
            subplot(numel(metrics), numel(conds), (i - 1) * numel(conds) + c); hold on
            S = T(T.V1_window == wins(w).name & T.condition == conds{c}, :);
            a = S.([metrics{i} '_A']); b = S.([metrics{i} '_B']); ok = isfinite(a) & isfinite(b);
            plot([1 2], [a(ok) b(ok)]', '-', 'Color', [0.8 0.8 0.8]);
            for k = 1:numel(animals)
                idx = ok & S.animal_id == animals(k);
                scatter([ones(sum(idx), 1); 2 * ones(sum(idx), 1)], [a(idx); b(idx)], 14, cols(k, :), 'filled', ...
                    'MarkerFaceAlpha', 0.7);
            end
            errorbar([0.75 2.25], [mean(a(ok)) mean(b(ok))], [std(a(ok)) std(b(ok))] / sqrt(sum(ok)), 'ko-', ...
                'MarkerFaceColor', 'k', 'LineWidth', 1.3, 'CapSize', 0);
            if ismember(metrics{i}, {'geomean', 'spearman'}), yline(0, '--'); end
            if strcmp(metrics{i}, 'agreement'), yline(0.5, '--'); end
            g = G(G.V1_window == wins(w).name & G.condition == conds{c} & G.metric == metrics{i}, :);
            t = sprintf('p=%s, %d/%d sess, %d/%d anim', p_str(g.p_delta), g.n_sessions_B_gt_A, g.n_sessions, ...
                g.n_animals_B_gt_A, g.n_animals);
            if i == 1, t = {strrep(conds{c}, '_', ' '), t}; end
            title(t, 'FontWeight', 'normal', 'FontSize', 8);
            set(gca, 'XTick', [1 2], 'XTickLabel', labels{c}, 'TickDir', 'out', 'Box', 'off', 'FontSize', 8);
            xlim([0.5 2.5]);
            if c == 1, ylabel(metric_labels{i}); end
        end
    end
    sgtitle({sprintf('Reactivation coherence: agreement vs magnitude (%s, HC 0-0.1 s)', wins(w).label), ...
        'original PLS-KDE log odds; spindle / SO from the dominant hemisphere; signrank of B - A across sessions'}, ...
        'FontSize', 11);
    savefig(f, fullfile(out_dir, [get(f, 'Name') '.fig']));
    exportgraphics(f, fullfile(out_dir, [get(f, 'Name') '.pdf']), 'ContentType', 'vector');
end
fprintf('Saved %s\n', out_dir);
end

function m = event_metrics(v, h)
if numel(v) < 5
    m = struct('geomean', NaN, 'agreement', NaN, 'spearman', NaN, 'abs_V1', NaN, 'abs_HC', NaN);
    return
end
m.geomean = mean(sign(h .* v) .* sqrt(abs(h .* v)));
nz = v ~= 0 & h ~= 0;
m.agreement = mean(sign(v(nz)) == sign(h(nz)));
m.spearman = corr(v, h, 'type', 'Spearman');
m.abs_V1 = mean(abs(v));
m.abs_HC = mean(abs(h));
end

function q = within_session_quartile(x)
% Quartile (1-4) within the session (as reactivation_AUC_coherence_PRE)
q = NaN(size(x));
pcts = prctile(x, 0:25:100);
for iBin = 1:4
    if iBin == 1
        in_q = x >= pcts(iBin) & x <= pcts(iBin + 1);
    else
        in_q = x > pcts(iBin) & x <= pcts(iBin + 1);
    end
    q(in_q) = iBin;
end
end

function s = p_str(p)
if p < 0.001
    s = sprintf('%.1e', p);
else
    s = sprintf('%.3f', p);
end
end
