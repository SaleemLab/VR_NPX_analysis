%% Per-animal reactivation coherence AUC (no bootstrap), POST sleep
% Same AUC as compare_session_AUC_PRE_POST, but events are pooled per animal
% (all its POST sessions) and condition percentiles are computed within animal.
% Non-duplicated ripple events (KDE_reactivation_ripples_PSTH_POST_unique), each
% for a pre (-0.2-0 s) and post (0-0.2 s) V1 window:
%   ripple_power : HPC->V1, within-animal ripple power Q1 (Low) vs Q4 (High)
%   spindle      : V1->HPC, hemisphere-matched within-animal spindle amplitude
%                  percentile Q1 (Low) vs Q4 (High)
%   SO_phase     : V1->HPC, hemisphere-matched SO peak vs trough
% Two-sided paired Wilcoxon signed-rank across animals (n = 6; minimum p = 0.031).
% Output: nature_neuroscience_revision_analysis\POST_animal
function compare_animal_AUC_POST()
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    out_dir = fullfile(rev_dir, 'POST_animal');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
    HC = (KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples + KDE_reactivation_ripples_PSTH.nan_mask)';
    V1 = (KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples + KDE_reactivation_ripples_PSTH.nan_mask)';
    clear KDE_reactivation_ripples_PSTH
    animal_id = ripple_info.subject_id(:);
    animals = unique(animal_id);
    nAnimal = numel(animals);

    % Within-animal percentiles
    ripple_q = within_group_quartile(ripple_info.ripple_power(:), animal_id);
    spindle_pct = NaN(size(ripple_info.spindle_amplitude));
    for a = animals'
        idx = animal_id == a;
        amp = ripple_info.spindle_amplitude(idx, :);
        spindle_pct(idx, :) = (tiedrank(amp) - 0.5) ./ size(amp, 1) * 100;
    end
    info_animal = ripple_info;
    info_animal.spindle_percentile = spindle_pct;   % hemisphere_conditions uses this field

    % Windows
    timebin = 0.01;
    bin_edges = -1:timebin:1;
    bin_centers = bin_edges(1:end-1) + timebin/2;
    bins_HPC = bin_centers > 0 & bin_centers < 0.1;
    v1_windows = struct( ...
        'name',  {'pre', 'post'}, ...
        'label', {'V1 -0.2-0 s (pre)', 'V1 0-0.2 s (post)'}, ...
        'bins',  {bin_centers > -0.2 & bin_centers < 0, bin_centers > 0 & bin_centers < 0.2});

    % Conditions: mask (events included), c1/c2 (Track 1 / Track 2 conditions)
    nEv = numel(animal_id);
    all_true = true(nEv, 1);
    Hs = hemisphere_conditions('spindle', info_animal);
    Hso = hemisphere_conditions('SO', ripple_info);
    sets = struct( ...
        'name',      {'ripple_power', 'spindle', 'SO_phase'}, ...
        'title',     {'Ripple power, within animal (HPC->V1)', 'Spindle percentile, within animal (V1->HPC)', 'SO phase (V1->HPC)'}, ...
        'direction', {'HPC->V1', 'V1->HPC', 'V1->HPC'}, ...
        'C', {struct('label', {'Low (Q1)', 'High (Q4)'}, 'mask', {ripple_q == 1, ripple_q == 4}, 'c1', all_true, 'c2', all_true), ...
              struct('label', {'Low (Q1)', 'High (Q4)'}, 'mask', all_true, 'c1', {Hs(1).c1, Hs(4).c1}, 'c2', {Hs(1).c2, Hs(4).c2}), ...
              struct('label', {Hso.label}, 'mask', all_true, 'c1', {Hso.c1}, 'c2', {Hso.c2})});

    mean_HPC = mean(HC(bins_HPC, :), 'omitnan');
    T = table(); S = table();
    for w = 1:numel(v1_windows)
        mean_V1 = mean(V1(v1_windows(w).bins, :), 'omitnan');
        for iSet = 1:numel(sets)
            C = sets(iSet).C;
            auc = NaN(nAnimal, 2); n_ev = NaN(nAnimal, 2);
            for iCond = 1:2
                for k = 1:nAnimal
                    in_a = animal_id == animals(k) & C(iCond).mask;
                    c1 = C(iCond).c1(in_a); c2 = C(iCond).c2(in_a);
                    if strcmp(sets(iSet).direction, 'HPC->V1')
                        r = reactivation_coherence_auc(mean_HPC(in_a), mean_V1(in_a), 0, c1, c2);
                    else
                        r = reactivation_coherence_auc(mean_V1(in_a), mean_HPC(in_a), 0, c1, c2);
                    end
                    auc(k, iCond) = r.AUC_mean; n_ev(k, iCond) = sum(in_a);
                end
                T = [T; table(animals, repmat(string(v1_windows(w).name), nAnimal, 1), repmat(string(sets(iSet).name), nAnimal, 1), ...
                    repmat(string(C(iCond).label), nAnimal, 1), repmat(string(sets(iSet).direction), nAnimal, 1), auc(:, iCond), n_ev(:, iCond), ...
                    'VariableNames', {'animal_id', 'V1_window', 'analysis', 'condition', 'direction', 'AUC', 'n_events'})]; %#ok<AGROW>
            end
            ok = all(isfinite(auc), 2);
            p = signrank(auc(ok, 1), auc(ok, 2));
            S = [S; table(string(sets(iSet).name), string(v1_windows(w).name), string(sprintf('%s vs %s', C.label)), ...
                sum(ok), mean(auc(ok, 1)), mean(auc(ok, 2)), sum(auc(ok, 2) > auc(ok, 1)), p, ...
                'VariableNames', {'analysis', 'V1_window', 'comparison', 'n_animals', 'mean_1', 'mean_2', 'n_animals_2_gt_1', 'p_signrank'})]; %#ok<AGROW>
        end
    end
    disp(S);
    save(fullfile(out_dir, 'animal_AUC_POST.mat'), 'T', 'S');
    writetable(T, fullfile(out_dir, 'animal_AUC_POST.csv'));
    writetable(S, fullfile(out_dir, 'animal_AUC_POST_stats.csv'));

    % --- Plot: rows = analysis, columns = V1 window; paired lines per animal ---
    col = [231, 41, 138] / 256;
    fig = figure('Name', 'Per-animal reactivation coherence AUC (POST)', 'Position', [100 100 700 1000]);
    tiledlayout(numel(sets), numel(v1_windows), 'TileSpacing', 'compact');
    for iSet = 1:numel(sets)
        for w = 1:numel(v1_windows)
            ax = nexttile; hold on;
            Tw = T(T.analysis == sets(iSet).name & T.V1_window == v1_windows(w).name, :);
            labels = {sets(iSet).C.label};
            Y = [Tw.AUC(Tw.condition == labels{1}) Tw.AUC(Tw.condition == labels{2})];
            plot([1 2], Y', '-', 'Color', [0 0 0 0.3], 'LineWidth', 1);
            for iCond = 1:2
                scatter(iCond * ones(nAnimal, 1), Y(:, iCond), 35, col, 'filled', 'MarkerFaceAlpha', 0.7);
                errorbar(iCond + (2*iCond - 3) * 0.25, mean(Y(:, iCond), 'omitnan'), std(Y(:, iCond), 'omitnan') / sqrt(sum(isfinite(Y(:, iCond)))), ...
                    'o', 'Color', 'k', 'MarkerFaceColor', col, 'LineWidth', 1.5, 'CapSize', 6);
            end
            yline(0, 'k--');
            set(ax, 'XTick', [1 2], 'XTickLabel', labels);
            xlim([0.5 2.5]); ylabel('Reactivation coherence AUC');
            st = S(S.analysis == sets(iSet).name & S.V1_window == v1_windows(w).name, :);
            title({sets(iSet).title, sprintf('%s: p = %.3g (%d/%d animals up, n = %d)', v1_windows(w).label, ...
                st.p_signrank, st.n_animals_2_gt_1, st.n_animals, st.n_animals)}, 'FontSize', 9, 'FontWeight', 'normal');
            set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 11);
        end
    end
    sgtitle('POST: per-animal reactivation coherence AUC (paired signed-rank across animals)', 'FontSize', 12);
    savefig(fig, fullfile(out_dir, 'animal_AUC_POST.fig'));
    exportgraphics(fig, fullfile(out_dir, 'animal_AUC_POST.png'), 'Resolution', 300);
    fprintf('Saved to %s\n', out_dir);
end

function q = within_group_quartile(metric, group_ids)
    % Quartile (1-4) of the metric within each group
    q = NaN(size(metric));
    for g = unique(group_ids)'
        idx = find(group_ids == g);
        pcts = prctile(metric(idx), 0:25:100);
        for iBin = 1:4
            if iBin == 1
                in_q = metric(idx) >= pcts(iBin) & metric(idx) <= pcts(iBin+1);
            else
                in_q = metric(idx) > pcts(iBin) & metric(idx) <= pcts(iBin+1);
            end
            q(idx(in_q)) = iBin;
        end
    end
end
