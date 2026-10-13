%% Mean |log-odds bias| of HPC and V1 for Low vs High conditions, PRE and POST sleep
% Non-duplicated ripple events (make_unique_ripple_data), capped z-scored log odds.
% HPC |bias|: mean over 0-0.1 s.  V1 |bias|: mean over pre (-0.2-0 s) or post
% (0-0.2 s) window.
% Conditions:
%   ripple power : within-session Q1 (Low) vs Q4 (High)
%   spindle      : within-session spindle percentile Q1 (Low) vs Q4 (High)
%   SO phase     : peak [-pi/2, pi/2] vs trough
% Spindle and SO use the hemisphere matched to the event's V1 bias sign in that
% window (V1 > 0 -> hemisphere 2, V1 < 0 -> hemisphere 1), as in the
% reactivation coherence AUC.
% Plots per-session means with paired Low vs High (Peak vs Trough) lines and a
% paired signed-rank test across sessions within each sleep (PRE 12, POST 22).
% Output: nature_neuroscience_revision_analysis\abs_bias
function plot_abs_bias_conditions()
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    out_dir = fullfile(rev_dir, 'abs_bias');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % --- Load non-duplicated data ---
    load(fullfile(rev_dir, 'ripple_info_PRE_unique.mat'), 'ripple_info_PRE');
    load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_PRE_unique.mat'), 'KDE_reactivation_ripples_PSTH_PRE');
    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');

    D(1).name = 'PRE';
    D(1).HC = (KDE_reactivation_ripples_PSTH_PRE.HPC_z_logodds_ripples + KDE_reactivation_ripples_PSTH_PRE.nan_mask)';
    D(1).V1 = (KDE_reactivation_ripples_PSTH_PRE.V1_z_logodds_ripples + KDE_reactivation_ripples_PSTH_PRE.nan_mask)';
    D(1).info = ripple_info_PRE;
    D(2).name = 'POST';
    D(2).HC = (KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples + KDE_reactivation_ripples_PSTH.nan_mask)';
    D(2).V1 = (KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples + KDE_reactivation_ripples_PSTH.nan_mask)';
    D(2).info = ripple_info;
    clear KDE_reactivation_ripples_PSTH KDE_reactivation_ripples_PSTH_PRE

    % --- Windows ---
    timebin = 0.01;
    bin_edges = -1:timebin:1;
    bin_centers = bin_edges(1:end-1) + timebin/2;
    bins_HPC = bin_centers > 0 & bin_centers < 0.1;
    v1_windows = struct( ...
        'name',  {'pre', 'post'}, ...
        'label', {'V1 -0.2-0 s (pre)', 'V1 0-0.2 s (post)'}, ...
        'bins',  {bin_centers > -0.2 & bin_centers < 0, bin_centers > 0 & bin_centers < 0.2});

    analyses = {'ripple_power', 'spindle', 'SO_phase'};
    analysis_titles = {'Ripple power', 'Spindle percentile', 'SO phase'};
    regions = {'HPC', 'V1'};

    % --- Mean |bias| per condition ---
    T = table();       % event-level summary
    Tsess = table();   % per-session means
    for iS = 1:2
        info = D(iS).info;
        sess = info.session_id(:);
        HPC_b = mean(D(iS).HC(bins_HPC, :), 'omitnan')';
        for w = 1:numel(v1_windows)
            V1_b = mean(D(iS).V1(v1_windows(w).bins, :), 'omitnan')';
            for iA = 1:numel(analyses)
                [masks, labels] = condition_masks(analyses{iA}, info, V1_b);
                for iCond = 1:2
                    for iR = 1:2
                        if iR == 1, vals = abs(HPC_b); else, vals = abs(V1_b); end
                        v = vals(masks{iCond});
                        T = [T; table(string(D(iS).name), string(v1_windows(w).name), string(analyses{iA}), string(labels{iCond}), ...
                            string(regions{iR}), mean(v, 'omitnan'), std(v, 'omitnan') / sqrt(sum(isfinite(v))), sum(isfinite(v)), ...
                            'VariableNames', {'sleep', 'V1_window', 'analysis', 'condition', 'region', 'mean_abs_bias', 'sem', 'n_events'})]; %#ok<AGROW>
                        for s = unique(sess)'
                            vs = vals(masks{iCond} & sess == s);
                            Tsess = [Tsess; table(s, string(D(iS).name), string(v1_windows(w).name), string(analyses{iA}), string(labels{iCond}), ...
                                string(regions{iR}), mean(vs, 'omitnan'), sum(isfinite(vs)), ...
                                'VariableNames', {'session_id', 'sleep', 'V1_window', 'analysis', 'condition', 'region', 'mean_abs_bias', 'n_events'})]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end

    % --- Paired signed-rank across sessions: condition A vs B within each sleep ---
    S = table();
    for w = 1:numel(v1_windows)
        for iA = 1:numel(analyses)
            for iR = 1:2
                for iS = 1:2
                    [Y, labels] = session_matrix(Tsess, D(iS).name, v1_windows(w).name, analyses{iA}, regions{iR});
                    ok = all(isfinite(Y), 2);
                    p = NaN;
                    if sum(ok) >= 3, p = signrank(Y(ok, 1), Y(ok, 2)); end
                    S = [S; table(string(D(iS).name), string(v1_windows(w).name), string(analyses{iA}), string(regions{iR}), ...
                        string(sprintf('%s vs %s', labels{1}, labels{2})), sum(ok), mean(Y(ok, 1)), mean(Y(ok, 2)), p, ...
                        'VariableNames', {'sleep', 'V1_window', 'analysis', 'region', 'comparison', 'n_sessions', 'mean_1', 'mean_2', 'p_signrank'})]; %#ok<AGROW>
                end
            end
        end
    end
    disp(S);

    save(fullfile(out_dir, 'abs_bias_conditions.mat'), 'T', 'Tsess', 'S');
    writetable(T, fullfile(out_dir, 'abs_bias_conditions.csv'));
    writetable(Tsess, fullfile(out_dir, 'abs_bias_conditions_per_session.csv'));
    writetable(S, fullfile(out_dir, 'abs_bias_conditions_stats.csv'));

    % --- Plot: one figure per V1 window; rows = region, columns = analysis.
    % Per-session means, paired lines between conditions within each sleep ---
    sleep_colors = {[0.5 0.5 0.5], [231, 41, 138] / 256};
    for w = 1:numel(v1_windows)
        wn = v1_windows(w).name;
        fig = figure('Name', ['Mean |bias| ' v1_windows(w).label], 'Position', [100 100 1300 800]);
        tiledlayout(2, numel(analyses), 'TileSpacing', 'compact');
        for iR = 1:2
            for iA = 1:numel(analyses)
                ax = nexttile; hold on;
                p_str = cell(1, 2);
                for iS = 1:2
                    [Y, labels] = session_matrix(Tsess, D(iS).name, wn, analyses{iA}, regions{iR});
                    xp = (iS - 1) * 2.5 + [1 2];
                    plot(xp, Y', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
                    for iCond = 1:2
                        scatter(xp(iCond) * ones(size(Y, 1), 1), Y(:, iCond), 25, sleep_colors{iS}, 'filled', 'MarkerFaceAlpha', 0.7);
                        m = mean(Y(:, iCond), 'omitnan'); se = std(Y(:, iCond), 'omitnan') / sqrt(sum(isfinite(Y(:, iCond))));
                        errorbar(xp(iCond) + (2*iCond - 3) * 0.3, m, se, 'o', 'Color', 'k', 'MarkerFaceColor', sleep_colors{iS}, 'LineWidth', 1.5, 'CapSize', 6);
                    end
                    st = S(S.sleep == D(iS).name & S.V1_window == wn & S.analysis == analyses{iA} & S.region == regions{iR}, :);
                    p_str{iS} = sprintf('%s p = %.3g (n = %d)', D(iS).name, st.p_signrank, st.n_sessions);
                end
                set(ax, 'XTick', [1 2 3.5 4.5], 'XTickLabel', ...
                    {['PRE ' labels{1}], ['PRE ' labels{2}], ['POST ' labels{1}], ['POST ' labels{2}]}, ...
                    'XTickLabelRotation', 30);
                xlim([0.4 5.1]);
                if iR == 1
                    ylabel('Mean |HPC bias| (0-0.1 s)');
                else
                    ylabel(sprintf('Mean |V1 bias| (%s)', v1_windows(w).label));
                end
                title({sprintf('%s: %s', regions{iR}, analysis_titles{iA}), p_str{1}, p_str{2}}, 'FontSize', 10, 'FontWeight', 'normal');
                set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 11);
            end
        end
        sgtitle(sprintf('Mean |log-odds bias| per session (paired signed-rank), %s', v1_windows(w).label), 'FontSize', 13);
        savefig(fig, fullfile(out_dir, ['abs_bias_conditions_' wn '_window.fig']));
        exportgraphics(fig, fullfile(out_dir, ['abs_bias_conditions_' wn '_window.png']), 'Resolution', 300);
    end
    fprintf('Saved to %s\n', out_dir);
end

function [Y, labels] = session_matrix(Tsess, sleep, window, analysis, region)
    % Sessions x 2 conditions matrix of per-session mean |bias|
    Ts = Tsess(Tsess.sleep == sleep & Tsess.V1_window == window & Tsess.analysis == analysis & Tsess.region == region, :);
    labels = cellstr(unique(Ts.condition, 'stable'))';
    sessions = unique(Ts.session_id);
    Y = NaN(numel(sessions), 2);
    for iCond = 1:2
        Tc = Ts(Ts.condition == labels{iCond}, :);
        [~, loc] = ismember(Tc.session_id, sessions);
        Y(loc, iCond) = Tc.mean_abs_bias;
    end
end

function [masks, labels] = condition_masks(analysis, info, V1_b)
    % Logical event masks for the two conditions of an analysis
    switch analysis
        case 'ripple_power'
            q = within_session_quartile(info.ripple_power(:), info.session_id(:));
            masks = {q == 1, q == 4};
            labels = {'Low', 'High'};
        case {'spindle', 'SO_phase'}
            % Feature from the hemisphere matched to the V1 bias sign
            if strcmp(analysis, 'spindle'), F = info.spindle_percentile; else, F = info.SO_phase; end
            f = NaN(size(V1_b));
            f(V1_b > 0) = F(V1_b > 0, 2);
            f(V1_b < 0) = F(V1_b < 0, 1);
            if strcmp(analysis, 'spindle')
                masks = {f <= 25, f > 75};
                labels = {'Low', 'High'};
            else
                masks = {f >= -pi/2 & f <= pi/2, (f >= -pi & f < -pi/2) | (f > pi/2 & f <= pi)};
                labels = {'Peak', 'Trough'};
            end
    end
end

function q = within_session_quartile(metric, session_ids)
    % Quartile (1-4) of the metric within each session
    q = NaN(size(metric));
    for s = unique(session_ids)'
        idx = find(session_ids == s);
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
