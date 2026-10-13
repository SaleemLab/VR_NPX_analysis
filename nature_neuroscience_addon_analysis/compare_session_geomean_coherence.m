%% Per-session geomean reactivation coherence, Low vs High and PRE vs POST
% Event coherence = sign(HPC*V1) * sqrt(|HPC*V1|) (as in predicting_coherence_with_GAMs),
% with HPC = mean z log odds over 0-0.1 s and V1 = mean z log odds over the pre
% (-0.2-0 s) or post (0-0.2 s) window. Non-duplicated ripple events
% (make_unique_ripple_data). Session means per condition:
%   all_ripples  : all events
%   ripple_power : within-session ripple power Q1 (Low) vs Q4 (High)
%   spindle      : spindle percentile <= 25 (Low) vs > 75 (High)
%   SO_phase     : SO peak [-pi/2, pi/2] vs trough
% Spindle and SO use the hemisphere matched to the event's V1 bias sign in that
% window (V1 > 0 -> hemisphere 2, V1 < 0 -> hemisphere 1).
% Two-sided paired Wilcoxon signed-rank tests across sessions:
%   condition A vs B within each sleep (all sessions of that sleep) and each
%   condition PRE vs POST (sessions with both).
% Output: nature_neuroscience_revision_analysis\PRE_vs_POST\session_geomean_coherence_*
function compare_session_geomean_coherence()
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    out_dir = fullfile(rev_dir, 'PRE_vs_POST');
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

    analyses = {'all_ripples', 'ripple_power', 'spindle', 'SO_phase'};
    titles = {'All Ripples', 'Ripple Power, within session', 'Spindle Percentile (V1-matched hemisphere)', 'SO Phase (V1-matched hemisphere)'};

    % --- Per-session mean geomean coherence ---
    Tsess = table();
    for iS = 1:2
        info = D(iS).info;
        sess = info.session_id(:);
        HPC_b = mean(D(iS).HC(bins_HPC, :), 'omitnan')';
        for w = 1:numel(v1_windows)
            V1_b = mean(D(iS).V1(v1_windows(w).bins, :), 'omitnan')';
            coh = sign(HPC_b .* V1_b) .* sqrt(abs(HPC_b .* V1_b));
            for iA = 1:numel(analyses)
                [masks, labels] = condition_masks(analyses{iA}, info, V1_b);
                for iCond = 1:numel(masks)
                    for s = unique(sess)'
                        v = coh(masks{iCond} & sess == s);
                        Tsess = [Tsess; table(s, string(D(iS).name), string(v1_windows(w).name), string(analyses{iA}), ...
                            string(labels{iCond}), mean(v, 'omitnan'), sum(isfinite(v)), ...
                            'VariableNames', {'session_id', 'sleep', 'V1_window', 'analysis', 'condition', 'geomean_coherence', 'n_events'})]; %#ok<AGROW>
                    end
                end
            end
        end
    end

    % --- Paired signed-rank tests ---
    S = table();
    for w = 1:numel(v1_windows)
        wn = v1_windows(w).name;
        for iA = 1:numel(analyses)
            [Ypre, labels] = session_matrix(Tsess, 'PRE', wn, analyses{iA});
            Ypost = session_matrix(Tsess, 'POST', wn, analyses{iA});
            if numel(labels) == 2
                S = [S; paired_test(analyses{iA}, wn, sprintf('PRE: %s vs %s', labels{:}), Ypre(:, 1), Ypre(:, 2))]; %#ok<AGROW>
                S = [S; paired_test(analyses{iA}, wn, sprintf('POST: %s vs %s', labels{:}), Ypost(:, 1), Ypost(:, 2))]; %#ok<AGROW>
            end
            for iCond = 1:numel(labels)
                S = [S; paired_test(analyses{iA}, wn, sprintf('%s: PRE vs POST', labels{iCond}), Ypre(:, iCond), Ypost(:, iCond))]; %#ok<AGROW>
            end
        end
    end
    disp(S);

    save(fullfile(out_dir, 'session_geomean_coherence.mat'), 'Tsess', 'S');
    writetable(Tsess(Tsess.n_events > 0, :), fullfile(out_dir, 'session_geomean_coherence.csv'));
    writetable(S, fullfile(out_dir, 'session_geomean_coherence_stats.csv'));

    for iA = 1:numel(analyses)
        plot_analysis(Tsess, S(S.analysis == analyses{iA}, :), analyses{iA}, titles{iA}, v1_windows, ...
            fullfile(out_dir, ['session_geomean_coherence_' analyses{iA}]));
    end
    fprintf('Saved to %s\n', out_dir);
end

function [masks, labels] = condition_masks(analysis, info, V1_b)
    % Logical event masks per condition
    switch analysis
        case 'all_ripples'
            masks = {true(size(V1_b))};
            labels = {'All'};
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

function [Y, labels] = session_matrix(Tsess, sleep, window, analysis)
    % Sessions (1:22 index) x conditions matrix of session means; NaN where missing
    Ts = Tsess(Tsess.sleep == sleep & Tsess.V1_window == window & Tsess.analysis == analysis, :);
    labels = cellstr(unique(Ts.condition, 'stable'))';
    Y = NaN(max(Tsess.session_id), numel(labels));
    for iCond = 1:numel(labels)
        Tc = Ts(Ts.condition == labels{iCond} & Ts.n_events > 0, :);
        Y(Tc.session_id, iCond) = Tc.geomean_coherence;
    end
end

function row = paired_test(analysis, window, comparison, x, y)
    ok = isfinite(x) & isfinite(y);
    p = NaN;
    if sum(ok) >= 3
        p = signrank(x(ok), y(ok));
    end
    row = table(string(analysis), string(window), string(comparison), sum(ok), mean(x(ok)), mean(y(ok)), p, ...
        'VariableNames', {'analysis', 'V1_window', 'comparison', 'n_sessions', 'mean_1', 'mean_2', 'p_signrank'});
end

function plot_analysis(Tsess, S, analysis, title_str, v1_windows, save_path)
    % One tile per V1 window; per-session lines connect conditions within each sleep
    % (or PRE-POST for all ripples)
    sleep_colors = {[0.5 0.5 0.5], [231, 41, 138] / 256};
    sleeps = {'PRE', 'POST'};
    nW = numel(v1_windows);
    fig = figure('Name', ['Geomean coherence ' title_str], 'Position', [100 100 420*nW 520]);
    tiledlayout(1, nW, 'TileSpacing', 'compact');
    for w = 1:nW
        wn = v1_windows(w).name;
        ax = nexttile; hold on;
        Y = cell(1, 2);
        for iS = 1:2
            [Y{iS}, labels] = session_matrix(Tsess, sleeps{iS}, wn, analysis);
        end
        nCond = numel(labels);
        if nCond == 2
            xpos = {[1 2], [3.5 4.5]};
            xlabels = {['PRE ' labels{1}], ['PRE ' labels{2}], ['POST ' labels{1}], ['POST ' labels{2}]};
            for iS = 1:2
                plot(xpos{iS}, Y{iS}', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
            end
        else
            xpos = {1, 2};
            xlabels = {'PRE', 'POST'};
            plot([1 2], [Y{1} Y{2}]', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
        end
        for iS = 1:2
            for iCond = 1:nCond
                y = Y{iS}(:, iCond); y = y(isfinite(y));
                x = xpos{iS}(iCond);
                scatter(x * ones(size(y)), y, 25, sleep_colors{iS}, 'filled', 'MarkerFaceAlpha', 0.7);
                errorbar(x + 0.3 * (2 * (iCond == nCond) - 1) * (nCond == 2) + 0.3 * (nCond == 1), mean(y), std(y) / sqrt(numel(y)), ...
                    'o', 'Color', 'k', 'MarkerFaceColor', sleep_colors{iS}, 'LineWidth', 1.5, 'CapSize', 6);
            end
        end
        yline(0, 'k--');
        set(ax, 'XTick', [xpos{:}], 'XTickLabel', xlabels, 'XTickLabelRotation', 30);
        xlim([min([xpos{:}]) - 0.6, max([xpos{:}]) + 0.6]);
        ylabel('Geomean reactivation coherence (session mean)');
        Sw = S(S.V1_window == wn, :);
        p_lines = arrayfun(@(i) sprintf('%s: p = %.3g (n = %d)', Sw.comparison(i), Sw.p_signrank(i), Sw.n_sessions(i)), ...
            1:height(Sw), 'UniformOutput', false);
        title([{v1_windows(w).label}, p_lines], 'FontSize', 9, 'FontWeight', 'normal');
        set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 11);
    end
    linkaxes(findobj(fig, 'Type', 'axes'), 'y');
    sgtitle(['Geomean coherence: ' title_str ' (paired signed-rank)'], 'FontSize', 12);
    savefig(fig, [save_path '.fig']);
    exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
end
