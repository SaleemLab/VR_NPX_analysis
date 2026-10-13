%% Per-session reactivation coherence AUC (no bootstrap), PRE vs POST sleep
% Non-duplicated ripple events (make_unique_ripple_data), capped log odds.
% Analyses (each for a pre (-0.2-0 s) and post (0-0.2 s) V1 window):
%   all_ripples     HPC->V1, all events
%   ripple_power    HPC->V1, within-session ripple power Q1 (Low) vs Q4 (High)
%   spindle         V1->HPC, hemisphere-matched within-session spindle percentile Q1 vs Q4
%   SO_phase        V1->HPC, hemisphere-matched SO peak vs trough
% Paired Wilcoxon signed-rank tests across sessions:
%   condition A vs B within each sleep (all sessions of that sleep), and each
%   condition PRE vs POST (sessions with both PRE and POST).
function compare_session_AUC_PRE_POST()
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
    fprintf('Loading PRE and POST data...\n');
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
    for iS = 1:2
        D(iS).session_id = D(iS).info.session_id(:);
        assert(numel(D(iS).session_id) == size(D(iS).HC, 2), '%s session ids do not match PSTH rows', D(iS).name);
        fprintf('%s: %d events\n', D(iS).name, numel(D(iS).session_id));
    end

    % --- Windows ---
    timebin = 0.01;
    bin_edges = -1:timebin:1;
    bin_centers = bin_edges(1:end-1) + timebin/2;
    bins_HPC = bin_centers > 0 & bin_centers < 0.1;
    v1_windows = struct( ...
        'name',  {'pre', 'post'}, ...
        'label', {'V1 -0.2-0 s (pre)', 'V1 0-0.2 s (post)'}, ...
        'bins',  {bin_centers > -0.2 & bin_centers < 0, bin_centers > 0 & bin_centers < 0.2});

    % AUC is computed for every session of each sleep (NaN where a sleep has no
    % data). Within-sleep tests use all sessions of that sleep (PRE 12, POST 22);
    % PRE vs POST tests use only sessions with both.
    sessions = union(unique(D(1).session_id), unique(D(2).session_id));
    nSess = numel(sessions);
    fprintf('PRE sessions: %d, POST sessions: %d, both: %d\n', numel(unique(D(1).session_id)), ...
        numel(unique(D(2).session_id)), numel(intersect(unique(D(1).session_id), unique(D(2).session_id))));

    sets = struct( ...
        'name',      {'all_ripples', 'ripple_power', 'spindle', 'SO_phase'}, ...
        'direction', {'HPC->V1', 'HPC->V1', 'V1->HPC', 'V1->HPC'}, ...
        'title',     {'All Ripples (HPC->V1)', 'Ripple Power, within session (HPC->V1)', ...
                      'Spindle Percentile, within session (V1->HPC)', 'SO Phase (V1->HPC)'});

    % --- Per-session AUC ---
    T = table();
    res = struct('session_id', sessions(:), 'window_labels', {{v1_windows.label}});
    for iSet = 1:numel(sets)
        for iS = 1:2
            C = make_conditions(sets(iSet).name, D(iS).info);
            mean_HPC = mean(D(iS).HC(bins_HPC, :), 'omitnan');
            for w = 1:numel(v1_windows)
                mean_V1 = mean(D(iS).V1(v1_windows(w).bins, :), 'omitnan');
                for iCond = 1:numel(C)
                    auc = NaN(nSess, 1); n_ev = NaN(nSess, 1);
                    for k = 1:nSess
                        in_sess = D(iS).session_id == sessions(k) & C(iCond).mask;
                        c1 = C(iCond).c1(in_sess); c2 = C(iCond).c2(in_sess);
                        if strcmp(sets(iSet).direction, 'HPC->V1')
                            r = reactivation_coherence_auc(mean_HPC(in_sess), mean_V1(in_sess), 0, c1, c2);
                        else
                            r = reactivation_coherence_auc(mean_V1(in_sess), mean_HPC(in_sess), 0, c1, c2);
                        end
                        auc(k) = r.AUC_mean; n_ev(k) = sum(in_sess);
                    end
                    res.(sets(iSet).name).labels{iCond} = C(iCond).label;
                    res.(sets(iSet).name).(D(iS).name).(v1_windows(w).name)(:, iCond) = auc;
                    res.(sets(iSet).name).(D(iS).name).n_events(:, iCond) = n_ev;
                    T = [T; table(res.session_id, repmat(string(D(iS).name), nSess, 1), repmat(string(sets(iSet).name), nSess, 1), ...
                        repmat(string(C(iCond).label), nSess, 1), repmat(string(sets(iSet).direction), nSess, 1), ...
                        repmat(string(v1_windows(w).name), nSess, 1), auc, n_ev, ...
                        'VariableNames', {'session_id', 'sleep', 'analysis', 'condition', 'direction', 'V1_window', 'AUC', 'n_events'})]; %#ok<AGROW>
                end
            end
        end
    end

    % --- Paired signed-rank tests ---
    S = table();
    for iSet = 1:numel(sets)
        R = res.(sets(iSet).name);
        nCond = numel(R.labels);
        for w = 1:numel(v1_windows)
            wn = v1_windows(w).name;
            % Condition A vs B within each sleep
            if nCond == 2
                for iS = 1:2
                    S = [S; paired_test(sets(iSet).name, wn, sprintf('%s: %s vs %s', D(iS).name, R.labels{1}, R.labels{2}), ...
                        R.(D(iS).name).(wn)(:, 1), R.(D(iS).name).(wn)(:, 2))]; %#ok<AGROW>
                end
            end
            % Each condition PRE vs POST
            for iCond = 1:nCond
                S = [S; paired_test(sets(iSet).name, wn, sprintf('%s: PRE vs POST', R.labels{iCond}), ...
                    R.PRE.(wn)(:, iCond), R.POST.(wn)(:, iCond))]; %#ok<AGROW>
            end
        end
    end
    disp(S);
    res.stats = S;

    save(fullfile(out_dir, 'session_AUC_PRE_vs_POST.mat'), 'res');
    T = T(T.n_events > 0, :);   % drop sessions without data for that sleep
    writetable(T, fullfile(out_dir, 'session_AUC_PRE_vs_POST.csv'));
    writetable(S, fullfile(out_dir, 'session_AUC_PRE_vs_POST_stats.csv'));

    for iSet = 1:numel(sets)
        plot_set(res.(sets(iSet).name), S(S.analysis == sets(iSet).name, :), v1_windows, sets(iSet).title, ...
            fullfile(out_dir, ['session_AUC_PRE_vs_POST_' sets(iSet).name]));
        % Within-sleep condition comparison (e.g. Low vs High), one figure per sleep
        if numel(res.(sets(iSet).name).labels) == 2
            for iS = 1:2
                plot_within_sleep(res.(sets(iSet).name), D(iS).name, S(S.analysis == sets(iSet).name, :), v1_windows, ...
                    sets(iSet).title, fullfile(out_dir, ['session_AUC_' D(iS).name '_' sets(iSet).name]));
            end
        end
    end
    fprintf('Saved to %s\n', out_dir);
end

function C = make_conditions(name, info)
    % mask: events included; c1/c2: Track 1 / Track 2 conditions (V1->HPC only)
    nEv = numel(info.session_id);
    all_true = true(nEv, 1);
    switch name
        case 'all_ripples'
            C = struct('label', 'All Ripples', 'mask', all_true, 'c1', all_true, 'c2', all_true);
        case 'ripple_power'
            q = within_session_quartile(info.ripple_power(:), info.session_id(:));
            C = struct('label', {'Low (Q1)', 'High (Q4)'}, 'mask', {q == 1, q == 4}, 'c1', all_true, 'c2', all_true);
        case 'spindle'
            H = hemisphere_conditions('spindle', info);
            C = struct('label', {'Low (Q1)', 'High (Q4)'}, 'mask', all_true, 'c1', {H(1).c1, H(4).c1}, 'c2', {H(1).c2, H(4).c2});
        case 'SO_phase'
            H = hemisphere_conditions('SO', info);
            C = struct('label', {H.label}, 'mask', all_true, 'c1', {H.c1}, 'c2', {H.c2});
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

function row = paired_test(analysis, window, comparison, x, y)
    ok = isfinite(x) & isfinite(y);
    p = NaN;
    if sum(ok) >= 3
        p = signrank(x(ok), y(ok));
    end
    row = table(string(analysis), string(window), string(comparison), sum(ok), mean(x(ok)), mean(y(ok)), p, ...
        'VariableNames', {'analysis', 'V1_window', 'comparison', 'n', 'mean_1', 'mean_2', 'p_signrank'});
end

function plot_within_sleep(R, sleep, S, v1_windows, title_str, save_path)
    % Condition A vs B within one sleep; one tile per V1 window
    nW = numel(v1_windows);
    if strcmp(sleep, 'PRE'), col = [0.5 0.5 0.5]; else, col = [231, 41, 138] / 256; end
    comparison = sprintf('%s: %s vs %s', sleep, R.labels{1}, R.labels{2});
    fig = figure('Name', [sleep ' per-session ' title_str], 'Position', [640 100 300*nW 420]);
    tiledlayout(1, nW, 'TileSpacing', 'compact');
    for w = 1:nW
        wn = v1_windows(w).name;
        Y = R.(sleep).(wn);
        ax = nexttile; hold on;
        plot([1 2], Y', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
        for iCond = 1:2
            scatter(iCond * ones(size(Y,1),1), Y(:, iCond), 30, col, 'filled', 'MarkerFaceAlpha', 0.7);
            m = mean(Y(:, iCond), 'omitnan'); se = std(Y(:, iCond), 'omitnan') / sqrt(sum(isfinite(Y(:, iCond))));
            errorbar(iCond + (2*iCond - 3) * 0.25, m, se, 'o', 'Color', 'k', 'MarkerFaceColor', col, 'LineWidth', 1.5, 'CapSize', 6);
        end
        yline(0, 'k--');
        set(ax, 'XTick', [1 2], 'XTickLabel', R.labels);
        xlim([0.5 2.5]); ylabel('Reactivation coherence AUC');
        st = S(S.V1_window == wn & S.comparison == comparison, :);
        title({v1_windows(w).label, sprintf('signed-rank p = %.3g (n = %d)', st.p_signrank, st.n)}, 'FontSize', 10, 'FontWeight', 'normal');
        set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 11);
    end
    linkaxes(findobj(fig, 'Type', 'axes'), 'y');
    sgtitle([sleep ': ' title_str], 'FontSize', 12);
    savefig(fig, [save_path '.fig']);
    exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
end

function plot_set(R, S, v1_windows, title_str, save_path)
    % One tile per V1 window; per-session lines connect conditions within each sleep
    nCond = numel(R.labels); nW = numel(v1_windows);
    pre_color = [0.5 0.5 0.5];
    post_color = [231, 41, 138] / 256;
    if nCond == 2
        xpos = {[1 2], [3.5 4.5]};          % PRE: A, B | POST: A, B
        xlabels = {['PRE ' R.labels{1}], ['PRE ' R.labels{2}], ['POST ' R.labels{1}], ['POST ' R.labels{2}]};
    else
        xpos = {1, 2};
        xlabels = {'PRE', 'POST'};
    end
    fig = figure('Name', ['Per-session ' title_str], 'Position', [640 100 380*nW 460]);
    tiledlayout(1, nW, 'TileSpacing', 'compact');
    for w = 1:nW
        wn = v1_windows(w).name;
        ax = nexttile; hold on;
        sleeps = {'PRE', 'POST'}; cols = {pre_color, post_color};
        for iS = 1:2
            Y = R.(sleeps{iS}).(wn);
            xp = xpos{iS};
            if nCond == 2
                plot(xp, Y', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
            end
            for iCond = 1:nCond
                scatter(xp(iCond) * ones(size(Y,1),1), Y(:, iCond), 25, cols{iS}, 'filled', 'MarkerFaceAlpha', 0.7);
                m = mean(Y(:, iCond), 'omitnan'); se = std(Y(:, iCond), 'omitnan') / sqrt(sum(isfinite(Y(:, iCond))));
                errorbar(xp(iCond) + 0.25, m, se, 'o', 'Color', 'k', 'MarkerFaceColor', cols{iS}, 'LineWidth', 1.5, 'CapSize', 6);
            end
        end
        if nCond == 1
            Y = [R.PRE.(wn) R.POST.(wn)];
            plot([1 2], Y', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
        end
        yline(0, 'k--');
        set(ax, 'XTick', [xpos{:}], 'XTickLabel', xlabels, 'XTickLabelRotation', 30);
        xlim([min([xpos{:}]) - 0.5, max([xpos{:}]) + 0.6]);
        ylabel('Reactivation coherence AUC');
        Sw = S(S.V1_window == wn, :);
        p_lines = arrayfun(@(i) sprintf('%s: p = %.3g', Sw.comparison(i), Sw.p_signrank(i)), 1:height(Sw), 'UniformOutput', false);
        title([{v1_windows(w).label}, p_lines], 'FontSize', 9, 'FontWeight', 'normal');
        set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 11);
    end
    linkaxes(findobj(fig, 'Type', 'axes'), 'y');
    sgtitle(['Per-session ' title_str ' (paired signed-rank)'], 'FontSize', 12);
    savefig(fig, [save_path '.fig']);
    exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
end
