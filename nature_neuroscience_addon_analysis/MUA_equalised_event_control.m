%% MUA-equalised event control: Low vs High reactivation coherence (POST)
% Firing-rate control for the condition comparisons (R1.2, R4.2/4/6): within each
% session, Low and High events are subsampled so their joint distribution of
% HC MUA and dominant-hemisphere V1 MUA is equal, then coherence is recomputed.
% Decoding is unchanged (decoder trained on RUN; each event's bias is independent
% of which events are selected).
%
% Event coherence = sign(HPC*V1) * sqrt(|HPC*V1|) (as compare_session_geomean_coherence),
% HPC = mean z log odds 0-0.1 s, V1 = mean z log odds over pre (-0.2-0 s) or post
% (0-0.2 s) window. Dominant hemisphere from the V1 bias sign in that window
% (V1 > 0 -> hemisphere 2 / V1 R, V1 < 0 -> hemisphere 1 / V1 L).
% MUA (z-scored, 10 ms bins, -2 to 2 s): HC MUA 0-0.1 s, V1 MUA of the dominant
% hemisphere over the same window as V1 bias.
% Conditions (within session):
%   ripple_power : ripple power percentile <= 25 (Low) vs > 75 (High)
%   spindle      : dominant-hemisphere spindle percentile <= 25 (Low) vs > 75 (High),
%                  re-ranked within session over unique events
%   SO_phase     : dominant-hemisphere SO peak [-pi/2, pi/2] vs trough
% Matching: within session, HC MUA and V1 MUA are each split into nQ quantile bins
% (edges from both conditions pooled); in every joint bin both conditions keep
% min(nA, nB) randomly drawn events. Repeated nDraw times (seeded).
% Outputs per session: unmatched and matched (mean over draws) geomean coherence,
% paired signed-rank tests across sessions; pooled AUC (V1 selector, HPC measure;
% reactivation_coherence_auc) per draw.
% Output: nature_neuroscience_revision_analysis\MUA_equalised_event_control\geomean
function MUA_equalised_event_control()
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    out_dir = fullfile(rev_dir, 'MUA_equalised_event_control', 'geomean');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    nDraw = 200;
    nQ = 5;          % quantile bins per MUA dimension (nQ x nQ joint bins)
    min_events = 10; % per condition per session after matching

    % --- Load non-duplicated POST data (rows aligned to ripple_info) ---
    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
    HC = (KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples + KDE_reactivation_ripples_PSTH.nan_mask)';
    V1 = (KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples + KDE_reactivation_ripples_PSTH.nan_mask)';
    clear KDE_reactivation_ripples_PSTH
    load(fullfile(rev_dir, 'ripple_MUA_HPC.mat'), 'ripple_MUA_HPC');
    load(fullfile(rev_dir, 'ripple_MUA_V1_L.mat'), 'ripple_MUA_V1_L');
    load(fullfile(rev_dir, 'ripple_MUA_V1_R.mat'), 'ripple_MUA_V1_R');
    n_ev = numel(ripple_info.session_id);
    assert(size(ripple_MUA_HPC, 1) == n_ev && size(HC, 2) == n_ev, 'MUA / PSTH rows do not match ripple_info');

    % --- Windows ---
    timebin = 0.01;
    bin_edges = -1:timebin:1;
    bin_centers = bin_edges(1:end-1) + timebin/2;
    mua_centers = -2 + timebin/2:timebin:2 - timebin/2;
    bins_HPC = bin_centers > 0 & bin_centers < 0.1;
    v1_windows = struct( ...
        'name',  {'pre', 'post'}, ...
        'label', {'V1 -0.2-0 s (pre)', 'V1 0-0.2 s (post)'}, ...
        'range', {[-0.2 0], [0 0.2]});

    sess = ripple_info.session_id(:);
    sessions = unique(sess)';
    HPC_b = mean(HC(bins_HPC, :), 'omitnan')';
    HC_MUA = double(mean(ripple_MUA_HPC(:, mua_centers > 0 & mua_centers < 0.1), 2));
    ripple_pct = within_session_percentile(ripple_info.ripple_power(:), sess);
    spindle_pct = within_session_percentile(ripple_info.spindle_percentile, sess);

    analyses = {'ripple_power', 'spindle', 'SO_phase'};
    titles = {'Ripple power (within-session percentile)', 'Spindle-band power (dominant hemisphere)', 'SO phase (dominant hemisphere)'};

    T = table();          % per session x window x analysis x condition
    P = table();          % pooled AUC per draw
    M = struct();         % example matched events (draw 1) for MUA distribution plots
    for w = 1:numel(v1_windows)
        wn = v1_windows(w).name;
        in_w = mua_centers > v1_windows(w).range(1) & mua_centers < v1_windows(w).range(2);
        in_wb = bin_centers > v1_windows(w).range(1) & bin_centers < v1_windows(w).range(2);
        V1_b = mean(V1(in_wb, :), 'omitnan')';
        coh = sign(HPC_b .* V1_b) .* sqrt(abs(HPC_b .* V1_b));
        hemi = NaN(n_ev, 1); hemi(V1_b < 0) = 1; hemi(V1_b > 0) = 2;
        V1_L = double(mean(ripple_MUA_V1_L(:, in_w), 2));
        V1_R = double(mean(ripple_MUA_V1_R(:, in_w), 2));
        V1_MUA = NaN(n_ev, 1);
        V1_MUA(hemi == 1) = V1_L(hemi == 1);
        V1_MUA(hemi == 2) = V1_R(hemi == 2);
        valid = isfinite(coh) & isfinite(HC_MUA) & isfinite(V1_MUA) & ~isnan(hemi);

        for iA = 1:numel(analyses)
            [masks, labels] = condition_masks(analyses{iA}, ripple_info, ripple_pct, spindle_pct, hemi);
            A = masks{1} & valid; B = masks{2} & valid;
            pooled_idx = cell(nDraw, 2);
            for s = sessions
                ia = find(A & sess == s); ib = find(B & sess == s);
                [draws_a, draws_b] = match_events(ia, ib, HC_MUA, V1_MUA, nQ, nDraw, s);
                n_match = numel(draws_a{1});
                coh_m = [NaN NaN]; hc_m = [NaN NaN]; v1_m = [NaN NaN];
                if n_match >= min_events
                    ca = cellfun(@(i) mean(coh(i)), draws_a); cb = cellfun(@(i) mean(coh(i)), draws_b);
                    coh_m = [mean(ca) mean(cb)];
                    hc_m = [mean(cellfun(@(i) mean(HC_MUA(i)), draws_a)) mean(cellfun(@(i) mean(HC_MUA(i)), draws_b))];
                    v1_m = [mean(cellfun(@(i) mean(V1_MUA(i)), draws_a)) mean(cellfun(@(i) mean(V1_MUA(i)), draws_b))];
                    for d = 1:nDraw
                        pooled_idx{d, 1} = [pooled_idx{d, 1}; draws_a{d}];
                        pooled_idx{d, 2} = [pooled_idx{d, 2}; draws_b{d}];
                    end
                end
                idx_c = {ia, ib};
                for c = 1:2
                    i = idx_c{c};
                    T = [T; table(s, string(wn), string(analyses{iA}), string(labels{c}), ...
                        mean(coh(i)), coh_m(c), numel(i), n_match * (n_match >= min_events), ...
                        mean(HC_MUA(i)), hc_m(c), mean(V1_MUA(i)), v1_m(c), ...
                        'VariableNames', {'session_id', 'V1_window', 'analysis', 'condition', ...
                        'coh_unmatched', 'coh_matched', 'n_unmatched', 'n_matched', ...
                        'HC_MUA_unmatched', 'HC_MUA_matched', 'V1_MUA_unmatched', 'V1_MUA_matched'})]; %#ok<AGROW>
                end
            end

            % Pooled AUC (V1 selector -> HPC measure): unmatched and per matched draw
            auc_un = [auc_point(V1_b(A), HPC_b(A)) auc_point(V1_b(B), HPC_b(B))];
            auc_m = NaN(nDraw, 2);
            parfor d = 1:nDraw
                auc_m(d, :) = [auc_point(V1_b(pooled_idx{d, 1}), HPC_b(pooled_idx{d, 1})), ...
                    auc_point(V1_b(pooled_idx{d, 2}), HPC_b(pooled_idx{d, 2}))];
            end
            P = [P; table(string(wn), string(analyses{iA}), string(labels{1}), string(labels{2}), ...
                auc_un(1), auc_un(2), {auc_m(:, 1)}, {auc_m(:, 2)}, ...
                mean(auc_m(:, 2) - auc_m(:, 1)), prctile(auc_m(:, 2) - auc_m(:, 1), 2.5), ...
                prctile(auc_m(:, 2) - auc_m(:, 1), 97.5), mean(auc_m(:, 2) - auc_m(:, 1) <= 0), ...
                sum(A), sum(B), numel(pooled_idx{1, 1}), ...
                'VariableNames', {'V1_window', 'analysis', 'cond_A', 'cond_B', ...
                'AUC_A_unmatched', 'AUC_B_unmatched', 'AUC_A_matched', 'AUC_B_matched', ...
                'AUC_diff_matched_mean', 'AUC_diff_CI_lo', 'AUC_diff_CI_hi', 'p_draw_diff_le0', ...
                'n_A_unmatched', 'n_B_unmatched', 'n_per_cond_matched'})]; %#ok<AGROW>

            key = sprintf('%s_%s', analyses{iA}, wn);
            M.(key) = struct('labels', {labels}, 'A', find(A), 'B', find(B), ...
                'A_m', pooled_idx{1, 1}, 'B_m', pooled_idx{1, 2}, 'HC_MUA', HC_MUA, 'V1_MUA', V1_MUA);
            fprintf('%s %s: %d/%d events retained per condition (pooled)\n', wn, analyses{iA}, ...
                numel(pooled_idx{1, 1}), min(sum(A), sum(B)));
        end
    end

    % --- Session-level paired tests (B vs A; one-sided B > A) ---
    S = table();
    for w = 1:numel(v1_windows)
        wn = v1_windows(w).name;
        for iA = 1:numel(analyses)
            for meas = ["unmatched", "matched"]
                [Y, labels] = session_matrix(T, wn, analyses{iA}, "coh_" + meas);
                S = [S; paired_test(analyses{iA}, wn, meas, labels, Y(:, 1), Y(:, 2))]; %#ok<AGROW>
            end
        end
    end
    disp(S);
    Pout = removevars(P, {'AUC_A_matched', 'AUC_B_matched'});
    disp(Pout);

    save(fullfile(out_dir, 'MUA_equalised_event_control.mat'), 'T', 'S', 'P', 'nDraw', 'nQ', 'min_events');
    writetable(T, fullfile(out_dir, 'MUA_equalised_session_coherence.csv'));
    writetable(S, fullfile(out_dir, 'MUA_equalised_session_stats.csv'));
    writetable(Pout, fullfile(out_dir, 'MUA_equalised_pooled_AUC.csv'));

    for iA = 1:numel(analyses)
        plot_session_coherence(T, S, analyses{iA}, titles{iA}, v1_windows, ...
            fullfile(out_dir, ['MUA_equalised_session_coherence_' analyses{iA}]));
        plot_MUA_distributions(M, analyses{iA}, titles{iA}, v1_windows, ...
            fullfile(out_dir, ['MUA_equalised_MUA_distributions_' analyses{iA}]));
    end
    plot_pooled_AUC(P, analyses, titles, v1_windows, fullfile(out_dir, 'MUA_equalised_pooled_AUC'));
    fprintf('Saved to %s\n', out_dir);
end

function pct = within_session_percentile(X, session_ids)
    % Percentile (0-100) of each column of X within each session
    pct = NaN(size(X));
    for s = unique(session_ids)'
        idx = session_ids == s;
        for c = 1:size(X, 2)
            x = X(idx, c);
            ok = isfinite(x);
            p = NaN(size(x));
            p(ok) = (tiedrank(x(ok)) - 0.5) ./ sum(ok) * 100;
            pct(idx, c) = p;
        end
    end
end

function [masks, labels] = condition_masks(analysis, info, ripple_pct, spindle_pct, hemi)
    % Logical event masks {A, B}; spindle / SO from the dominant hemisphere
    n = numel(hemi);
    switch analysis
        case 'ripple_power'
            masks = {ripple_pct <= 25, ripple_pct > 75};
            labels = {'Low', 'High'};
        case 'spindle'
            f = dominant_feature(spindle_pct, hemi, n);
            masks = {f <= 25, f > 75};
            labels = {'Low', 'High'};
        case 'SO_phase'
            f = dominant_feature(info.SO_phase, hemi, n);
            masks = {f >= -pi/2 & f <= pi/2, (f >= -pi & f < -pi/2) | (f > pi/2 & f <= pi)};
            labels = {'Peak', 'Trough'};
    end
end

function f = dominant_feature(F, hemi, n)
    f = NaN(n, 1);
    f(hemi == 1) = F(hemi == 1, 1);
    f(hemi == 2) = F(hemi == 2, 2);
end

function [draws_a, draws_b] = match_events(ia, ib, x1, x2, nQ, nDraw, seed)
    % Equalise the joint (x1, x2) distribution of event sets ia and ib by
    % subsampling within nQ x nQ quantile bins (edges from both sets pooled)
    draws_a = repmat({zeros(0, 1)}, nDraw, 1); draws_b = draws_a;
    if isempty(ia) || isempty(ib)
        return
    end
    pooled = [ia; ib];
    ca = joint_bin(x1, x2, pooled, ia, nQ);
    cb = joint_bin(x1, x2, pooled, ib, nQ);
    cells = unique([ca; cb])';
    for d = 1:nDraw
        s = RandStream('philox4x32_10', 'Seed', seed * 100000 + d);
        sa = []; sb = [];
        for c = cells
            ea = ia(ca == c); eb = ib(cb == c);
            k = min(numel(ea), numel(eb));
            if k > 0
                sa = [sa; ea(randperm(s, numel(ea), k))]; %#ok<AGROW>
                sb = [sb; eb(randperm(s, numel(eb), k))]; %#ok<AGROW>
            end
        end
        draws_a{d} = sa; draws_b{d} = sb;
    end
end

function cid = joint_bin(x1, x2, pooled, idx, nQ)
    q1 = quantile_bin(x1(pooled), x1(idx), nQ);
    q2 = quantile_bin(x2(pooled), x2(idx), nQ);
    cid = (q1 - 1) * nQ + q2;
end

function q = quantile_bin(ref, x, nQ)
    edges = unique(prctile(ref, linspace(0, 100, nQ + 1)));
    edges(1) = -Inf; edges(end) = Inf;
    q = discretize(x, edges);
end

function a = auc_point(sel, meas)
    r = reactivation_coherence_auc(sel, meas, 0);
    a = r.AUC_mean;
end

function [Y, labels] = session_matrix(T, window, analysis, var)
    Ts = T(T.V1_window == window & T.analysis == analysis, :);
    labels = cellstr(unique(Ts.condition, 'stable'))';
    Y = NaN(max(T.session_id), numel(labels));
    for c = 1:numel(labels)
        Tc = Ts(Ts.condition == labels{c}, :);
        Y(Tc.session_id, c) = Tc.(var);
    end
end

function row = paired_test(analysis, window, meas, labels, a, b)
    ok = isfinite(a) & isfinite(b);
    p2 = NaN; p1 = NaN;
    if sum(ok) >= 3
        p2 = signrank(b(ok), a(ok));
        p1 = signrank(b(ok), a(ok), 'tail', 'right');
    end
    row = table(string(analysis), string(window), string(meas), string(sprintf('%s vs %s', labels{2}, labels{1})), ...
        sum(ok), mean(a(ok)), mean(b(ok)), mean(b(ok) - a(ok)), p2, p1, ...
        'VariableNames', {'analysis', 'V1_window', 'events', 'comparison', 'n_sessions', ...
        'mean_A', 'mean_B', 'mean_diff', 'p_signrank_two', 'p_signrank_one'});
end

function plot_session_coherence(T, S, analysis, title_str, v1_windows, save_path)
    % Rows: V1 window; columns: unmatched vs MUA-matched; lines = sessions
    cond_colors = {[0.6 0.6 0.6], [231, 41, 138] / 256};
    meas = ["unmatched", "matched"];
    fig = figure('Name', title_str, 'Position', [100 100 760 380 * numel(v1_windows)]);
    tiledlayout(numel(v1_windows), 2, 'TileSpacing', 'compact');
    for w = 1:numel(v1_windows)
        wn = v1_windows(w).name;
        for m = 1:2
            ax = nexttile; hold on;
            [Y, labels] = session_matrix(T, wn, analysis, "coh_" + meas(m));
            ok = all(isfinite(Y), 2);
            plot([1 2], Y(ok, :)', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
            for c = 1:2
                y = Y(ok, c);
                scatter(c * ones(size(y)), y, 25, cond_colors{c}, 'filled', 'MarkerFaceAlpha', 0.7);
                errorbar(c + 0.25 * (2 * c - 3), mean(y), std(y) / sqrt(numel(y)), 'o', 'Color', 'k', ...
                    'MarkerFaceColor', cond_colors{c}, 'LineWidth', 1.5, 'CapSize', 6);
            end
            yline(0, 'k--');
            Sr = S(S.analysis == analysis & S.V1_window == wn & S.events == meas(m), :);
            title({sprintf('%s, %s', v1_windows(w).label, meas(m)), ...
                sprintf('%s: p = %.3g (one-sided %.3g), n = %d', Sr.comparison, Sr.p_signrank_two, Sr.p_signrank_one, Sr.n_sessions)}, ...
                'FontSize', 9, 'FontWeight', 'normal');
            set(ax, 'XTick', [1 2], 'XTickLabel', labels, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 11);
            xlim([0.4 2.6]);
            ylabel('Geomean coherence (session mean)');
        end
    end
    linkaxes(findobj(fig, 'Type', 'axes'), 'y');
    sgtitle([title_str ': unmatched vs MUA-equalised events'], 'FontSize', 12);
    savefig(fig, [save_path '.fig']);
    exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
end

function plot_MUA_distributions(M, analysis, title_str, v1_windows, save_path)
    % HC and dominant V1 MUA, condition A vs B, before and after matching (draw 1)
    cond_colors = {[0.5 0.5 0.5], [231, 41, 138] / 256};
    fig = figure('Name', [title_str ' MUA'], 'Position', [100 100 1100 360 * numel(v1_windows)]);
    tiledlayout(numel(v1_windows), 4, 'TileSpacing', 'compact');
    vars = {'HC_MUA', 'V1_MUA'}; var_labels = {'HC MUA 0-0.1 s (z)', 'Dominant V1 MUA (z)'};
    for w = 1:numel(v1_windows)
        D = M.(sprintf('%s_%s', analysis, v1_windows(w).name));
        for v = 1:2
            x = D.(vars{v});
            edges = linspace(prctile(x([D.A; D.B]), 0.5), prctile(x([D.A; D.B]), 99.5), 50);
            sets = {{D.A, D.B}, {D.A_m, D.B_m}}; set_names = {'unmatched', 'matched'};
            for k = 1:2
                ax = nexttile; hold on;
                for c = 1:2
                    histogram(x(sets{k}{c}), edges, 'Normalization', 'probability', 'FaceColor', cond_colors{c}, ...
                        'EdgeColor', 'none', 'FaceAlpha', 0.5);
                end
                [~, pks] = kstest2(x(sets{k}{1}), x(sets{k}{2}));
                title({sprintf('%s, %s', v1_windows(w).label, set_names{k}), sprintf('n = %d / %d, KS p = %.2g', ...
                    numel(sets{k}{1}), numel(sets{k}{2}), pks)}, 'FontSize', 9, 'FontWeight', 'normal');
                xlabel(var_labels{v}); ylabel('Probability');
                if w == 1 && v == 1 && k == 1, legend(D.labels, 'Box', 'off'); end
                set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
            end
        end
    end
    sgtitle([title_str ': MUA before and after equalisation (pooled sessions, draw 1)'], 'FontSize', 12);
    savefig(fig, [save_path '.fig']);
    exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
end

function plot_pooled_AUC(P, analyses, titles, v1_windows, save_path)
    % Pooled V1->HPC coherence AUC: unmatched (diamond) and matched draws (violin-like scatter)
    cond_colors = {[0.5 0.5 0.5], [231, 41, 138] / 256};
    fig = figure('Name', 'Pooled AUC MUA-equalised', 'Position', [100 100 360 * numel(analyses) 340 * numel(v1_windows)]);
    tiledlayout(numel(v1_windows), numel(analyses), 'TileSpacing', 'compact');
    for w = 1:numel(v1_windows)
        for iA = 1:numel(analyses)
            ax = nexttile; hold on;
            Pr = P(P.V1_window == v1_windows(w).name & P.analysis == analyses{iA}, :);
            vals = {Pr.AUC_A_matched{1}, Pr.AUC_B_matched{1}};
            un = [Pr.AUC_A_unmatched Pr.AUC_B_unmatched];
            for c = 1:2
                swarm = c + 0.15 * (rand(size(vals{c})) - 0.5);
                scatter(swarm, vals{c}, 6, cond_colors{c}, 'filled', 'MarkerFaceAlpha', 0.3);
                plot(c + [-0.2 0.2], mean(vals{c}) * [1 1], 'k-', 'LineWidth', 2);
                scatter(c + 0.3, un(c), 50, cond_colors{c}, 'd', 'filled', 'MarkerEdgeColor', 'k');
            end
            yline(0, 'k--');
            set(ax, 'XTick', [1 2], 'XTickLabel', {Pr.cond_A, Pr.cond_B}, 'TickDir', 'out', 'Box', 'off', ...
                'Color', 'none', 'FontSize', 10);
            xlim([0.5 2.6]);
            ylabel('V1\rightarrowHPC coherence AUC');
            title({sprintf('%s, %s', titles{iA}, v1_windows(w).label), ...
                sprintf('matched diff %.3f [%.3f %.3f], P(diff<=0) = %.3f', Pr.AUC_diff_matched_mean, ...
                Pr.AUC_diff_CI_lo, Pr.AUC_diff_CI_hi, Pr.p_draw_diff_le0)}, 'FontSize', 8, 'FontWeight', 'normal');
        end
    end
    sgtitle('Pooled AUC: matched draws (dots, bar = mean) vs unmatched (diamond)', 'FontSize', 11);
    savefig(fig, [save_path '.fig']);
    exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
end
