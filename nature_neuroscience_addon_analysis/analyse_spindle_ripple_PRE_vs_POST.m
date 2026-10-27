%% Spindle-band power (dominant hemisphere) and ripple power distribution: PRE vs POST sleep
% Control for the PRE vs POST reactivation comparison: are ripples in POST
% accompanied by different spindle-band power, or drawn from a different ripple
% power distribution, than in PRE?
% Non-duplicated ripple events (ripple_info_PRE_unique / ripple_info, rows aligned
% to the *_unique KDE PSTH files).
% Spindle: z-scored spindle-band Hilbert amplitude at ripple peak (z-scored within
% each recording, so PRE vs POST compares power relative to that sleep's own mean);
% POST uses spindle_amplitude_peak (untruncated, see below).
% Dominant hemisphere from the V1 bias sign in each V1 window, as in
% MUA_equalised_geomean (V1 > 0 -> hemisphere 2 / V1 R, V1 < 0 -> hemisphere 1).
% Per session and sleep, for each V1 window (pre -0.2-0 s, post 0-0.2 s):
%   spindle_dom          : mean dominant-hemisphere spindle z
%   spindle_nondom       : mean non-dominant-hemisphere spindle z
%   spindle_dom_minus_nondom : mean (dominant - non-dominant) per event
%   p_spindle_dom_high   : fraction of events with dominant spindle z > 1
% Ripple power (peak z, averaged over merged bilateral pairs), per session:
%   ripple_power_mean / median / p90, n_events
%   SWS_duration_min     : cumulative SWS (behavioural_state_merged_all_PRE/POST.mat)
%   ripple_rate          : non-duplicated ripple events per minute of SWS
%   ripple_rate_first30  : same, restricted to the first 30 min of cumulative SWS
%                          (all SWS if the session has less than 30 min)
%   ripple_rate_matched  : paired sessions only, first min(30, PRE SWS, POST SWS) min
%                          of cumulative SWS in both sleeps (rate declines across SWS)
% Dominant spindle z, P(dominant spindle z > 1), ripple power and ripple rate are also
% computed for the first-30 and equal-SWS epochs (suffix _first30 / _matched), plotted in
% spindle_by_SWS_epoch, ripple_power_by_SWS_epoch and ripple_rate_by_SWS_epoch.
% PRE vs POST: paired signed-rank and exact sign-flip permutation (12 sessions with
% both), unpaired rank-sum (all sessions); all two-sided.
% Output: nature_neuroscience_revision_analysis\spindle_ripple_PRE_vs_POST
function analyse_spindle_ripple_PRE_vs_POST()
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    out_dir = fullfile(rev_dir, 'spindle_ripple_PRE_vs_POST');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % --- Load non-duplicated data (V1 bias only needed for the dominant hemisphere) ---
    load(fullfile(rev_dir, 'ripple_info_PRE_unique.mat'), 'ripple_info_PRE');
    load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_PRE_unique.mat'), 'KDE_reactivation_ripples_PSTH_PRE');
    D(1).name = "PRE";
    D(1).info = ripple_info_PRE;
    D(1).V1 = (KDE_reactivation_ripples_PSTH_PRE.V1_z_logodds_ripples + KDE_reactivation_ripples_PSTH_PRE.nan_mask)';
    clear KDE_reactivation_ripples_PSTH_PRE
    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
    D(2).name = "POST";
    D(2).info = ripple_info;
    D(2).V1 = (KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples + KDE_reactivation_ripples_PSTH.nan_mask)';
    clear KDE_reactivation_ripples_PSTH
    for iS = 1:2
        assert(numel(D(iS).info.session_id) == size(D(iS).V1, 2), '%s ripple_info / PSTH rows differ', D(iS).name);
    end

    % --- Windows ---
    timebin = 0.01;
    bin_edges = -1:timebin:1;
    bin_centers = bin_edges(1:end-1) + timebin/2;
    v1_windows = struct( ...
        'name',  {'pre', 'post'}, ...
        'label', {'V1 -0.2-0 s (pre)', 'V1 0-0.2 s (post)'}, ...
        'bins',  {bin_centers > -0.2 & bin_centers < 0, bin_centers > 0 & bin_centers < 0.2});
    spindle_high_z = 1;
    sws_cap_min = 30;    % ripple rate over the first 30 min of cumulative SWS

    % --- Per-event dominant / non-dominant spindle z ---
    % POST ripple_info.spindle_amplitude has events below the 10th within-session
    % percentile set to NaN; spindle_amplitude_peak is the same quantity without that
    % truncation, so it is used to keep the POST distribution comparable to PRE.
    for iS = 1:2
        amp = D(iS).info.spindle_amplitude;
        if isfield(D(iS).info, 'spindle_amplitude_peak')
            ok = isfinite(amp);
            assert(max(abs(amp(ok) - D(iS).info.spindle_amplitude_peak(ok))) < 1e-9, ...
                '%s spindle_amplitude_peak differs from spindle_amplitude', D(iS).name);
            amp = D(iS).info.spindle_amplitude_peak;
        end
        n = size(amp, 1);
        for w = 1:numel(v1_windows)
            V1_b = mean(D(iS).V1(v1_windows(w).bins, :), 'omitnan')';
            hemi = NaN(n, 1); hemi(V1_b < 0) = 1; hemi(V1_b > 0) = 2;
            dom = NaN(n, 1); nondom = NaN(n, 1);
            dom(hemi == 1) = amp(hemi == 1, 1);    nondom(hemi == 1) = amp(hemi == 1, 2);
            dom(hemi == 2) = amp(hemi == 2, 2);    nondom(hemi == 2) = amp(hemi == 2, 1);
            D(iS).spindle_dom.(v1_windows(w).name) = dom;
            D(iS).spindle_nondom.(v1_windows(w).name) = nondom;
        end
    end

    % --- Event selection by cumulative SWS epoch (all / first30 / matched, see sws_epoch_masks) ---
    epochs = {'all', 'first30', 'matched'};
    epoch_labels = {'All SWS', sprintf('First %d min SWS', sws_cap_min), 'Equal SWS (PRE = POST)'};
    E = sws_epoch_masks(D(1).info, D(2).info, analysis_folder, sws_cap_min);
    for iS = 1:2
        for f = {'sws_session_id', 'sws_bouts', 'sws_duration', 'peaktime', 'sel', 'used_min'}
            D(iS).(f{1}) = E(iS).(f{1});
        end
    end
    paired_ids = intersect(D(1).sws_session_id, D(2).sws_session_id);

    % --- Per-session metrics (suffix _first30 / _matched for the restricted epochs) ---
    T = table();
    for iS = 1:2
        info = D(iS).info;
        sess = info.session_id(:);
        for k = 1:numel(D(iS).sws_session_id)
            s = D(iS).sws_session_id(k);
            row = table(s, info.subject_id(find(sess == s, 1)), D(iS).name, ...
                'VariableNames', {'session_id', 'subject_id', 'sleep'});
            for e = 1:numel(epochs)
                in_e = sess == s & D(iS).sel.(epochs{e});
                sfx = ''; if e > 1, sfx = ['_' epochs{e}]; end
                used = D(iS).used_min.(epochs{e})(k);
                rp = info.ripple_power(in_e);
                row.(['n_events' sfx]) = sum(in_e);
                row.(['SWS_duration' sfx '_min']) = used;
                row.(['ripple_rate' sfx]) = sum(in_e) / used;
                row.(['ripple_power_mean' sfx]) = mean(rp);
                row.(['ripple_power_median' sfx]) = median(rp);
                row.(['ripple_power_p90' sfx]) = quantile_or_nan(rp, 0.9);
                for w = 1:numel(v1_windows)
                    wn = v1_windows(w).name;
                    dom = D(iS).spindle_dom.(wn)(in_e); nondom = D(iS).spindle_nondom.(wn)(in_e);
                    ok = isfinite(dom) & isfinite(nondom);
                    row.(['spindle_dom_' wn sfx]) = mean(dom(ok));
                    row.(['p_spindle_dom_high_' wn sfx]) = mean(dom(ok) > spindle_high_z);
                    if e == 1
                        row.(['spindle_nondom_' wn]) = mean(nondom(ok));
                        row.(['spindle_dom_minus_nondom_' wn]) = mean(dom(ok) - nondom(ok));
                        row.(['n_spindle_' wn]) = sum(ok);
                    end
                end
            end
            T = [T; row]; %#ok<AGROW>
        end
    end
    writetable(T, fullfile(out_dir, 'spindle_ripple_per_session.csv'));

    % --- PRE vs POST tests (two-sided) ---
    metrics = {'spindle_dom_pre', 'spindle_dom_post', 'spindle_nondom_pre', 'spindle_nondom_post', ...
        'spindle_dom_minus_nondom_pre', 'spindle_dom_minus_nondom_post', 'p_spindle_dom_high_pre', 'p_spindle_dom_high_post', ...
        'ripple_power_mean', 'ripple_power_median', 'ripple_power_p90', 'n_events', 'SWS_duration_min', 'ripple_rate'};
    metric_labels = {'Dominant spindle z (V1 pre)', 'Dominant spindle z (V1 post)', ...
        'Non-dominant spindle z (V1 pre)', 'Non-dominant spindle z (V1 post)', ...
        'Dominant - non-dominant z (V1 pre)', 'Dominant - non-dominant z (V1 post)', ...
        sprintf('P(dominant spindle z > %g) (V1 pre)', spindle_high_z), sprintf('P(dominant spindle z > %g) (V1 post)', spindle_high_z), ...
        'Ripple power, mean (peak z)', 'Ripple power, median (peak z)', 'Ripple power, 90th pct (peak z)', 'Number of ripple events', ...
        'SWS duration used (min)', 'Ripple rate (events / min SWS)'};
    % Spindle, ripple power and ripple rate restricted to the first-30 / equal-SWS epochs
    epoch_base = {'spindle_dom_pre', 'spindle_dom_post', 'p_spindle_dom_high_pre', 'p_spindle_dom_high_post', ...
        'ripple_power_mean', 'ripple_power_median', 'ripple_power_p90', 'SWS_duration_min', 'n_events', 'ripple_rate'};
    for e = 2:numel(epochs)
        for b = 1:numel(epoch_base)
            metric_labels{end+1} = metric_labels{strcmp(metrics, epoch_base{b})}; %#ok<AGROW>
            metrics{end+1} = epoch_metric(epoch_base{b}, epochs{e}); %#ok<AGROW>
        end
    end
    S = table();
    for m = 1:numel(metrics)
        [x, y] = pre_post_matrix(T, metrics{m});
        S = [S; prepost_tests(metrics{m}, x, y)]; %#ok<AGROW>
    end
    disp(S);
    writetable(S, fullfile(out_dir, 'spindle_ripple_stats.csv'));

    % --- Event-level distributions (matched sessions) and per-session histograms ---
    pre_ids = unique(D(1).info.session_id(:));
    dists = struct( ...
        'name',  {'spindle_dom_pre', 'spindle_dom_post', 'ripple_power'}, ...
        'label', {'Dominant spindle z (V1 pre)', 'Dominant spindle z (V1 post)', 'Ripple power (peak z)'}, ...
        'edges', {linspace(-1.5, 5, 40), linspace(-1.5, 5, 40), linspace(2, 25, 47)});
    H = struct();
    for k = 1:numel(dists)
        for iS = 1:2
            v = event_values(D(iS), dists(k).name);
            sess = D(iS).info.session_id(:);
            matched = ismember(sess, pre_ids) & isfinite(v);
            H.(dists(k).name).(D(iS).name).values = v(matched);
            H.(dists(k).name).(D(iS).name).session = sess(matched);
            hs = NaN(numel(pre_ids), numel(dists(k).edges) - 1);
            for i = 1:numel(pre_ids)
                vs = v(sess == pre_ids(i) & isfinite(v));
                hs(i, :) = histcounts(min(max(vs, dists(k).edges(1)), dists(k).edges(end)), dists(k).edges, ...
                    'Normalization', 'probability');
            end
            H.(dists(k).name).(D(iS).name).session_hist = hs;
        end
        % Per-session two-sample KS (PRE vs POST events of the same session)
        ks = NaN(numel(pre_ids), 2);
        for i = 1:numel(pre_ids)
            a = H.(dists(k).name).PRE.values(H.(dists(k).name).PRE.session == pre_ids(i));
            b = H.(dists(k).name).POST.values(H.(dists(k).name).POST.session == pre_ids(i));
            [~, ks(i, 2), ks(i, 1)] = kstest2(a, b);
        end
        H.(dists(k).name).ks_D = ks(:, 1); H.(dists(k).name).ks_p = ks(:, 2);
        fprintf('%s: per-session KS PRE vs POST, %d/%d sessions p < 0.05, median D = %.3f\n', ...
            dists(k).name, sum(ks(:, 2) < 0.05), numel(pre_ids), median(ks(:, 1)));
    end
    save(fullfile(out_dir, 'spindle_ripple_PRE_vs_POST.mat'), 'T', 'S', 'H', 'pre_ids', 'spindle_high_z');

    % --- Plot 1: spindle power (incl. non-dominant), all SWS, matched sessions ---
    fig = figure('Position', [50 50 1600 750]);
    tiledlayout(2, 4, 'TileSpacing', 'compact');
    for m = 1:8
        ax = nexttile;
        [x, y] = pre_post_matrix(T, metrics{m});
        n = plot_matched_pre_post(ax, x, y, [1 2]);
        st = S(S.metric == metrics{m}, :);
        set(ax, 'XTick', [1 2], 'XTickLabel', {'PRE', 'POST'}, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
        xlim(ax, [0.4 2.6]);
        ylabel(ax, metric_labels{m}, 'FontSize', 9);
        title(ax, {sprintf('perm p = %.3g, signrank p = %.3g', st.p_perm_two_sided, st.p_signrank_two_sided), ...
            sprintf('n = %d matched; %d/%d POST > PRE', n, st.n_POST_gt_PRE, n)}, 'FontSize', 8, 'FontWeight', 'normal');
    end
    sgtitle('Spindle-band power at ripples (dominant hemisphere): PRE vs POST, matched sessions', 'FontSize', 12);
    savefig(fig, fullfile(out_dir, 'spindle_metrics_matched.fig'));
    exportgraphics(fig, fullfile(out_dir, 'spindle_metrics_matched.png'), 'Resolution', 300);

    % --- Plot 1b: spindle power / ripple power / ripple rate by SWS epoch (rows), matched sessions ---
    % (+ pooled CDF column for spindle and ripple power)
    sleep_colors = {[0.5 0.5 0.5], [231, 41, 138] / 256};
    efigs = struct( ...
        'file',  {'spindle_by_SWS_epoch', 'ripple_power_by_SWS_epoch', 'ripple_rate_by_SWS_epoch'}, ...
        'title', {'Dominant-hemisphere spindle z at ripples by SWS epoch: PRE vs POST, matched sessions', ...
                  'Ripple power by SWS epoch: PRE vs POST, matched sessions', ...
                  'Ripple rate by SWS epoch: PRE vs POST, matched sessions'}, ...
        'base',  {epoch_base(1:4), epoch_base(5:7), epoch_base(8:10)}, ...
        'dist',  {'spindle_dom_post', 'ripple_power', ''}, ...
        'dist_label', {'Dominant spindle z (V1 post)', 'Ripple power (peak z)', ''}, ...
        'xlim',  {[-1.5 5], [5 25], []});
    for iF = 1:numel(efigs)
        has_dist = ~isempty(efigs(iF).dist);
        nCol = numel(efigs(iF).base) + has_dist;
        fig = figure('Position', [50 50 380 * nCol 1000]);
        tiledlayout(numel(epochs), nCol, 'TileSpacing', 'compact');
        for e = 1:numel(epochs)
            for b = 1:numel(efigs(iF).base)
                mname = epoch_metric(efigs(iF).base{b}, epochs{e});
                ax = nexttile;
                [x, y] = pre_post_matrix(T, mname);
                n = plot_matched_pre_post(ax, x, y, [1 2]);
                st = S(S.metric == mname, :);
                set(ax, 'XTick', [1 2], 'XTickLabel', {'PRE', 'POST'}, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
                xlim(ax, [0.4 2.6]);
                ylabel(ax, metric_labels{strcmp(metrics, efigs(iF).base{b})}, 'FontSize', 9);
                title(ax, {epoch_labels{e}, sprintf('perm p = %.3g, signrank p = %.3g', st.p_perm_two_sided, ...
                    st.p_signrank_two_sided), sprintf('n = %d matched; %d/%d POST > PRE', n, st.n_POST_gt_PRE, n)}, ...
                    'FontSize', 8, 'FontWeight', 'normal');
            end
            if ~has_dist
                continue
            end
            ax = nexttile; hold(ax, 'on');
            v = cell(1, 2);
            for iS = 1:2
                vv = event_values(D(iS), efigs(iF).dist);
                keep = D(iS).sel.(epochs{e}) & ismember(D(iS).info.session_id(:), paired_ids) & isfinite(vv);
                v{iS} = vv(keep);
                [f, xs] = ecdf(v{iS});
                plot(ax, xs, f, '-', 'Color', sleep_colors{iS}, 'LineWidth', 1.5, ...
                    'DisplayName', sprintf('%s (n = %d)', D(iS).name, numel(v{iS})));
            end
            [~, ~, ks_d] = kstest2(v{1}, v{2});
            xlim(ax, efigs(iF).xlim); xlabel(ax, efigs(iF).dist_label); ylabel(ax, 'Cumulative fraction');
            legend(ax, 'Location', 'southeast', 'Box', 'off');
            title(ax, {sprintf('%s: pooled ripples', epoch_labels{e}), sprintf('KS D = %.3f', ks_d)}, ...
                'FontSize', 8, 'FontWeight', 'normal');
            set(ax, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
        end
        sgtitle(efigs(iF).title, 'FontSize', 12);
        savefig(fig, fullfile(out_dir, [efigs(iF).file '.fig']));
        exportgraphics(fig, fullfile(out_dir, [efigs(iF).file '.png']), 'Resolution', 300);
    end

    % --- Plot 2: distributions (top: mean +/- SEM of session histograms; bottom: pooled CDF) ---
    sleep_colors = {[0.5 0.5 0.5], [231, 41, 138] / 256};
    fig = figure('Position', [100 100 1400 700]);
    tiledlayout(2, numel(dists), 'TileSpacing', 'compact');
    for k = 1:numel(dists)
        nexttile(k); hold on;
        e = dists(k).edges; centers = e(1:end-1) + diff(e) / 2;
        for iS = 1:2
            M = H.(dists(k).name).(D(iS).name).session_hist;
            mu = mean(M, 1); se = std(M, [], 1) / sqrt(size(M, 1));
            fill([centers fliplr(centers)], [mu - se fliplr(mu + se)], sleep_colors{iS}, 'FaceAlpha', 0.2, ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(centers, mu, '-', 'Color', sleep_colors{iS}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%s (n = %d sessions)', D(iS).name, size(M, 1)));
        end
        xlabel(dists(k).label); ylabel('Fraction of ripples');
        legend('Location', 'best', 'Box', 'off');
        title({dists(k).label, sprintf('per-session KS: %d/%d p < 0.05, median D = %.3f', ...
            sum(H.(dists(k).name).ks_p < 0.05), numel(pre_ids), median(H.(dists(k).name).ks_D))}, ...
            'FontSize', 9, 'FontWeight', 'normal');
        set(gca, 'TickDir', 'out', 'Box', 'off');

        nexttile(k + numel(dists)); hold on;
        for iS = 1:2
            v = H.(dists(k).name).(D(iS).name).values;
            [f, xs] = ecdf(v);
            plot(xs, f, '-', 'Color', sleep_colors{iS}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%s (n = %d ripples)', D(iS).name, numel(v)));
        end
        xlim(e([1 end])); xlabel(dists(k).label); ylabel('Cumulative fraction');
        legend('Location', 'southeast', 'Box', 'off');
        title('Pooled ripples, matched sessions', 'FontSize', 9, 'FontWeight', 'normal');
        set(gca, 'TickDir', 'out', 'Box', 'off');
    end
    sgtitle('Distributions at ripple times: PRE vs POST (matched sessions)', 'FontSize', 12);
    savefig(fig, fullfile(out_dir, 'spindle_ripple_distributions.fig'));
    exportgraphics(fig, fullfile(out_dir, 'spindle_ripple_distributions.png'), 'Resolution', 300);
    fprintf('Saved to %s\n', out_dir);
end

function name = epoch_metric(base, epoch)
% Per-session metric name for an SWS epoch: base (all), base_first30, base_matched;
% SWS_duration_min -> SWS_duration_first30_min
if strcmp(epoch, 'all')
    name = base;
elseif endsWith(base, '_min')
    name = [extractBefore(base, strlength(base) - 3) '_' epoch '_min'];
else
    name = [base '_' epoch];
end
end

function q = quantile_or_nan(x, p)
% prctile that returns NaN for an empty input
if isempty(x), q = NaN; else, q = prctile(x, 100 * p); end
end

function v = event_values(Ds, name)
% Per-event values for the distribution plots
switch name
    case 'ripple_power'
        v = Ds.info.ripple_power(:);
    case 'spindle_dom_pre'
        v = Ds.spindle_dom.pre;
    case 'spindle_dom_post'
        v = Ds.spindle_dom.post;
end
end

function [x, y] = pre_post_matrix(T, metric)
% PRE (x) and POST (y) values indexed by session_id; NaN where missing
n = max(T.session_id);
x = NaN(n, 1); y = NaN(n, 1);
Tp = T(T.sleep == "PRE", :);  x(Tp.session_id) = Tp.(metric);
Tq = T(T.sleep == "POST", :); y(Tq.session_id) = Tq.(metric);
end

function row = prepost_tests(metric, x, y)
% Paired signed-rank and permutation (sessions with both), unpaired rank-sum (all); two-sided
ok = isfinite(x) & isfinite(y);
d = y(ok) - x(ok);
ps2 = signrank(y(ok), x(ok));
pp2 = paired_perm_test(d, 'right');
xa = x(isfinite(x)); ya = y(isfinite(y));
pr2 = ranksum(ya, xa);
row = table(string(metric), sum(ok), mean(x(ok)), mean(y(ok)), mean(d), sum(d > 0), ps2, pp2, ...
    numel(xa), numel(ya), mean(xa), mean(ya), pr2, ...
    'VariableNames', {'metric', 'n_paired', 'mean_PRE', 'mean_POST', 'mean_diff', 'n_POST_gt_PRE', ...
    'p_signrank_two_sided', 'p_perm_two_sided', 'n_pre_all', 'n_post_all', 'mean_PRE_all', 'mean_POST_all', ...
    'p_ranksum_two_sided'});
end
