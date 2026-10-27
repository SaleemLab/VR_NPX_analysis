%% V1->HPC reactivation coherence by spindle power and SO phase, inside and outside ripples, POST (R1)
% Does V1-HPC track-bias coherence depend on the dominant hemisphere's spindle
% power and SO phase only around ripples, or also in ripple-free windows?
%
% Windows: 200 ms, non-overlapping, starting every 0.2 s from -1 to 0.8 s around
% each non-duplicated POST SWS ripple onset (centres -0.9 to 0.9 s). Windows inside
% +-0.2 s of onset (-0.2-0, 0-0.2 s): every event. Others: ripple-free (no ripple of
% either probe overlaps) and closer to their own event's onset than to any other
% ripple onset (bilateral copies within 0.05 s count as one), as the |V1 bias| time
% course (nonripple_V1_abs_bias_POST).
% Ripple windows (-0.2-0, 0-0.2 s): the static AUC of the paper
%   (reactivation_AUC_coherence_PRE, POST_unique): z log odds + nan_mask, HPC bias
%   0-0.1 s, V1 bias over the window; spindle percentile (Q1 vs Q4; bottom 10% NaN)
%   and SO phase at ripple onset from ripple_info (hemisphere_conditions).
% Ripple-free windows: z log odds without nan_mask, V1 and HPC bias both the mean
%   over the same window; oscillations at that time (extract_SO_spindle_traces_POST,
%   reference V1 shank per hemisphere): SO phase at the window centre, peak
%   [-pi/2, pi/2] vs trough; spindle amplitude mean over the window, percentile
%   within session and hemisphere over all windows of the session, Low <= 25,
%   High > 75.
% Hemisphere matched to the V1 bias, as hemisphere_conditions: Track 1 events
% (V1 bias >= th) use hemisphere 2, Track 2 (V1 bias <= -th) hemisphere 1.
% AUC (V1 -> HPC, reactivation_coherence_auc): per window centre, pooled across
% sessions, nBoot event resamples (shared by both conditions, so the difference is
% paired); shuffled control resamples HPC only. p = fraction of resamples with
% AUC(High) - AUC(Low) <= 0 (Trough - Peak).
% Output: nature_neuroscience_revision_analysis\nonripple_oscillation_AUC
function nonripple_oscillation_AUC_POST()
    rev_dir = addon_setup();
    out_dir = fullfile(rev_dir, 'nonripple_oscillation_AUC');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    win_len = 0.2;
    win_starts = -1:win_len:(1 - win_len);
    nBoot = 1000;

    % --- Load non-duplicated POST data (rows aligned to ripple_info) ---
    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
    V1 = KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples;    % events x 200, no nan_mask
    HC = KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples;
    nan_mask = KDE_reactivation_ripples_PSTH.nan_mask;
    clear KDE_reactivation_ripples_PSTH
    load(fullfile(rev_dir, 'ripples_all_POST.mat'), 'ripples_all');
    load(fullfile(rev_dir, 'SO_spindle_traces_POST.mat'), 'traces', 'dt');
    sws1 = ripples_all(1).SWS_index == 1; sws2 = ripples_all(2).SWS_index == 1;
    onset_sws = [ripples_all(1).onset(sws1); ripples_all(2).onset(sws2)];
    sess_sws = [ripples_all(1).session_count(sws1); ripples_all(2).session_count(sws2)];
    assert(isequal(sess_sws(ripple_info.event_id(:)), ripple_info.session_count(:)), 'event_id mapping mismatch');
    t0 = onset_sws(ripple_info.event_id(:));
    sess = ripple_info.session_id(:);
    assert(isequal(sess, ripple_info.session_count(:)), 'session_id and session_count differ');

    % --- Window bias (events x windows) ---
    timebin = 0.01;
    bin_edges = -1:timebin:1;
    bin_centers = bin_edges(1:end-1) + timebin/2;
    nW = numel(win_starts);
    Wb = zeros(numel(bin_centers), nW);
    for w = 1:nW
        Wb(:, w) = bin_centers > win_starts(w) & bin_centers < win_starts(w) + win_len;
    end
    V1_w = nanmean_window(V1, Wb);
    HC_w = nanmean_window(HC, Wb);
    % Static AUC inputs of the ripple windows (as reactivation_AUC_coherence_PRE)
    V1_static = nanmean_window(V1 + nan_mask, Wb);
    HPC_static = nanmean_window(HC + nan_mask, double(bin_centers > 0 & bin_centers < 0.1)');
    clear V1 HC nan_mask
    t_rel = round((win_starts + win_len / 2) * 1000) / 1000;
    peri_w = win_starts >= -0.2 - 1e-9 & win_starts + win_len <= 0.2 + 1e-9;

    % --- Windows and their dominant-hemisphere oscillation features ---
    X = table();
    for s = unique(sess)'
        ev = find(sess == s);
        [on, off] = session_ripples(ripples_all, s);
        [E, K] = ndgrid(ev, 1:nW);
        ws = t0(E) + win_starts(K); we = ws + win_len; c = ws + win_len / 2;
        ip = discretize(we, [-Inf; on; Inf]) - 1;            % last ripple with onset < window end
        prev_off = -Inf(size(ws)); prev_off(ip > 0) = off(ip(ip > 0));
        next_on = Inf(size(ws)); hn = ip < numel(on); next_on(hn) = on(ip(hn) + 1);
        free = ws >= prev_off & next_on >= we;
        jn = discretize(c, [-Inf; on; Inf]) - 1;             % nearest ripple onset to the window centre
        on_prev = -Inf(size(c)); on_prev(jn > 0) = on(jn(jn > 0));
        on_next = Inf(size(c)); hn = jn < numel(on); on_next(hn) = on(jn(hn) + 1);
        nearest_on = on_next; use_prev = c - on_prev <= on_next - c; nearest_on(use_prev) = on_prev(use_prev);
        own = abs(nearest_on - t0(E)) <= 0.05;
        keep = peri_w(K) | (free & own);
        lin = sub2ind(size(V1_w), E(keep), K(keep));

        % SO phase at the window centre, spindle amplitude averaged over the window
        tr = traces(s);
        kc = round((c(keep) - tr.t(1)) / dt) + 1;
        k0 = round((ws(keep) - tr.t(1)) / dt) + 1;
        n_in = round(win_len / dt);
        ok = kc >= 1 & kc <= numel(tr.t) & k0 >= 1 & k0 + n_in - 1 <= numel(tr.t);
        SO = NaN(numel(kc), 2); SP = NaN(numel(kc), 2);
        SO(ok, :) = double(tr.SO_phase(kc(ok), :));
        idx = k0(ok) + (0:n_in - 1);                         % windows x samples
        for h = 1:2
            a = double(tr.spindle_amp(:, h));
            SP(ok, h) = mean(a(idx), 2);
        end
        SP_pct = NaN(size(SP));
        for h = 1:2
            g = isfinite(SP(:, h));
            SP_pct(g, h) = (tiedrank(SP(g, h)) - 0.5) / sum(g) * 100;
        end
        X = [X; table(repmat(s, sum(keep(:)), 1), t_rel(K(keep))', peri_w(K(keep))', V1_w(lin), HC_w(lin), ...
            SO, SP_pct, 'VariableNames', {'session', 't', 'peri', 'V1_bias', 'HC_bias', 'SO_phase', 'spindle_pct'})]; %#ok<AGROW>
    end
    X = X(isfinite(X.V1_bias) & isfinite(X.HC_bias) & all(isfinite(X.SO_phase), 2), :);

    % --- Conditions (columns: hemisphere 1, 2) ---
    S = AUC_low_high_style();
    peak = X.SO_phase >= -pi/2 & X.SO_phase <= pi/2;
    trough = (X.SO_phase >= -pi & X.SO_phase < -pi/2) | (X.SO_phase > pi/2 & X.SO_phase <= pi);
    an = struct('name', {'spindle', 'SO_phase'}, 'title', {S.spindle.title, S.SO_phase.title}, ...
        'labels', {S.spindle.labels, S.SO_phase.labels}, 'colors', {S.spindle.colors, S.SO_phase.colors}, ...
        'cond', {{X.spindle_pct <= 25, X.spindle_pct > 75}, {peak, trough}}, ...
        'ripple_cond', {hemisphere_conditions('spindle', ripple_info), hemisphere_conditions('SO', ripple_info)}, ...
        'ripple_idx', {[1 4], [1 2]});   % ripple windows: Q1 vs Q4, Peak vs Trough at ripple onset

    % --- Bootstrapped AUC per window centre ---
    ts = unique(X.t)';
    R = table();
    for iA = 1:numel(an)
        for j = 1:numel(ts)
            in_t = find(X.t == ts(j));
            is_peri = any(X.peri(in_t));
            if is_peri   % static AUC of the paper, all events
                sel = V1_static(:, abs(t_rel - ts(j)) < 1e-9)'; meas = HPC_static';
            else
                sel = X.V1_bias(in_t)'; meas = X.HC_bias(in_t)';
            end
            I = boot_indices(numel(sel), nBoot);
            a = NaN(nBoot, 2); as = a; n_ev = [0 0];
            for c = 1:2
                if is_peri
                    H = an(iA).ripple_cond(an(iA).ripple_idx(c));
                    c1 = H.c1(:)'; c2 = H.c2(:)';
                else
                    c1 = an(iA).cond{c}(in_t, 2)'; c2 = an(iA).cond{c}(in_t, 1)';
                end
                [a(:, c), as(:, c)] = boot_auc(sel, meas, c1, c2, I);
                n_ev(c) = sum((sel > 0 & c1) | (sel < 0 & c2));
            end
            d = a(:, 2) - a(:, 1);
            ci = @(x) reshape(prctile(x, [2.5 97.5], 1), 1, 2, 2);   % 1 x [lo hi] x condition
            R = [R; table(string(an(iA).name), ts(j), is_peri, numel(sel), n_ev, ...
                mean(a, 1, 'omitnan'), ci(a), mean(as, 1, 'omitnan'), ci(as), ...
                mean(d, 'omitnan'), prctile(d, [2.5 97.5]), mean(d <= 0), ...
                'VariableNames', {'analysis', 't', 'peri', 'n_windows', 'n_events', 'AUC', 'AUC_CI', ...
                'AUC_shuffled', 'AUC_shuffled_CI', 'diff', 'diff_CI', 'p_diff'})]; %#ok<AGROW>
            fprintf('%s t = %+.1f: n = %d, %s %.3f vs %s %.3f, p = %.3f\n', an(iA).name, ts(j), numel(sel), ...
                an(iA).labels{1}, mean(a(:, 1), 'omitnan'), an(iA).labels{2}, mean(a(:, 2), 'omitnan'), mean(d <= 0));
        end
    end

    save(fullfile(out_dir, 'nonripple_oscillation_AUC.mat'), 'R', 'X', 'nBoot', 'win_len', '-v7.3');
    Rcsv = table(R.analysis, R.t, R.peri, R.n_windows, ...
        'VariableNames', {'analysis', 't', 'peri', 'n_windows'});
    for c = 1:2   % condition 1 = Low / Peak, 2 = High / Trough
        Rcsv.(sprintf('n_events_%d', c)) = R.n_events(:, c);
        Rcsv.(sprintf('AUC_%d', c)) = R.AUC(:, c);
        Rcsv.(sprintf('AUC_CI_lo_%d', c)) = R.AUC_CI(:, 1, c);
        Rcsv.(sprintf('AUC_CI_hi_%d', c)) = R.AUC_CI(:, 2, c);
        Rcsv.(sprintf('AUC_shuffled_%d', c)) = R.AUC_shuffled(:, c);
        Rcsv.(sprintf('AUC_shuffled_CI_hi_%d', c)) = R.AUC_shuffled_CI(:, 2, c);
    end
    Rcsv.diff = R.diff; Rcsv.diff_CI_lo = R.diff_CI(:, 1); Rcsv.diff_CI_hi = R.diff_CI(:, 2); Rcsv.p_diff = R.p_diff;
    writetable(Rcsv, fullfile(out_dir, 'nonripple_oscillation_AUC.csv'));
    plot_AUC_time(R, an, win_len, fullfile(out_dir, 'nonripple_oscillation_AUC_time'));
    fprintf('Saved to %s\n', out_dir);
end

function X = nanmean_window(P, Wb)
    % Mean of P (events x bins) over the bins of each window column, ignoring NaN
    ok = isfinite(P);
    P(~ok) = 0;
    X = (P * Wb) ./ (double(ok) * Wb);
end

function [on, off] = session_ripples(ripples_all, s)
    % Sorted ripple onsets of both probes (bilateral copies with equal onset merged)
    % and the latest offset among ripples starting up to each onset
    on = [ripples_all(1).onset(ripples_all(1).session_count == s); ripples_all(2).onset(ripples_all(2).session_count == s)];
    off = [ripples_all(1).offset(ripples_all(1).session_count == s); ripples_all(2).offset(ripples_all(2).session_count == s)];
    [on, ~, g] = unique(on);
    off = cummax(accumarray(g, off, [], @max));
end

function I = boot_indices(n, nBoot)
    % nBoot x n event resamples, seeded per iteration as reactivation_coherence_auc
    I = zeros(nBoot, n, 'int32');
    for b = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', b);
        I(b, :) = randi(s, n, n, 1)';
    end
end

function [a, as] = boot_auc(sel, meas, c1, c2, I)
    % Vectorised reactivation_coherence_auc: real (sel, meas, conditions resampled
    % together) and shuffled (only meas resampled) AUC per resample, in chunks
    nBoot = size(I, 1); n = numel(sel);
    a = NaN(nBoot, 1); as = a;
    if n < 10, return; end
    th = prctile(abs(sel), 0:10:100); th = th(1:end-1);
    for b0 = 1:100:nBoot
        rows = b0:min(b0 + 99, nBoot);
        Ib = I(rows, :);
        S = sel(Ib); M = meas(Ib);
        ok = ~isnan(M); M(~ok) = 0;
        d = track_diff(S, M, ok & c1(Ib), ok & c2(Ib), th);
        ds = track_diff(repmat(sel, numel(rows), 1), M, ok & c1, ok & c2, th);
        a(rows) = trapz(th, d, 2) / (max(th) - min(th));
        as(rows) = trapz(th, ds, 2) / (max(th) - min(th));
    end
end

function d = track_diff(S, M, w1, w2, th)
    % mean(M | S >= th(j) & w1) - mean(M | S <= -th(j) & w2) for every threshold and
    % resample (rows); events are binned by the number of thresholds they pass, then
    % reverse cumulative sums over bins give the nested threshold sets
    [nBoot, n] = size(M); nT = numel(th);
    jp = zeros(nBoot, n); jn = jp;
    for j = 1:nT
        jp = jp + (S >= th(j)); jn = jn + (S <= -th(j));
    end
    row = repmat((1:nBoot)', 1, n);
    pass = @(jj, x) flip(cumsum(flip(reshape(accumarray(row(:) + nBoot * jj(:), x(:), [nBoot * (nT + 1) 1]), ...
        nBoot, nT + 1), 2), 2), 2);
    s1 = pass(jp, M .* w1); n1 = pass(jp, double(w1 & true(nBoot, 1)));
    s2 = pass(jn, M .* w2); n2 = pass(jn, double(w2 & true(nBoot, 1)));
    d = s1(:, 2:end) ./ n1(:, 2:end) - s2(:, 2:end) ./ n2(:, 2:end);
end

function save_figure(fig, save_path)
    savefig(fig, [save_path '.fig']);
    exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
    try
        exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
    catch err
        warning('nonripple:pdf', 'PDF not written (open in a viewer?): %s\n%s', [save_path '.pdf'], err.message);
    end
end

function plot_AUC_time(R, an, win_len, save_path)
    % Per window, as AUC_low_high_panel: Shuffle (mean of the two shuffled means,
    % CI pooled as min lower / max upper), then Low and High (Peak and Trough) bars
    % with bootstrap 95% CI. '*' above a bar: its CI does not overlap the shuffled CI.
    % Bracket: one-sided bootstrap p of High > Low (Trough > Peak) when < 0.05.
    fig = figure('Name', 'Oscillation AUC in and out of ripples', 'Position', [80 80 1400 400 * numel(an)]);
    tiledlayout(numel(an), 1, 'TileSpacing', 'compact');
    group_x = [1 2 3]; group_gap = 4;
    for iA = 1:numel(an)
        ax = nexttile; hold(ax, 'on');
        Ra = R(R.analysis == an(iA).name, :);
        nG = height(Ra);
        lims = [min([Ra.AUC_CI(:); Ra.AUC_shuffled_CI(:)]) max([Ra.AUC_CI(:); Ra.AUC_shuffled_CI(:)])];
        span = diff(lims);
        xt = NaN(1, nG); xtl = cell(1, nG); h = gobjects(1, 3);
        for g = 1:nG
            x = group_x + (g - 1) * group_gap;
            if Ra.peri(g)
                patch(ax, x([1 3 3 1]) + [-0.6 0.6 0.6 -0.6], lims(1) - span * [0.3 0.3 -1.5 -1.5], ...
                    [231 41 138] / 256, 'FaceAlpha', 0.08, 'EdgeColor', 'none');
            end
            ci_s = [min(Ra.AUC_shuffled_CI(g, 1, :)) max(Ra.AUC_shuffled_CI(g, 2, :))];
            m_s = mean(Ra.AUC_shuffled(g, :));
            h(1) = errorbar(ax, x(1), m_s, m_s - ci_s(1), ci_s(2) - m_s, 'ko', 'MarkerFaceColor', 'k', ...
                'MarkerSize', 4, 'LineWidth', 1.2, 'CapSize', 5);
            for c = 1:2
                m = Ra.AUC(g, c); lo = Ra.AUC_CI(g, 1, c); hi = Ra.AUC_CI(g, 2, c);
                h(c + 1) = bar(ax, x(c + 1), m, 0.7, 'FaceColor', an(iA).colors(c, :), 'EdgeColor', 'none');
                errorbar(ax, x(c + 1), m, m - lo, hi - m, 'k', 'LineWidth', 1.2, 'CapSize', 5);
                if lo > ci_s(2) || hi < ci_s(1)
                    text(ax, x(c + 1), hi + 0.01 * span, '*', 'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'bottom', 'FontSize', 14);
                end
            end
            st = stars(Ra.p_diff(g));
            if ~isempty(st)
                yb = max(Ra.AUC_CI(g, 2, :)) + 0.12 * span;
                plot(ax, x([2 2 3 3]), yb + [-0.02 0 0 -0.02] * span, 'k', 'LineWidth', 1);
                text(ax, mean(x(2:3)), yb, st, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
            end
            text(ax, x(2), lims(1) - 0.2 * span, sprintf('n = %d', Ra.n_windows(g)), 'HorizontalAlignment', 'center', ...
                'FontSize', 7, 'Color', [0.4 0.4 0.4]);
            t0 = Ra.t(g) - win_len / 2;
            xt(g) = x(2); xtl{g} = sprintf('%g to %g', round(t0, 2), round(t0 + win_len, 2));
        end
        yline(ax, 0, 'k--');
        set(ax, 'XTick', xt, 'XTickLabel', xtl, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
        xlim(ax, [0.3, xt(end) + 1.7]); ylim(ax, [lims(1) - 0.3 * span, lims(2) + 0.35 * span]);
        legend(ax, h, [{'Shuffle'}, an(iA).labels], 'Box', 'off', 'Location', 'northwest');
        ylabel(ax, 'Reactivation coherence AUC (V1\rightarrowHPC)');
        xlabel(ax, 'Window relative to ripple onset (s)');
        title(ax, an(iA).title, 'FontWeight', 'normal');
    end
    sgtitle({sprintf(['POST: reactivation coherence AUC in %d ms non-overlapping windows; pink: ripple windows ' ...
        '(static AUC, oscillations at ripple onset), others ripple-free (oscillations in the window)'], round(1000 * win_len)), ...
        'pooled sessions, bootstrap 95% CI; * CI outside shuffled CI; bracket: one-sided bootstrap p (High > Low, Trough > Peak)'}, ...
        'FontSize', 10);
    save_figure(fig, save_path);
end

function s = stars(p)
    if p < 0.001, s = '***'; elseif p < 0.01, s = '**'; elseif p < 0.05, s = '*'; else, s = ''; end
end
