%% |V1 track bias| in ripple-free 200 ms windows vs peri-ripple windows, POST (R4.1)
% Is V1 track bias specific to ripples, or as strong in ripple-free periods with
% the same V1 activity?
%
% Bias: original KDE z-scored log odds PSTH (-1 to 1 s around each non-duplicated
% POST SWS ripple, 20 ms bins, 10 ms step), WITHOUT nan_mask. Window bias = mean
% over the bins in a 200 ms window; |bias| is its absolute value.
%   Peri-ripple windows: 0-0.2 s (post) and -0.2-0 s (pre) of every event.
%   Ripple-free windows: 200 ms windows starting every 0.1 s from -1 to 0.8 s,
%     overlapping no ripple of either probe (SWS or awake, ripples_all_POST
%     onset-offset), at least min_gap s from the nearest ripple edge. Windows of
%     neighbouring events cover the same time, so within a session, windows
%     overlapping an earlier kept window (in absolute time) are dropped.
%   Distance = gap between the window edge and the nearest ripple edge; signed
%     negative when the nearest ripple follows the window (pre-ripple), positive
%     when it precedes it (post-ripple).
% MUA (z-scored, 10 ms bins, -2 to 2 s; ripple_MUA_*): mean over the same window.
% Equalisation: within session, peri-ripple and ripple-free windows are subsampled
% so their joint distribution of V1 L and V1 R MUA is equal (nQ x nQ quantile
% bins, edges from both sets pooled; min(nA, nB) windows kept per bin, nDraw
% seeded draws). Matching on both hemispheres does not depend on the bias sign.
% Distance curve: per session and distance bin, mean |bias| and V1 MUA, and the
% matched difference |bias|(ripple-free bin) - |bias|(peri-ripple post), each bin
% equalised separately against the post windows.
% Time course: |bias| vs window centre relative to ripple onset, 100 ms
% non-overlapping windows (centres -0.95 to 0.95 s). Windows inside +-0.2 s of
% onset: every event. Others:
% ripple-free windows (no ripple overlap, no min_gap) whose nearest ripple onset
% is the event's own (within 0.05 s), so each stretch of time is counted for one
% ripple only. Plotted as mean +- SE across sessions, and as session means and
% animal means without shade. MUA-matched version: within session, the pooled
% peri-ripple windows are the reference (kept whole) and the ripple-free windows
% at each centre are subsampled to the reference V1 L x V1 R MUA distribution
% (nQ x nQ cells, reference quantile edges, nDraw draws). MUA matching check: V1
% MUA distributions at -0.4 to -0.3 s and 0.3 to 0.4 s, before / after matching,
% against the reference.
% Output: nature_neuroscience_revision_analysis\nonripple_V1_bias
function nonripple_V1_abs_bias_POST()
    rev_dir = addon_setup();
    out_dir = fullfile(rev_dir, 'nonripple_V1_bias');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    win_len = 0.2;
    win_starts = -1:0.1:(1 - win_len);
    min_gap = 0.1;             % s, ripple-free window to nearest ripple edge
    nDraw = 100;
    nQ = 5;                    % quantile bins per MUA dimension
    min_windows = 10;          % per set per session after matching
    min_windows_dist = 5;      % per distance bin per session after matching
    min_windows_tc = 20;       % per time-course point per session after matching
    tc_len = 0.1;              % s, non-overlapping time-course windows
    tc_name = 'nonripple_V1_abs_bias_time_course_100ms';
    dist_edges = 0.1:0.2:0.9;  % |distance| bins (s); 0.2 s wide so each holds one 0.1 s grid step after de-overlapping

    % --- Load non-duplicated POST data (rows aligned to ripple_info) ---
    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
    V1 = KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples;    % events x 200, no nan_mask
    HC = KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples;
    clear KDE_reactivation_ripples_PSTH
    load(fullfile(rev_dir, 'ripple_MUA_HPC.mat'), 'ripple_MUA_HPC');
    load(fullfile(rev_dir, 'ripple_MUA_V1_L.mat'), 'ripple_MUA_V1_L');
    load(fullfile(rev_dir, 'ripple_MUA_V1_R.mat'), 'ripple_MUA_V1_R');
    load(fullfile(rev_dir, 'ripples_all_POST.mat'), 'ripples_all');
    n_ev = numel(ripple_info.session_id);
    assert(size(V1, 1) == n_ev && size(ripple_MUA_V1_L, 1) == n_ev, 'PSTH / MUA rows do not match ripple_info');

    % Onset of each event: PSTH rows are probe 1 SWS ripples, then probe 2 (event_id)
    sws1 = ripples_all(1).SWS_index == 1; sws2 = ripples_all(2).SWS_index == 1;
    onset_sws = [ripples_all(1).onset(sws1); ripples_all(2).onset(sws2)];
    sess_sws = [ripples_all(1).session_count(sws1); ripples_all(2).session_count(sws2)];
    assert(isequal(sess_sws(ripple_info.event_id(:)), ripple_info.session_count(:)), 'event_id mapping mismatch');
    t0 = onset_sws(ripple_info.event_id(:));
    sess = ripple_info.session_id(:);
    sessions = unique(sess)';

    timebin = 0.01;
    bin_edges = -1:timebin:1;
    bin_centers = bin_edges(1:end-1) + timebin/2;
    mua_centers = -2 + timebin/2:timebin:2 - timebin/2;

    % --- Time course: 100 ms non-overlapping windows ---
    C = time_course_windows(tc_len, -1:tc_len:(1 - tc_len), V1, ripple_MUA_V1_L, ripple_MUA_V1_R, ...
        bin_centers, mua_centers, t0, sess, ripple_info.subject_id(:), ripples_all);
    TC = time_course(C, sessions, tc_len, nQ, nDraw, min_windows_tc);
    clear C
    fprintf('Time course, %d ms windows:\n', round(1000 * tc_len));
    disp(array2table([TC.t; sum(TC.n_windows, 1); mean(TC.session_mean, 1, 'omitnan'); ...
        sum(TC.n_matched, 1); sum(isfinite(TC.session_mean_matched), 1); ...
        mean(TC.session_mean_matched, 1, 'omitnan'); mean(TC.MUA_raw, 1, 'omitnan'); ...
        mean(TC.MUA_matched, 1, 'omitnan')]', 'VariableNames', {'t', 'n_windows', 'absV1', ...
        'n_matched', 'n_sessions_matched', 'absV1_matched', 'V1_MUA', 'V1_MUA_matched'}));
    fprintf('Peri-ripple reference V1 MUA (mean of sessions): %.3f\n', mean(TC.MUA_reference));

    % --- Window means (events x windows) ---
    nW = numel(win_starts);
    Wb = zeros(numel(bin_centers), nW); Wm = zeros(numel(mua_centers), nW);
    for w = 1:nW
        Wb(:, w) = bin_centers > win_starts(w) & bin_centers < win_starts(w) + win_len;
        Wm(:, w) = mua_centers > win_starts(w) & mua_centers < win_starts(w) + win_len;
    end
    V1_w = nanmean_window(V1, Wb);
    HC_w = nanmean_window(HC, Wb);
    MUA_L = double(ripple_MUA_V1_L) * (Wm ./ sum(Wm, 1));
    MUA_R = double(ripple_MUA_V1_R) * (Wm ./ sum(Wm, 1));
    MUA_HC = double(ripple_MUA_HPC) * (Wm ./ sum(Wm, 1));
    clear V1 HC ripple_MUA_V1_L ripple_MUA_V1_R ripple_MUA_HPC
    i_pre = find(abs(win_starts + 0.2) < 1e-9); i_post = find(abs(win_starts) < 1e-9);

    % --- Window table: peri-ripple (pre, post) and ripple-free windows ---
    W = table();
    for s = sessions
        ev = find(sess == s);
        [on, off] = session_ripples(ripples_all, s);

        for k = [i_pre i_post]
            if k == i_pre, type = "ripple_pre"; else, type = "ripple_post"; end
            W = [W; window_rows(s, ev, k, win_starts(k), type, t0(ev) + win_starts(k), NaN(numel(ev), 1), ...
                V1_w, HC_w, MUA_L, MUA_R, MUA_HC)]; %#ok<AGROW>
        end

        % Ripple-free candidates: all events x windows
        [E, K] = ndgrid(ev, 1:nW);
        ws = t0(E) + win_starts(K); we = ws + win_len;
        ip = discretize(we, [-Inf; on; Inf]) - 1;            % last ripple with onset < window end
        prev_off = -Inf(size(ws)); prev_off(ip > 0) = off(ip(ip > 0));
        next_on = Inf(size(ws)); has_next = ip < numel(on); next_on(has_next) = on(ip(has_next) + 1);
        since_prev = ws - prev_off; to_next = next_on - we;
        free = since_prev >= min_gap & to_next >= min_gap;
        signed_d = since_prev; pre = to_next < since_prev; signed_d(pre) = -to_next(pre);
        E = E(free); K = K(free); ws = ws(free); signed_d = signed_d(free);
        % Drop windows overlapping an earlier kept window (same time from another event)
        [ws, o] = sort(ws); E = E(o); K = K(o); signed_d = signed_d(o);
        keep = false(size(ws)); last_end = -Inf;
        for i = 1:numel(ws)
            if ws(i) >= last_end - 1e-9
                keep(i) = true; last_end = ws(i) + win_len;
            end
        end
        lin = sub2ind(size(V1_w), E(keep), K(keep));
        W = [W; table(repmat(s, sum(keep), 1), E(keep), win_starts(K(keep))', repmat("nonripple", sum(keep), 1), ...
            ws(keep), signed_d(keep), V1_w(lin), HC_w(lin), MUA_L(lin), MUA_R(lin), MUA_HC(lin), ...
            'VariableNames', window_vars())]; %#ok<AGROW>
    end
    W.abs_V1 = abs(W.V1_bias);
    W.abs_HC = abs(W.HC_bias);
    W.V1_MUA = (W.MUA_V1_L + W.MUA_V1_R) / 2;
    W = W(isfinite(W.V1_bias) & isfinite(W.MUA_V1_L) & isfinite(W.MUA_V1_R), :);
    fprintf('Windows: %d post, %d pre, %d ripple-free\n', sum(W.type == "ripple_post"), ...
        sum(W.type == "ripple_pre"), sum(W.type == "nonripple"));

    % --- Peri-ripple vs ripple-free: per session, unmatched and V1 MUA-equalised ---
    refs = ["ripple_post", "ripple_pre"];
    T = table();      % per session x reference
    M = struct();     % matched windows of draw 1 (pooled distributions)
    for r = refs
        M.(char(r)) = struct('ref', [], 'free', []);
        for s = sessions
            ia = find(W.session == s & W.type == r);
            ib = find(W.session == s & W.type == "nonripple");
            [da, db] = match_sets(ia, ib, W.MUA_V1_L, W.MUA_V1_R, nQ, nDraw, s);
            n_m = numel(da{1});
            m_ref = NaN; m_free = NaN; mua_ref = NaN; mua_free = NaN;
            if n_m >= min_windows
                m_ref = mean(cellfun(@(i) mean(W.abs_V1(i)), da));
                m_free = mean(cellfun(@(i) mean(W.abs_V1(i)), db));
                mua_ref = mean(cellfun(@(i) mean(W.V1_MUA(i)), da));
                mua_free = mean(cellfun(@(i) mean(W.V1_MUA(i)), db));
                M.(char(r)).ref = [M.(char(r)).ref; da{1}];
                M.(char(r)).free = [M.(char(r)).free; db{1}];
            end
            T = [T; table(s, r, numel(ia), numel(ib), n_m * (n_m >= min_windows), ...
                mean(W.abs_V1(ia)), mean(W.abs_V1(ib)), m_ref, m_free, ...
                mean(W.V1_MUA(ia)), mean(W.V1_MUA(ib)), mua_ref, mua_free, ...
                mean(W.abs_HC(ia), 'omitnan'), mean(W.abs_HC(ib), 'omitnan'), ...
                'VariableNames', {'session_id', 'reference', 'n_ref', 'n_free', 'n_matched', ...
                'absV1_ref', 'absV1_free', 'absV1_ref_matched', 'absV1_free_matched', ...
                'V1MUA_ref', 'V1MUA_free', 'V1MUA_ref_matched', 'V1MUA_free_matched', ...
                'absHC_ref', 'absHC_free'})]; %#ok<AGROW>
        end
    end
    S = table();
    for r = refs
        Tr = T(T.reference == r, :);
        S = [S; paired_row(r, "unmatched", Tr.absV1_ref, Tr.absV1_free)]; %#ok<AGROW>
        S = [S; paired_row(r, "V1 MUA equalised", Tr.absV1_ref_matched, Tr.absV1_free_matched)]; %#ok<AGROW>
    end
    disp(S);

    % --- |bias| as a function of distance from the nearest ripple ---
    free_rows = find(W.type == "nonripple");
    F = W(free_rows, :);
    dist_centers = (dist_edges(1:end-1) + dist_edges(2:end)) / 2;
    sd_centers = [-fliplr(dist_centers) dist_centers];          % pre-ripple (<0), post-ripple (>0)
    nD = numel(sd_centers);
    db_idx = discretize(abs(F.signed_dist), dist_edges);
    db_idx(F.signed_dist < 0) = numel(dist_centers) + 1 - db_idx(F.signed_dist < 0);
    db_idx(F.signed_dist >= 0) = numel(dist_centers) + db_idx(F.signed_dist >= 0);
    n_sess = numel(sessions);
    D = struct('absV1', NaN(n_sess, nD), 'V1MUA', NaN(n_sess, nD), 'absHC', NaN(n_sess, nD), ...
        'n', zeros(n_sess, nD), 'dAbsV1_matched', NaN(n_sess, nD), 'n_matched', zeros(n_sess, nD), ...
        'absV1_post', NaN(n_sess, 1), 'absV1_pre', NaN(n_sess, 1), 'V1MUA_post', NaN(n_sess, 1), ...
        'V1MUA_pre', NaN(n_sess, 1), 'centers', sd_centers);
    for si = 1:n_sess
        s = sessions(si);
        ia = find(W.session == s & W.type == "ripple_post");
        ipre = W.session == s & W.type == "ripple_pre";
        D.absV1_post(si) = mean(W.abs_V1(ia)); D.V1MUA_post(si) = mean(W.V1_MUA(ia));
        D.absV1_pre(si) = mean(W.abs_V1(ipre)); D.V1MUA_pre(si) = mean(W.V1_MUA(ipre));
        for d = 1:nD
            jb = F.session == s & db_idx == d;
            D.n(si, d) = sum(jb);
            if D.n(si, d) == 0, continue; end
            D.absV1(si, d) = mean(F.abs_V1(jb));
            D.V1MUA(si, d) = mean(F.V1_MUA(jb));
            D.absHC(si, d) = mean(F.abs_HC(jb), 'omitnan');
            ib = free_rows(jb);
            [da, dbm] = match_sets(ia, ib, W.MUA_V1_L, W.MUA_V1_R, nQ, nDraw, s * 100 + d);
            D.n_matched(si, d) = numel(da{1});
            if D.n_matched(si, d) >= min_windows_dist
                D.dAbsV1_matched(si, d) = mean(cellfun(@(i, j) mean(W.abs_V1(j)) - mean(W.abs_V1(i)), da, dbm));
            end
        end
    end
    SD = table();
    for d = 1:nD
        x = D.dAbsV1_matched(:, d); ok = isfinite(x);
        p = NaN; if sum(ok) >= 3, p = signrank(x(ok)); end
        y = D.absV1(:, d) - D.absV1_post; oky = isfinite(y);
        p_un = NaN; if sum(oky) >= 3, p_un = signrank(y(oky)); end
        SD = [SD; table(sd_centers(d), sum(D.n(:, d)), sum(oky), mean(D.absV1(oky, d)), mean(y(oky)), p_un, ...
            sum(ok), mean(x(ok)), std(x(ok)) / sqrt(sum(ok)), p, ...
            'VariableNames', {'signed_distance', 'n_windows', 'n_sessions', 'absV1', 'dAbsV1_vs_post', 'p_unmatched', ...
            'n_sessions_matched', 'dAbsV1_matched', 'SE_matched', 'p_matched'})]; %#ok<AGROW>
    end
    disp(SD);

    vars = [{'session_id', 'subject_id'}, compose('t_%+.2f', TC.t)];
    writetable(array2table([TC.sessions' TC.session_subject TC.session_mean], 'VariableNames', vars), ...
        fullfile(out_dir, [tc_name '_per_session.csv']));
    writetable(array2table([TC.sessions' TC.session_subject TC.session_mean_matched], 'VariableNames', vars), ...
        fullfile(out_dir, [tc_name '_MUA_matched_per_session.csv']));
    save(fullfile(out_dir, 'nonripple_V1_abs_bias.mat'), 'W', 'T', 'S', 'D', 'SD', 'M', 'TC', 'nDraw', 'nQ', ...
        'min_gap', 'win_len', 'min_windows', 'min_windows_dist', 'dist_edges', '-v7.3');
    writetable(T, fullfile(out_dir, 'nonripple_V1_abs_bias_per_session.csv'));
    writetable(S, fullfile(out_dir, 'nonripple_V1_abs_bias_stats.csv'));
    writetable(SD, fullfile(out_dir, 'nonripple_V1_abs_bias_distance_stats.csv'));

    plot_distributions(W, M, refs, fullfile(out_dir, 'nonripple_V1_abs_bias_distributions'));
    plot_session_means(T, S, refs, fullfile(out_dir, 'nonripple_V1_abs_bias_session_means'));
    plot_distance(D, SD, min_gap, fullfile(out_dir, 'nonripple_V1_abs_bias_distance'));
    for matched = [false true]
        plot_time_course(TC, fullfile(out_dir, tc_name), matched);
    end
    plot_MUA_matching_check(TC, [-0.35 0.35], fullfile(out_dir, [tc_name '_MUA_matching_check']));
    fprintf('Saved to %s\n', out_dir);
end

function [on, off] = session_ripples(ripples_all, s)
    % Sorted ripple onsets of both probes (bilateral copies with equal onset merged)
    % and the latest offset among ripples starting up to each onset
    on = [ripples_all(1).onset(ripples_all(1).session_count == s); ripples_all(2).onset(ripples_all(2).session_count == s)];
    off = [ripples_all(1).offset(ripples_all(1).session_count == s); ripples_all(2).offset(ripples_all(2).session_count == s)];
    [on, ~, g] = unique(on);
    off = cummax(accumarray(g, off, [], @max));
end

function C = time_course_windows(tc_len, tc_starts, V1, MUA_L, MUA_R, bin_centers, mua_centers, t0, sess, subj, ripples_all)
    % |V1 bias| of windows of length tc_len starting at tc_starts (s from onset).
    % Windows inside +-0.2 s of onset from every event; the others must be
    % ripple-free and closer to their own event's onset than to any other ripple
    % onset (bilateral copies within 0.05 s count as one)
    nW = numel(tc_starts);
    Wb = zeros(numel(bin_centers), nW); Wm = zeros(numel(mua_centers), nW);
    for w = 1:nW
        Wb(:, w) = bin_centers > tc_starts(w) & bin_centers < tc_starts(w) + tc_len;
        Wm(:, w) = mua_centers > tc_starts(w) & mua_centers < tc_starts(w) + tc_len;
    end
    V1_w = nanmean_window(V1, Wb);
    MUA_Lw = double(MUA_L) * (Wm ./ sum(Wm, 1));
    MUA_Rw = double(MUA_R) * (Wm ./ sum(Wm, 1));
    t_rel = round((tc_starts + tc_len / 2) * 1000) / 1000;
    peri_w = tc_starts >= -0.2 - 1e-9 & tc_starts + tc_len <= 0.2 + 1e-9;
    C = table();
    for s = unique(sess)'
        ev = find(sess == s);
        [on, off] = session_ripples(ripples_all, s);
        [E, K] = ndgrid(ev, 1:nW);
        ws = t0(E) + tc_starts(K); we = ws + tc_len; c = ws + tc_len / 2;
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
        C = [C; table(repmat(s, sum(keep(:)), 1), subj(E(keep)), t_rel(K(keep))', peri_w(K(keep))', ...
            abs(V1_w(lin)), MUA_Lw(lin), MUA_Rw(lin), (MUA_Lw(lin) + MUA_Rw(lin)) / 2, ...
            'VariableNames', {'session', 'subject', 't', 'peri', 'abs_V1', 'MUA_V1_L', 'MUA_V1_R', 'V1_MUA'})]; %#ok<AGROW>
    end
    C = C(isfinite(C.abs_V1), :);
end

function TC = time_course(C, sessions, win_len, nQ, nDraw, min_matched)
    % Session and animal mean |V1 bias| at each window centre: raw, and with the
    % ripple-free windows of each centre subsampled so their V1 L x V1 R MUA
    % distribution matches the session's pooled peri-ripple windows (kept whole)
    t = unique(C.t)';
    n_sess = numel(sessions); nT = numel(t);
    peri = arrayfun(@(x) any(C.peri(C.t == x)), t);
    TC = struct('win_len', win_len, 't', t, 'peri', peri, 'sessions', sessions, 'session_mean', NaN(n_sess, nT), ...
        'n_windows', zeros(n_sess, nT), 'session_subject', NaN(n_sess, 1), 'subjects', [], 'animal_mean', [], ...
        'session_mean_matched', NaN(n_sess, nT), 'n_matched', zeros(n_sess, nT), 'animal_mean_matched', [], ...
        'MUA_raw', NaN(n_sess, nT), 'MUA_matched', NaN(n_sess, nT), 'MUA_reference', NaN(n_sess, 1), ...
        'nQ', nQ, 'nDraw', nDraw, 'min_matched', min_matched, 'C', C, 'draw1', {cell(n_sess, nT)});
    for si = 1:n_sess
        rows_s = find(C.session == sessions(si));
        Cs = C(rows_s, :);
        TC.session_subject(si) = mode(Cs.subject);
        ir = find(Cs.peri);
        TC.MUA_reference(si) = mean(Cs.V1_MUA(ir));
        for j = 1:nT
            ij = find(Cs.t == t(j));
            TC.n_windows(si, j) = numel(ij);
            if isempty(ij), continue; end
            TC.session_mean(si, j) = mean(Cs.abs_V1(ij));
            TC.MUA_raw(si, j) = mean(Cs.V1_MUA(ij));
            if peri(j)
                TC.session_mean_matched(si, j) = TC.session_mean(si, j);
                TC.MUA_matched(si, j) = TC.MUA_raw(si, j);
                continue
            end
            draws = match_to_reference(ir, ij, Cs.MUA_V1_L, Cs.MUA_V1_R, nQ, nDraw, si * 1000 + j);
            TC.n_matched(si, j) = numel(draws{1});
            if TC.n_matched(si, j) >= min_matched
                TC.draw1{si, j} = rows_s(draws{1});   % rows of TC.C, for the MUA distribution check
                TC.session_mean_matched(si, j) = mean(cellfun(@(i) mean(Cs.abs_V1(i)), draws));
                TC.MUA_matched(si, j) = mean(cellfun(@(i) mean(Cs.V1_MUA(i)), draws));
            end
        end
    end
    TC.subjects = unique(TC.session_subject)';
    TC.animal_mean = NaN(numel(TC.subjects), nT); TC.animal_mean_matched = TC.animal_mean;
    for a = 1:numel(TC.subjects)
        in_a = TC.session_subject == TC.subjects(a);
        TC.animal_mean(a, :) = mean(TC.session_mean(in_a, :), 1, 'omitnan');
        TC.animal_mean_matched(a, :) = mean(TC.session_mean_matched(in_a, :), 1, 'omitnan');
    end
end

function draws = match_to_reference(ir, ib, x1, x2, nQ, nDraw, seed)
    % Subsample windows ib so their joint (x1, x2) distribution matches the
    % reference windows ir, which are kept whole: nQ x nQ cells with quantile
    % edges from the reference; per cell round(p_ref * N) windows are drawn, with
    % N as large as the scarcest cell allows
    draws = repmat({zeros(0, 1)}, nDraw, 1);
    if isempty(ir) || isempty(ib)
        return
    end
    cr = joint_bin(x1, x2, ir, ir, nQ);
    cb = joint_bin(x1, x2, ir, ib, nQ);
    cells = unique(cr)';
    p = arrayfun(@(c) mean(cr == c), cells);
    nb = arrayfun(@(c) sum(cb == c), cells);
    N = floor(min(nb ./ p));
    k = round(p * N);
    if N == 0
        return
    end
    for d = 1:nDraw
        st = RandStream('philox4x32_10', 'Seed', seed * 1000 + d);
        sb = zeros(0, 1);
        for c = 1:numel(cells)
            eb = ib(cb == cells(c));
            sb = [sb; eb(randperm(st, numel(eb), k(c)))]; %#ok<AGROW>
        end
        draws{d} = sb;
    end
end

function plot_time_course(TC, save_path, matched)
    % |V1 bias| vs window centre relative to ripple onset, three versions:
    % mean +- SE across sessions; session means; animal means (no shade).
    % matched: ripple-free points V1 MUA-matched to the peri-ripple windows
    c_peri = [231, 41, 138] / 256; c_out = [0.35 0.35 0.35];
    if matched
        f_sess = 'session_mean_matched'; f_anim = 'animal_mean_matched'; save_path = [save_path '_MUA_matched'];
        grey_str = 'grey: ripple-free, V1 MUA matched to pink windows (dashed: unmatched)';
    else
        f_sess = 'session_mean'; f_anim = 'animal_mean';
        grey_str = 'grey: ripple-free';
    end
    versions = {'sessions_SE', 'sessions', 'animals'};
    for v = 1:numel(versions)
        fig = figure('Name', ['|V1 bias| time course, ' versions{v}], 'Position', [100 100 560 440]);
        ax = axes(fig); hold(ax, 'on');
        switch versions{v}
            case 'sessions_SE'
                Y = TC.(f_sess);
                se = std(Y, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(Y), 1));
                if matched
                    plot(ax, TC.t, mean(TC.session_mean, 1, 'omitnan'), '--', 'Color', [c_out 0.6], 'LineWidth', 1);
                end
                segmented_line(ax, TC.t, TC.peri, mean(Y, 1, 'omitnan'), se, c_peri, c_out, 2, 1, true);
                ttl = sprintf('Mean \\pm SE across sessions (n = %d)', size(Y, 1));
            case {'sessions', 'animals'}
                if strcmp(versions{v}, 'sessions')
                    Y = TC.(f_sess); unit = 'session';
                else
                    Y = TC.(f_anim); unit = 'animal';
                end
                if matched
                    grey_str = 'grey: ripple-free, V1 MUA matched to pink windows';
                end
                for i = 1:size(Y, 1)
                    segmented_line(ax, TC.t, TC.peri, Y(i, :), [], c_peri, c_out, 0.75, 0.35, false);
                end
                segmented_line(ax, TC.t, TC.peri, mean(Y, 1, 'omitnan'), [], c_peri, c_out, 2.5, 1, true);
                ttl = sprintf('Thin: each %s; bold: mean of %ss (n = %d)', unit, unit, size(Y, 1));
        end
        xline(ax, 0, 'k:');
        xlim(ax, [-1 1]);
        xlabel(ax, 'Window centre relative to ripple onset (s)');
        len_ms = round(1000 * TC.win_len);
        ylabel(ax, sprintf('Mean |V1 bias| (%d ms window)', len_ms));
        if len_ms == 100, step_str = 'non-overlapping'; else, step_str = 'every 100 ms'; end
        title(ax, {sprintf('POST |V1 bias|, %d ms windows %s', len_ms, step_str), ttl, ...
            'pink: windows inside \pm0.2 s of onset (all ripples)', grey_str}, 'FontSize', 9, 'FontWeight', 'normal');
        set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
        save_figure(fig, [save_path '_' versions{v}]);
    end
end

function plot_MUA_matching_check(TC, t_show, save_path)
    % V1 MUA of the ripple-free windows at the centres t_show, before and after
    % matching (draw 1), against the peri-ripple reference windows, pooled across
    % sessions. Reference windows are weighted so each session contributes as many
    % windows as it does to the matched set (matching is done within session).
    c_peri = [231, 41, 138] / 256; c_out = [0.35 0.35 0.35];
    C = TC.C;
    edges = -2:0.1:4;
    vars = {'MUA_V1_L', 'MUA_V1_R', 'V1_MUA'};
    var_labels = {'V1 L MUA (z)', 'V1 R MUA (z)', 'V1 MUA, mean L and R (z)'};
    fig = figure('Name', 'MUA matching check', 'Position', [80 80 1200 330 * numel(t_show)]);
    tiledlayout(numel(t_show), numel(vars), 'TileSpacing', 'compact');
    for r = 1:numel(t_show)
        j = find(abs(TC.t - t_show(r)) < 1e-9);
        i_un = find(C.t == TC.t(j) & ismember(C.session, TC.sessions(TC.n_matched(:, j) >= TC.min_matched)));
        i_m = vertcat(TC.draw1{:, j});
        i_ref = find(C.peri & ismember(C.session, C.session(i_m)));
        w_ref = zeros(size(i_ref));
        for s = unique(C.session(i_m))'
            in_s = C.session(i_ref) == s;
            w_ref(in_s) = sum(C.session(i_m) == s) / sum(in_s);
        end
        for v = 1:numel(vars)
            ax = nexttile; hold(ax, 'on');
            x = C.(vars{v});
            stairs_prob(ax, x(i_ref), w_ref, edges, c_peri, '-', 1.5);
            stairs_prob(ax, x(i_un), ones(size(i_un)), edges, c_out, '--', 1.2);
            stairs_prob(ax, x(i_m), ones(size(i_m)), edges, c_out, '-', 1.5);
            legend(ax, {sprintf('Pink windows (\\pm0.2 s), mean %.3f', sum(w_ref .* x(i_ref)) / sum(w_ref)), ...
                sprintf('Ripple-free, unmatched (n = %d), mean %.3f', numel(i_un), mean(x(i_un))), ...
                sprintf('Ripple-free, matched (n = %d), mean %.3f', numel(i_m), mean(x(i_m)))}, ...
                'Box', 'off', 'Location', 'northeast', 'FontSize', 8);
            xlabel(ax, var_labels{v}); ylabel(ax, 'Probability');
            title(ax, sprintf('Window %.2f to %.2f s', TC.t(j) - TC.win_len / 2, TC.t(j) + TC.win_len / 2), ...
                'FontWeight', 'normal');
            set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
            xlim(ax, [-1.5 3]);
            yl = ylim(ax); ylim(ax, [0 yl(2) * 1.6]);   % room for the legend above the curves
        end
    end
    sgtitle(sprintf(['POST V1 MUA, %d ms windows: ripple-free before / after matching (draw 1) vs peri-ripple ' ...
        'reference (pooled sessions, reference weighted to matched session counts)'], round(1000 * TC.win_len)), 'FontSize', 11);
    save_figure(fig, save_path);
end

function stairs_prob(ax, x, w, edges, c, ls, lw)
    % Weighted histogram as probability, drawn as a stairs outline
    b = discretize(x, edges);
    ok = ~isnan(b);
    p = accumarray(b(ok), w(ok), [numel(edges) - 1, 1]) / sum(w(ok));
    stairs(ax, edges, [p; p(end)], ls, 'Color', c, 'LineWidth', lw);
end

function segmented_line(ax, t, peri, mu, se, c_peri, c_out, lw, alpha, markers)
    % One connected line (and shade): peri-ripple stretch pink, rest grey
    i_p = find(peri);
    segs = {1:i_p(1), i_p(1):i_p(end), i_p(end):numel(t)};
    cols = {c_out, c_peri, c_out};
    for k = 1:numel(segs)
        i = segs{k};
        if ~isempty(se)
            fill(ax, [t(i) fliplr(t(i))], [mu(i) + se(i) fliplr(mu(i) - se(i))], cols{k}, 'FaceAlpha', 0.2, ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        plot(ax, t(i), mu(i), '-', 'Color', [cols{k} alpha], 'LineWidth', lw, 'HandleVisibility', 'off');
    end
    if markers
        scatter(ax, t(~peri), mu(~peri), 22, c_out, 'filled');
        scatter(ax, t(peri), mu(peri), 22, c_peri, 'filled');
    end
end

function X = nanmean_window(P, Wb)
    % Mean of P (events x bins) over the bins of each window column, ignoring NaN
    ok = isfinite(P);
    P(~ok) = 0;
    X = (P * Wb) ./ (double(ok) * Wb);
end

function v = window_vars()
    v = {'session', 'event', 'window_start', 'type', 'window_start_abs', 'signed_dist', ...
        'V1_bias', 'HC_bias', 'MUA_V1_L', 'MUA_V1_R', 'MUA_HC'};
end

function R = window_rows(s, ev, k, ws_rel, type, ws, sd, V1_w, HC_w, MUA_L, MUA_R, MUA_HC)
    n = numel(ev);
    R = table(repmat(s, n, 1), ev, repmat(ws_rel, n, 1), repmat(type, n, 1), ws, sd, V1_w(ev, k), HC_w(ev, k), ...
        MUA_L(ev, k), MUA_R(ev, k), MUA_HC(ev, k), 'VariableNames', window_vars());
end

function [draws_a, draws_b] = match_sets(ia, ib, x1, x2, nQ, nDraw, seed)
    % Equalise the joint (x1, x2) distribution of window sets ia and ib by
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
        st = RandStream('philox4x32_10', 'Seed', seed * 100000 + d);
        sa = []; sb = [];
        for c = cells
            ea = ia(ca == c); eb = ib(cb == c);
            k = min(numel(ea), numel(eb));
            if k > 0
                sa = [sa; ea(randperm(st, numel(ea), k))]; %#ok<AGROW>
                sb = [sb; eb(randperm(st, numel(eb), k))]; %#ok<AGROW>
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

function row = paired_test(a, b)
    ok = isfinite(a) & isfinite(b);
    p = NaN; if sum(ok) >= 3, p = signrank(a(ok), b(ok)); end
    row = [sum(ok), mean(a(ok)), std(a(ok)) / sqrt(sum(ok)), mean(b(ok)), std(b(ok)) / sqrt(sum(ok)), p];
end

function row = paired_row(ref, meas, a, b)
    v = paired_test(a, b);
    row = table(ref, string(meas), v(1), v(2), v(3), v(4), v(5), v(4) / v(2), v(6), ...
        'VariableNames', {'reference', 'measure', 'n_sessions', 'absV1_ref', 'SE_ref', ...
        'absV1_free', 'SE_free', 'ratio_free_to_ref', 'p_signrank'});
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

function plot_distributions(W, M, refs, save_path)
    % Pooled |V1 bias| and V1 MUA distributions: unmatched and equalised (draw 1)
    ref_labels = struct('ripple_post', 'Ripple 0-0.2 s', 'ripple_pre', 'Ripple -0.2-0 s');
    c_ref = [231, 41, 138] / 256; c_free = [0.4 0.4 0.4];
    edges_b = 0:0.05:3; edges_m = -1.5:0.1:4;
    fig = figure('Name', '|V1 bias| distributions', 'Position', [80 80 1400 720]);
    tiledlayout(2, 4, 'TileSpacing', 'compact');
    for ri = 1:numel(refs)
        r = char(refs(ri));
        ia_all = find(W.type == r); ib_all = find(W.type == "nonripple");
        sets = {ia_all, ib_all; M.(r).ref, M.(r).free};
        set_names = ["unmatched", "V1 MUA equalised"];
        for m = 1:2
            ax = nexttile; hold on;
            histogram(ax, W.abs_V1(sets{m, 1}), edges_b, 'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
                'EdgeColor', c_ref, 'LineWidth', 1.5);
            histogram(ax, W.abs_V1(sets{m, 2}), edges_b, 'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
                'EdgeColor', c_free, 'LineWidth', 1.5);
            xline(median(W.abs_V1(sets{m, 1})), '--', 'Color', c_ref, 'HandleVisibility', 'off');
            xline(median(W.abs_V1(sets{m, 2})), '--', 'Color', c_free, 'HandleVisibility', 'off');
            legend({sprintf('%s (n = %d)', ref_labels.(r), numel(sets{m, 1})), ...
                sprintf('Ripple-free (n = %d)', numel(sets{m, 2}))}, 'Box', 'off', 'Location', 'northeast');
            title(sprintf('|V1 bias|, %s', set_names(m)), 'FontWeight', 'normal');
            xlabel('|V1 bias| (z log odds, 200 ms mean)'); ylabel('Probability');
            set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);

            ax = nexttile; hold on;
            histogram(ax, W.V1_MUA(sets{m, 1}), edges_m, 'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
                'EdgeColor', c_ref, 'LineWidth', 1.5);
            histogram(ax, W.V1_MUA(sets{m, 2}), edges_m, 'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
                'EdgeColor', c_free, 'LineWidth', 1.5);
            title(sprintf('V1 MUA (mean L, R), %s', set_names(m)), 'FontWeight', 'normal');
            xlabel('V1 MUA (z)'); ylabel('Probability');
            set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
        end
    end
    sgtitle('POST: |V1 track bias| in peri-ripple vs ripple-free 200 ms windows (pooled sessions; equalised = draw 1)', 'FontSize', 12);
    save_figure(fig, save_path);
end

function plot_session_means(T, S, refs, save_path)
    % Per-session mean |V1 bias|, paired lines peri-ripple vs ripple-free
    ref_labels = struct('ripple_post', 'Ripple 0-0.2 s', 'ripple_pre', 'Ripple -0.2-0 s');
    cols = {[231, 41, 138] / 256, [0.5 0.5 0.5]};
    fig = figure('Name', '|V1 bias| session means', 'Position', [100 100 1000 420]);
    tiledlayout(1, 4, 'TileSpacing', 'compact');
    meas = {'absV1_ref', 'absV1_free'; 'absV1_ref_matched', 'absV1_free_matched'};
    meas_names = ["unmatched", "V1 MUA equalised"];
    for ri = 1:numel(refs)
        Tr = T(T.reference == refs(ri), :);
        for m = 1:2
            ax = nexttile; hold on;
            Y = [Tr.(meas{m, 1}) Tr.(meas{m, 2})];
            Y = Y(all(isfinite(Y), 2), :);
            plot([1 2], Y', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
            for c = 1:2
                scatter(c * ones(size(Y, 1), 1), Y(:, c), 25, cols{c}, 'filled', 'MarkerFaceAlpha', 0.7);
                errorbar(c + 0.25 * (2 * c - 3), mean(Y(:, c)), std(Y(:, c)) / sqrt(size(Y, 1)), 'o', 'Color', 'k', ...
                    'MarkerFaceColor', cols{c}, 'LineWidth', 1.5, 'CapSize', 6);
            end
            st = S(S.reference == refs(ri) & S.measure == meas_names(m), :);
            title({ref_labels.(char(refs(ri))), char(meas_names(m)), sprintf('ratio %.2f, p = %.3g, n = %d', ...
                st.ratio_free_to_ref, st.p_signrank, st.n_sessions)}, 'FontSize', 9, 'FontWeight', 'normal');
            set(ax, 'XTick', [1 2], 'XTickLabel', {'Ripple', 'Ripple-free'}, 'TickDir', 'out', 'Box', 'off', ...
                'Color', 'none', 'FontSize', 10);
            xlim([0.4 2.6]);
            ylabel('Mean |V1 bias| (session)');
        end
    end
    sgtitle('POST: per-session mean |V1 bias|, peri-ripple vs ripple-free (paired signed-rank)', 'FontSize', 11);
    save_figure(fig, save_path);
end

function plot_distance(D, SD, min_gap, save_path)
    % |V1 bias| and V1 MUA vs distance from the nearest ripple, mean +- SE across sessions
    x = D.centers;
    c = [0.3 0.3 0.3]; c_r = [231, 41, 138] / 256;
    fig = figure('Name', '|V1 bias| vs distance from ripple', 'Position', [100 100 1350 380]);
    tiledlayout(1, 4, 'TileSpacing', 'compact');
    panels = {D.absV1, 'Mean |V1 bias|'; D.V1MUA, 'V1 MUA (z, mean L, R)'; D.absHC, 'Mean |HPC bias|'};
    for p = 1:size(panels, 1)
        ax = nexttile; hold on;
        Y = panels{p, 1};
        plot_curve(ax, x, Y, c);
        if p == 1
            errorbar(ax, [-0.05 0.05], [mean(D.absV1_pre) mean(D.absV1_post)], ...
                [std(D.absV1_pre) std(D.absV1_post)] / sqrt(numel(D.absV1_post)), 'o', 'Color', c_r, ...
                'MarkerFaceColor', c_r, 'LineWidth', 1.5);
        elseif p == 2
            errorbar(ax, [-0.05 0.05], [mean(D.V1MUA_pre) mean(D.V1MUA_post)], ...
                [std(D.V1MUA_pre) std(D.V1MUA_post)] / sqrt(numel(D.V1MUA_post)), 'o', 'Color', c_r, ...
                'MarkerFaceColor', c_r, 'LineWidth', 1.5);
        end
        decorate(ax, min_gap);
        ylabel(panels{p, 2});
        title(panels{p, 2}, 'FontWeight', 'normal');
    end
    ax = nexttile; hold on;
    plot_curve(ax, x, D.dAbsV1_matched, c);
    yline(ax, 0, 'k--');
    sig = SD.p_matched < 0.05;
    m = mean(D.dAbsV1_matched, 1, 'omitnan');
    scatter(ax, x(sig), m(sig), 40, 'k', '*');
    decorate(ax, min_gap);
    ylabel('\Delta|V1 bias| (ripple-free - ripple 0-0.2 s)');
    title({'V1 MUA equalised per distance bin', '* signed-rank p < 0.05'}, 'FontWeight', 'normal');
    sgtitle(sprintf(['POST: ripple-free 200 ms windows vs distance from nearest ripple ' ...
        '(<0 before ripple, >0 after; pink = peri-ripple windows); mean \\pm SE across sessions']), 'FontSize', 11);
    save_figure(fig, save_path);
end

function plot_curve(ax, x, Y, c)
    n = sum(isfinite(Y), 1);
    mu = mean(Y, 1, 'omitnan'); se = std(Y, 0, 1, 'omitnan') ./ sqrt(n);
    for side = [-1 1]
        k = sign(x) == side & n >= 3;
        fill(ax, [x(k) fliplr(x(k))], [mu(k) + se(k) fliplr(mu(k) - se(k))], c, 'FaceAlpha', 0.2, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(ax, x(k), mu(k), '-o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4, 'LineWidth', 1.5);
    end
end

function decorate(ax, min_gap)
    xline(ax, 0, 'k:');
    set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
    xlabel(ax, 'Distance from ripple (s)');
    text(ax, 0.02, 0.02, sprintf('min gap %.1f s', min_gap), 'Units', 'normalized', 'FontSize', 8, 'Color', [0.4 0.4 0.4]);
end
