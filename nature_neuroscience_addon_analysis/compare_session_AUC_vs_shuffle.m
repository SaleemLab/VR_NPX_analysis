%% Per-session reactivation coherence AUC vs the session's own shuffle (95th prctile)
% For each session, sleep (PRE, POST), V1 window (pre -0.2-0 s, post 0-0.2 s)
% and condition, the actual AUC (no bootstrap) is compared with the 95th
% percentile of 1000 within-session shuffles. Shuffle: the measured-region bias
% is permuted across events while the selector and its Track 1/2 conditions stay
% fixed (same null as the shuffled control in reactivation_coherence_auc).
% Paired Wilcoxon signed-rank across sessions: actual AUC vs shuffle 95th.
% Non-duplicated ripple events (make_unique_ripple_data).
% Conditions:
%   all_ripples  : HPC->V1, all events
%   ripple_power : HPC->V1, within-session ripple power Q1 (Low) / Q4 (High)
%   spindle      : V1->HPC, hemisphere-matched within-session spindle percentile Q1 / Q4
%   SO_phase     : V1->HPC, hemisphere-matched SO peak / trough
% Output: nature_neuroscience_revision_analysis\session_AUC_vs_shuffle
function compare_session_AUC_vs_shuffle()
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    out_dir = fullfile(rev_dir, 'session_AUC_vs_shuffle');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    nPerm = 1000;

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
    directions = {'HPC->V1', 'HPC->V1', 'V1->HPC', 'V1->HPC'};
    titles = {'All Ripples (HPC->V1)', 'Ripple Power, within session (HPC->V1)', ...
              'Spindle Percentile, within session (V1->HPC)', 'SO Phase (V1->HPC)'};

    % --- Actual AUC and shuffle 95th per session ---
    T = table();
    for iS = 1:2
        sess = D(iS).info.session_id(:);
        sessions = unique(sess);
        mean_HPC = mean(D(iS).HC(bins_HPC, :), 'omitnan');
        for w = 1:numel(v1_windows)
            mean_V1 = mean(D(iS).V1(v1_windows(w).bins, :), 'omitnan');
            for iA = 1:numel(analyses)
                C = make_conditions(analyses{iA}, D(iS).info);
                for iCond = 1:numel(C)
                    for k = 1:numel(sessions)
                        in_s = sess == sessions(k) & C(iCond).mask;
                        c1 = C(iCond).c1(in_s); c2 = C(iCond).c2(in_s);
                        if strcmp(directions{iA}, 'HPC->V1')
                            sel = mean_HPC(in_s); meas = mean_V1(in_s);
                        else
                            sel = mean_V1(in_s); meas = mean_HPC(in_s);
                        end
                        r = reactivation_coherence_auc(sel, meas, 0, c1, c2);
                        shuf = shuffle_auc(sel, meas, c1, c2, nPerm);
                        T = [T; table(sessions(k), string(D(iS).name), string(v1_windows(w).name), string(analyses{iA}), ...
                            string(C(iCond).label), string(directions{iA}), r.AUC_mean, prctile(shuf, 95), mean(shuf, 'omitnan'), ...
                            mean(shuf >= r.AUC_mean), sum(in_s), ...
                            'VariableNames', {'session_id', 'sleep', 'V1_window', 'analysis', 'condition', 'direction', ...
                            'AUC', 'shuffle_95th', 'shuffle_mean', 'p_perm', 'n_events'})]; %#ok<AGROW>
                    end
                    fprintf('%s %s %s %s done\n', D(iS).name, v1_windows(w).name, analyses{iA}, C(iCond).label);
                end
            end
        end
    end

    % --- Signed-rank: actual vs shuffle 95th across sessions (two- and one-sided) ---
    writetable(T, fullfile(out_dir, 'session_AUC_vs_shuffle.csv'));
    S = session_vs_shuffle_stats_plot(T, 'AUC', 'shuffle_95th', 'Shuffle 95th', out_dir, 'session_AUC_vs_shuffle_');
    save(fullfile(out_dir, 'session_AUC_vs_shuffle.mat'), 'T', 'S');
    fprintf('Saved to %s\n', out_dir);
end

function shuf = shuffle_auc(sel, meas, c1, c2, nPerm)
    % AUC with meas permuted across events (selector and conditions fixed)
    sel = sel(:)'; meas = meas(:)'; c1 = logical(c1(:)'); c2 = logical(c2(:)');
    n = numel(sel);
    shuf = NaN(nPerm, 1);
    if n < 10
        return
    end
    th = prctile(abs(sel), 0:10:100); th = th(1:end-1);
    parfor iPerm = 1:nPerm
        s = RandStream('philox4x32_10', 'Seed', iPerm);
        mp = meas(randperm(s, n));
        d = NaN(1, numel(th));
        for j = 1:numel(th)
            t1 = sel >= th(j) & c1;
            t2 = sel <= -th(j) & c2;
            if any(t1) && any(t2)
                d(j) = mean(mp(t1), 'omitnan') - mean(mp(t2), 'omitnan');
            end
        end
        shuf(iPerm) = trapz(th, d) / (max(th) - min(th));
    end
end

function C = make_conditions(name, info)
    % mask: events included; c1/c2: Track 1 / Track 2 conditions (V1->HPC only)
    nEv = numel(info.session_id);
    all_true = true(nEv, 1);
    switch name
        case 'all_ripples'
            C = struct('label', 'All', 'mask', all_true, 'c1', all_true, 'c2', all_true);
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
