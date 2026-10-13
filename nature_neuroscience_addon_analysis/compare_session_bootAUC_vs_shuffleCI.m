%% Per-session bootstrapped reactivation coherence AUC vs the session's shuffle 95% CI
% For each session, sleep (PRE, POST), V1 window (pre -0.2-0 s, post 0-0.2 s)
% and condition, events of the session are bootstrapped 1000 times
% (reactivation_coherence_auc): the bootstrapped mean AUC is compared with the
% upper bound (97.5th prctile) of the session's shuffled-AUC 95% CI (selector
% and its Track 1/2 conditions un-resampled, measured bias resampled).
% Paired Wilcoxon signed-rank across sessions: bootstrapped mean vs shuffle CI upper.
% Non-duplicated ripple events (make_unique_ripple_data).
% Conditions:
%   all_ripples  : HPC->V1, all events
%   ripple_power : HPC->V1, within-session ripple power Q1 (Low) / Q4 (High)
%   spindle      : V1->HPC, hemisphere-matched within-session spindle percentile Q1 / Q4
%   SO_phase     : V1->HPC, hemisphere-matched SO peak / trough
%   SO_trough_bilateral : V1->HPC, hemisphere-matched SO trough with the other
%                  hemisphere in peak (Unilateral Trough) or trough (Bilateral Trough)
% analyses_to_run (optional): cell array of analysis names; only these are
%   (re)computed and merged into the existing saved T (other rows kept as is).
% Output: nature_neuroscience_revision_analysis\session_AUC_vs_shuffle\session_bootAUC_vs_shuffleCI_*
function compare_session_bootAUC_vs_shuffleCI(analyses_to_run)
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
    nBoot = 1000;

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

    analyses = {'all_ripples', 'ripple_power', 'spindle', 'SO_phase', 'SO_trough_bilateral'};
    directions = {'HPC->V1', 'HPC->V1', 'V1->HPC', 'V1->HPC', 'V1->HPC'};

    % --- Analyses to compute; keep previously saved rows of the others ---
    T = table();
    if nargin >= 1 && ~isempty(analyses_to_run)
        keep = ismember(analyses, analyses_to_run);
        analyses = analyses(keep); directions = directions(keep);
        saved = load(fullfile(out_dir, 'session_bootAUC_vs_shuffleCI.mat'), 'T');
        T = saved.T(~ismember(saved.T.analysis, string(analyses)), :);
    end

    % --- Bootstrapped AUC and shuffle CI per session ---
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
                            r = reactivation_coherence_auc(mean_HPC(in_s), mean_V1(in_s), nBoot, c1, c2);
                        else
                            r = reactivation_coherence_auc(mean_V1(in_s), mean_HPC(in_s), nBoot, c1, c2);
                        end
                        T = [T; table(sessions(k), string(D(iS).name), string(v1_windows(w).name), string(analyses{iA}), ...
                            string(C(iCond).label), string(directions{iA}), r.AUC_mean, r.AUC_CI(1), r.AUC_CI(2), ...
                            r.AUC_mean_shuffled, r.AUC_CI_shuffled(1), r.AUC_CI_shuffled(2), sum(in_s), ...
                            'VariableNames', {'session_id', 'sleep', 'V1_window', 'analysis', 'condition', 'direction', ...
                            'AUC_boot_mean', 'AUC_CI_lo', 'AUC_CI_hi', 'shuffle_mean', 'shuffle_CI_lo', 'shuffle_CI_hi', 'n_events'})]; %#ok<AGROW>
                    end
                    fprintf('%s %s %s %s done\n', D(iS).name, v1_windows(w).name, analyses{iA}, C(iCond).label);
                end
            end
        end
    end

    % --- Signed-rank: bootstrapped mean vs shuffle CI upper across sessions (two- and one-sided) ---
    writetable(T, fullfile(out_dir, 'session_bootAUC_vs_shuffleCI.csv'));
    S = session_vs_shuffle_stats_plot(T, 'AUC_boot_mean', 'shuffle_CI_hi', 'Shuffle CI upper', out_dir, 'session_bootAUC_vs_shuffleCI_');
    save(fullfile(out_dir, 'session_bootAUC_vs_shuffleCI.mat'), 'T', 'S');
    fprintf('Saved to %s\n', out_dir);
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
        case 'SO_trough_bilateral'
            H = hemisphere_conditions('SO_trough_bilateral', info);
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
