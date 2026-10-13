%% Per-session bootstrapped geomean reactivation coherence: Low vs High and PRE vs POST
% Event coherence = sign(HPC*V1) * sqrt(|HPC*V1|), HPC = mean z log odds 0-0.1 s,
% V1 = mean z log odds over the pre (-0.2-0 s) or post (0-0.2 s) window.
% For each session and condition, events are bootstrapped 1000 times (seeded) and
% the mean coherence is taken; the bootstrapped mean and 95% CI are saved.
% Conditions (as compare_session_geomean_coherence):
%   all_ripples  : all events
%   ripple_power : within-session ripple power Q1 (Low) vs Q4 (High)
%   spindle      : spindle percentile <= 25 (Low) vs > 75 (High)
%   SO_phase     : SO peak [-pi/2, pi/2] vs trough
%   SO_trough_bilateral : matched-hemisphere SO trough with the other hemisphere
%                  in peak (Unilateral Trough) vs also in trough (Bilateral Trough)
% Spindle and SO use the hemisphere matched to the event's V1 bias sign.
% Paired signed-rank tests and plots via session_paired_stats_plot (two- and
% one-sided p saved; plots show one-sided). Other PRE vs POST tests:
% plot_session_bootGeomean_PRE_vs_POST.
% Non-duplicated ripple events. Output: nature_neuroscience_revision_analysis\session_bootGeomean
function compare_session_bootGeomean()
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    out_dir = fullfile(rev_dir, 'session_bootGeomean');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    nBoot = 1000;

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

    timebin = 0.01;
    bin_edges = -1:timebin:1;
    bin_centers = bin_edges(1:end-1) + timebin/2;
    bins_HPC = bin_centers > 0 & bin_centers < 0.1;
    v1_windows = struct('name', {'pre', 'post'}, ...
        'bins', {bin_centers > -0.2 & bin_centers < 0, bin_centers > 0 & bin_centers < 0.2});
    analyses = {'all_ripples', 'ripple_power', 'spindle', 'SO_phase', 'SO_trough_bilateral'};

    T = table();
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
                        v = v(isfinite(v));
                        [bm, ci] = boot_mean(v, nBoot);
                        T = [T; table(s, string(D(iS).name), string(v1_windows(w).name), string(analyses{iA}), ...
                            string(labels{iCond}), bm, ci(1), ci(2), mean(v), numel(v), ...
                            'VariableNames', {'session_id', 'sleep', 'V1_window', 'analysis', 'condition', ...
                            'geomean_boot_mean', 'geomean_CI_lo', 'geomean_CI_hi', 'geomean_mean', 'n_events'})]; %#ok<AGROW>
                    end
                end
            end
        end
    end
    save(fullfile(out_dir, 'session_bootGeomean.mat'), 'T');
    writetable(T, fullfile(out_dir, 'session_bootGeomean.csv'));

    S = session_paired_stats_plot(T, 'geomean_boot_mean', 'Geomean coherence (bootstrapped session mean)', ...
        out_dir, 'session_bootGeomean_');
    save(fullfile(out_dir, 'session_bootGeomean.mat'), 'T', 'S');
    fprintf('Saved to %s\n', out_dir);
end

function [bm, ci] = boot_mean(v, nBoot)
    % Bootstrapped mean and 95% CI of v (seeded resampling)
    bm = NaN; ci = [NaN NaN];
    n = numel(v);
    if n < 10
        return
    end
    m = NaN(nBoot, 1);
    for iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        m(iBoot) = mean(v(randi(s, n, n, 1)));
    end
    bm = mean(m);
    ci = prctile(m, [2.5 97.5]);
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
        case 'SO_trough_bilateral'
            % Matched hemisphere in trough; other hemisphere in peak (unilateral) or trough (bilateral)
            F = info.SO_phase;
            f = NaN(size(V1_b)); g = NaN(size(V1_b));
            f(V1_b > 0) = F(V1_b > 0, 2); g(V1_b > 0) = F(V1_b > 0, 1);
            f(V1_b < 0) = F(V1_b < 0, 1); g(V1_b < 0) = F(V1_b < 0, 2);
            is_trough = @(x) (x >= -pi & x < -pi/2) | (x > pi/2 & x <= pi);
            is_peak = @(x) x >= -pi/2 & x <= pi/2;
            masks = {is_trough(f) & is_peak(g), is_trough(f) & is_trough(g)};
            labels = {'Unilateral Trough', 'Bilateral Trough'};
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
