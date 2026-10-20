%% main_02_session_level_stats: per-session and per-animal statistics, PRE vs POST
% Reactivation coherence computed per session (or animal) instead of pooled events:
% condition comparisons within each sleep, PRE vs POST, and each session against
% its own shuffle; plus the |log-odds bias| of each condition.
%
% Code: 02_session_level_stats\ (helpers in shared\)
% Inputs (main_00): ripple_info.mat, ripple_info_PRE_unique.mat,
%   KDE_reactivation_ripples_PSTH_PRE_unique.mat, KDE_reactivation_ripples_PSTH_POST_unique.mat
% Outputs (nature_neuroscience_revision_analysis):
%   PRE_vs_POST, session_AUC_vs_shuffle, session_bootAUC_PRE_vs_POST(_permutation,
%   _ranksum, _difference), session_bootGeomean, session_bootGeomean_PRE_vs_POST(
%   _permutation, _ranksum, _difference), POST_animal, abs_bias,
%   POST_original_bootAUC_per_session, POST_original_bootAUC_per_animal
rev_dir = addon_setup();
run_analysis = true;   % false: only rerun the plot_session_* tests / figures on saved results
test_methods = {'signrank', 'permutation', 'ranksum', 'difference'};

%% 1. Per-session AUC and geomean coherence (no bootstrap), PRE vs POST -> PRE_vs_POST
if run_analysis
    compare_session_AUC_PRE_POST();
    compare_session_geomean_coherence();
end

%% 2. Per-session AUC vs the session's own shuffle -> session_AUC_vs_shuffle
if run_analysis
    compare_session_AUC_vs_shuffle();       % actual AUC vs shuffle 95th prctile
    compare_session_bootAUC_vs_shuffleCI(); % bootstrapped AUC vs shuffle 95% CI upper bound
end
plot_session_AUC_vs_shuffle_mean();         % actual AUC vs shuffle mean

%% 3. Per-session bootstrapped AUC: Low vs High and PRE vs POST -> session_bootAUC_PRE_vs_POST*
for k = 1:numel(test_methods)
    plot_session_bootAUC_PRE_vs_POST(test_methods{k});
end

%% 4. Per-session bootstrapped geomean coherence -> session_bootGeomean, session_bootGeomean_PRE_vs_POST*
if run_analysis
    compare_session_bootGeomean();
end
for k = 1:numel(test_methods)
    plot_session_bootGeomean_PRE_vs_POST(test_methods{k});
end

%% 5. Per-animal AUC, POST -> POST_animal
if run_analysis
    compare_animal_AUC_POST();
end

%% 6. Mean |log-odds bias| of HPC and V1 per condition -> abs_bias
if run_analysis
    plot_abs_bias_conditions();
end

%% 7. Original POST decoder: per-session and per-animal bootstrapped AUC (95% CI) and
% shuffle, static and temporal (tables), plus the summary of temporal traces
% (all ripples, ripple power, spindle, SO phase, dominant SO trough uni/bilateral)
%   -> POST_original_bootAUC_per_session, POST_original_bootAUC_per_animal
if run_analysis
    boot_AUC_per_unit_POST_original('session');
    boot_AUC_per_unit_POST_original('animal');
end
