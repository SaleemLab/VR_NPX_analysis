%% main_01_reactivation_AUC: pooled reactivation coherence AUC, PRE and POST sleep
% HPC->V1 (all ripples, ripple power) and V1->HPC (spindle, SO phase, dominant SO
% trough with non-dominant peak vs trough) AUC for the pre (-0.2-0 s) and post
% (0-0.2 s) V1 windows plus the temporal AUC, per sleep (event bootstrap), and the
% same for the original POST decoder with a hierarchical (session -> event)
% bootstrap (R1.6).
%
% Code: 01_reactivation_AUC\ (helpers in shared\)
% Inputs (main_00): ripple_info.mat, ripple_info_PRE_unique.mat,
%   KDE_reactivation_ripples_PSTH_PRE_unique.mat, KDE_reactivation_ripples_PSTH_POST_unique.mat
% Outputs (nature_neuroscience_revision_analysis):
%   PRE, PRE_equal_n, POST_matched, POST_matched_equal_n, POST   (AUC_*.mat, figures,
%       AUC_low_vs_high.fig/png/pdf)
%   POST_original_AUC_hierarchical
%   awake   (as the sleep folders, plus KDE_temporal_AUC_ripple_percentile_awake; inputs
%       ripple_info_awake_unique.mat, KDE_reactivation_ripples_PSTH_awake_unique.mat)
rev_dir = addon_setup();
run_analysis = true;   % false: only redraw the low vs high bar plots from saved results

%% 1. AUC per sleep (folder = sleep, '_equal_n' when PRE and POST events per session are equalised)
% 'POST_matched': POST events of the 12 sessions with PRE sleep; 'POST': all 22 sessions
runs = struct('sleep',   {'PRE', 'PRE', 'POST_matched', 'POST_matched', 'POST'}, ...
              'equal_n', {false, true,  false,          true,           false});
for k = 1:numel(runs)
    if run_analysis
        reactivation_AUC_coherence_PRE(runs(k).sleep, runs(k).equal_n);
    end
    folder = runs(k).sleep;
    if runs(k).equal_n, folder = [folder '_equal_n']; end
    plot_AUC_low_vs_high(folder);
end

%% 2. Original POST decoder, hierarchical (session -> event) bootstrap
if run_analysis
    hierarchical_AUC_POST_original();
end

%% 3. Awake (non-SWS POST) ripples: same sections, plus temporal ripple power AUC
if run_analysis
    reactivation_AUC_coherence_PRE('awake');
end
plot_AUC_low_vs_high('awake');
