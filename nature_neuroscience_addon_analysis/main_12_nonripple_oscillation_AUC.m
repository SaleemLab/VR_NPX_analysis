%% main_12_nonripple_oscillation_AUC: oscillation-gated V1-HPC coherence outside ripples (R1), POST
% Reviewer 1 asked whether reactivation occurs in the absence of ripples at the
% 'good' SO phase with coincident spindles. V1->HPC reactivation coherence AUC,
% spindle Low vs High and SO Peak vs Trough (dominant hemisphere from the V1 bias
% sign), in 200 ms non-overlapping windows from -1 to 1 s around ripple onset:
% windows within +-0.2 s of onset vs ripple-free windows further away.
%
% Inputs (main_00): ripple_info.mat, KDE_reactivation_ripples_PSTH_POST_unique.mat,
%   ripples_all_POST.mat; analysis folder slow_waves_all_POST.mat; per-session
%   extracted_LFP.mat (step 1 only)
% Outputs (nature_neuroscience_revision_analysis):
%   SO_spindle_traces_POST.mat (step 1: continuous SO phase / spindle amplitude,
%   10 ms, reference V1 shank per hemisphere)
%   nonripple_oscillation_AUC\ nonripple_oscillation_AUC.mat / .csv,
%   nonripple_oscillation_AUC_time (.fig/.png/.pdf)
rev_dir = addon_setup();
step.traces = false;   % 1: re-extract SO / spindle traces from LFP (~30 min, network drive)

%% 1. Continuous SO phase and spindle amplitude -> SO_spindle_traces_POST.mat
if step.traces
    extract_SO_spindle_traces_POST();
end

%% 2. AUC per 200 ms window, Low vs High spindle and SO Peak vs Trough -> nonripple_oscillation_AUC
nonripple_oscillation_AUC_POST();
