%% main_04_MUA_equalised_event_control: firing-rate control (R1.2, R4.2/4/6), POST
% Within each session, Low and High (Peak and Trough) events are subsampled so their
% joint distribution of HC MUA and dominant-hemisphere V1 MUA is equal, then the
% reactivation coherence is recomputed: per-session geomean coherence, and pooled
% AUC (ripple power, spindle, SO phase, dominant SO trough with non-dominant peak
% vs trough; static and temporal).
%
% Code: 04_MUA_equalised_event_control\ (helpers in shared\)
% Inputs (main_00): ripple_info.mat, KDE_reactivation_ripples_PSTH_POST_unique.mat,
%   ripple_MUA_HPC.mat, ripple_MUA_V1_L.mat, ripple_MUA_V1_R.mat
% Outputs (nature_neuroscience_revision_analysis\MUA_equalised_event_control):
%   geomean\  MUA_equalised_event_control.mat, session coherence / MUA distribution figures
%   AUC\      MUA_equalised_AUC_static.*, MUA_equalised_AUC_temporal_*.*,
%             MUA_equalised_AUC_low_vs_high.fig/png/pdf
%   MUA_PSTH\ dominant V1 MUA PSTH Low vs High (per session, mean +- SE), session
%             window means (-0.2-0 s, 0-0.2 s) and signed-rank tests
rev_dir = addon_setup();
run_analysis = true;   % false: only redraw the low vs high bar plot from saved results

%% 1. Per-session geomean coherence, MUA-equalised -> geomean
if run_analysis
    MUA_equalised_event_control();
end

%% 2. Pooled AUC, MUA-equalised (static + temporal) -> AUC
% MUA_equalised_AUC({'SO_trough_uni_bi'}) reruns one analysis and keeps the others
if run_analysis
    MUA_equalised_AUC();
end

%% 3. Low vs high bar plot -> AUC\MUA_equalised_AUC_low_vs_high
plot_AUC_low_vs_high_controls('MUA_equalised');

%% 4. Dominant V1 MUA PSTH, Low vs High (Peak vs Trough) -> MUA_PSTH
if run_analysis
    MUA_PSTH_low_vs_high();
end
