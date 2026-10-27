%% main_04_MUA_equalised_event_control: firing-rate control (R1.2, R4.2/4/6), POST
% Within each session, Low and High (Peak and Trough) events are subsampled so their
% joint MUA distribution is equal (MUA_match_prepare / MUA_match_draw), then the
% reactivation coherence is recomputed. Two versions (MUA_equalised_variant), each
% in its own folder under nature_neuroscience_revision_analysis:
%   'dominant'         : HC MUA x dominant-hemisphere V1 MUA
%                        -> MUA_equalised_event_control
%   'with_nondominant' : HC MUA x dominant x non-dominant-hemisphere V1 MUA
%                        -> MUA_equalised_event_control_with_nondominant
%
% Code: 04_MUA_equalised_event_control\ (helpers in shared\)
% Inputs (main_00): ripple_info.mat, KDE_reactivation_ripples_PSTH_POST_unique.mat,
%   ripple_MUA_HPC.mat, ripple_MUA_V1_L.mat, ripple_MUA_V1_R.mat
% Outputs, in each version's folder:
%   geomean\  per-session geomean coherence, unmatched vs matched, signed-rank
%             across sessions; MUA distributions before / after matching
%   AUC\      pooled AUC, static (MUA_equalised_AUC_static.*) and temporal
%             (MUA_equalised_AUC_temporal_*.*), low vs high bar and temporal plots
%   MUA_PSTH\ V1 MUA PSTH Low vs High, all events (dominant hemisphere in
%             MUA_equalised_event_control, non-dominant in ..._with_nondominant)
rev_dir = addon_setup();
run_analysis = true;   % false: only redraw the low vs high plots from saved results
variants = {'dominant', 'with_nondominant'};

%% 1. Geomean: per-session geomean coherence, MUA-equalised -> geomean
if run_analysis
    for iv = 1:numel(variants)
        MUA_equalised_geomean(variants{iv});
    end
end

%% 2. AUC: pooled, MUA-equalised (static + temporal) -> AUC
% MUA_equalised_AUC('dominant', {'SO_trough_uni_bi'}) reruns one analysis and keeps the others
if run_analysis
    for iv = 1:numel(variants)
        MUA_equalised_AUC(variants{iv});
    end
end

%% 3. AUC low vs high bar and temporal plots -> AUC\MUA_equalised[_temporal]_AUC_low_vs_high
for iv = 1:numel(variants)
    v = MUA_equalised_variant(variants{iv});
    plot_AUC_low_vs_high_controls(v.key);
    plot_temporal_AUC_low_vs_high(v.key);
end

%% 4. MUA: V1 MUA PSTH, Low vs High (Peak vs Trough), all events -> MUA_PSTH
if run_analysis
    MUA_PSTH_low_vs_high('dominant');
    MUA_PSTH_low_vs_high('non_dominant');
end
