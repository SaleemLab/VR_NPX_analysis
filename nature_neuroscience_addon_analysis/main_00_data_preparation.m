%% main_00_data_preparation: shared input files of all revision analyses
% Builds the files in nature_neuroscience_revision_analysis that the other
% main_*.m scripts read. These steps were run once; most overwrite the shared
% inputs (some modify files in place), so every step is off by default. Turn on
% only the steps to rebuild, in this order.
%
% Code: 00_data_preparation\ (helpers in shared\)
%
% Files from the original (pre-revision) pipeline, not regenerated here:
%   ripple_info.mat (non-duplicated POST ripple events, event_id), ripples_all_POST.mat,
%   ripple_MUA_HPC.mat, ripple_MUA_V1_L.mat, ripple_MUA_V1_R.mat,
%   behavioural_state_merged_all_PRE/POST.mat, slow_waves_all_PRE/POST.mat (analysis
%   folder), V1-HPC sleep reactivation\KDE_reactivation_ripples_PSTH.mat, and the
%   per-session session_clusters_*.mat / KDE_reactivation*_ripples.mat.
%
% Outputs (nature_neuroscience_revision_analysis):
%   session_clusters_all_RUN1/PRE/POST.mat, ripples_all_PRE.mat, ripple_info_PRE.mat,
%   KDE_reactivation_ripples_PSTH_PRE.mat, ripple_info_awake.mat,
%   KDE_reactivation_ripples_PSTH_awake.mat, KDE_reactivation_ripples_PSTH.mat (capped),
%   KDE_reactivation_ripples_PSTH_POST_unique.mat, KDE_reactivation_ripples_PSTH_PRE_unique.mat,
%   ripple_info_PRE_unique.mat, ripple_info_awake_unique.mat,
%   KDE_reactivation_ripples_PSTH_awake_unique.mat, spindle onset / peak fields of ripple_info.mat,
%   RUN1_decoding_input.mat
rev_dir = addon_setup();
step.session_clusters   = false;  % 1
step.PRE_ripple_info    = false;  % 2
step.PRE_PSTH           = false;  % 3
step.awake_and_copy     = false;  % 4
step.cap_POST_PSTH      = false;  % 5
step.unique_events      = false;  % 6
step.unique_awake_events = false; % 6b
step.spindle_onset_peak = false;  % 7
step.RUN1_decoding      = false;  % 8

%% 1. Spatial-cell summaries of RUN1 and PRE sleep (saved in the analysis folder)
if step.session_clusters
    extract_session_clusters_all_RUN1();
    extract_session_clusters_all_PRE();
end

%% 2. PRE ripple info: SO phase / spindle per ripple, session and subject ids as POST
if step.PRE_ripple_info
    extract_pre_data();             % ripple_info_PRE.mat (+ copies ripples_all_PRE.mat)
    fix_pre_ids();                  % map PRE session / subject ids onto the 22 POST sessions
    add_spindle_percentile_PRE();   % within-session spindle percentile (was all NaN)
end

%% 3. PRE KDE ripple PSTH (log odds capped as POST)
if step.PRE_PSTH
    extract_KDE_reactivation_ripples_PSTH_PRE();
end

%% 4. Awake ripples, and copy of the POST inputs into the revision folder
if step.awake_and_copy
    extract_awake_data();
end

%% 5. Cap +/-Inf POST log odds (revision-folder copy only)
if step.cap_POST_PSTH
    cap_KDE_reactivation_ripples_PSTH_POST();
end

%% 6. Non-duplicated (bilateral-merged) PRE and POST ripple events
if step.unique_events
    make_unique_ripple_data();
end
if step.unique_awake_events
    make_unique_awake_ripple_data();   % also caps +/-Inf awake log odds
end

%% 7. Spindle-band amplitude at ripple onset and peak (POST ripple_info, backup kept)
if step.spindle_onset_peak
    add_spindle_onset_peak_to_ripple_info();
end

%% 8. RUN1 spikes, cells and behaviour for the RUN / replay decoders (main_07, main_08)
if step.RUN1_decoding
    prepare_RUN1_decoding_input();
end
