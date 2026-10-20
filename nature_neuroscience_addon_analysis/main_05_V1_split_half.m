%% main_05_V1_split_half: split-half V1 control (R1.3), POST
% V1 neurons of each session are split into two matched halves A and B, each with
% its own PLS-KDE decoder. One half sets the dominant hemisphere (whose spindle / SO
% defines the condition), the other selects Track 1 / Track 2 events for the
% V1->HPC coherence AUC (A->B, B->A; Full->Full = all neurons, same code).
% AUC with event bootstrap and with hierarchical (session -> event) bootstrap (R1.6).
%
% Code: 05_V1_split_half\ (helpers in shared\: PLS_KDE_train, PLS_KDE_decode,
%   cap_log_odds_PSTH, plot_AUC_low_vs_high_controls)
% Inputs (main_00): ripple_info.mat, ripples_all_POST.mat, KDE_reactivation_ripples_PSTH.mat,
%   KDE_reactivation_ripples_PSTH_POST_unique.mat; per-session spikes and
%   extracted_ripple_events.mat of the original pipeline
% Outputs:
%   split_dir\POST_session_XX.mat (per-session decoders; local disk, large)
%   nature_neuroscience_revision_analysis\
%     KDE_reactivation_ripples_PSTH_POST_unique_V1_split_half.mat
%     split_half_V1_AUC_POST\               (split_half_AUC_static_*, _temporal_*, _low_vs_high_*)
%     split_half_V1_AUC_POST_hierarchical\  (hier_AUC_static_*, _temporal_*, _low_vs_high_*)
rev_dir = addon_setup();
run_decoding = false;  % true: retrain the split-half decoders (slow; skips sessions already saved)
run_analysis = true;   % false: only redraw the low vs high bar plots from saved results
overwrite    = false;  % true: redo sessions that already have a decoder output
split_dir = 'C:\Data\PLS_split_half';   % per-session decoder outputs

%% 1. Split-half decoding per session (halves A, B and all neurons)
if run_decoding
    run_V1_split_half_decoding_POST('include_full', true, 'out_dir', split_dir, 'overwrite', overwrite);
end

%% 2. Merge sessions into the ripple PSTH file (non-duplicated events)
if run_decoding
    assemble_V1_split_half_PSTH_POST(split_dir);
end

%% 3. Split-half AUC, event bootstrap -> split_half_V1_AUC_POST
if run_analysis
    split_half_V1_AUC_POST();
end
plot_AUC_low_vs_high_controls('split_half');

%% 4. Split-half AUC, hierarchical (session -> event) bootstrap -> split_half_V1_AUC_POST_hierarchical
if run_analysis
    split_half_V1_AUC_POST_hierarchical();
end
plot_AUC_low_vs_high_controls('split_half_hierarchical');
