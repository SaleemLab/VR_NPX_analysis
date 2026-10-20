%% main_06_low_rate_decoder: decoders without high-rate neurons (interneuron control), POST
% PLS-KDE decoders of HPC and V1 built only from spatial cells with mean RUN1 rate
% < 10 Hz, applied to the POST ripple windows; the reactivation coherence AUC is
% then recomputed with these log odds (as main_01 'POST').
% The Bayesian decoders with the same cell filter are in main_07 (cell_filter
% 'rate_lt10').
%
% Code: 06_low_rate_decoder\ (helpers in shared\: PLS_KDE_train, PLS_KDE_decode,
%   cap_log_odds_PSTH, plot_AUC_low_vs_high; 01_reactivation_AUC\reactivation_AUC_coherence_PRE)
% Inputs (main_00): ripple_info.mat, ripple_info_PRE_unique.mat, ripples_all_POST.mat,
%   KDE_reactivation_ripples_PSTH.mat, KDE_reactivation_ripples_PSTH_POST_unique.mat;
%   per-session spikes and extracted_ripple_events.mat of the original pipeline
% Outputs (nature_neuroscience_revision_analysis):
%   PLS_KDE_rate_lt10\POST_session_XX_<region>.mat (per-session decoders)
%   KDE_reactivation_ripples_PSTH_POST_unique_rate_lt10.mat
%   POST_rate_lt10\ (AUC_*.mat, figures, AUC_low_vs_high)
%   POST_rate_lt10_bootAUC_per_animal\ (per-animal bootstrapped AUC, temporal summaries)
rev_dir = addon_setup();
run_decoding = false;  % true: retrain the decoders (slow; skips sessions already saved)
run_analysis = true;   % false: only redraw the low vs high bar plot from saved results
overwrite    = false;  % true: redo sessions that already have a decoder output

%% 1. PLS-KDE decoding per session and region, neurons < 10 Hz -> PLS_KDE_rate_lt10
if run_decoding
    run_PLS_KDE_rate_lt10_POST('max_rate', 10, 'overwrite', overwrite);
end

%% 2. Merge sessions into the ripple PSTH file (non-duplicated events)
if run_decoding
    assemble_PLS_KDE_rate_lt10_PSTH_POST();
end

%% 3. Reactivation coherence AUC with the < 10 Hz decoders -> POST_rate_lt10
if run_analysis
    reactivation_AUC_coherence_PRE('POST_rate_lt10');
end
plot_AUC_low_vs_high('POST_rate_lt10');

%% 4. Per-animal bootstrapped temporal AUC with the < 10 Hz decoders -> POST_rate_lt10_bootAUC_per_animal
% (summary_temporal and summary_temporal_diff; as main_02 section 7 for the original decoder)
if run_analysis
    boot_AUC_per_unit_POST_original('animal', [], false, 'rate_lt10');
end
