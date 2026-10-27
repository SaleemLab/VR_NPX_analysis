%% main_10_pairwise_coupling: pairwise V1-HC neuron coupling across POST ripples
% R1 summary ("interactions between V1 and CA1 neurons ... pairwise neuronal
% correlations") and R4.2 ("is activity in the two regions more synchronised").
% For every V1 x HPC pair (spatial cells, mean RUN1 rate < 10 Hz), the correlation across
% ripples of the V1 count (-0.2-0 s or 0-0.2 s) and the HPC count (0-0.1 s) is related to
% the pair's shared track preference (z_track_diff of Extended Data Fig. 6, product
% z_V1 * z_HC), against a null that permutes preferences within each region.
% Gain control (leave-one-out sorted-unit MUA, partialled out of both cells):
%   'total_V1'       : V1_L + V1_R MUA and HPC MUA (keeps the hemisphere difference)
%   'own_hemisphere' : MUA of the V1 cell's own hemisphere and HPC MUA (removes it)
%
% Code: extract_pairwise_counts_POST, pairwise_V1_HC_coupling_POST, plot_pairwise_V1_HC_coupling,
%   plot_pairwise_summary, pairwise_coupling_oscillation_POST (pairwise_partial_corr_matrix),
%   plot_pairwise_oscillation_effects
% Inputs: ripple_info.mat, ripples_all_POST.mat, KDE_reactivation_ripples_PSTH_POST_unique.mat;
%   per-session session_clusters_RUN1.mat / session_clusters_Sleep.mat of the original pipeline
% Outputs (nature_neuroscience_revision_analysis\pairwise_V1_HC_coupling):
%   counts\session_XX.mat (per-event single-cell and MUA counts)
%   gain_total_V1\, gain_own_hemisphere\ (coupling_POST.mat, session_stats.csv,
%   group_stats.csv, figures)
rev_dir = addon_setup();
run_counts   = false;  % true: re-extract counts from the per-session spike files
run_analysis = true;   % false: only redraw figures from saved coupling results
overwrite    = false;  % true: redo sessions that already have a counts file

%% 1. Per-event spike counts and sorted-unit MUA -> counts\
if run_counts
    extract_pairwise_counts_POST('max_rate', 10, 'overwrite', overwrite);
end

%% 2. Pair coupling and preference statistics, both gain controls
for gain_control = {'total_V1', 'own_hemisphere'}
    if run_analysis
        pairwise_V1_HC_coupling_POST(gain_control{1}, 'nPerm', 1000, 'min_active', 10);
    end
    plot_pairwise_V1_HC_coupling(gain_control{1});
    % Revision summary figures (partial correlation; needs main_02 per-session bootAUC)
    plot_pairwise_summary(gain_control{1});
end

%% 3. Oscillation conditions (total V1 gain): ripple power Q1/Q4; own-hemisphere spindle
% Q1/Q4, SO peak/trough, bilateral/unilateral trough; with binomial spike thinning
% -> gain_total_V1\oscillation_effects
if run_analysis
    pairwise_coupling_oscillation_POST('nPerm', 1000, 'nThin', 20, 'nPermThin', 50, 'min_active', 10);
end
plot_pairwise_oscillation_effects();
