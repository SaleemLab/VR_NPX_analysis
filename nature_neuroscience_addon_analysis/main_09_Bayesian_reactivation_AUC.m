%% main_09_Bayesian_reactivation_AUC: reactivation coherence AUC with Bayesian log odds, POST
% Multinomial Bayesian track log odds (decoders of main_07) on the PLS-KDE PSTH grid
% (20 ms windows, 10 ms step, -1 to 1 s from ripple onset); windows without spikes
% have a flat posterior (log odds 0) and are kept; neighbouring ripples are masked with
% the KDE nan_mask. Then the AUC as main_01 'POST' and the low vs high bar plot.
%
% Inputs: Bayesian_replay_decoding(_rate_lt10)\ session files (main_07), ripple_info.mat,
%   ripples_all_POST.mat, KDE_reactivation_ripples_PSTH_POST_unique.mat
% Outputs (nature_neuroscience_revision_analysis):
%   Bayes_reactivation_ripples_PSTH_POST_unique(_rate_lt10).mat
%   POST_Bayes\, POST_Bayes_rate_lt10\ (AUC_*.mat, figures, AUC_low_vs_high)
rev_dir = addon_setup();
run_decoding = true;   % false: only redo the AUC from the saved PSTH files
cell_filters = {'rate_lt10', 'all'};
tags = {'POST_Bayes_rate_lt10', 'POST_Bayes'};

for k = 1:numel(cell_filters)
    if run_decoding
        run_Bayesian_logodds_PSTH_POST('cell_filter', cell_filters{k});
    end
    reactivation_AUC_coherence_PRE(tags{k});
    plot_AUC_low_vs_high(tags{k});
end
