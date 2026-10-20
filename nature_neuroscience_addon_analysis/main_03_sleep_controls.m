%% main_03_sleep_controls: PRE vs POST sleep controls
% Do POST ripples differ from PRE in interhemispheric SO synchrony, in
% spindle-band power or in the ripple power distribution?
%
% Code: 03_sleep_controls\ (helpers in shared\)
% Inputs (main_00): ripple_info.mat, ripple_info_PRE_unique.mat,
%   KDE_reactivation_ripples_PSTH_PRE_unique.mat, KDE_reactivation_ripples_PSTH_POST_unique.mat,
%   ripples_all_PRE/POST.mat, behavioural_state_merged_all_PRE/POST.mat (SWS epochs)
% Outputs (nature_neuroscience_revision_analysis):
%   SO_bilateral_synchrony_PRE_vs_POST, spindle_ripple_PRE_vs_POST
rev_dir = addon_setup();

%% 1. Interhemispheric SO synchrony at ripple times
analyse_SO_bilateral_synchrony_PRE_vs_POST();

%% 2. Dominant-hemisphere spindle-band power and ripple power distribution
analyse_spindle_ripple_PRE_vs_POST();
