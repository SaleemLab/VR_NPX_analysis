%% main_11_nonripple_V1_bias: |V1 track bias| outside ripples (R4.1), POST
% Is V1 track bias a ripple-specific feature or a general property of V1 during
% sleep? Compares |V1 bias| in ripple-free 200 ms windows with peri-ripple windows
% (0-0.2 s and -0.2-0 s), before and after equalising V1 L / V1 R MUA, and shows
% |V1 bias| as a function of distance from the nearest ripple, and as a time
% course of 100 ms non-overlapping windows (raw and V1 MUA-matched to the
% peri-ripple windows).
% Bias from the original KDE log-odds ripple PSTH (+-1 s around each ripple)
% without nan_mask, so ripple-free windows lie within 1 s of a ripple.
%
% Inputs (main_00): ripple_info.mat, KDE_reactivation_ripples_PSTH_POST_unique.mat,
%   ripple_MUA_HPC.mat, ripple_MUA_V1_L.mat, ripple_MUA_V1_R.mat, ripples_all_POST.mat
% Outputs (nature_neuroscience_revision_analysis\nonripple_V1_bias):
%   nonripple_V1_abs_bias.mat (window table W, session table T, stats S / SD,
%   time course TC), *_distributions, *_session_means, *_distance,
%   *_time_course_100ms[_MUA_matched]_{sessions_SE,sessions,animals},
%   *_time_course_100ms_MUA_matching_check (.fig/.png/.pdf), stats .csv
rev_dir = addon_setup();

nonripple_V1_abs_bias_POST();
