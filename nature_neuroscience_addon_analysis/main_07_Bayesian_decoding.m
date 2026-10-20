%% main_07_Bayesian_decoding: Bayesian position decoding, RUN1 (CV) and POST replay
% Naive Bayes decoding over the joint state space of both tracks (10 cm bins):
% 10-fold cross-validated decoding of RUN1 laps (Poisson and multinomial), and
% multinomial replay decoding of POST ripples (20 ms bins, -0.5 to 0.5 s) with
% weighted-correlation scores vs shuffles, then the replay follow-up analyses.
% Each is run with all spatial cells and with cells < 10 Hz ('rate_lt10').
%
% Code: 07_Bayesian_decoding\ (helpers in shared\: RUN1_session_behaviour)
% Inputs (main_00): RUN1_decoding_input.mat, ripple_info.mat, ripples_all_POST.mat;
%   session_clusters_Sleep.mat of each session (original pipeline)
% Outputs (nature_neuroscience_revision_analysis):
%   Bayesian_RUN_decoding, Bayesian_RUN_decoding_rate_lt10
%   Bayesian_replay_decoding, Bayesian_replay_decoding_rate_lt10 (each with analysis\)
rev_dir = addon_setup();
run_decoding = true;   % false: only rebuild the summaries / figures from saved sessions
overwrite    = false;  % true: redo sessions that already have an output file
cell_filters = {'all', 'rate_lt10'};

%% 1. RUN1 decoding, 10-fold CV -> Bayesian_RUN_decoding(_rate_lt10)
for k = 1:numel(cell_filters)
    run_Bayesian_RUN_decoding_CV('cell_filter', cell_filters{k}, 'overwrite', overwrite, ...
        'summary_only', ~run_decoding);
end

%% 2. POST replay decoding -> Bayesian_replay_decoding(_rate_lt10)
for k = 1:numel(cell_filters)
    run_Bayesian_replay_decoding_POST('cell_filter', cell_filters{k}, 'overwrite', overwrite, ...
        'summary_only', ~run_decoding);
end

%% 3. Replay follow-up: pre vs post window, posteriors, firing by place-field peak -> <folder>\analysis
analyse_Bayesian_replay_POST('Bayesian_replay_decoding');
analyse_Bayesian_replay_POST('Bayesian_replay_decoding_rate_lt10');
