%% main_08_Lseq_decoding: sequence-likelihood (L_seq; Huh et al. 2026) decoding
% Track 1 / Track 2 sequence templates from RUN1 running spikes: 10-fold
% cross-validated scoring of held-out laps and pseudo-events (RUN1), then scoring
% of POST ripples on the same events and windows as the Bayesian replay decoding.
%
% Code: 08_Lseq_decoding\ (helpers in shared\: RUN1_session_behaviour)
% Inputs (main_00): RUN1_decoding_input.mat; Bayesian_replay_decoding\ (main_07,
%   events and windows); session_clusters_Sleep.mat of each session (original pipeline)
% Outputs (nature_neuroscience_revision_analysis):
%   Lseq_RUN_decoding (incl. <animal>\session_<date>_templates_* figures)
%   Lseq_replay_decoding (incl. <animal>\session_<date>_examples figures)
rev_dir = addon_setup();
run_decoding = true;   % false: only rebuild the summaries / figures from saved sessions
overwrite    = false;  % true: redo sessions that already have an output file

%% 1. RUN1 decoding, 10-fold CV -> Lseq_RUN_decoding
run_Lseq_RUN_decoding_CV('overwrite', overwrite, 'summary_only', ~run_decoding);

%% 2. Template matrices of every session -> Lseq_RUN_decoding\<animal>
files = dir(fullfile(rev_dir, 'Lseq_RUN_decoding', '*', 'session_*.mat'));
for k = 1:numel(files)
    [~, subject] = fileparts(files(k).folder);
    session = regexp(files(k).name, '^session_(\d+)\.mat$', 'tokens', 'once');
    if ~isempty(session)
        plot_Lseq_templates_session(subject, session{1});
    end
end

%% 3. POST replay scoring -> Lseq_replay_decoding (needs main_07 section 2)
run_Lseq_replay_decoding_POST('overwrite', overwrite, 'summary_only', ~run_decoding);

%% 4. Example events (also drawn per session by section 3) -> Lseq_replay_decoding\<animal>
plot_Lseq_replay_examples();
