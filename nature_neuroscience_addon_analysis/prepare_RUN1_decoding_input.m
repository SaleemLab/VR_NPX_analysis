function prepare_RUN1_decoding_input()
% Minimal per-session RUN1 input for the RUN decoding analyses
% (run_Bayesian_RUN_decoding_CV.m, run_Lseq_RUN_decoding_CV.m).
% Spikes, cell info and behaviour come from session_clusters_all_RUN1.mat (spatial
% cells: peak_percentile > 0.95 & odd_even_stability > 0.95 on either track).
% Lap start/end times are not in that summary, so they are read from each session's
% session_clusters_RUN1.mat (same session order: experiment_info sessions
% [4 5 6 17 ...], as in extract_session_clusters_all_RUN1.m).
% Output: nature_neuroscience_revision_analysis\RUN1_decoding_input.mat (struct array D)

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
if exist('D:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'P:\corticohippocampal_replay';
end
revision_folder = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');

fprintf('Loading session_clusters_all_RUN1.mat ...\n');
S = load(fullfile(revision_folder, 'session_clusters_all_RUN1.mat'), 'session_clusters_all');
S = S.session_clusters_all;

SUBJECTS = {'M24016', 'M24017', 'M24018', 'M24062', 'M24064', 'M24065'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS, 'bilateral');
experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);

D = struct();
for nsession = 1:numel(S.spike_id)
    options = experiment_info(nsession).session(1).probe(1);
    assert(strcmp(options.SUBJECT, strtrim(S.subject(nsession, :))) && strcmp(options.SESSION, strtrim(S.session(nsession, :))), ...
        'Session %d: summary and experiment_info disagree', nsession);
    fprintf('Session %d: %s %s\n', nsession, options.SUBJECT, options.SESSION);

    L = load(fullfile(options.ANALYSIS_DATAPATH, '..', 'session_clusters_RUN1.mat'), 'session_clusters');
    sc = L.session_clusters;
    lap_start = unwrap_cell(sc.start_time_all); lap_end = unwrap_cell(sc.end_time_all);
    lap_track = unwrap_cell(sc.track_ID_all);
    track_summary = unwrap_cell(S.track_ID_all{nsession});
    assert(isequal(lap_track(:), track_summary(:)), 'Session %d: lap track IDs differ from summary', nsession);

    D(nsession).subject = options.SUBJECT;
    D(nsession).session = options.SESSION;
    D(nsession).session_day = S.session_day(nsession);
    D(nsession).cell_id = S.spatial_cell_id{nsession}(:);
    D(nsession).region = string(S.region{nsession}(:));
    D(nsession).cell_type = S.cell_type{nsession}(:);
    D(nsession).interneuron = S.interneurons_all{nsession}(:);
    D(nsession).spike_times = S.spike_times{nsession}(:);
    D(nsession).spike_id = S.spike_id{nsession}(:);
    D(nsession).tvec = unwrap_cell(S.tvec{nsession}); D(nsession).tvec = D(nsession).tvec(:);
    D(nsession).position = unwrap_cell(S.position{nsession}); D(nsession).position = D(nsession).position(:);
    D(nsession).speed = unwrap_cell(S.speed{nsession}); D(nsession).speed = D(nsession).speed(:);
    D(nsession).lap_times = [lap_start(:) lap_end(:)];
    D(nsession).lap_track = lap_track(:);
    fprintf('  %d cells (HPC %d, V1_L %d, V1_R %d), %d laps (T1 %d, T2 %d)\n', numel(D(nsession).cell_id), ...
        sum(contains(D(nsession).region, 'HPC')), sum(D(nsession).region == "V1_L"), sum(D(nsession).region == "V1_R"), ...
        numel(lap_track), sum(lap_track == 1), sum(lap_track == 2));
end

save(fullfile(revision_folder, 'RUN1_decoding_input.mat'), 'D', '-v7.3');
fprintf('Saved RUN1_decoding_input.mat\n');
end

function x = unwrap_cell(x)
if iscell(x), x = x{1}; end
end
