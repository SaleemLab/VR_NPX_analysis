function extract_pairwise_counts_POST(varargin)
% Per-event spike counts of single V1 and HPC neurons and of the sorted-unit MUA of
% each region, around the non-duplicated POST ripples (ripple_info rows), for the
% pairwise V1-HC coupling analysis (pairwise_V1_HC_coupling_POST).
%   neurons : spatial cells of the decoder (peak_percentile > 0.95 & odd_even_stability
%             > 0.95 on either track) with spikes in RUN1 and POST sleep and mean RUN1
%             rate < max_rate (as run_PLS_KDE_rate_lt10_POST)
%   track preference : z_track_diff of extract_V1_ripple_PSTH (Extended Data Fig. 6),
%             (mean Track 1 rate - mean Track 2 rate) / SD of the pooled lap x position
%             raw rate maps (spatial_response_raw, speed > 1 cm/s); > 0 = prefers Track 1 (L)
%   windows (s from ripple onset, [a b)): V1 pre [-0.2 0], V1 post [0 0.2], HC [0 0.1]
%   MUA : summed spikes of all sorted clusters of V1_L, V1_R and HPC (any cell type or
%         rate) in the same windows; the coupling analysis subtracts each cell's own
%         spikes (leave-one-out)
% Options:
%   'sessions' : indices into the 22 POST sessions (default all)
%   'max_rate' : firing-rate cut-off in Hz (default 10)
%   'overwrite': redo existing outputs (default false)
%   'out_name' : output subfolder (default 'counts'; e.g. 'counts_all_rates' with max_rate Inf)
%   'sleep'    : 'POST' (default; ripple_info.mat, session_clusters_Sleep.mat) or 'PRE'
%                (ripple_info_PRE_unique.mat, ripples_all_PRE.mat, session_clusters_Sleep_PRE.mat;
%                sessions without PRE sleep are skipped). Session numbers are the 22 POST
%                sessions in both cases.
% Output: nature_neuroscience_revision_analysis\pairwise_V1_HC_coupling\<out_name>\session_XX.mat

p = inputParser;
addParameter(p, 'sessions', 1:22, @isnumeric);
addParameter(p, 'max_rate', 10, @isnumeric);
addParameter(p, 'overwrite', false, @islogical);
addParameter(p, 'out_name', 'counts', @ischar);
addParameter(p, 'sleep', 'POST', @(s) ismember(s, {'POST', 'PRE'}));
parse(p, varargin{:});
P = p.Results;

rev_dir = addon_setup();
out_dir = fullfile(rev_dir, 'pairwise_V1_HC_coupling', P.out_name);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

windows = struct('name', {'V1_pre', 'V1_post', 'HC'}, 'range', {[-0.2 0], [0 0.2], [0 0.1]});

% Onsets of the non-duplicated ripples (rows of ripple_info / ripple_info_PRE_unique)
if strcmp(P.sleep, 'POST')
    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    load(fullfile(rev_dir, 'ripples_all_POST.mat'), 'ripples_all');
    session_field = 'session_count';
    stim_name = 'SleepChronic'; sleep_file = 'session_clusters_Sleep.mat';
else
    load(fullfile(rev_dir, 'ripple_info_PRE_unique.mat'), 'ripple_info_PRE');
    ripple_info = ripple_info_PRE; clear ripple_info_PRE
    load(fullfile(rev_dir, 'ripples_all_PRE.mat'), 'ripples_all');
    session_field = 'mapped_session_id';                 % PRE sessions mapped to the 22 POST sessions
    stim_name = 'SleepChronic_PRE'; sleep_file = 'session_clusters_Sleep_PRE.mat';
end
onset_all = [ripples_all(1).SWS_onset; ripples_all(2).SWS_onset];
session_all = [ripples_all(1).(session_field)(ripples_all(1).SWS_index == 1); ...
    ripples_all(2).(session_field)(ripples_all(2).SWS_index == 1)];
assert(isequal(session_all(ripple_info.event_id(:)), ripple_info.session_id(:)), ...
    'ripples_all %s SWS rows do not match the ripple info table', P.sleep);
clear ripples_all

SUBJECTS = {'M24016', 'M24017', 'M24018', 'M24062', 'M24064', 'M24065'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS, 'bilateral');
experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);

for nsession = P.sessions
    f_out = fullfile(out_dir, sprintf('session_%02d.mat', nsession));
    if exist(f_out, 'file') && ~P.overwrite
        fprintf('Session %d exists, skipping\n', nsession);
        continue
    end
    rows = find(ripple_info.session_id(:) == nsession);
    if isempty(rows)
        fprintf('Session %d: no %s ripples, skipping\n', nsession, P.sleep);
        continue
    end
    is_sleep = strcmp(experiment_info(nsession).StimulusName, stim_name);
    assert(sum(is_sleep) == 1, 'Session %d: expected one %s stimulus', nsession, stim_name);
    session_info = experiment_info(nsession).session(is_sleep);
    options = session_info.probe(1);
    fprintf('\n=== Session %d (%s): %s %s ===\n', nsession, P.sleep, options.SUBJECT, options.SESSION);
    [sc_RUN, sc_sleep, T_run] = load_session(options, sleep_file);

    C = struct();
    C.session_count = nsession; C.SUBJECT = options.SUBJECT; C.SESSION = options.SESSION;
    C.sleep = P.sleep; C.options = P; C.windows = windows;
    C.ripple_info_rows = rows;
    C.subject_id = ripple_info.subject_id(rows(1));
    C.onset = onset_all(ripple_info.event_id(rows));
    n_ev = numel(rows);

    % Single neurons and their counts
    for region = {'V1', 'HPC'}
        N = select_neurons(sc_RUN, sc_sleep, region{1}, T_run, P.max_rate);
        wins = windows(contains({windows.name}, 'V1'));
        if strcmp(region{1}, 'HPC'), wins = windows(strcmp({windows.name}, 'HC')); end
        for w = 1:numel(wins)
            N.counts.(wins(w).name) = zeros(n_ev, numel(N.cluster_id));
            for c = 1:numel(N.cluster_id)
                t = sort(sc_sleep.spike_times{sc_sleep.cluster_id == N.cluster_id(c)});
                N.counts.(wins(w).name)(:, c) = window_counts(t, C.onset, wins(w).range);
            end
        end
        C.(region{1}) = N;
        fprintf('  %s: %d spatial cells, %d < %g Hz\n', region{1}, N.n_spatial, numel(N.cluster_id), P.max_rate);
    end

    % Sorted-unit MUA per region (all clusters)
    for region = {'V1_L', 'V1_R', 'HPC'}
        is_reg = contains(sc_sleep.region, region{1});
        t = sort(vertcat(sc_sleep.spike_times{is_reg}));
        C.MUA.(region{1}).n_clusters = sum(is_reg);
        for w = 1:numel(windows)
            C.MUA.(region{1}).(windows(w).name) = window_counts(t, C.onset, windows(w).range);
        end
    end
    fprintf('  MUA clusters: V1_L %d, V1_R %d, HPC %d; %d events\n', C.MUA.V1_L.n_clusters, ...
        C.MUA.V1_R.n_clusters, C.MUA.HPC.n_clusters, n_ev);
    save(f_out, 'C', '-v7.3');
    fprintf('  saved %s\n', f_out);
end
end

function [sc_RUN, sc_sleep, T_run] = load_session(options, sleep_file)
f_run = fullfile(options.ANALYSIS_DATAPATH, '..', 'session_clusters_RUN.mat');
if exist(fullfile(options.ANALYSIS_DATAPATH, '..', 'session_clusters_RUN1.mat'), 'file')
    f_run = fullfile(options.ANALYSIS_DATAPATH, '..', 'session_clusters_RUN1.mat');
end
S = load(f_run, 'session_clusters'); sc_RUN = S.session_clusters;
S = load(fullfile(options.ANALYSIS_DATAPATH, '..', sleep_file), 'session_clusters'); sc_sleep = S.session_clusters;
T_run = sc_RUN.tvec{1}(end) - sc_RUN.tvec{1}(1);
end

function N = select_neurons(sc_RUN, sc_sleep, region, T_run, max_rate)
% Spatial cells of the region with spikes in RUN and sleep (as run_PLS_KDE_rate_lt10_POST),
% mean RUN1 rate < max_rate; track preference z_track_diff as extract_V1_ripple_PSTH
spatial = (sc_RUN.peak_percentile(:, 1) > 0.95 & sc_RUN.odd_even_stability(:, 1) > 0.95) | ...
    (sc_RUN.peak_percentile(:, 2) > 0.95 & sc_RUN.odd_even_stability(:, 2) > 0.95);
ids = sc_RUN.cluster_id(spatial & contains(sc_RUN.region, region));
has_run = arrayfun(@(c) ~isempty(sc_RUN.spike_times{sc_RUN.cluster_id == c}), ids);
has_sleep = ismember(ids, sc_sleep.cluster_id);
has_sleep(has_sleep) = arrayfun(@(c) ~isempty(sc_sleep.spike_times{sc_sleep.cluster_id == c}), ids(has_sleep));
ids = ids(has_run & has_sleep);
N.n_spatial = numel(ids);
rate_RUN = arrayfun(@(c) numel(sc_RUN.spike_times{sc_RUN.cluster_id == c}), ids) / T_run;
keep = rate_RUN < max_rate;
N.cluster_id = ids(keep);
N.rate_RUN = rate_RUN(keep);
N.region = arrayfun(@(c) string(sc_RUN.region(sc_RUN.cluster_id == c)), N.cluster_id);
N.z_track_diff = NaN(numel(N.cluster_id), 1);
for c = 1:numel(N.cluster_id)
    idx = find(sc_RUN.cluster_id == N.cluster_id(c));
    r1 = sc_RUN.spatial_response_raw{idx, 1};
    r2 = sc_RUN.spatial_response_raw{idx, 2};
    FR_dist = reshape([r1; r2], 1, []);
    N.z_track_diff(c) = (mean(mean(r1)) - mean(FR_dist)) ./ std(FR_dist) ...
        - (mean(mean(r2)) - mean(FR_dist)) ./ std(FR_dist);
end
end

function n = window_counts(t, onset, range)
% Spikes of sorted times t in [onset + range(1), onset + range(2))
n = n_below(t, onset + range(2)) - n_below(t, onset + range(1));
end

function n = n_below(t, x)
% Number of elements of sorted column t that are < x
[u, ~, ic] = unique(x(:));
cs = cumsum(histcounts(t, [-Inf; u; Inf]));
n = reshape(cs(ic), size(x));
end
