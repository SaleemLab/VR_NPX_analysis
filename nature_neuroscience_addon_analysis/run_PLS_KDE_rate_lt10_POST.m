function run_PLS_KDE_rate_lt10_POST(varargin)
% PLS-KDE track decoding of POST sleep ripple windows using only neurons with mean RUN1
% firing rate < 10 Hz (control for interneuron contribution). Same decoder as the full
% decoder of run_V1_split_half_decoding_POST (PLS_KDE_train on RUN, PLS_KDE_decode on
% [onset-1 onset+1] of every ripple of each probe), applied to HPC and to V1:
%   neurons: spatial cells (peak_percentile > 0.95 & odd_even_stability > 0.95 on
%   either track), region HPC (HPC_L + HPC_R) or V1 (V1_L + V1_R), with spikes in RUN
%   and POST sleep, and mean rate over the whole RUN1 session < max_rate (10 Hz).
% Options:
%   'sessions' : indices into the 22 POST sessions (default all)
%   'regions'  : {'HPC','V1'} (default both; one output file per region and session)
%   'max_rate' : firing-rate cut-off in Hz (default 10)
%   'nShuffle', 'nBootBW' : passed to PLS_KDE_train (default 1000, 1000)
%   'out_dir'  : default nature_neuroscience_revision_analysis\PLS_KDE_rate_lt10
%   'overwrite': redo existing outputs (default false)
% Output per session and region: POST_session_XX_<region>.mat (struct out: neurons,
% decoder, decoded ripple windows per probe); assemble with
% assemble_PLS_KDE_rate_lt10_PSTH_POST.

p = inputParser;
addParameter(p, 'sessions', 1:22, @isnumeric);
addParameter(p, 'regions', {'HPC', 'V1'}, @iscell);
addParameter(p, 'max_rate', 10, @isnumeric);
addParameter(p, 'nShuffle', 1000, @isnumeric);
addParameter(p, 'nBootBW', 1000, @isnumeric);
addParameter(p, 'out_dir', '', @ischar);
addParameter(p, 'overwrite', false, @islogical);
parse(p, varargin{:});
P = p.Results;

if exist('D:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'P:\corticohippocampal_replay';
end
if isempty(P.out_dir)
    P.out_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis', 'PLS_KDE_rate_lt10');
end
if ~exist(P.out_dir, 'dir'), mkdir(P.out_dir); end

SUBJECTS = {'M24016', 'M24017', 'M24018', 'M24062', 'M24064', 'M24065'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS, 'bilateral');
experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);

for nsession = P.sessions
    files = cellfun(@(r) fullfile(P.out_dir, sprintf('POST_session_%02d_%s.mat', nsession, r)), P.regions, 'UniformOutput', false);
    if all(cellfun(@(f) exist(f, 'file') > 0, files)) && ~P.overwrite
        fprintf('Session %d exists, skipping\n', nsession);
        continue
    end
    is_post = strcmp(experiment_info(nsession).StimulusName, 'SleepChronic');
    assert(sum(is_post) == 1, 'Session %d: expected one SleepChronic stimulus', nsession);
    session_info = experiment_info(nsession).session(is_post);
    options = session_info.probe(1);
    fprintf('\n=== Session %d: %s %s ===\n', nsession, options.SUBJECT, options.SESSION);
    [run, sc_RUN, sc_sleep, ripple_windows] = load_session(session_info);
    for nr = 1:numel(P.regions)
        if exist(files{nr}, 'file') && ~P.overwrite, continue; end
        out = struct();
        out.session_count = nsession; out.SUBJECT = options.SUBJECT; out.SESSION = options.SESSION;
        out.region = P.regions{nr}; out.options = P;
        out.neurons = select_neurons(sc_RUN, sc_sleep, run, P.regions{nr}, P.max_rate);
        ids = out.neurons.cluster_id(out.neurons.keep);
        fprintf('  %s: %d spatial neurons, %d with mean RUN rate < %g Hz\n', P.regions{nr}, numel(out.neurons.keep), numel(ids), P.max_rate);
        spikes_RUN = renumber(flatten_spikes(sc_RUN, ids), ids);
        spikes_sleep = renumber(flatten_spikes(sc_sleep, ids), ids);
        tic
        dec = PLS_KDE_train(run.tvec, run.position, run.speed, run.lap_times, run.track_ID, spikes_RUN, ...
            'nShuffle', P.nShuffle, 'nBootBW', P.nBootBW, 'nNeurons', numel(ids));
        fprintf('  trained: bandwidth %.1f, CV AUC %.3f, CV hit rate %.3f, %d good bins (%.0f s)\n', ...
            dec.bandwidth, dec.bandwidth_AUC, mean(dec.cv_hit_rate), dec.n_good_bins, toc);
        out.cluster_id = ids;
        out.dec = dec;
        for nprobe = 1:numel(ripple_windows)
            out.R(nprobe) = PLS_KDE_decode(dec, spikes_sleep, ripple_windows{nprobe});
        end
        save(files{nr}, 'out', '-v7.3');
        fprintf('  saved %s\n', files{nr});
    end
end
end

function [run, sc_RUN, sc_sleep, ripple_windows] = load_session(session_info)
% RUN behaviour, RUN and POST sleep clusters, and ripple windows per probe (as
% run_V1_split_half_decoding_POST)
options = session_info.probe(1);
f_run = fullfile(options.ANALYSIS_DATAPATH, '..', 'session_clusters_RUN.mat');
if exist(fullfile(options.ANALYSIS_DATAPATH, '..', 'session_clusters_RUN1.mat'), 'file')
    f_run = fullfile(options.ANALYSIS_DATAPATH, '..', 'session_clusters_RUN1.mat');
end
S = load(f_run, 'session_clusters'); sc_RUN = S.session_clusters;
S = load(fullfile(options.ANALYSIS_DATAPATH, '..', 'session_clusters_Sleep.mat'), 'session_clusters'); sc_sleep = S.session_clusters;
S = load(fullfile(options.ANALYSIS_DATAPATH, 'extracted_ripple_events.mat'), 'ripples'); ripples = S.ripples;
run.tvec = sc_RUN.tvec{1};
run.position = sc_RUN.position{1};
run.speed = sc_RUN.speed{1};
run.lap_times = [sc_RUN.start_time_all{1} sc_RUN.end_time_all{1}'];
run.track_ID = sc_RUN.track_ID_all{1};
ripple_windows = cell(1, numel(ripples));
for nprobe = 1:numel(ripples)
    ripple_probe_no = session_info.probe(nprobe).probe_id + 1;
    ripple_windows{ripple_probe_no} = [ripples(nprobe).onset - 1 ripples(nprobe).onset + 1];
end
end

function N = select_neurons(sc_RUN, sc_sleep, run, region, max_rate)
% Spatial cells of the region with spikes in RUN and sleep (decoder order: sorted
% cluster id, as the original analysis); keep = mean RUN1 rate < max_rate
spatial = (sc_RUN.peak_percentile(:, 1) > 0.95 & sc_RUN.odd_even_stability(:, 1) > 0.95) | ...
    (sc_RUN.peak_percentile(:, 2) > 0.95 & sc_RUN.odd_even_stability(:, 2) > 0.95);
ids = sc_RUN.cluster_id(spatial & contains(sc_RUN.region, region));
spikes_RUN = flatten_spikes(sc_RUN, ids);
spikes_sleep = flatten_spikes(sc_sleep, ids);
N.cluster_id = intersect(unique(spikes_RUN(:, 2)), unique(spikes_sleep(:, 2)));
T_run = run.tvec(end) - run.tvec(1);
N.rate_RUN = arrayfun(@(c) sum(spikes_RUN(:, 2) == c), N.cluster_id) / T_run;
N.keep = N.rate_RUN < max_rate;
end

function spikes = flatten_spikes(sc, ids)
% [time cluster_id] for the given clusters, sorted by time
keep = ismember(sc.cluster_id, ids);
t = vertcat(sc.spike_times{keep});
id = vertcat(sc.spike_id{keep});
[t, order] = sort(t);
spikes = [t id(order)];
end

function spikes = renumber(spikes, ids)
% Keep the given clusters and renumber them 1..numel(ids)
[keep, new_id] = ismember(spikes(:, 2), ids);
spikes = [spikes(keep, 1) new_id(keep)];
end
