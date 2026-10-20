function run_Bayesian_logodds_PSTH_POST(varargin)
% Bayesian (multinomial) track log odds of POST sleep ripples on the PLS-KDE PSTH grid,
% so the reactivation coherence AUC (reactivation_AUC_coherence_PRE('POST_Bayes...'))
% can be computed with Bayesian instead of PLS-KDE log odds.
% Decoder per session and region (HPC, V1) = the one of run_Bayesian_replay_decoding_POST
% (RUN1 rate maps of all laps, 10 cm bins of both tracks, 0.1 Hz floor, same decoder
% cells), read from its session files. Grid as PLS_KDE_decode: 20 ms windows stepped by
% 10 ms, window k = [onset - 1 + 0.01 (k - 1), + 0.02), k = 1..200.
%   log odds = log(sum P(Track 1) / sum P(Track 2)) of the joint posterior.
% Windows with no spikes have a flat posterior, i.e. log odds exactly 0; they are kept
% (all windows used). Bins after the next / before the previous ripple of the same probe
% are excluded downstream with the existing KDE nan_mask (same mask as the KDE PSTH).
% Rows: all SWS ripples, probe 1 then probe 2, sessions in order (as
% KDE_reactivation_ripples_PSTH.mat); log odds z-scored per session and probe over all
% windows of all its SWS ripples; then rows ripple_info.event_id (non-duplicated events).
% Options: 'cell_filter' ('rate_lt10' (default) or 'all'), 'sessions' (1:22).
% Output: nature_neuroscience_revision_analysis\
%   Bayes_reactivation_ripples_PSTH_POST_unique[_rate_lt10].mat (struct
%   Bayes_reactivation_ripples_PSTH: HPC/V1 _logodds_ripples, _z_logodds_ripples,
%   _n_spikes_ripples (spikes per window), nan_mask, event_id, n_cells)

p = inputParser;
addParameter(p, 'cell_filter', 'rate_lt10', @(x) ismember(x, {'all', 'rate_lt10'}));
addParameter(p, 'sessions', 1:22, @isnumeric);
parse(p, varargin{:});
P = p.Results;
P.nbins = 200; P.step = 0.01; P.bin_size = 0.02; P.nPos = 14;
P.regions = {'HPC', 'V1'};

rev_dir = get_revision_folder();
dec_dir = fullfile(rev_dir, 'Bayesian_replay_decoding');
suffix = '';
if ~strcmp(P.cell_filter, 'all'), suffix = ['_' P.cell_filter]; dec_dir = [dec_dir suffix]; end
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))

load(fullfile(rev_dir, 'ripples_all_POST.mat'), 'ripples_all');
load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
K = load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
K = K.KDE_reactivation_ripples_PSTH;

files = dir(fullfile(dec_dir, 'M*', 'session_*.mat'));
rows = struct();
for r = 1:2
    rows.(P.regions{r}) = struct('lo', {cell(2, max(P.sessions))}, 'n', {cell(2, max(P.sessions))});
end
n_cells = NaN(max(P.sessions), 2);
for k = 1:numel(files)
    tic
    S = load(fullfile(files(k).folder, files(k).name), 'out'); out = S.out;
    nsession = out.session_index;
    if ~ismember(nsession, P.sessions), continue; end
    if ~isfield(out, 'decoder_cell'), out.decoder_cell = true(numel(out.cell_id), 1); end   % older 'all' outputs
    [t, cell_idx] = load_sleep_spikes(session_experiment_info(nsession), out.cell_id);
    for r = 1:2
        cells = region_cells(out.region_labels, P.regions{r}) & out.decoder_cell(:);
        assert(isequal(out.cell_id(cells), out.region(r).cell_id), 'Decoder cells differ in session %d', nsession);
        nCr = sum(cells); n_cells(nsession, r) = nCr;
        new_idx = zeros(numel(out.cell_id), 1); new_idx(cells) = 1:nCr;
        in = cells(cell_idx);
        tr = t(in); cr = new_idx(cell_idx(in));
        mr = max(out.rate_maps(cells, :), out.params.rate_floor);
        logp = log(mr ./ sum(mr, 1));                         % nCells x 28 (multinomial)
        for nprobe = 1:2
            idx = find(ripples_all(nprobe).session_count == nsession & ripples_all(nprobe).SWS_index == 1);
            onsets = ripples_all(nprobe).onset(idx);
            [lo, ns] = decode_events(tr, cr, nCr, logp, onsets(:), P);
            rows.(P.regions{r}).lo{nprobe, nsession} = lo;
            rows.(P.regions{r}).n{nprobe, nsession} = ns;
        end
    end
    fprintf('Session %2d: %3d HPC, %3d V1 cells (%.0f s)\n', nsession, n_cells(nsession, :), toc);
end

assert(all(isfinite(n_cells(P.sessions, 1))), 'Missing decoder sessions %s', mat2str(P.sessions(~isfinite(n_cells(P.sessions, 1)))));
n_sws = sum(ripples_all(1).SWS_index == 1) + sum(ripples_all(2).SWS_index == 1);
PSTH = struct();
for r = 1:2
    R = rows.(P.regions{r});
    Z = R.lo;
    for nprobe = 1:2
        for nsession = P.sessions
            x = R.lo{nprobe, nsession};
            Z{nprobe, nsession} = (x - mean(x(:))) / std(x(:));
        end
    end
    full_lo = vertcat(R.lo{1, :}, R.lo{2, :});
    full_z = vertcat(Z{1, :}, Z{2, :});
    full_n = vertcat(R.n{1, :}, R.n{2, :});
    assert(size(full_lo, 1) == n_sws, 'Row count %d differs from %d SWS ripples', size(full_lo, 1), n_sws);
    PSTH.([P.regions{r} '_logodds_ripples']) = full_lo(ripple_info.event_id, :);
    PSTH.([P.regions{r} '_z_logodds_ripples']) = full_z(ripple_info.event_id, :);
    PSTH.([P.regions{r} '_n_spikes_ripples']) = full_n(ripple_info.event_id, :);
end
PSTH.nan_mask = K.nan_mask;
PSTH.event_id = ripple_info.event_id;
PSTH.n_cells = n_cells;
PSTH.cell_filter = P.cell_filter;

% Summary and check against PLS-KDE (event means of z log odds, nan_mask applied)
bc = (-1:0.01:1); bc = bc(1:end-1) + 0.005;
wins = {bc > 0 & bc < 0.1, bc > -0.2 & bc < 0, bc > 0 & bc < 0.2}; wn = {'0-0.1 s', '-0.2-0 s', '0-0.2 s'};
for r = 1:2
    reg = P.regions{r};
    keep = isfinite(PSTH.nan_mask);
    fprintf('%s: %.1f%% of unmasked windows have no spikes (log odds 0)\n', reg, ...
        100 * mean(PSTH.([reg '_n_spikes_ripples'])(keep) == 0));
    for w = 1:numel(wins)
        x = mean(K.([reg '_z_logodds_ripples'])(:, wins{w}) + PSTH.nan_mask(:, wins{w}), 2, 'omitnan');
        y = mean(PSTH.([reg '_z_logodds_ripples'])(:, wins{w}) + PSTH.nan_mask(:, wins{w}), 2, 'omitnan');
        ok = isfinite(x) & isfinite(y);
        fprintf('  %s KDE vs Bayesian (%s, %d events): r = %.3f, sign agreement %.1f%%\n', reg, wn{w}, sum(ok), ...
            corr(x(ok), y(ok)), 100 * mean(sign(x(ok)) == sign(y(ok))));
    end
end

Bayes_reactivation_ripples_PSTH = PSTH;
out_file = fullfile(rev_dir, ['Bayes_reactivation_ripples_PSTH_POST_unique' suffix '.mat']);
save(out_file, 'Bayes_reactivation_ripples_PSTH', '-v7.3');
fprintf('Saved %s\n', out_file);
end

function [lo, ns] = decode_events(tr, cr, nCr, logp, onsets, P)
% Log odds and spike count of each 20 ms window (10 ms step) from onset - 1 s, per event
nE = numel(onsets);
lo = NaN(nE, P.nbins); ns = NaN(nE, P.nbins);
nfine = P.nbins + round(P.bin_size / P.step) - 1;    % 10 ms bins covered by the windows
first = first_ge(tr, onsets - 1);
last = first_ge(tr, onsets - 1 + nfine * P.step) - 1;
for e = 1:nE
    s = (first(e):last(e))';
    b = floor(round((tr(s) - (onsets(e) - 1)) / P.step, 6)) + 1;
    ok = b >= 1 & b <= nfine;
    c = accumarray([cr(s(ok)) b(ok)], 1, [nCr nfine]);
    counts = c(:, 1:P.nbins) + c(:, 2:P.nbins + 1);
    L = logp' * counts;                                % 28 x nbins; 0 for empty windows
    lo(e, :) = logsumexp(L(1:P.nPos, :)) - logsumexp(L(P.nPos + 1:end, :));
    ns(e, :) = sum(counts, 1);
end
end

function y = logsumexp(L)
m = max(L, [], 1);
y = m + log(sum(exp(L - m), 1));
end

function idx = first_ge(t, x)
% Index of the first element of sorted t that is >= x (numel(t)+1 if none)
idx = zeros(size(x));
n = numel(t);
for k = 1:numel(x)
    lo = 1; hi = n + 1;
    while lo < hi
        mid = floor((lo + hi) / 2);
        if t(mid) < x(k), lo = mid + 1; else, hi = mid; end
    end
    idx(k) = lo;
end
end

function cells = region_cells(region, name)
switch name
    case 'HPC', cells = contains(region, 'HPC');
    case 'V1', cells = contains(region, 'V1');
end
cells = cells(:);
end

function [t, cell_idx] = load_sleep_spikes(info, cell_id)
% POST sleep (SleepChronic) spikes of the given clusters: times and index into cell_id
is_post = strcmp(info.StimulusName, 'SleepChronic');
assert(sum(is_post) == 1, 'Expected one SleepChronic stimulus');
options = info.session(is_post).probe(1);
S = load(fullfile(options.ANALYSIS_DATAPATH, '..', 'session_clusters_Sleep.mat'), 'session_clusters');
sc = S.session_clusters;
keep = ismember(sc.cluster_id, cell_id);
t = vertcat(sc.spike_times{keep}); id = vertcat(sc.spike_id{keep});
[~, cell_idx] = ismember(id, cell_id);
[t, o] = sort(t); cell_idx = cell_idx(o);
end

function info = session_experiment_info(nsession)
% experiment_info of the 22 bilateral sessions (index = session number)
persistent INFO
if isempty(INFO)
    SUBJECTS = {'M24016', 'M24017', 'M24018', 'M24062', 'M24064', 'M24065'};
    INFO = subject_session_stimuli_mapping(SUBJECTS, 'bilateral');
    INFO = INFO([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
end
info = INFO(nsession);
end

function folder = get_revision_folder()
if exist('D:\corticohippocampal_replay', 'dir') > 0
    folder = 'D:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    folder = 'P:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
end
end
