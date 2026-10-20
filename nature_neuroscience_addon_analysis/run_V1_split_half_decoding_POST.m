function run_V1_split_half_decoding_POST(varargin)
% Split-half V1 PLS-KDE decoding of POST sleep ripple windows (R1.3 control).
% Per session, V1 neurons used in the original decoder (spatially reliable on
% either track, present in RUN and POST sleep; as main_reactivation_analysis.m)
% are split into two disjoint halves A and B matched for track preference,
% selectivity and firing rate: within each preferred-track group (higher mean RUN
% rate on Track 1 or Track 2 laps), neurons are paired by |selectivity index| and
% POST sleep firing rate and one of each pair goes to each half (seeded coin
% flip; see split_neurons). Each half gets its own decoder (PLS_KDE_train on RUN) and
% decodes the same ripple windows [onset-1 onset+1] of both ripple probes
% (PLS_KDE_decode). One half can then select the dominant hemisphere and the
% other measure coherence.
% Options:
%   'sessions'     : indices into the 22 POST sessions (default all)
%   'nSplits'      : random splits per session (default 1; seed 1000*session + split)
%   'include_full' : also decode with all neurons (default false)
%   'nShuffle', 'nBootBW' : passed to PLS_KDE_train (default 1000, 1000)
%   'out_dir'      : output folder (default nature_neuroscience_revision_analysis\
%                    V1_split_half_decoding)
%   'overwrite'    : redo sessions with an existing output file (default false)
% Output per session: POST_session_XX.mat (struct out), decoder figures, and a
% bin-by-bin comparison of log odds with the original KDE_reactivation_V1_ripples.

p = inputParser;
addParameter(p, 'sessions', 1:22, @isnumeric);
addParameter(p, 'nSplits', 1, @isnumeric);
addParameter(p, 'include_full', false, @islogical);
addParameter(p, 'nShuffle', 1000, @isnumeric);
addParameter(p, 'nBootBW', 1000, @isnumeric);
addParameter(p, 'out_dir', '', @ischar);
addParameter(p, 'overwrite', false, @islogical);
parse(p, varargin{:});
P = p.Results;

out_dir = P.out_dir;
if isempty(out_dir)
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    out_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis', 'V1_split_half_decoding');
end
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

SUBJECTS = {'M24016', 'M24017', 'M24018', 'M24062', 'M24064', 'M24065'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS, 'bilateral');
experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);

for nsession = P.sessions
    out_file = fullfile(out_dir, sprintf('POST_session_%02d.mat', nsession));
    if exist(out_file, 'file') && ~P.overwrite
        fprintf('Session %d exists, skipping\n', nsession);
        continue
    end
    is_post = strcmp(experiment_info(nsession).StimulusName, 'SleepChronic');
    assert(sum(is_post) == 1, 'Session %d: expected one SleepChronic stimulus', nsession);
    session_info = experiment_info(nsession).session(is_post);
    options = session_info.probe(1);
    fprintf('\n=== Session %d: %s %s ===\n', nsession, options.SUBJECT, options.SESSION);

    [run, V1_RUN, V1_sleep, N, ripple_windows] = load_session(session_info);

    out = struct();
    out.session_count = nsession;
    out.SUBJECT = options.SUBJECT;
    out.SESSION = options.SESSION;
    out.neurons = N;
    out.options = P;

    names = {}; D = {};
    if P.include_full
        out.full = decode_subset(run, V1_RUN, V1_sleep, N.cluster_id, ripple_windows, P);
        names{end+1} = 'Full'; D{end+1} = out.full;
    end
    for split = 1:P.nSplits
        seed = 1000 * nsession + split;
        [ids_A, ids_B] = split_neurons(N, seed);
        print_split(split, N, ids_A, ids_B);
        out.split(split).seed = seed;
        out.split(split).A = decode_subset(run, V1_RUN, V1_sleep, ids_A, ripple_windows, P);
        out.split(split).B = decode_subset(run, V1_RUN, V1_sleep, ids_B, ripple_windows, P);
        names = [names, {sprintf('Half A (split %d)', split), sprintf('Half B (split %d)', split)}]; %#ok<AGROW>
        D = [D, {out.split(split).A, out.split(split).B}]; %#ok<AGROW>
    end
    save(out_file, 'out', '-v7.3');
    fprintf('Saved %s\n', out_file);

    % Figures and checks (failures here do not stop the batch)
    tag = sprintf('POST_session_%02d', nsession);
    title_str = sprintf('%s %s %s', tag, options.SUBJECT, options.SESSION);
    try
        plot_split_balance(N, out.split(1).A.cluster_id, out.split(1).B.cluster_id, title_str, fullfile(out_dir, tag));
        plot_decoders(D, names, title_str, fullfile(out_dir, tag));
        compare_with_original(D, names, options, title_str, fullfile(out_dir, tag));
        if P.include_full
            compare_V1_split_half_event_means(out_file);
        end
    catch err
        warning('Session %d plotting/comparison failed: %s', nsession, err.message);
    end
    close all
end
end

function [run, V1_RUN, V1_sleep, N, ripple_windows] = load_session(session_info)
% RUN behaviour, V1 spikes (RUN and POST sleep), neuron info and ripple windows
% per probe, with the same neuron selection as main_reactivation_analysis.m
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

spatial = (sc_RUN.peak_percentile(:, 1) > 0.95 & sc_RUN.odd_even_stability(:, 1) > 0.95) | ...
    (sc_RUN.peak_percentile(:, 2) > 0.95 & sc_RUN.odd_even_stability(:, 2) > 0.95);
V1_ids = sc_RUN.cluster_id(spatial & contains(sc_RUN.region, 'V1'));

V1_RUN = flatten_spikes(sc_RUN, V1_ids);
V1_sleep = flatten_spikes(sc_sleep, V1_ids);
% Neurons with spikes in both RUN and sleep, sorted by cluster id (decoder order)
N.cluster_id = intersect(unique(V1_RUN(:, 2)), unique(V1_sleep(:, 2)));
nN = numel(N.cluster_id);
N.hemisphere = strings(nN, 1);
N.rate_T1 = NaN(nN, 1); N.rate_T2 = NaN(nN, 1);
laps = {run.lap_times(run.track_ID == 1, :), run.lap_times(run.track_ID == 2, :)};
dur = cellfun(@(l) sum(l(:, 2) - l(:, 1)), laps);
for i = 1:nN
    N.hemisphere(i) = sc_RUN.region(sc_RUN.cluster_id == N.cluster_id(i));
    t = V1_RUN(V1_RUN(:, 2) == N.cluster_id(i), 1);
    N.rate_T1(i) = sum(InIntervals(t, laps{1})) / dur(1);
    N.rate_T2(i) = sum(InIntervals(t, laps{2})) / dur(2);
end
N.preferred_track = 1 + (N.rate_T2 > N.rate_T1);
N.SI = (N.rate_T1 - N.rate_T2) ./ (N.rate_T1 + N.rate_T2);   % track selectivity index
N.SI(~isfinite(N.SI)) = 0;
% Mean firing rate over the POST sleep session
T_sleep = V1_sleep(end, 1) - V1_sleep(1, 1);
N.rate_sleep = arrayfun(@(c) sum(V1_sleep(:, 2) == c), N.cluster_id) / T_sleep;

ripple_windows = cell(1, numel(ripples));
for nprobe = 1:numel(ripples)
    ripple_probe_no = session_info.probe(nprobe).probe_id + 1;
    ripple_windows{ripple_probe_no} = [ripples(nprobe).onset - 1 ripples(nprobe).onset + 1];
end
end

function spikes = flatten_spikes(sc, ids)
% [time cluster_id] for the given clusters, sorted by time
keep = ismember(sc.cluster_id, ids);
t = vertcat(sc.spike_times{keep});
id = vertcat(sc.spike_id{keep});
[t, order] = sort(t);
spikes = [t id(order)];
end

function [ids_A, ids_B] = split_neurons(N, seed)
% Halves matched for track preference, selectivity and sleep firing rate: within
% each preferred-track group, neurons are ranked by |SI| and by log sleep rate
% (ranks scaled 0-1), the closest pair in that rank space is taken greedily and
% removed until none are left, and a seeded coin flip sends one neuron of each
% pair to A and the other to B. An odd leftover goes to A (Track 1) or B (Track 2).
s = RandStream('mrg32k3a', 'Seed', seed);
ids_A = []; ids_B = [];
for track = 1:2
    in = find(N.preferred_track == track);
    n = numel(in);
    if n == 0, continue; end
    F = [rank01(abs(N.SI(in))) rank01(log(N.rate_sleep(in) + eps))];
    D = abs(F(:, 1) - F(:, 1)') + abs(F(:, 2) - F(:, 2)');
    D(logical(eye(n))) = Inf;
    left = true(n, 1);
    while sum(left) >= 2
        Dl = D; Dl(~left, :) = Inf; Dl(:, ~left) = Inf;
        [~, k] = min(Dl(:));
        [i, j] = ind2sub([n n], k);
        if rand(s) < 0.5, [i, j] = deal(j, i); end
        ids_A = [ids_A; N.cluster_id(in(i))]; %#ok<AGROW>
        ids_B = [ids_B; N.cluster_id(in(j))]; %#ok<AGROW>
        left([i j]) = false;
    end
    if any(left)
        if track == 1, ids_A = [ids_A; N.cluster_id(in(left))]; else, ids_B = [ids_B; N.cluster_id(in(left))]; end %#ok<AGROW>
    end
end
ids_A = sort(ids_A); ids_B = sort(ids_B);
end

function r = rank01(x)
r = (tiedrank(x) - 1) / max(numel(x) - 1, 1);
end

function print_split(split, N, ids_A, ids_B)
ids = {ids_A, ids_B}; lab = 'AB';
for h = 1:2
    in = ismember(N.cluster_id, ids{h});
    fprintf('Split %d half %s: %d neurons (T1-pref %d, T2-pref %d; V1_L %d, V1_R %d) | median |SI| %.2f, sleep rate %.2f Hz (sum %.0f), RUN rate %.2f Hz\n', ...
        split, lab(h), sum(in), sum(in & N.preferred_track == 1), sum(in & N.preferred_track == 2), ...
        sum(in & N.hemisphere == "V1_L"), sum(in & N.hemisphere == "V1_R"), median(abs(N.SI(in))), ...
        median(N.rate_sleep(in)), sum(N.rate_sleep(in)), median((N.rate_T1(in) + N.rate_T2(in)) / 2));
end
A = ismember(N.cluster_id, ids_A); B = ismember(N.cluster_id, ids_B);
fprintf('  KS p: SI %.2f, |SI| %.2f, sleep rate %.2f, RUN rate %.2f\n', ks_p(N.SI(A), N.SI(B)), ...
    ks_p(abs(N.SI(A)), abs(N.SI(B))), ks_p(N.rate_sleep(A), N.rate_sleep(B)), ...
    ks_p(N.rate_T1(A) + N.rate_T2(A), N.rate_T1(B) + N.rate_T2(B)));
end

function p = ks_p(a, b)
[~, p] = kstest2(a, b);
end

function plot_split_balance(N, ids_A, ids_B, title_str, save_path)
% Cumulative distributions of signed SI, |SI|, RUN and sleep firing rate per half
A = ismember(N.cluster_id, ids_A); B = ismember(N.cluster_id, ids_B);
vals = {N.SI, abs(N.SI), (N.rate_T1 + N.rate_T2) / 2, N.rate_sleep};
names = {'Track selectivity index (T1-T2)/(T1+T2)', '|Selectivity index|', 'RUN firing rate (Hz)', 'POST sleep firing rate (Hz)'};
fig = figure('Name', [title_str ' split balance'], 'Position', [100 100 1500 360]);
tiledlayout(1, 4, 'TileSpacing', 'compact');
for k = 1:4
    nexttile; hold on;
    [fa, xa] = ecdf(vals{k}(A)); [fb, xb] = ecdf(vals{k}(B));
    stairs(xa, fa, 'LineWidth', 1.5, 'Color', [0.2 0.6 0.8]);
    stairs(xb, fb, 'LineWidth', 1.5, 'Color', [231 41 138] / 256);
    if k >= 3, set(gca, 'XScale', 'log'); end
    xlabel(names{k}); ylabel('Cumulative fraction');
    title(sprintf('A median %.2f, B median %.2f, KS p = %.2f', median(vals{k}(A)), median(vals{k}(B)), ...
        ks_p(vals{k}(A), vals{k}(B))), 'FontWeight', 'normal');
    if k == 1
        legend({sprintf('Half A (n=%d: T1 %d, T2 %d)', sum(A), sum(A & N.preferred_track == 1), sum(A & N.preferred_track == 2)), ...
            sprintf('Half B (n=%d: T1 %d, T2 %d)', sum(B), sum(B & N.preferred_track == 1), sum(B & N.preferred_track == 2))}, ...
            'Box', 'off', 'Location', 'southeast');
    end
    set(gca, 'TickDir', 'out', 'Box', 'off');
end
sgtitle([title_str ': split balance (track preference x selectivity x sleep rate matched pairs)'], 'Interpreter', 'none');
savefig(fig, [save_path '_split_balance.fig']);
exportgraphics(fig, [save_path '_split_balance.png'], 'Resolution', 200);
exportgraphics(fig, [save_path '_split_balance.pdf'], 'ContentType', 'vector');
end

function D = decode_subset(run, V1_RUN, V1_sleep, ids, ripple_windows, P)
% Train on RUN with the given neurons (renumbered 1..n in cluster id order) and
% decode the ripple windows of each probe
ids = sort(ids(:));
spikes_RUN = renumber(V1_RUN, ids);
spikes_sleep = renumber(V1_sleep, ids);
tic
dec = PLS_KDE_train(run.tvec, run.position, run.speed, run.lap_times, run.track_ID, spikes_RUN, ...
    'nShuffle', P.nShuffle, 'nBootBW', P.nBootBW, 'nNeurons', numel(ids));
fprintf('  trained on %d neurons: bandwidth %.1f, CV AUC %.3f, CV hit rate %.3f, %d good bins (%.0f s)\n', ...
    numel(ids), dec.bandwidth, dec.bandwidth_AUC, mean(dec.cv_hit_rate), dec.n_good_bins, toc);
D.cluster_id = ids;
D.dec = dec;
for nprobe = 1:numel(ripple_windows)
    tic
    D.R(nprobe) = PLS_KDE_decode(dec, spikes_sleep, ripple_windows{nprobe});
    fprintf('  decoded probe %d ripples (%d events, %.0f s)\n', nprobe, size(ripple_windows{nprobe}, 1), toc);
end
end

function spikes = renumber(spikes, ids)
% Keep the given clusters and renumber them 1..numel(ids)
[keep, new_id] = ismember(spikes(:, 2), ids);
spikes = [spikes(keep, 1) new_id(keep)];
end

function plot_decoders(D, names, title_str, save_path)
% Per decoder: PLS projection of RUN training bins (Track 1 blue, Track 2 red)
% and cross-validated classification accuracy by position vs neuron-ID shuffle
% (grey: shuffle mean and 95% CI; shaded columns: good position bins)
col = [44, 123, 182; 215, 25, 28] / 256;
nD = numel(D);
fig = figure('Name', [title_str ' PLS decoders'], 'Position', [100 100 1100 360 * nD]);
tiledlayout(nD, 2, 'TileSpacing', 'compact');
for d = 1:nD
    dec = D{d}.dec;
    nexttile; hold on;
    P = {dec.X_T1_proj, dec.X_T2_proj};
    for t = 1:2
        k = 1:size(P{t}, 1);
        if numel(k) > 4000, k = round(linspace(1, numel(k), 4000)); end
        scatter3(P{t}(k, 1), P{t}(k, 2), P{t}(k, 3), 6, col(t, :), 'filled', 'MarkerFaceAlpha', 0.15);
    end
    view(40, 20); grid on;
    xlabel('PLS 1'); ylabel('PLS 2'); zlabel('PLS 3');
    title(sprintf('%s (%d neurons): PLS projection, RUN', names{d}, dec.nNeurons), 'FontWeight', 'normal');
    legend({'Track 1', 'Track 2'}, 'Box', 'off', 'Location', 'northeast');

    nexttile; hold on;
    x = dec.pos_bins(:)';
    good = ismember(x, dec.good_position);
    yl = [0 1.05];
    for g = find(good)
        patch(x(g) + [-2.5 2.5 2.5 -2.5], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.95 0.85], 'EdgeColor', 'none');
    end
    ci = dec.acc_shuffle_CI;
    fill([x fliplr(x)], [ci(1, :) fliplr(ci(2, :))], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    plot(x, dec.acc_shuffle_mean, 'k', 'LineWidth', 1);
    acc = {dec.acc_T1, dec.acc_T2}; h = gobjects(1, 2);
    for t = 1:2
        m = mean(acc{t}, 1, 'omitnan');
        se = std(acc{t}, 0, 1, 'omitnan') / sqrt(size(acc{t}, 1));
        fill([x fliplr(x)], [m + se fliplr(m - se)], col(t, :), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        h(t) = plot(x, m, 'Color', col(t, :), 'LineWidth', 1.5);
    end
    ylim(yl); xlim([min(x) - 2.5, max(x) + 2.5]);
    xlabel('Position (cm)'); ylabel('Proportion correct (10-fold CV)');
    title(sprintf('%s: %d/%d good bins, CV hit %.2f, KDE bw %.1f (AUC %.2f)', names{d}, sum(good), numel(x), ...
        mean(dec.cv_hit_rate), dec.bandwidth, dec.bandwidth_AUC), 'FontWeight', 'normal');
    legend(h, {'Track 1', 'Track 2'}, 'Box', 'off', 'Location', 'southeast');
    set(gca, 'TickDir', 'out', 'Box', 'off');
end
sgtitle([title_str ': PLS projection and position-wise classification (green: good bins)'], 'Interpreter', 'none');
savefig(fig, [save_path '_PLS_decoders.fig']);
exportgraphics(fig, [save_path '_PLS_decoders.png'], 'Resolution', 200);
exportgraphics(fig, [save_path '_PLS_decoders.pdf'], 'ContentType', 'vector');
print(fig, [save_path '_PLS_decoders.svg'], '-dsvg', '-vector');
end

function compare_with_original(D, names, options, title_str, save_path)
% Bin-by-bin log odds (and z log odds per probe) of each decoder vs the original
% KDE_reactivation_V1_ripples output of this session, and between decoders
f_old = fullfile(options.ANALYSIS_DATAPATH, 'KDE_reactivation_V1_ripples.mat');
S = load(f_old, 'KDE_reactivation_V1_ripples'); old = S.KDE_reactivation_V1_ripples;
nD = numel(D);
L = cell(1, nD + 1); Lz = L;
for nprobe = 1:numel(old)
    lo = log(old(nprobe).event_T1_probability(:) ./ old(nprobe).event_T2_probability(:));
    L{1} = [L{1}; lo]; Lz{1} = [Lz{1}; zfinite(lo)];
    for d = 1:nD
        R = D{d}.R(nprobe);
        assert(numel(R.T1_probability) == numel(lo), 'Bin count differs from original (probe %d)', nprobe);
        lo = log(R.T1_probability(:) ./ R.T2_probability(:));
        L{d+1} = [L{d+1}; lo]; Lz{d+1} = [Lz{d+1}; zfinite(lo)];
    end
end
labels = [{'Original'}, names];
nL = numel(labels);
r = NaN(nL); rz = NaN(nL);
for i = 1:nL
    for j = 1:nL
        ok = isfinite(L{i}) & isfinite(L{j});
        r(i, j) = corr(L{i}(ok), L{j}(ok));
        rz(i, j) = corr(Lz{i}(ok), Lz{j}(ok));
    end
end
T = array2table(r, 'VariableNames', labels, 'RowNames', labels);
fprintf('\nPearson r of bin log odds (finite bins, both probes):\n'); disp(T);
writetable(T, [save_path '_logodds_correlation.csv'], 'WriteRowNames', true);
save([save_path '_logodds_comparison.mat'], 'r', 'rz', 'labels');

pairs = [ones(nL - 1, 1) (2:nL)'];         % original vs each decoder
if nL >= 3, pairs(end+1, :) = [nL - 1 nL]; end % half A vs half B (last split)
fig = figure('Name', [title_str ' log odds comparison'], 'Position', [100 100 380 * size(pairs, 1) 380]);
tiledlayout(1, size(pairs, 1), 'TileSpacing', 'compact');
for k = 1:size(pairs, 1)
    i = pairs(k, 1); j = pairs(k, 2);
    ok = find(isfinite(L{i}) & isfinite(L{j}));
    if numel(ok) > 20000, ok = ok(round(linspace(1, numel(ok), 20000))); end
    nexttile; hold on;
    scatter(L{i}(ok), L{j}(ok), 3, 'k', 'filled', 'MarkerFaceAlpha', 0.05);
    lim = prctile([L{i}(ok); L{j}(ok)], [0.5 99.5]);
    plot(lim, lim, 'r--'); xlim(lim); ylim(lim); axis square;
    xlabel([labels{i} ' log odds']); ylabel([labels{j} ' log odds']);
    title(sprintf('r = %.3f (z per probe r = %.3f)', r(i, j), rz(i, j)), 'FontWeight', 'normal');
    set(gca, 'TickDir', 'out', 'Box', 'off');
end
sgtitle([title_str ': bin log odds, ripple windows'], 'Interpreter', 'none');
savefig(fig, [save_path '_logodds_comparison.fig']);
exportgraphics(fig, [save_path '_logodds_comparison.png'], 'Resolution', 200);
exportgraphics(fig, [save_path '_logodds_comparison.pdf'], 'ContentType', 'vector');
print(fig, [save_path '_logodds_comparison.svg'], '-dsvg', '-vector');
end

function z = zfinite(x)
fin = isfinite(x);
z = (x - mean(x(fin))) / std(x(fin));
end
