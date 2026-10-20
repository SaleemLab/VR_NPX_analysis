function plot_Lseq_replay_examples(varargin)
% Example POST ripple events from run_Lseq_replay_decoding_POST, one column per event and
% four raster rows: HPC cells ordered by RUN place-field peak on Track 1, HPC ordered on
% Track 2, V1 ordered on Track 1, V1 ordered on Track 2 (only cells active in the shown
% window; y = rank of peak position). Each row is titled with that region's
% L_seq z and two-sided p for that track's (raw) template (post window).
% Columns: the n_sig HPC events with the largest |z| among significant events (p < 0.05
% on either track template), then n_nonsig scored non-significant HPC events drawn at
% random (seeded) for comparison.
% Red dashed: post window (onset to end); blue dotted: start of the pre window.
% Options: 'sessions' (all with Lseq output), 'n_sig' (8), 'n_nonsig' (4), 'seed' (2026).
% Output: Lseq_replay_decoding\<animal>\session_<date>_examples.(fig/png/pdf)

p = inputParser;
addParameter(p, 'sessions', [], @isnumeric);
addParameter(p, 'n_sig', 8, @isnumeric);
addParameter(p, 'n_nonsig', 4, @isnumeric);
addParameter(p, 'seed', 2026, @isnumeric);
parse(p, varargin{:});
P = p.Results;

rev = get_revision_folder();
in_dir = fullfile(rev, 'Lseq_replay_decoding');
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
SUBJECTS = {'M24016', 'M24017', 'M24018', 'M24062', 'M24064', 'M24065'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS, 'bilateral');
experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);

files = dir(fullfile(in_dir, 'M*', 'session_*.mat'));
for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name), 'out'); out = S.out;
    if ~isempty(P.sessions) && ~ismember(out.session_index, P.sessions), continue; end
    fprintf('Session %d %s %s\n', out.session_index, out.subject, out.session);
    info = experiment_info(out.session_index);
    is_post = strcmp(info.StimulusName, 'SleepChronic');
    options = info.session(is_post).probe(1);
    Sl = load(fullfile(options.ANALYSIS_DATAPATH, '..', 'session_clusters_Sleep.mat'), 'session_clusters');
    sc = Sl.session_clusters;
    keep = ismember(sc.cluster_id, out.cell_id);
    t = vertcat(sc.spike_times{keep}); id = vertcat(sc.spike_id{keep});
    [~, ci] = ismember(id, out.cell_id); [t, o] = sort(t); ci = ci(o);
    plot_session(out, t, ci, P, fullfile(files(k).folder, strrep(files(k).name, '.mat', '_examples')));
end
end

function plot_session(out, t, ci, P, save_path)
X = out.region(1).post;
pmin = min(X.p(:, 1:2), [], 2);
sig = find(X.scored & pmin < 0.05);
nonsig = find(X.scored & pmin >= 0.05);
[~, o] = sort(max(abs(X.z(sig, 1:2)), [], 2), 'descend');
ev = sig(o(1:min(P.n_sig, numel(o))));
rs = RandStream('mrg32k3a', 'Seed', P.seed + out.session_index);
ev_ns = nonsig(randperm(rs, numel(nonsig), min(P.n_nonsig, numel(nonsig))));
ev = [ev(:); ev_ns(:)]; is_sig = [true(min(P.n_sig, numel(sig)), 1); false(numel(ev_ns), 1)];
nE = numel(ev);
if nE == 0, return; end
nPos = 14; tb = out.time_bin;
regs = {'HPC', 'V1'};
fig = figure('Name', 'Lseq examples', 'Position', [10 10 230 * nE 1150], 'Color', 'w');
tl = tiledlayout(4, nE, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:nE
    e = ev(k);
    onset = out.events.onset(e); dur = out.events.n_post_bins(e) * tb;
    win = [-dur - 0.05, dur + 0.05];
    row = 0;
    for r = 1:2
        rc = find(region_cells(out.region_labels, regs{r}));
        s = t >= onset + win(1) & t < onset + win(2) & ismember(ci, rc);
        active = unique(ci(s));
        Y = out.region(r).post;
        for track = 1:2
            row = row + 1;
            nexttile(tl, (row - 1) * nE + k); hold on;
            [~, pk] = max(out.rate_maps(active, (track - 1) * nPos + (1:nPos)), [], 2);
            [~, ord] = sort(pk);
            yrow = zeros(numel(out.cell_id), 1); yrow(active(ord)) = 1:numel(active);
            col = [44 123 182; 215 25 28] / 256;
            scatter(t(s) - onset, yrow(ci(s)), 10, col(track, :), '|', 'LineWidth', 1.1);
            xline([0 dur], 'r--'); xline(-dur, 'b:');
            xlim(win); ylim([0 numel(active) + 1]);
            zt = Y.z(e, track); pt = Y.p(e, track);
            if ~isfinite(zt), ttl = sprintf('%s T%d: not scored', regs{r}, track);
            else, ttl = sprintf('%s T%d z %.1f p %.3f', regs{r}, track, zt, pt); end
            if pt < 0.05, fw = 'bold'; else, fw = 'normal'; end
            if row == 1
                lab = {'non-sig', 'SIG'};
                ttl = {sprintf('ev %d %s', e, lab{is_sig(k) + 1}), ttl};
            end
            title(ttl, 'FontWeight', fw, 'FontSize', 7);
            if k == 1, ylabel(sprintf('%s cells, T%d peak order', regs{r}, track), 'FontSize', 8); end
            if row == 4, xlabel('Time from onset (s)'); end
            set(gca, 'TickDir', 'out', 'FontSize', 7);
        end
    end
end
title(tl, sprintf(['%s %s (session %d): L_{seq} track templates (cell-ID shuffle, two-sided). Columns: %d most significant ' ...
    'HPC events, then %d random non-significant. Rows: HPC by T1 / T2 peak, V1 by T1 / T2 peak (bold: p < 0.05)'], ...
    out.subject, out.session, out.session_index, sum(is_sig), sum(~is_sig)), 'FontSize', 9);
save_figure(fig, save_path);
close(fig);
end

function cells = region_cells(region, name)
switch name
    case 'HPC', cells = contains(region, 'HPC');
    case 'V1', cells = contains(region, 'V1');
end
end

function save_figure(fig, save_path)
for attempt = 1:10
    try
        savefig(fig, [save_path '.fig']);
        exportgraphics(fig, [save_path '.png'], 'Resolution', 150);
        exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
        return
    catch err
        if attempt == 10, rethrow(err); end
        pause(3);
    end
end
end

function folder = get_revision_folder()
if exist('D:\corticohippocampal_replay', 'dir') > 0
    folder = 'D:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    folder = 'P:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
end
end
