function run_Bayesian_RUN_decoding_CV(varargin)
% 10-fold cross-validated Bayesian decoding of RUN1 position on both tracks
% (10 cm bins, 14 per track, 28 joint states; posterior of each time bin sums to 1
% over both tracks), with a Poisson ('poisson') and a multinomial naive Bayes
% ('multinomial', invariant to global rate scaling) decoder; see
% bayes_decode_two_track.m.
% Input: RUN1_decoding_input.mat (prepare_RUN1_decoding_input.m): all spatial cells
% (peak_percentile > 0.95 & odd_even_stability > 0.95 on either track) of each region.
% Laps of each track are split into 10 folds; rate maps (occupancy-normalised, running
% samples only) come from the other 9 folds of both tracks and the held-out laps are
% decoded in running time bins.
% Options:
%   'sessions'        : indices into the 22 sessions (default all)
%   'time_bin'        : decoding bin in s (default 0.1)
%   'speed_threshold' : running threshold in cm/s (default 5)
%   'rate_floor'      : minimum rate map value in Hz (default 0.1)
%   'regions'         : default {'HPC','V1','V1_L','V1_R'}
%   'methods'         : default {'poisson','multinomial'}
%   'cell_filter'     : 'all' (default) or 'rate_lt10' (cells with mean RUN1 firing rate
%                       < 10 Hz, i.e. excluding putative interneurons)
%   'out_dir'         : default nature_neuroscience_revision_analysis\Bayesian_RUN_decoding
%                       (with _<cell_filter> appended unless 'all')
%   'overwrite'       : redo sessions with an existing output (default false)
%   'summary_only'    : only rebuild animal / all-animal summaries (default false)
% Output: out_dir\<method>\<animal>\session_<date>_<bin>ms.mat (+ figures), per-animal
% and all-animal confusion matrices, and out_dir\metrics_<bin>ms.csv.

p = inputParser;
addParameter(p, 'sessions', 1:22, @isnumeric);
addParameter(p, 'time_bin', 0.1, @isnumeric);
addParameter(p, 'speed_threshold', 5, @isnumeric);
addParameter(p, 'rate_floor', 0.1, @isnumeric);
addParameter(p, 'regions', {'HPC', 'V1', 'V1_L', 'V1_R'}, @iscell);
addParameter(p, 'methods', {'poisson', 'multinomial'}, @iscell);
addParameter(p, 'out_dir', '', @ischar);
addParameter(p, 'overwrite', false, @islogical);
addParameter(p, 'summary_only', false, @islogical);
addParameter(p, 'cell_filter', 'all', @(x) ismember(x, {'all', 'rate_lt10'}));
parse(p, varargin{:});
P = p.Results;
P.pos_edges = 0:10:140;
P.nFolds = 10;

revision_folder = get_revision_folder();
if isempty(P.out_dir)
    P.out_dir = fullfile(revision_folder, 'Bayesian_RUN_decoding');
    if ~strcmp(P.cell_filter, 'all'), P.out_dir = [P.out_dir '_' P.cell_filter]; end
end
tag = sprintf('%dms', round(P.time_bin * 1000));

if ~P.summary_only
    S = load(fullfile(revision_folder, 'RUN1_decoding_input.mat'), 'D'); D = S.D;
    for nsession = P.sessions
        decode_session(D(nsession), nsession, P, tag);
    end
end
summarise_all(P, tag);
end

function decode_session(Ds, nsession, P, tag)
if all(arrayfun(@(m) exist(files_for(P, Ds, m, tag), 'file') > 0, 1:numel(P.methods))) && ~P.overwrite
    fprintf('Session %d exists, skipping\n', nsession);
    return
end
fprintf('\n=== Session %d: %s %s ===\n', nsession, Ds.subject, Ds.session);
[beh, spk, folds] = RUN1_session_behaviour(Ds, nsession, P);
nPos = numel(P.pos_edges) - 1;
nCells = numel(Ds.cell_id);
nLaps = size(Ds.lap_times, 1);

% Training statistics per lap: occupancy (running samples) and running spike counts
r = beh.running;
occ = accumarray([beh.lap(r) beh.pos_bin(r)], beh.dt(r), [nLaps nPos]);
r = spk.running;
cnt = accumarray([spk.cell(r) spk.pos_bin(r) spk.lap(r)], 1, [nCells nPos nLaps]);

% Test bins: consecutive time bins within each lap, kept if mean speed > threshold
% and position defined; true position is the sample nearest the bin centre
B = test_bins(Ds, beh, spk, P, nCells);
fprintf('  %d laps (T1 %d, T2 %d), %d running test bins of %s\n', nLaps, sum(Ds.lap_track == 1), ...
    sum(Ds.lap_track == 2), numel(B.lap), tag);

% Decoder cells (mean firing rate over the whole RUN1 session)
run_rate = accumarray(spk.cell, 1, [nCells 1]) / (Ds.tvec(end) - Ds.tvec(1));
switch P.cell_filter
    case 'all', use = true(nCells, 1);
    case 'rate_lt10', use = run_rate < 10;
end

for m = 1:numel(P.methods)
    out = struct();
    out.subject = Ds.subject; out.session = Ds.session; out.session_index = nsession;
    out.method = P.methods{m}; out.params = P;
    out.lap_times = Ds.lap_times; out.lap_track = Ds.lap_track; out.lap_fold = folds;
    out.bins = rmfield(B, 'counts');
    for nr = 1:numel(P.regions)
        cells = region_cells(Ds.region, P.regions{nr}) & use;
        R = struct('name', P.regions{nr}, 'cell_id', Ds.cell_id(cells), 'n_cells', sum(cells));
        post = NaN(2 * nPos, numel(B.lap));
        for f = 1:P.nFolds
            train = folds ~= f;
            maps = zeros(sum(cells), 2 * nPos);
            for track = 1:2
                laps = train & Ds.lap_track == track;
                o = sum(occ(laps, :), 1);
                c = sum(cnt(cells, :, laps), 3);
                rate = c ./ o; rate(:, o == 0) = 0;
                maps(:, (track - 1) * nPos + (1:nPos)) = rate;
            end
            test = folds(B.lap) == f;
            post(:, test) = bayes_decode_two_track(B.counts(cells, test), maps, P.time_bin, P.methods{m}, P.rate_floor);
        end
        R.posterior = single(post);
        R.n_spikes = full(sum(B.counts(cells, :), 1));
        R.metrics = decoding_metrics(post, B.state, R.n_spikes > 0, nPos, P.pos_edges);
        out.region(nr) = R;
        fprintf('  %-12s %-4s (%3d cells): track acc %.2f (T1 %.2f, T2 %.2f) | median error own %.1f cm, other-track template %.1f cm\n', ...
            P.methods{m}, P.regions{nr}, R.n_cells, R.metrics.track_acc, R.metrics.track_acc_by_track(1), ...
            R.metrics.track_acc_by_track(2), R.metrics.median_error_own, R.metrics.median_error_other);
    end
    folder = fileparts(files_for(P, Ds, m, tag));
    if ~exist(folder, 'dir'), mkdir(folder); end
    save(files_for(P, Ds, m, tag), 'out', '-v7.3');
    title_str = sprintf('%s %s (session %d): %s, %s bins, 10-fold CV', Ds.subject, Ds.session, nsession, ...
        P.methods{m}, tag);
    try
        plot_confusion_summary({out.region.metrics}, {out.region.name}, [out.region.n_cells], title_str, ...
            strrep(files_for(P, Ds, m, tag), '.mat', ''));
    catch err
        warning('Session %d plotting failed: %s', nsession, err.message);
    end
    close all
end
end

function f = files_for(P, Ds, m, tag)
f = fullfile(P.out_dir, P.methods{m}, Ds.subject, sprintf('session_%s_%s.mat', Ds.session, tag));
end

function cells = region_cells(region, name)
switch name
    case 'HPC', cells = contains(region, 'HPC');
    case 'V1', cells = contains(region, 'V1');
    otherwise, cells = region == string(name);
end
end

function B = test_bins(Ds, beh, spk, P, nCells)
lap = []; t0 = []; pos = []; spd = []; cnt_rows = []; cnt_cols = [];
nb = 0;
for k = 1:size(Ds.lap_times, 1)
    edges = Ds.lap_times(k, 1):P.time_bin:Ds.lap_times(k, 2);
    if numel(edges) < 2, continue; end
    n = numel(edges) - 1;
    centres = edges(1:n) + P.time_bin / 2;
    in = beh.lap == k;
    sb = discretize(beh.tvec(in), edges);
    v = beh.speed(in);
    ok = ~isnan(sb) & ~isnan(v);
    sp = accumarray(sb(ok), v(ok), [n 1], @mean, NaN);
    sample = interp1(beh.tvec, (1:numel(beh.tvec))', centres', 'nearest');
    pb = beh.pos_bin(sample);
    keep = sp > P.speed_threshold & ~isnan(pb);
    idx = zeros(n, 1); idx(keep) = nb + (1:sum(keep));
    s = spk.lap == k;
    kb = discretize(spk.t(s), edges);
    c = spk.cell(s);
    ok = ~isnan(kb); kb = kb(ok); c = c(ok);
    ok = idx(kb) > 0;
    cnt_rows = [cnt_rows; c(ok)]; cnt_cols = [cnt_cols; idx(kb(ok))]; %#ok<AGROW>
    lap = [lap; k * ones(sum(keep), 1)]; t0 = [t0; edges(keep)']; %#ok<AGROW>
    pos = [pos; pb(keep)]; spd = [spd; sp(keep)]; %#ok<AGROW>
    nb = nb + sum(keep);
end
B.lap = lap; B.t_start = t0; B.pos_bin = pos; B.speed = spd;
B.track = Ds.lap_track(lap);
B.state = (B.track - 1) * (numel(P.pos_edges) - 1) + pos;
B.counts = sparse(cnt_rows, cnt_cols, 1, nCells, nb);
end

function M = decoding_metrics(post, state, valid, nPos, pos_edges)
% Confusion matrices and summary measures over bins with >= 1 spike.
% joint_*      : 28 x 28, rows = true state, columns = decoded state (rows sum to 1)
% template_*   : {data track, template track} 14 x 14, posterior restricted to the
%                template track and renormalised (= decoding with that template alone)
% track_mass   : 2 x 2, mean posterior summed over each track, rows = true track
bin_cm = pos_edges(2) - pos_edges(1);
post = post(:, valid); state = state(valid);
true_track = ceil(state / nPos); true_pos = state - (true_track - 1) * nPos;
nS = 2 * nPos;
M.n_bins = [sum(true_track == 1) sum(true_track == 2)];
onehot = sparse(1:numel(state), state, 1, numel(state), nS);
n_state = full(sum(onehot, 1))';
M.joint_soft = (full(onehot' * post')) ./ n_state;
[~, map] = max(post, [], 1);
M.joint_MAP = full(onehot' * sparse(1:numel(map), map, 1, numel(map), nS)) ./ n_state;
M.template_soft = cell(2); M.template_MAP = cell(2);
M.track_mass = NaN(2);
M.median_error_by_template = NaN(2);
M.track_mass_by_pos = NaN(2, nPos);
for a = 1:2
    sel = true_track == a;
    oh = sparse(1:sum(sel), true_pos(sel), 1, sum(sel), nPos);
    n_pos = full(sum(oh, 1))';
    for b = 1:2
        rows = (b - 1) * nPos + (1:nPos);
        pb = post(rows, sel);
        M.track_mass(a, b) = mean(sum(pb, 1));
        pb = pb ./ sum(pb, 1);
        M.template_soft{a, b} = full(oh' * pb') ./ n_pos;
        [~, mb] = max(pb, [], 1);
        M.template_MAP{a, b} = full(oh' * sparse(1:numel(mb), mb, 1, numel(mb), nPos)) ./ n_pos;
        M.median_error_by_template(a, b) = median(abs(mb(:) - true_pos(sel))) * bin_cm;
    end
    mass = sum(post((a - 1) * nPos + (1:nPos), sel), 1);
    M.track_mass_by_pos(a, :) = accumarray(true_pos(sel), mass(:), [nPos 1], @mean, NaN)';
    M.track_acc_by_track(a) = mean(mass > 0.5);
end
mass_true = sum(post .* (ceil((1:nS)' / nPos) == true_track'), 1);
M.track_acc = mean(mass_true > 0.5);
M.median_error_own = mean(diag(M.median_error_by_template));                                  % mean of T1|T1 and T2|T2
M.median_error_other = mean([M.median_error_by_template(1, 2) M.median_error_by_template(2, 1)]);  % mean of T1|T2 and T2|T1
map_track = ceil(map / nPos); map_pos = map - (map_track - 1) * nPos;
correct = map_track(:) == true_track;
M.joint_track_correct = mean(correct);
M.joint_median_error = median(abs(map_pos(correct)' - true_pos(correct))) * bin_cm;
end

function plot_confusion_summary(Ms, names, n_cells, title_str, save_path)
% Rows = regions. Columns: joint 28 x 28 confusion (mean posterior), the four
% template-specific 14 x 14 confusions (data track | template track), and the
% posterior on the true track by position.
col = [44, 123, 182; 215, 25, 28] / 256;
nR = numel(Ms);
fig = figure('Name', title_str, 'Position', [50 50 1900 330 * nR], 'Color', 'w');
tl = tiledlayout(nR, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
for nr = 1:nR
    M = Ms{nr};
    nPos = size(M.template_soft{1, 1}, 1);
    nexttile;
    imagesc(M.joint_soft); axis square xy; hold on;
    plot([0.5 2 * nPos + 0.5], [nPos nPos] + 0.5, 'w-', 'LineWidth', 1);
    plot([nPos nPos] + 0.5, [0.5 2 * nPos + 0.5], 'w-', 'LineWidth', 1);
    set(gca, 'XTick', [nPos / 2 1.5 * nPos] + 0.5, 'XTickLabel', {'T1', 'T2'}, ...
        'YTick', [nPos / 2 1.5 * nPos] + 0.5, 'YTickLabel', {'T1', 'T2'}, 'TickDir', 'out');
    clim([0 max(prctile(M.joint_soft(:), 99), eps)]); colorbar;
    xlabel('Decoded state'); ylabel('True state');
    title(sprintf('%s (n=%d): joint, track acc %.2f', names{nr}, n_cells(nr), M.track_acc), 'FontWeight', 'normal', 'Interpreter', 'none');
    cl = max(cellfun(@(x) max(x(:)), M.template_soft(:)));
    for a = 1:2
        for b = 1:2
            nexttile;
            imagesc(M.template_soft{a, b}, [0 cl]); axis square xy; colorbar;
            set(gca, 'XTick', [0.5 nPos / 2 + 0.5 nPos + 0.5], 'XTickLabel', {'0', '70', '140'}, ...
                'YTick', [0.5 nPos / 2 + 0.5 nPos + 0.5], 'YTickLabel', {'0', '70', '140'}, 'TickDir', 'out');
            xlabel(sprintf('Decoded position, T%d template (cm)', b)); ylabel(sprintf('True position, T%d (cm)', a));
            title(sprintf('T%d data | T%d template: err %.0f cm', a, b, M.median_error_by_template(a, b)), 'FontWeight', 'normal');
        end
    end
    nexttile; hold on;
    x = (1:nPos) * 10 - 5;
    for a = 1:2
        plot(x, M.track_mass_by_pos(a, :), 'Color', col(a, :), 'LineWidth', 1.5);
    end
    plot([0 140], [0.5 0.5], 'k--');
    ylim([0 1]); xlim([0 140]); axis square;
    xlabel('True position (cm)'); ylabel('Posterior on true track');
    legend({sprintf('T1 (%.2f)', M.track_acc_by_track(1)), sprintf('T2 (%.2f)', M.track_acc_by_track(2))}, 'Box', 'off', 'Location', 'southeast');
    set(gca, 'TickDir', 'out', 'Box', 'off');
end
colormap(fig, flipud(gray));
title(tl, title_str, 'Interpreter', 'none');
save_figure(fig, save_path);
end

function summarise_all(P, tag)
% Per animal: mean confusion over that animal's sessions; all animals: mean over
% animal means. Also a CSV of per-session metrics and a per-session scatter.
rows = {};
for m = 1:numel(P.methods)
    files = dir(fullfile(P.out_dir, P.methods{m}, 'M*', sprintf('session_*_%s.mat', tag)));
    if isempty(files), continue; end
    sessions = struct('subject', {}, 'metrics', {}, 'names', {}, 'n_cells', {});
    for k = 1:numel(files)
        S = load(fullfile(files(k).folder, files(k).name), 'out'); out = S.out;
        sessions(k).subject = out.subject;
        sessions(k).metrics = {out.region.metrics};
        sessions(k).names = {out.region.name};
        sessions(k).n_cells = [out.region.n_cells];
        for nr = 1:numel(out.region)
            M = out.region(nr).metrics;
            rows(end+1, :) = {P.methods{m}, out.subject, out.session, out.region(nr).name, out.region(nr).n_cells, ...
                sum(M.n_bins), M.track_acc, M.track_acc_by_track(1), M.track_acc_by_track(2), ...
                M.track_mass(1, 1), M.track_mass(2, 2), M.median_error_by_template(1, 1), M.median_error_by_template(1, 2), ...
                M.median_error_by_template(2, 1), M.median_error_by_template(2, 2), M.joint_track_correct, M.joint_median_error}; %#ok<AGROW>
        end
    end
    names = sessions(1).names;
    animals = unique({sessions.subject});
    animal_means = cell(1, numel(animals));
    for na = 1:numel(animals)
        in = strcmp({sessions.subject}, animals{na});
        animal_means{na} = average_metrics({sessions(in).metrics});
        plot_confusion_summary(animal_means{na}, names, round(mean(vertcat(sessions(in).n_cells), 1)), ...
            sprintf('%s: %s, %d sessions, %s bins (mean n cells)', animals{na}, P.methods{m}, sum(in), tag), ...
            fullfile(P.out_dir, P.methods{m}, animals{na}, sprintf('%s_summary_%s', animals{na}, tag)));
        close all
    end
    plot_confusion_summary(average_metrics(animal_means), names, round(mean(vertcat(sessions.n_cells), 1)), ...
        sprintf('All animals (%d animals, %d sessions): %s, %s bins (mean of animal means)', numel(animals), ...
        numel(sessions), P.methods{m}, tag), fullfile(P.out_dir, P.methods{m}, sprintf('all_animals_summary_%s', tag)));
    close all
end
if isempty(rows), return; end
T = cell2table(rows, 'VariableNames', {'method', 'subject', 'session', 'region', 'n_cells', 'n_bins', ...
    'track_acc', 'track_acc_T1', 'track_acc_T2', 'track_mass_T1', 'track_mass_T2', 'err_T1data_T1tmpl', ...
    'err_T1data_T2tmpl', 'err_T2data_T1tmpl', 'err_T2data_T2tmpl', 'joint_track_correct', 'joint_median_error'});
writetable(T, fullfile(P.out_dir, sprintf('metrics_%s.csv', tag)));
plot_session_metrics(T, P, tag);
end

function A = average_metrics(list)
% Element-wise mean of metric structs (cell of per-region cells)
nR = numel(list{1});
A = cell(1, nR);
for nr = 1:nR
    Ms = cellfun(@(x) x{nr}, list, 'UniformOutput', false);
    A{nr} = Ms{1};
    for fn = fieldnames(A{nr})'
        v = cellfun(@(x) x.(fn{1}), Ms, 'UniformOutput', false);
        if iscell(v{1})
            for c = 1:numel(v{1})
                vc = cellfun(@(x) x{c}, v, 'UniformOutput', false);
                A{nr}.(fn{1}){c} = mean(cat(3, vc{:}), 3, 'omitnan');
            end
        else
            A{nr}.(fn{1}) = mean(cat(3, v{:}), 3, 'omitnan');
        end
    end
end
end

function plot_session_metrics(T, P, tag)
% Per session: track accuracy and median position error with own vs other-track
% template, per region and method (colours = animals)
animals = unique(T.subject);
cmap = lines(numel(animals));
fig = figure('Name', 'Bayesian RUN decoding sessions', 'Position', [50 50 1500 700], 'Color', 'w');
tl = tiledlayout(numel(P.methods), 3, 'TileSpacing', 'compact');
for m = 1:numel(P.methods)
    Tm = T(strcmp(T.method, P.methods{m}), :);
    measures = {{'track_acc'}, {'err_T1data_T1tmpl', 'err_T2data_T2tmpl'}, {'err_T1data_T2tmpl', 'err_T2data_T1tmpl'}};
    labels = {'Track accuracy (posterior > 0.5)', 'Median error, own template (cm)', 'Median error, other-track template (cm)'};
    for k = 1:3
        nexttile; hold on;
        for nr = 1:numel(P.regions)
            in = strcmp(Tm.region, P.regions{nr});
            y = mean(Tm{in, measures{k}}, 2);
            [~, ia] = ismember(Tm.subject(in), animals);
            scatter(nr + 0.15 * randn(sum(in), 1), y, 20, cmap(ia, :), 'filled', 'MarkerFaceAlpha', 0.7);
            plot(nr + [-0.3 0.3], median(y) * [1 1], 'k-', 'LineWidth', 2);
        end
        if k == 1, plot([0.5 numel(P.regions) + 0.5], [0.5 0.5], 'k--'); ylim([0.3 1]); else, ylim([0 70]); end
        set(gca, 'XTick', 1:numel(P.regions), 'XTickLabel', P.regions, 'TickDir', 'out', 'TickLabelInterpreter', 'none');
        ylabel(labels{k}); title(P.methods{m}, 'FontWeight', 'normal');
    end
end
title(tl, sprintf('Bayesian RUN decoding (10-fold CV, %s bins): sessions (dots, colour = animal), bars = median', tag));
save_figure(fig, fullfile(P.out_dir, sprintf('session_metrics_%s', tag)));
close(fig);
end

function save_figure(fig, save_path)
savefig(fig, [save_path '.fig']);
exportgraphics(fig, [save_path '.png'], 'Resolution', 150);
exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
end

function folder = get_revision_folder()
if exist('D:\corticohippocampal_replay', 'dir') > 0
    folder = 'D:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    folder = 'P:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
end
end
