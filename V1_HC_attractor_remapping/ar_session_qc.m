function Q = ar_session_qc()
% Per session and region: number of cells, fraction of cells whose rate-map peak
% (either track) lies in the first 20 cm, cross-validated in-lap track accuracy
% (running windows, lobo3 and lobo5) at the track start (0-20 cm) and on the
% rest of the track, and onset-epoch (0-0.5 s) accuracy over all laps.
% Answers whether a session has enough cells covering the track start.
% Output: results\qc\session_qc.csv and session_qc.(fig/png/pdf)
P = ar_setup();
S = load(fullfile(P.inputs, 'RUN1_decoding_input.mat'), 'D'); D = S.D;
Lt = load(fullfile(P.results, 'data', 'lap_tables.mat'), 'L'); L = Lt.L;
decs = {'bayes_multinomial', 'pls'};
rows = {};
for s = 1:numel(D)
    f = fullfile(P.results, 'decoding', sprintf('session_%02d.mat', s));
    if ~exist(f, 'file'), continue; end
    X = load(f, 'out'); out = X.out; G = out.G; T = L(s).laps;
    peak_pos = rate_map_peaks(D(s), T, P);
    run = G.lap > 0 & G.speed > P.bayes_speed;
    true1 = false(size(G.lap)); true1(G.lap > 0) = T.track(G.lap(G.lap > 0)) == 1;
    start = run & G.pos < 20; rest = run & G.pos >= 20;
    onset = G.lap > 0 & G.t_from_on >= 0.05 & G.t_from_on <= 0.45;
    for r = 1:numel(out.region)
        R = out.region(r);
        cells = ismember(D(s).cell_id, R.cell_id);
        frac_start = mean(peak_pos(cells) < 20);
        for d = 1:numel(decs)
            if isempty(R.(decs{d}).LO)
                acc = NaN(1, 6);
            else
                lo3 = double(R.(decs{d}).LO); lo5 = double(R.(decs{d}).LO5);
                acc = [track_acc(lo3, true1, start), track_acc(lo3, true1, rest), track_acc(lo3, true1, onset), ...
                    track_acc(lo5, true1, start), track_acc(lo5, true1, rest), track_acc(lo5, true1, onset)];
            end
            rows(end + 1, :) = [{s, L(s).animal_index, string(L(s).subject), string(R.name), string(decs{d}), R.n_cells, frac_start}, num2cell(acc)]; %#ok<AGROW>
        end
    end
end
Q = cell2table(rows, 'VariableNames', {'session', 'animal', 'subject', 'region', 'decoder', 'n_cells', 'frac_peak_0_20cm', ...
    'acc_start_lobo3', 'acc_rest_lobo3', 'acc_onset_lobo3', 'acc_start_lobo5', 'acc_rest_lobo5', 'acc_onset_lobo5'});
out_dir = fullfile(P.results, 'qc');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
writetable(Q, fullfile(out_dir, 'session_qc.csv'));

regions = {'HPC', 'V1', 'V1_L', 'V1_R'};
acol = lines(6);
fig = figure('Position', [50 50 1600 700], 'Color', 'w');
tl = tiledlayout(2, numel(regions), 'TileSpacing', 'compact');
for d = 1:numel(decs)
    for r = 1:numel(regions)
        x = Q(Q.region == regions{r} & Q.decoder == decs{d}, :);
        nexttile(tl); hold on;
        vals = [x.acc_start_lobo3 x.acc_rest_lobo3 x.acc_onset_lobo3 x.acc_onset_lobo5];
        for k = 1:size(vals, 2)
            scatter(k + 0.15 * randn(height(x), 1), vals(:, k), 20, acol(x.animal, :), 'filled', 'MarkerFaceAlpha', 0.7);
            plot(k + [-0.3 0.3], median(vals(:, k), 'omitnan') * [1 1], 'k-', 'LineWidth', 2);
        end
        yline(0.5, 'k:');
        set(gca, 'XTick', 1:4, 'XTickLabel', {'0-20 cm', '20-140 cm', 'onset 0-0.5 s', 'onset (lobo5)'}, 'TickDir', 'out', 'Box', 'off');
        ylim([0.3 1]); ylabel('Track accuracy (CV)');
        title(sprintf('%s, %s (median %d cells)', regions{r}, decs{d}, round(median(x.n_cells))), 'Interpreter', 'none', 'FontWeight', 'normal');
    end
end
title(tl, 'Session QC: cross-validated track accuracy at the track start (dots = sessions, colour = animal)');
ar_save_figure(fig, fullfile(out_dir, 'session_qc'));
close(fig);
end

function a = track_acc(lo, true1, sel)
sel = sel & isfinite(lo);
a = mean((lo(sel) > 0) == true1(sel) & lo(sel) ~= 0);
end

function peak = rate_map_peaks(Ds, T, P)
% position (cm, bin centre) of the highest rate over both tracks' rate maps
tvec = Ds.tvec; nLaps = height(T); nPos = numel(P.pos_edges) - 1;
lap = zeros(size(tvec));
for k = 1:nLaps, lap(tvec >= T.t_on(k) & tvec < T.t_off(k)) = k; end
pb = discretize(Ds.position, P.pos_edges);
run = lap > 0 & ~isnan(pb) & Ds.speed > P.bayes_speed;
tr = zeros(size(tvec)); tr(lap > 0) = T.track(lap(lap > 0));
dt = [diff(tvec); median(diff(tvec))];
occ = accumarray([tr(run) pb(run)], dt(run), [2 nPos]);
[~, ci] = ismember(Ds.spike_id, Ds.cell_id);
ok = ci > 0 & Ds.spike_times >= tvec(1) & Ds.spike_times <= tvec(end);
smp = interp1(tvec, (1:numel(tvec))', Ds.spike_times(ok), 'nearest');
ci = ci(ok); keep = run(smp);
cnt = accumarray([ci(keep) tr(smp(keep)) pb(smp(keep))], 1, [numel(Ds.cell_id) 2 nPos]);
rate = cnt ./ reshape(occ, [1 2 nPos]);
rate = reshape(permute(rate, [1 3 2]), numel(Ds.cell_id), []);   % cells x (pos T1, pos T2)
[~, k] = max(rate, [], 2);
k = mod(k - 1, nPos) + 1;
peak = (k - 0.5) * (P.pos_edges(2) - P.pos_edges(1));
end
