function plot_Lseq_templates_session(subject, session, varargin)
% Template matrices of one session from run_Lseq_RUN_decoding_CV output (whole-session
% templates saved in out.region(nr).diagnostics). One figure per region:
%   rows    : cells sorted by Track 1 rate-map peak (top) and Track 2 peak (bottom)
%   col 1-4 : W_T1, W_T2, W_shared = (W_T1 + W_T2)/2, W_specific = (W_T1 - W_T2)/2
%             W(a,b) > 0 = cell a fires before cell b, so a forward sequence in the
%             sorting order is red above the diagonal. Pairs with < min_pairs ordered
%             spike pairs in either track template are grey.
%   col 5   : rate-map peak on T1 vs T2 per cell (top); W_T1 and W_T2 against the
%             peak distance on their own track (bottom)
%   col 6   : W_shared against the mean peak distance (dT1 + dT2)/2 (top) and
%             W_specific against half the difference (dT1 - dT2)/2 (bottom); mean +- s.e.m.
%             per distance (10 cm bins)
% 'scale', 'P' draws the matrices as probabilities instead, as Huh et al. 2026 Fig. 1f:
%   P = 1 / (1 + exp(-W)) on [0 1], cell 1 at the bottom, pixel (column j, row i) =
%   P(row cell i fires after column cell j), so a forward sequence is high above the
%   diagonal (upper left). For shared and specific this is the probability implied by
%   that component's W alone.
% Example: plot_Lseq_templates_session('M24017', '20240603', 'scale', 'P')
% Options: 'regions' (all in the file), 'min_pairs' (as used by the decoding), 'scale'
% ('W' or 'P'), 'out_dir'.

p = inputParser;
addParameter(p, 'regions', {}, @iscell);
addParameter(p, 'scale', 'W', @(x) any(strcmp(x, {'W', 'P'})));
addParameter(p, 'min_pairs', [], @isnumeric);
addParameter(p, 'out_dir', '', @ischar);
parse(p, varargin{:});
P = p.Results;
if isempty(P.out_dir), P.out_dir = fullfile(get_revision_folder(), 'Lseq_RUN_decoding'); end

S = load(fullfile(P.out_dir, subject, sprintf('session_%s.mat', session)), 'out');
out = S.out;
if isempty(P.min_pairs), P.min_pairs = out.params.min_pairs; end
regions = {out.region.name};
if ~isempty(P.regions), regions = intersect(regions, P.regions, 'stable'); end

for name = regions
    R = out.region(strcmp({out.region.name}, name{1}));
    d = R.diagnostics;
    ok = d.N_T1 >= P.min_pairs & d.N_T2 >= P.min_pairs;
    ok(logical(eye(size(ok)))) = false;
    W = {double(d.W_T1), double(d.W_T2), double(d.W_T1 + d.W_T2) / 2, double(d.W_T1 - d.W_T2) / 2};
    for k = 1:4, W{k}(~ok) = NaN; end
    wmax = prctile(abs(W{1}(ok)), 99);
    peak = d.peak_bin;
    titles = {'W Track 1', 'W Track 2', 'W shared (T1+T2)/2', 'W specific (T1-T2)/2'};

    fig = figure('Name', name{1}, 'Position', [30 30 2100 720], 'Color', 'w');
    tl = tiledlayout(2, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
    for row = 1:2
        [~, order] = sort(peak(:, row));   % NaN peaks (silent on that track) last
        n_sorted = sum(~isnan(peak(:, row)));
        for k = 1:4
            nexttile(6 * (row - 1) + k);
            Wk = W{k}(order, order);
            hold on;
            if strcmp(P.scale, 'W')
                imagesc(Wk, 'AlphaData', ~isnan(Wk), [-wmax wmax]);
                set(gca, 'YDir', 'reverse'); colormap(gca, diverging_map());
                cb_label = 'W (log odds a before b)'; yl = 'cell a'; xl = 'cell b';
            else
                Pk = 1 ./ (1 + exp(-Wk'));   % (i,j) = P(j before i) = P(row i after column j)
                imagesc(Pk, 'AlphaData', ~isnan(Pk), [0 1]);
                set(gca, 'YDir', 'normal'); colormap(gca, parula);
                cb_label = 'P(row cell fires after column cell)'; yl = 'cell'; xl = 'cell';
            end
            set(gca, 'Color', [0.8 0.8 0.8], 'TickDir', 'out'); axis square tight;
            plot([0.5 n_sorted + 0.5], [n_sorted + 0.5 n_sorted + 0.5], 'k-', [n_sorted + 0.5 n_sorted + 0.5], [0.5 n_sorted + 0.5], 'k-');
            if k == 4, cb = colorbar; cb.Label.String = cb_label; end
            title(strrep(titles{k}, 'W ', [P.scale ' ']), 'FontWeight', 'normal');
            if k == 1
                ylabel(sprintf('%s (sorted by T%d peak)', yl, row));
            end
            xlabel(xl);
        end
    end

    % Rate-map peaks on the two tracks
    nexttile(5); hold on;
    both = all(~isnan(peak), 2);
    jit = 0.15 * randn(sum(both), 2);
    scatter(peak(both, 1) + jit(:, 1), peak(both, 2) + jit(:, 2), 12, 'k', 'filled', 'MarkerFaceAlpha', 0.4);
    lim = [0.5 numel(out.params.pos_edges) - 0.5];
    plot(lim, lim, 'r--'); axis square; xlim(lim); ylim(lim);
    xlabel('T1 peak (10 cm bin)'); ylabel('T2 peak (10 cm bin)');
    title(sprintf('peaks: r = %.2f, order agreement %.2f', corr(peak(both, 1), peak(both, 2)), d.peak_order_agreement), 'FontWeight', 'normal');
    set(gca, 'TickDir', 'out');

    % W against peak distances (pair a,b: delta = peak(b) - peak(a))
    U = triu(ok, 1) & both & both';
    [a, b] = find(U);
    dT1 = peak(b, 1) - peak(a, 1); dT2 = peak(b, 2) - peak(a, 2);
    idx = sub2ind(size(U), a, b);

    nexttile(11); hold on;
    h1 = profile_plot(dT1, W{1}(idx), [0.8 0.2 0.2]);
    h2 = profile_plot(dT2, W{2}(idx), [0.2 0.4 0.8]);
    legend([h1 h2], {'W_{T1} vs \Delta peak T1', 'W_{T2} vs \Delta peak T2'}, 'Location', 'northwest', 'Box', 'off');
    xlabel('\Delta peak, b - a (10 cm bins)'); ylabel('mean W');
    title('own-track templates', 'FontWeight', 'normal'); set(gca, 'TickDir', 'out');

    nexttile(6); hold on;
    profile_plot((dT1 + dT2) / 2, W{3}(idx), [0.45 0.45 0.45]);
    xlabel('(\Delta T1 + \Delta T2) / 2 (10 cm bins)'); ylabel('mean W shared');
    title(sprintf('shared vs mean distance (r = %.2f)', corr((dT1 + dT2) / 2, W{3}(idx))), 'FontWeight', 'normal');
    set(gca, 'TickDir', 'out');

    nexttile(12); hold on;
    profile_plot((dT1 - dT2) / 2, W{4}(idx), [0.1 0.1 0.1]);
    xlabel('(\Delta T1 - \Delta T2) / 2 (10 cm bins)'); ylabel('mean W specific');
    title(sprintf('specific vs distance difference (r = %.2f)', corr((dT1 - dT2) / 2, W{4}(idx))), 'FontWeight', 'normal');
    set(gca, 'TickDir', 'out');

    title(tl, sprintf('%s %s %s (%d cells, %d pairs with >= %d spike pairs on both tracks): whole-session L_{seq} templates', ...
        subject, session, name{1}, R.n_cells, sum(triu(ok, 1), 'all'), P.min_pairs), 'Interpreter', 'tex');
    suffix = ''; if strcmp(P.scale, 'P'), suffix = '_P'; end
    save_figure(fig, fullfile(P.out_dir, subject, sprintf('session_%s_templates_%s%s', session, name{1}, suffix)));
    close(fig);
end
end

function h = profile_plot(x, y, col)
% Mean +- s.e.m. of y per value of x (x on a 0.5-bin grid), bins with >= 10 pairs
x = round(2 * x) / 2;
xs = unique(x);
m = NaN(size(xs)); e = m;
for k = 1:numel(xs)
    v = y(x == xs(k));
    if numel(v) < 10, continue; end
    m(k) = mean(v); e(k) = std(v) / sqrt(numel(v));
end
keep = ~isnan(m); xs = xs(keep); m = m(keep); e = e(keep);
plot(xs([1 end]), [0 0], 'k:');
plot([0 0], [-1 1] * max(abs(m) + e, [], 'omitnan'), 'k:');
fill([xs; flipud(xs)], [m - e; flipud(m + e)], col, 'EdgeColor', 'none', 'FaceAlpha', 0.25, 'HandleVisibility', 'off');
h = plot(xs, m, '-', 'Color', col, 'LineWidth', 1.5);
end

function cmap = diverging_map()
n = 128;
up = [linspace(0.13, 1, n)' linspace(0.4, 1, n)' linspace(0.67, 1, n)'];
dn = [linspace(1, 0.7, n)' linspace(1, 0.09, n)' linspace(1, 0.17, n)'];
cmap = [up; dn];
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
