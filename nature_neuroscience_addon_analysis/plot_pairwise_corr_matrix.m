function plot_pairwise_corr_matrix(cell_set, window, example_session)
% Visualisation of the V1-HC pair correlation matrix with V1 and HC clusters sorted by track
% preference (z_track_diff; Track L first). Shared preference predicts positive pair
% correlation for agreeing pairs (L-L, R-R blocks) and negative for disagreeing pairs
% (L-R, R-L blocks). Pairs and r from pairwise_V1_HC_coupling_POST('total_V1') (POST ripples,
% HC 0-0.1 s); r is shown minus the session mean over pairs (raw r has a positive offset).
%   left  : example session, every V1 x HC pair (strips: preference, blue Track L, red Track R)
%   right : all sessions, clusters binned into within-session preference quartiles, mean r
%           per bin averaged across sessions
% Rows: raw correlation, partial correlation (V1 and HC MUA removed).
% cell_set: 'all_cells' (default), 'lt10', 'rate_0p5_5'; window: 'V1_post' (default) or
% 'V1_pre'; example_session: default = session with the largest partial effect
% (rho_minus_null_partial in session_summary.csv)
% Output: pairwise_V1_HC_coupling\gain_total_V1[_<cell_set>]\pair_corr_matrix_<window>.fig/.pdf

if nargin < 1 || isempty(cell_set), cell_set = 'all_cells'; end
if nargin < 2 || isempty(window), window = 'V1_post'; end
rev_dir = addon_setup();
base = fullfile(rev_dir, 'pairwise_V1_HC_coupling');
[~, out_name, cells_txt] = pairwise_cell_set_folders('total_V1', cell_set);
out_dir = fullfile(base, out_name);
load(fullfile(out_dir, 'coupling_POST.mat'), 'S');
SS = readtable(fullfile(out_dir, 'session_summary.csv'));
if nargin < 3 || isempty(example_session)
    [~, i] = max(SS.(['rho_minus_null_partial_' window]));
    example_session = SS.session_id(i);
end
win_label = struct('V1_pre', 'V1 -0.2-0 s', 'V1_post', 'V1 0-0.2 s').(window);
meas = {'r_raw', 'r_partial'};
meas_lab = {'raw', 'partial (MUA removed)'};
nb = 4;

B = NaN(nb, nb, numel(S), 2);
for s = 1:numel(S)
    W = S(s).(window);
    zv = W.z_V1(:); zh = W.z_HC(:);
    bv = ceil(nb * tiedrank(-zv) / numel(zv));            % bin 1 = strongest Track L
    bh = ceil(nb * tiedrank(-zh) / numel(zh));
    for m = 1:2
        R = W.(meas{m}); R = R - mean(R(isfinite(R)));
        for i = 1:nb
            for j = 1:nb
                B(i, j, s, m) = mean(R(bv == i, bh == j), 'all', 'omitnan');
            end
        end
    end
    if S(s).session_count == example_session, ex = W; end
end
Bm = squeeze(mean(B, 3, 'omitnan'));

corr_map = diverging([0.37 0.23 0.60], [0.90 0.38 0.00], 256);     % purple (-) to orange (+)
pref_map = diverging([0.878 0.188 0.184], [0.173 0.373 0.859], 256); % Track R red to Track L blue
[~, ov] = sort(ex.z_V1(:), 'descend'); [~, oh] = sort(ex.z_HC(:), 'descend');
nL_v = sum(ex.z_V1 > 0); nL_h = sum(ex.z_HC > 0);
z_max = max(abs([ex.z_V1(:); ex.z_HC(:)]));

f = figure('Name', ['pair_corr_matrix_' window], 'Position', [60 40 1100 950], 'Color', 'w');
row_y = [0.52 0.06]; h = 0.33;
for m = 1:2
    % example session, all pairs
    R = ex.(meas{m}); R = R - mean(R(isfinite(R)));
    R = R(ov, oh);
    ax = axes('Position', [0.10 row_y(m) 0.36 h]);
    imagesc(ax, R); colormap(ax, corr_map);
    cl = prctile(abs(R(:)), 98); set(ax, 'CLim', [-cl cl]);
    hold(ax, 'on');
    xline(ax, nL_h + 0.5, 'k-', 'LineWidth', 1); yline(ax, nL_v + 0.5, 'k-', 'LineWidth', 1);
    set(ax, 'XTick', [], 'YTick', [], 'Box', 'on', 'FontSize', 9);
    xlabel(ax, sprintf('HC clusters (n = %d), Track L \\rightarrow R preferring', numel(oh)));
    ylabel(ax, sprintf('V1 clusters (n = %d), Track L \\rightarrow R', numel(ov)));
    title(ax, {sprintf('Example session %d: %s', example_session, meas_lab{m}), ''}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    cb = colorbar(ax, 'Position', [0.475 row_y(m) 0.01 h]);
    ylabel(cb, 'pair r - session mean');
    % preference strips
    as = axes('Position', [0.085 row_y(m) 0.012 h]);
    imagesc(as, ex.z_V1(ov)); colormap(as, pref_map); set(as, 'CLim', [-z_max z_max]); axis(as, 'off');
    as = axes('Position', [0.10 row_y(m) + h + 0.004 0.36 0.012]);
    imagesc(as, ex.z_HC(oh)'); colormap(as, pref_map); set(as, 'CLim', [-z_max z_max]); axis(as, 'off');

    % all sessions, preference quartiles
    ax = axes('Position', [0.60 row_y(m) 0.30 h]);
    M = Bm(:, :, m);
    imagesc(ax, M); colormap(ax, corr_map);
    cl = max(abs(M(:))); set(ax, 'CLim', [-cl cl]);
    for i = 1:nb
        for j = 1:nb
            text(ax, j, i, sprintf('%.3f', M(i, j)), 'HorizontalAlignment', 'center', 'FontSize', 8);
        end
    end
    tl = [{'Track L'}, repmat({''}, 1, nb - 2), {'Track R'}];
    set(ax, 'XTick', 1:nb, 'XTickLabel', tl, 'YTick', 1:nb, 'YTickLabel', tl, 'TickDir', 'out', 'FontSize', 9);
    axis(ax, 'square');
    xlabel(ax, 'HC cluster preference (quartile)'); ylabel(ax, 'V1 cluster preference (quartile)');
    title(ax, {sprintf('All sessions (n = %d): %s', numel(S), meas_lab{m}), 'mean r per quartile bin - session mean'}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    cb = colorbar(ax, 'Position', [0.92 row_y(m) 0.01 h]);
    ylabel(cb, 'pair r - session mean');
end
sgtitle(sprintf('V1-HC pair correlation sorted by track preference (POST ripples, %s, HC 0-0.1 s; %s)', ...
    win_label, cells_txt), 'FontSize', 11);
savefig(f, fullfile(out_dir, ['pair_corr_matrix_' window '.fig']));
exportgraphics(f, fullfile(out_dir, ['pair_corr_matrix_' window '.pdf']), 'ContentType', 'vector');
agree = @(M) mean(M([1 nb], [1 nb]) .* [1 -1; -1 1], 'all');
fprintf('%s example session %d; quartile corners (agree - disagree)/2: raw %.4f, partial %.4f\n', window, ...
    example_session, agree(Bm(:, :, 1)), agree(Bm(:, :, 2)));
fprintf('Saved %s\n', fullfile(out_dir, ['pair_corr_matrix_' window '.pdf']));
end

function map = diverging(c_low, c_high, n)
% Linear colormap c_low -> white -> c_high
t = linspace(0, 1, n / 2)';
map = [c_low + t * (1 - c_low); 1 - t * (1 - c_high)];
end
