function plot_temporal_AUC_low_vs_high(which)
% Temporal version of plot_AUC_low_vs_high: one figure per analysis folder with
% four panels (ripple power Q1 vs Q4, spindle Q1 vs Q4, SO peak vs trough, dominant
% SO trough with non-dominant peak (unilateral) vs trough (bilateral)), each the
% two contrasting conditions' temporal reactivation coherence AUC (sliding 0.1 s
% V1 window, HPC 0-0.1 s) with 95% bootstrap CI, over the shuffle (grey; mean of
% the two conditions' shuffled means, CI pooled as min lower / max upper, as in
% AUC_low_high_panel). Colours and labels: AUC_low_high_style.
% x = V1 window centre; analyses saved with the window starting at each bin
% (V1->HPC of the hierarchical / split-half runs) are shifted by half a window.
% Bars above the traces: bins where the condition's 95% CI does not overlap the
% shuffled 95% CI (pale: first, full: second condition; the '*' criterion of the
% bar plot) and, where the bootstrap difference was saved, bins with two-sided
% p < 0.05 for second vs first condition (black).
% which: folder(s) of reactivation_AUC_coherence_PRE (e.g. 'PRE', 'POST',
%   'POST_matched_equal_n', 'POST_Bayes', 'awake'), 'POST_original_AUC_hierarchical',
%   'MUA_equalised', 'MUA_equalised_with_nondominant', 'split_half' or
%   'split_half_hierarchical' (one figure per dominance half -> selector half).
%   Default: all.
% Saves <prefix>temporal_AUC_low_vs_high[_<direction>].fig/.png/.pdf in each folder.
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    if nargin < 1
        which = {'PRE', 'PRE_equal_n', 'POST_matched', 'POST_matched_equal_n', 'POST', 'POST_rate_lt10', ...
            'POST_Bayes', 'POST_Bayes_rate_lt10', 'POST_Bayes_raw_rate_lt10', 'awake', ...
            'POST_original_AUC_hierarchical', 'MUA_equalised', 'MUA_equalised_with_nondominant', ...
            'split_half', 'split_half_hierarchical'};
    end
    which = cellstr(which);
    names = {'ripple_power', 'spindle', 'SO_phase', 'SO_trough_uni_bi'};
    for k = 1:numel(which)
        switch which{k}
            case 'POST_original_AUC_hierarchical'
                in_dir = fullfile(rev_dir, which{k});
                % hierarchical: HPC->V1 window centred, V1->HPC starting at each bin
                files = strcat('hier_AUC_temporal_', names);
                conds = {[1 4], [1 4], [1 2], [2 1]};   % trough_uni_bi saved as {bilateral, unilateral}
                D = load_set(in_dir, files, conds, @(L) L.AUC, hier_fields());
                plot_figure(D, names, 'POST original decoder, hierarchical (session -> event) bootstrap', ...
                    fullfile(in_dir, 'hier_temporal_AUC_low_vs_high'));
            case {'MUA_equalised', 'MUA_equalised_with_nondominant'}
                if strcmp(which{k}, 'MUA_equalised')
                    v = MUA_equalised_variant('dominant');
                else
                    v = MUA_equalised_variant('with_nondominant');
                end
                in_dir = fullfile(rev_dir, v.folder, 'AUC');
                files = strcat('MUA_equalised_AUC_temporal_', names);
                f = struct('m', 'matched_mean', 'ci', 'matched_CI', 'ms', 'shuffled_mean', 'cis', 'shuffled_CI', 'p', 'diff_p');
                D = load_set(in_dir, files, repmat({[1 2]}, 1, 4), @(L) L.AUC, f);
                if all(cellfun(@isempty, D)), continue; end
                plot_figure(D, names, ['POST MUA-equalised (' v.label ')'], ...
                    fullfile(in_dir, 'MUA_equalised_temporal_AUC_low_vs_high'));
            case {'split_half', 'split_half_hierarchical'}
                if strcmp(which{k}, 'split_half')
                    in_dir = fullfile(rev_dir, 'split_half_V1_AUC_POST'); prefix = 'split_half';
                    fig_title = 'POST split-half V1';
                else
                    in_dir = fullfile(rev_dir, 'split_half_V1_AUC_POST_hierarchical'); prefix = 'hier';
                    fig_title = 'POST split-half V1, hierarchical bootstrap';
                end
                files = strcat([prefix '_AUC_temporal_'], names);
                conds = {[1 4], [1 4], [1 2], [2 1]};
                first = find(cellfun(@(f) exist(fullfile(in_dir, [f '.mat']), 'file') > 0, files), 1);
                dirs = getfield(load(fullfile(in_dir, [files{first} '.mat']), 'dirs'), 'dirs'); %#ok<GFLD>
                for d = 1:numel(dirs)
                    dn = matlab.lang.makeValidName(dirs(d).name);
                    D = load_set(in_dir, files, conds, @(L) L.AUC.(dn), hier_fields());
                    plot_figure(D, names, sprintf('%s: %s (dominance half -> selector half)', fig_title, dirs(d).name), ...
                        fullfile(in_dir, [prefix '_temporal_AUC_low_vs_high_' strrep(dirs(d).name, '->', '_to_')]));
                end
            otherwise   % reactivation_AUC_coherence_PRE folder (event bootstrap)
                in_dir = fullfile(rev_dir, which{k});
                files = strcat('KDE_temporal_AUC_', {'ripple_percentile_', 'spindle_percentile_', 'SO_phase_', ...
                    'SO_trough_bilateral_'}, which{k});
                conds = {[1 4], [1 4], [1 2], [1 2]};   % SO_trough_bilateral saved as {unilateral, bilateral}
                f = struct('m', 'mean', 'ci', 'ci', 'ms', 'shifted_mean', 'cis', 'shifted_ci', 'p', '');
                D = load_set(in_dir, files, conds, @(L) L.AUC, f);
                plot_figure(D, names, strrep(which{k}, '_', ' '), fullfile(in_dir, 'temporal_AUC_low_vs_high'));
        end
    end
end

function f = hier_fields()
    f = struct('m', 'mean', 'ci', 'ci', 'ms', 'shuffled_mean', 'cis', 'shuffled_ci', 'p', 'diff_p');
end

function D = load_set(in_dir, files, conds, get_AUC, f)
% D{i}: t (window centres), m / lo / hi (nT x 2, conditions in conds order),
% ms / slo / shi (pooled shuffle, nT x 1), p (two-sided, nT x 1; [] if not saved)
    D = cell(1, numel(files));
    for i = 1:numel(files)
        fn = fullfile(in_dir, [files{i} '.mat']);
        if ~exist(fn, 'file')
            warning('plot_temporal_AUC_low_vs_high:missing', 'Missing %s', fn);
            continue
        end
        A = get_AUC(load(fn));
        c = conds{i};
        fprintf('%s: %s vs %s\n', files{i}, A.labels{c(1)}, A.labels{c(2)});
        if isfield(A, 'window_alignment')
            centred = strcmp(A.window_alignment, 'centred');
        else
            centred = strcmp(A.direction, 'HPC->V1');   % hierarchical_AUC_POST_original
        end
        t = A.time_bins(:) + ~centred * A.win_size / 2;
        keep = t >= -0.5 - 1e-9 & t <= 0.5 + 1e-9 & ~any(isnan(A.(f.m)(:, c)), 2);
        R.t = t(keep);
        R.m = A.(f.m)(keep, c);
        R.lo = A.(f.ci)(keep, c, 1); R.hi = A.(f.ci)(keep, c, 2);
        R.ms = mean(A.(f.ms)(keep, c), 2);
        R.slo = min(A.(f.cis)(keep, c, 1), [], 2); R.shi = max(A.(f.cis)(keep, c, 2), [], 2);
        R.p = [];
        if ~isempty(f.p)
            p = A.(f.p)(keep);
            R.p = min(1, 2 * min(p, 1 - p));   % saved one-sided, last - first condition
        end
        D{i} = R;
    end
end

function plot_figure(D, names, fig_title, save_path)
    style = AUC_low_high_style();
    ok = ~cellfun(@isempty, D);
    % shared y axis; room above the data for the significance bars
    lo = min(cellfun(@(R) min([R.lo(:); R.slo]), D(ok)));
    hi = max(cellfun(@(R) max([R.hi(:); R.shi]), D(ok)));
    y_lo = floor((lo - 0.02) * 20) / 20;
    bar_y = hi + 0.03 + (0:2) * 0.025;   % rows: first, second condition vs shuffle, second vs first
    y_lim = [y_lo, bar_y(3) + 0.03];

    fig = figure('Name', fig_title, 'Position', [200 200 1330 440]);
    tiledlayout(1, numel(names), 'TileSpacing', 'compact');
    for i = 1:numel(names)
        ax = nexttile; hold(ax, 'on');
        st = style.(names{i});
        if ~ok(i)
            title(ax, [cellstr(st.title), {'(not saved)'}], 'FontWeight', 'normal'); axis(ax, 'off');
            continue
        end
        R = D{i};
        t = R.t';
        h = gobjects(1, 3);
        fill(ax, [t fliplr(t)], [R.slo' fliplr(R.shi')], [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.35);
        h(3) = plot(ax, t, R.ms, 'k', 'LineWidth', 1);
        for c = 1:2
            fill(ax, [t fliplr(t)], [R.lo(:, c)' fliplr(R.hi(:, c)')], st.colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        end
        for c = 1:2
            h(c) = plot(ax, t, R.m(:, c), 'Color', st.colors(c, :) * (1 - 0.15 * (c == 1)), 'LineWidth', 2);
            sig = R.lo(:, c) > R.shi | R.hi(:, c) < R.slo;
            sig_bar(ax, R.t, sig, bar_y(c), st.colors(c, :));
        end
        if ~isempty(R.p)
            sig_bar(ax, R.t, R.p < 0.05, bar_y(3), [0 0 0]);
        end
        xline(ax, 0, 'k--', 'LineWidth', 1.5);
        yline(ax, 0, 'k:');
        xlim(ax, [-0.5 0.5]); ylim(ax, y_lim);
        set(ax, 'XTick', [-0.5 -0.2 0 0.2 0.5], 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12);
        xlabel(ax, 'V1 window centre (s)');
        if i == 1, ylabel(ax, 'Reactivation coherence AUC'); end
        title(ax, st.title, 'FontWeight', 'normal');
        legend(ax, h, [st.labels, {'Shuffle'}], 'Location', 'southoutside', 'Orientation', 'horizontal', ...
            'Box', 'off', 'FontSize', 10);
    end
    sgtitle({fig_title, ['bars: 95% CI outside shuffle 95% CI (pale: first, full: second condition); ' ...
        'black: second vs first, two-sided bootstrap p < 0.05']}, 'FontSize', 12);
    savefig(fig, [save_path '.fig']);
    exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
    exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
    close(fig);
end

function sig_bar(ax, t, sig, y, color)
% Horizontal bar over each run of significant bins (bin width = time step)
    if ~any(sig), return; end
    dt = median(diff(t));
    d = diff([0; sig(:); 0]);
    on = find(d == 1); off = find(d == -1) - 1;
    for r = 1:numel(on)
        plot(ax, [t(on(r)) - dt/2, t(off(r)) + dt/2], [y y], '-', 'Color', color, 'LineWidth', 4);
    end
end
