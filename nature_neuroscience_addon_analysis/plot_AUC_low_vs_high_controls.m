function plot_AUC_low_vs_high_controls(which)
% Low vs high AUC bar plots (plot_AUC_low_vs_high style and colours) for the
% POST control analyses, from their saved static AUC .mat files:
%   'MUA_equalised'           : MUA_equalised_AUC, MUA-matched AUC and its shuffle
%                               (one row; ripple power, spindle, SO phase,
%                               dominant SO trough with non-dominant peak vs trough)
%   'split_half'              : split_half_V1_AUC_POST (event bootstrap)
%   'split_half_hierarchical' : split_half_V1_AUC_POST_hierarchical (session ->
%                               event bootstrap)
%   Split-half: one figure per dominance half -> selector half (Full->Full,
%   A->B, B->A, and Avg = mean of A->B and B->A where saved; file suffix e.g.
%   _A_to_B, _Avg_AB_BA), shared y axis; panels = ripple power
%   Q1 vs Q4 (if run), spindle Q1 vs Q4, SO peak vs trough, dominant SO trough with
%   non-dominant peak (unilateral) vs trough (bilateral), dominant SO peak with
%   non-dominant peak (bilateral) vs trough (unilateral).
% Bracket: second vs first condition, two-sided bootstrap p (2 x min of the
% saved one-sided p and its complement).
% which: cell array / string of the above (default all). Saves
% <prefix>_AUC_low_vs_high[_<direction>].fig/.png/.pdf in each analysis folder.
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    if nargin < 1, which = {'MUA_equalised', 'split_half', 'split_half_hierarchical'}; end
    which = cellstr(which);
    if ismember('MUA_equalised', which)
        plot_MUA_equalised(fullfile(rev_dir, 'MUA_equalised_event_control', 'AUC'));
    end
    if ismember('split_half', which)
        plot_split_half(fullfile(rev_dir, 'split_half_V1_AUC_POST'), 'split_half_AUC_static_', ...
            'split_half', 'POST split-half V1 (dominance half -> selector half)');
    end
    if ismember('split_half_hierarchical', which)
        plot_split_half(fullfile(rev_dir, 'split_half_V1_AUC_POST_hierarchical'), 'hier_AUC_static_', ...
            'hier', 'POST split-half V1, hierarchical bootstrap (dominance half -> selector half)');
    end
end

function plot_MUA_equalised(in_dir)
    L = load(fullfile(in_dir, 'MUA_equalised_AUC_static.mat'), 'res');
    style = AUC_low_high_style();
    names = {'ripple_power', 'spindle', 'SO_phase', 'SO_trough_uni_bi'};
    names = names(isfield(L.res, names));
    windows = {'pre', 'post'};
    R = repmat({empty_R()}, 1, numel(names)); P = cell(1, numel(names));
    for iA = 1:numel(names)
        for w = 1:numel(windows)
            r = L.res.(names{iA}).(windows{w});
            R{iA}(w) = struct('m', r.matched_mean, 'ci', r.matched_CI, 'ms', r.shuffled_mean, 'cis', r.shuffled_CI);
            P{iA}(w) = two_sided(r.diff_p);
        end
    end
    fig = figure('Name', 'MUA-equalised low vs high', 'Position', [200 200 1330 470]);
    tiledlayout(1, numel(names), 'TileSpacing', 'compact');
    [y_lim, yb] = shared_ylim(R);
    for iA = 1:numel(names)
        st = style.(names{iA});
        ax = nexttile;
        AUC_low_high_panel(ax, R{iA}, st.labels, st.colors, {'V1 -0.2-0 s', 'V1 0-0.2 s'}, y_lim, P{iA}, yb);
        title(ax, st.title, 'FontWeight', 'normal');
    end
    sgtitle({'POST MUA-equalised (within-session HC and V1 MUA matched)', ...
        'bracket: two-sided bootstrap p, second vs first condition'}, 'FontSize', 13);
    save_fig(fig, fullfile(in_dir, 'MUA_equalised_AUC_low_vs_high'));
end

function plot_split_half(in_dir, file_prefix, out_prefix, fig_title)
    style = AUC_low_high_style();
    % analysis, saved condition indices shown as (first, second) bar
    sets = struct('name', {'ripple_power', 'spindle', 'SO_phase', 'SO_trough_uni_bi', 'SO_peak_uni_bi'}, ...
        'conds', {[1 4], [1 4], [1 2], [2 1], [1 2]});   % trough_uni_bi saved as {bilateral, unilateral}
    % ripple power only where it was run (split_half_V1_AUC_POST)
    sets = sets(arrayfun(@(s) exist(fullfile(in_dir, [file_prefix s.name '.mat']), 'file') > 0, sets));
    windows = {'pre', 'post'};
    L = cell(1, numel(sets));
    for iA = 1:numel(sets)
        L{iA} = load(fullfile(in_dir, [file_prefix sets(iA).name '.mat']), 'res', 'dirs');
    end
    dirs = L{1}.dirs;
    R = repmat({empty_R()}, numel(dirs), numel(sets)); P = cell(numel(dirs), numel(sets));
    for d = 1:numel(dirs)
        dn = matlab.lang.makeValidName(dirs(d).name);
        for iA = 1:numel(sets)
            c = sets(iA).conds;
            for w = 1:numel(windows)
                r = L{iA}.res.(dn).(windows{w});
                R{d, iA}(w) = struct('m', r.AUC_mean(c), 'ci', r.AUC_CI(:, c), ...
                    'ms', r.AUC_shuffled_mean(c), 'cis', r.AUC_shuffled_CI(:, c));
                P{d, iA}(w) = two_sided(r.diff_p);   % saved: last - first condition
            end
        end
    end
    % One figure per direction, same y axis for all
    [y_lim, yb] = shared_ylim(R);
    for d = 1:numel(dirs)
        % width scales with the panel count so each panel keeps the 4-panel size
        fig = figure('Name', [out_prefix ' low vs high ' dirs(d).name], 'Position', [200 200 1330 / 4 * numel(sets) 470]);
        tiledlayout(1, numel(sets), 'TileSpacing', 'compact');
        for iA = 1:numel(sets)
            st = style.(sets(iA).name);
            ax = nexttile;
            AUC_low_high_panel(ax, R{d, iA}, st.labels, st.colors, {'V1 -0.2-0 s', 'V1 0-0.2 s'}, y_lim, P{d, iA}, yb);
            title(ax, st.title, 'FontWeight', 'normal');
        end
        if strcmp(dirs(d).name, 'Avg')   % 2-fold cross-validated mean of A->B and B->A
            dir_label = 'mean of A->B and B->A (per bootstrap resample)'; suffix = 'Avg_AB_BA';
        else
            dir_label = dirs(d).name; suffix = strrep(dirs(d).name, '->', '_to_');
        end
        sgtitle({sprintf('%s: %s', fig_title, dir_label), ...
            'bracket: two-sided bootstrap p, second vs first condition'}, 'FontSize', 13);
        save_fig(fig, fullfile(in_dir, [out_prefix '_AUC_low_vs_high_' suffix]));
    end
end

function R = empty_R()
    R = struct('m', {}, 'ci', {}, 'ms', {}, 'cis', {});
end

function p2 = two_sided(p)
    p2 = min(1, 2 * min(p, 1 - p));
end

function [y_lim, yb] = shared_ylim(R)
    % Brackets on one row 0.12 above the highest condition CI of the figure
    % (clears the '*' / 'n.s.' markers); top leaves room for p and window labels
    hi = max(cellfun(@(x) max([x.ci], [], 'all'), R(:)));
    yb = hi + 0.12;
    y_lim = [-0.2 max(0.44, yb + 0.14)];
end

function save_fig(fig, save_path)
    savefig(fig, [save_path '.fig']);
    exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
    exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
end
