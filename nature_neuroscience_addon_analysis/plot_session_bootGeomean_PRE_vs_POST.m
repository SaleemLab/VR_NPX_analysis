%% Per-session bootstrapped geomean coherence: Low vs High and PRE vs POST
% Uses the per-session bootstrapped mean geomean coherence from compare_session_bootGeomean.
% Tests across sessions, test_method (as plot_session_bootAUC_PRE_vs_POST):
%   'signrank' (default) - paired Wilcoxon signed-rank
%   'permutation'        - exact paired sign-flip permutation test on the mean difference
%   'ranksum'            - unpaired Wilcoxon rank-sum using all sessions in each group
%   'difference'         - main test: condition x sleep interaction, per-session
%                          (High - Low etc.) POST vs PRE (session_difference_stats_plot)
% Paired tests also save matched-session PRE vs POST plots.
% Output: nature_neuroscience_revision_analysis\session_bootGeomean_PRE_vs_POST (signrank)
%         ...\session_bootGeomean_PRE_vs_POST_permutation, ..._ranksum, ..._difference
function plot_session_bootGeomean_PRE_vs_POST(test_method)
    if nargin < 1
        test_method = 'signrank';
    end
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    load(fullfile(rev_dir, 'session_bootGeomean', 'session_bootGeomean.mat'), 'T');
    out_dir = fullfile(rev_dir, 'session_bootGeomean_PRE_vs_POST');
    if ~strcmp(test_method, 'signrank')
        out_dir = [out_dir '_' test_method];
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    if strcmp(test_method, 'difference')
        S = session_difference_stats_plot(T, 'geomean_boot_mean', 'Geomean', out_dir, 'session_bootGeomean_diff_');
    else
        S = session_paired_stats_plot(T, 'geomean_boot_mean', 'Geomean coherence (bootstrapped session mean)', ...
            out_dir, 'session_bootGeomean_', test_method);
    end
    save(fullfile(out_dir, 'session_bootGeomean_stats.mat'), 'S');
    fprintf('Saved to %s\n', out_dir);
end
