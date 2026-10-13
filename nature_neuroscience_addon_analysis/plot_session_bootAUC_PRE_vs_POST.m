%% Per-session bootstrapped reactivation coherence AUC: Low vs High and PRE vs POST
% Uses the per-session bootstrapped mean AUC from compare_session_bootAUC_vs_shuffleCI.
% Tests across sessions (session_paired_stats_plot), test_method:
%   'signrank' (default) - paired Wilcoxon signed-rank
%   'permutation'        - exact paired sign-flip permutation test on the mean difference
%   'ranksum'            - unpaired Wilcoxon rank-sum using all sessions in each group
%   'difference'         - main test: condition x sleep interaction, per-session
%                          (High - Low etc.) POST vs PRE (session_difference_stats_plot)
% Comparisons:
%   within sleep (PRE 12, POST 22 sessions): High vs Low, Trough vs Peak
%   PRE vs POST (paired tests: 12 sessions with both; rank-sum: all PRE vs all POST)
% Both two-sided and one-sided p are saved; plots show the one-sided p
% (High > Low, Trough > Peak, Unilateral > Bilateral Trough; POST > PRE for
% High/Trough/Unilateral/All, POST < PRE for Low/Peak/Bilateral).
% Output: nature_neuroscience_revision_analysis\session_bootAUC_PRE_vs_POST (signrank)
%         ...\session_bootAUC_PRE_vs_POST_permutation, ..._ranksum, ..._difference
function plot_session_bootAUC_PRE_vs_POST(test_method)
    if nargin < 1
        test_method = 'signrank';
    end
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    load(fullfile(rev_dir, 'session_AUC_vs_shuffle', 'session_bootAUC_vs_shuffleCI.mat'), 'T');
    out_dir = fullfile(rev_dir, 'session_bootAUC_PRE_vs_POST');
    if ~strcmp(test_method, 'signrank')
        out_dir = [out_dir '_' test_method];
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    if strcmp(test_method, 'difference')
        S = session_difference_stats_plot(T, 'AUC_boot_mean', 'AUC', out_dir, 'session_bootAUC_diff_');
    else
        S = session_paired_stats_plot(T, 'AUC_boot_mean', 'Reactivation coherence AUC (bootstrapped)', ...
            out_dir, 'session_bootAUC_', test_method);
    end
    save(fullfile(out_dir, 'session_bootAUC_stats.mat'), 'S');
    fprintf('Saved to %s\n', out_dir);
end
