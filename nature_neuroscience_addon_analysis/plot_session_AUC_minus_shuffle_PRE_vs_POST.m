%% Per-session reactivation coherence AUC minus the session's shuffle reference: PRE vs POST
% Uses the results of compare_session_AUC_vs_shuffle (1000 within-session
% shuffles per session/condition). Per session/condition the value is
% AUC - shuffle reference, so PRE vs POST compares chance-corrected AUC.
% reference (default 'mean'):
%   'mean' - AUC - shuffle mean
%   '95th' - AUC - shuffle 95th percentile (> 0 means above the session's threshold)
% Tests across sessions via session_paired_stats_plot (test_method as there:
% 'signrank' (default), 'permutation', 'ranksum'); paired PRE vs POST uses the
% 12 sessions with both. test_method 'difference': condition x sleep interaction,
% per-session (High - Low, Trough - Peak) POST vs PRE (session_difference_stats_plot).
% Output: nature_neuroscience_revision_analysis\session_AUC_minus_shuffle_PRE_vs_POST (mean, signrank)
%         ...\session_AUC_minus_shuffle95_PRE_vs_POST (95th), each with
%         _permutation / _ranksum / _difference suffix
function plot_session_AUC_minus_shuffle_PRE_vs_POST(test_method, reference)
    if nargin < 1
        test_method = 'signrank';
    end
    if nargin < 2
        reference = 'mean';
    end
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    load(fullfile(rev_dir, 'session_AUC_vs_shuffle', 'session_AUC_vs_shuffle.mat'), 'T');
    switch reference
        case 'mean'
            T.AUC_minus_shuffle = T.AUC - T.shuffle_mean;
            tag = ''; y_label = 'AUC - shuffle mean';
        case '95th'
            T.AUC_minus_shuffle = T.AUC - T.shuffle_95th;
            tag = '95'; y_label = 'AUC - shuffle 95th percentile';
    end
    out_dir = fullfile(rev_dir, ['session_AUC_minus_shuffle' tag '_PRE_vs_POST']);
    if ~strcmp(test_method, 'signrank')
        out_dir = [out_dir '_' test_method];
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    prefix = ['session_AUC_minus_shuffle' tag '_'];
    if strcmp(test_method, 'difference')
        prefix = [prefix 'diff_'];
        S = session_difference_stats_plot(T, 'AUC_minus_shuffle', y_label, out_dir, prefix);
    else
        S = session_paired_stats_plot(T, 'AUC_minus_shuffle', y_label, out_dir, prefix, test_method);
    end
    save(fullfile(out_dir, [prefix 'stats.mat']), 'S');
    writetable(T, fullfile(out_dir, [prefix 'values.csv']));
    fprintf('Saved to %s\n', out_dir);
end
