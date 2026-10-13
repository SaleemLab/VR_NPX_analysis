%% Per-session reactivation coherence AUC vs the session's shuffle mean
% Uses the results of compare_session_AUC_vs_shuffle (1000 within-session
% shuffles per session/condition). Paired Wilcoxon signed-rank across sessions:
% actual AUC vs shuffle mean; two-sided and one-sided (actual > shuffle mean) p
% saved, plots show one-sided. Rows = sleep (PRE, POST), columns = V1 window.
% Output: nature_neuroscience_revision_analysis\session_AUC_vs_shuffle\session_AUC_vs_shuffle_mean_*
function plot_session_AUC_vs_shuffle_mean()
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    out_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis', 'session_AUC_vs_shuffle');
    load(fullfile(out_dir, 'session_AUC_vs_shuffle.mat'), 'T');
    session_vs_shuffle_stats_plot(T, 'AUC', 'shuffle_mean', 'Shuffle mean', out_dir, 'session_AUC_vs_shuffle_mean_');
    fprintf('Saved to %s\n', out_dir);
end
