function cap_KDE_reactivation_ripples_PSTH_POST()
% Cap +/-Inf log odds in POST ripple PSTH and save to the revision folder.
% Only the copy in nature_neuroscience_revision_analysis is overwritten;
% the original in 'V1-HPC sleep reactivation' is left untouched.

if exist('D:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'P:\corticohippocampal_replay';
end

load(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'KDE_reactivation_ripples_PSTH.mat'), 'KDE_reactivation_ripples_PSTH');

KDE_reactivation_ripples_PSTH = cap_log_odds_PSTH(KDE_reactivation_ripples_PSTH);

out_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
save(fullfile(out_dir, 'KDE_reactivation_ripples_PSTH.mat'), 'KDE_reactivation_ripples_PSTH', '-v7.3');
fprintf('Saved capped KDE_reactivation_ripples_PSTH.mat to %s\n', out_dir);
end
