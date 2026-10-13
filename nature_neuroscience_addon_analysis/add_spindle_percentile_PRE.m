function add_spindle_percentile_PRE()
% Fill ripple_info_PRE.spindle_percentile (previously all NaN) with the
% within-session tiedrank percentile of spindle amplitude, per shank (as in POST).
% Only the copy in nature_neuroscience_revision_analysis is overwritten.

if exist('D:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'P:\corticohippocampal_replay';
end

f_info = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis', 'ripple_info_PRE.mat');
load(f_info, 'ripple_info_PRE');

spindle_pct = NaN(size(ripple_info_PRE.spindle_amplitude));
unique_sessions = unique(ripple_info_PRE.session_id);
for i = 1:length(unique_sessions)
    idx = ripple_info_PRE.session_id == unique_sessions(i);
    amp = ripple_info_PRE.spindle_amplitude(idx, :);
    spindle_pct(idx, :) = (tiedrank(amp) - 0.5) ./ size(amp, 1) * 100;
end
ripple_info_PRE.spindle_percentile = spindle_pct;

save(f_info, 'ripple_info_PRE', '-v7.3');
fprintf('Saved spindle_percentile to %s (NaN frac %.3f)\n', f_info, mean(isnan(spindle_pct(:))));
end
