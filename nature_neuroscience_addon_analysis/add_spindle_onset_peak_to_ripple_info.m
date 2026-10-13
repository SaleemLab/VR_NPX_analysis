function add_spindle_onset_peak_to_ripple_info()
% Add spindle-band amplitude at ripple onset and at ripple peak (raw and
% within-session percentile) to the non-duplicated POST ripple_info in
% nature_neuroscience_revision_analysis, as in reactivation_coherence_KDE_glme.m
% (onset) and corticohippocampal_reactivation_mean_KDE_diff_V1_detection.m (peak):
% reference V1 shank per hemisphere from the best SO channel, percentile ranked
% within each probe-session over all ripples, then SWS ripples, then unique
% events (ripple_info.event_id). Columns = V1 hemisphere 1, 2.
% Adds: spindle_amplitude_onset, spindle_percentile_onset,
%       spindle_amplitude_peak,  spindle_percentile_peak
% Existing fields are left unchanged; a backup is saved before overwriting.

if exist('D:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'P:\corticohippocampal_replay';
end
rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
f_info = fullfile(rev_dir, 'ripple_info.mat');

load(fullfile(analysis_folder, 'slow_waves_all_POST.mat'), 'slow_waves_all');
load(fullfile(rev_dir, 'ripples_all_POST.mat'), 'ripples_all');
load(f_info, 'ripple_info');

% Reference V1 shank per session and hemisphere (best SO channel)
nSessions = max(ripples_all(1).session_count);
cortex_ref_shank = NaN(nSessions, 2);
for nsession = 1:nSessions
    for probe_no = 1:2
        cortex_ref_shank(nsession, probe_no) = find(slow_waves_all(probe_no).shank_id{nsession} == ...
            slow_waves_all(probe_no).shank{nsession}(slow_waves_all(probe_no).channel{nsession} == slow_waves_all(probe_no).best_channel(nsession)) ...
            & slow_waves_all(probe_no).probe_hemisphere{nsession} == probe_no);
    end
end

n_sws = sum(ripples_all(1).SWS_index == 1) + sum(ripples_all(2).SWS_index == 1);
assert(max(ripple_info.event_id) <= n_sws, 'ripple_info.event_id exceeds SWS ripple count');

fields = {'onset', 'spindle_amplitude_ripple_onset'; 'peak', 'spindle_amplitude_ripple_peaktime'};
for k = 1:size(fields, 1)
    amp = cell(1, 2); pct = cell(1, 2);
    for probe_no = 1:2
        for nsession = 1:nSessions
            temp = ripples_all(probe_no).(fields{k, 2}){nsession}(cortex_ref_shank(nsession, :), :);
            amp{probe_no} = [amp{probe_no} temp];
            pct{probe_no} = [pct{probe_no} (tiedrank(temp')' - 0.5) / length(temp) * 100];
        end
    end
    sws1 = ripples_all(1).SWS_index == 1; sws2 = ripples_all(2).SWS_index == 1;
    amp = [amp{1}(:, sws1) amp{2}(:, sws2)]';
    pct = [pct{1}(:, sws1) pct{2}(:, sws2)]';
    ripple_info.(['spindle_amplitude_' fields{k, 1}]) = amp(ripple_info.event_id, :);
    ripple_info.(['spindle_percentile_' fields{k, 1}]) = pct(ripple_info.event_id, :);
    fprintf('%s: max |diff| vs existing spindle_percentile = %.4g\n', fields{k, 1}, ...
        max(abs(ripple_info.(['spindle_percentile_' fields{k, 1}])(:) - ripple_info.spindle_percentile(:))));
end

copyfile(f_info, fullfile(rev_dir, 'ripple_info_before_spindle_onset_peak.mat'));
save(f_info, 'ripple_info', '-v7.3');
fprintf('Saved %s\n', f_info);
end
