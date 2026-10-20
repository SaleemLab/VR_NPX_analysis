function make_unique_awake_ripple_data()
% Save awake (non-SWS POST) ripple data restricted to non-duplicated
% (bilateral-merged) ripple events, in nature_neuroscience_revision_analysis:
%   ripple_info_awake_unique.mat
%   KDE_reactivation_ripples_PSTH_awake_unique.mat (rows = ripple_info_awake_unique.event_id)
% Ripples from both hemispheres within 0.05 s are merged with
% merge_bilateral_ripple_events (first one kept), as for PRE and POST.
% spindle_percentile = within-session percentile among awake ripples.
% +/-Inf log odds are capped with cap_log_odds_PSTH (percentiles over all awake
% rows, as for the awake fields of the capped POST PSTH).

if exist('D:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'P:\corticohippocampal_replay';
end
rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
max_diff = 0.05;

load(fullfile(rev_dir, 'ripples_all_POST.mat'), 'ripples_all');
load(fullfile(rev_dir, 'ripple_info_awake.mat'), 'ripple_info_awake');
load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_awake.mat'), 'KDE_reactivation_ripples_PSTH_awake');

% Row order of ripple_info_awake / awake PSTH: probe 1 non-SWS ripples, then probe 2
aw1 = ripples_all(1).SWS_index == 0;
aw2 = ripples_all(2).SWS_index == 0;
peaktimes = [ripples_all(1).peaktimes(aw1); ripples_all(2).peaktimes(aw2)];
hemisphere_id = [ones(sum(aw1), 1); 2 * ones(sum(aw2), 1)];
assert(numel(peaktimes) == numel(ripple_info_awake.session_id), 'ripples_all and ripple_info_awake rows differ');
assert(isequal(ripple_info_awake.event_id(:), [find(aw1); find(aw2)]), 'ripple_info_awake.event_id does not match non-SWS ripples');

event_ids_first = [];
event_ids_second = [];
sessions = unique(ripple_info_awake.session_id);
for i = 1:numel(sessions)
    rows = find(ripple_info_awake.session_id == sessions(i));
    [f, s] = merge_bilateral_ripple_events(hemisphere_id(rows), peaktimes(rows), max_diff);
    event_ids_first = [event_ids_first; rows(f(:))]; %#ok<AGROW>
    event_ids_second = [event_ids_second; rows(s(:))]; %#ok<AGROW>
end

% Same convention as POST ripple_info / ripple_info_PRE_unique
% (features from the first event; ripple power averaged over the pair)
info = ripple_info_awake;

% spindle_percentile of ripple_info_awake is all NaN: within-session tiedrank
% percentile of spindle amplitude per shank over all awake ripples, as
% add_spindle_percentile_PRE (0-100, no truncation)
info.spindle_percentile = NaN(size(info.spindle_amplitude));
for i = 1:numel(sessions)
    idx = info.session_id == sessions(i);
    amp = info.spindle_amplitude(idx, :);
    info.spindle_percentile(idx, :) = (tiedrank(amp) - 0.5) ./ size(amp, 1) * 100;
end

ripple_info_awake = struct();
ripple_info_awake.event_id = event_ids_first';                 % row index into full awake PSTH
ripple_info_awake.event_id_second = event_ids_second';
ripple_info_awake.ripples_all_index = info.event_id(event_ids_first);
ripple_info_awake.hemisphere_id = hemisphere_id(event_ids_first);
ripple_info_awake.session_id = info.session_id(event_ids_first);
ripple_info_awake.ripple_power = mean([info.ripple_power(event_ids_first) info.ripple_power(event_ids_second)], 2);
ripple_info_awake.SO_phase = info.SO_phase(event_ids_first, :);
ripple_info_awake.spindle_amplitude = info.spindle_amplitude(event_ids_first, :);
ripple_info_awake.spindle_percentile = info.spindle_percentile(event_ids_first, :);
save(fullfile(rev_dir, 'ripple_info_awake_unique.mat'), 'ripple_info_awake', '-v7.3');

KDE_reactivation_ripples_PSTH_awake = cap_log_odds_PSTH(KDE_reactivation_ripples_PSTH_awake);
fields = fieldnames(KDE_reactivation_ripples_PSTH_awake);
for i = 1:numel(fields)
    KDE_reactivation_ripples_PSTH_awake.(fields{i}) = KDE_reactivation_ripples_PSTH_awake.(fields{i})(event_ids_first, :);
end
KDE_reactivation_ripples_PSTH_awake.event_id = event_ids_first';
save(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_awake_unique.mat'), 'KDE_reactivation_ripples_PSTH_awake', '-v7.3');
fprintf('awake: %d of %d events kept (%d bilateral pairs merged)\n', numel(event_ids_first), numel(info.session_id), ...
    sum(event_ids_first ~= event_ids_second));
end
