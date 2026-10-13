function make_unique_ripple_data()
% Save PRE and POST ripple data restricted to non-duplicated (bilateral-merged)
% ripple events, in nature_neuroscience_revision_analysis:
%   KDE_reactivation_ripples_PSTH_POST_unique.mat  (rows = ripple_info.event_id)
%   KDE_reactivation_ripples_PSTH_PRE_unique.mat   (rows = ripple_info_PRE_unique.event_id)
%   ripple_info_PRE_unique.mat
% Ripples from both hemispheres within 0.05 s are merged with
% merge_bilateral_ripple_events (first one kept), as for POST ripple_info.

if exist('D:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'P:\corticohippocampal_replay';
end
rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
max_diff = 0.05;

%% POST: use existing ripple_info.event_id (index into the 64,194 SWS PSTH rows)
load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH.mat'), 'KDE_reactivation_ripples_PSTH');
KDE_reactivation_ripples_PSTH = select_sws_rows(KDE_reactivation_ripples_PSTH, ripple_info.event_id(:));
save(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH', '-v7.3');
fprintf('POST: saved %d unique events\n', numel(ripple_info.event_id));
clear KDE_reactivation_ripples_PSTH ripple_info

%% PRE: merge bilateral events within each session
load(fullfile(rev_dir, 'ripples_all_PRE.mat'), 'ripples_all');
load(fullfile(rev_dir, 'ripple_info_PRE.mat'), 'ripple_info_PRE');
load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_PRE.mat'), 'KDE_reactivation_ripples_PSTH_PRE');

% Row order of ripple_info_PRE / PRE PSTH: probe 1 SWS ripples, then probe 2
sws1 = ripples_all(1).SWS_index == 1;
sws2 = ripples_all(2).SWS_index == 1;
peaktimes = [ripples_all(1).peaktimes(sws1); ripples_all(2).peaktimes(sws2)];
hemisphere_id = [ones(sum(sws1), 1); 2 * ones(sum(sws2), 1)];
assert(numel(peaktimes) == numel(ripple_info_PRE.session_id), 'PRE ripples_all and ripple_info_PRE rows differ');

event_ids_first = [];
event_ids_second = [];
sessions = unique(ripple_info_PRE.session_id);
for i = 1:numel(sessions)
    rows = find(ripple_info_PRE.session_id == sessions(i));
    [f, s] = merge_bilateral_ripple_events(hemisphere_id(rows), peaktimes(rows), max_diff);
    event_ids_first = [event_ids_first; rows(f(:))]; %#ok<AGROW>
    event_ids_second = [event_ids_second; rows(s(:))]; %#ok<AGROW>
end

% ripple_info_PRE_unique: same convention as POST ripple_info
% (features from the first event; ripple power averaged over the pair)
info = ripple_info_PRE;
ripple_info_PRE = struct();
ripple_info_PRE.event_id = event_ids_first';                 % row index into full PRE PSTH
ripple_info_PRE.event_id_second = event_ids_second';
ripple_info_PRE.ripples_all_index = info.event_id(event_ids_first);
ripple_info_PRE.hemisphere_id = hemisphere_id(event_ids_first);
ripple_info_PRE.session_id = info.session_id(event_ids_first);
ripple_info_PRE.subject_id = info.subject_id(event_ids_first);
ripple_info_PRE.ripple_power = mean([info.ripple_power(event_ids_first) info.ripple_power(event_ids_second)], 2);
ripple_info_PRE.SO_phase = info.SO_phase(event_ids_first, :);
ripple_info_PRE.spindle_amplitude = info.spindle_amplitude(event_ids_first, :);
ripple_info_PRE.spindle_percentile = info.spindle_percentile(event_ids_first, :);
save(fullfile(rev_dir, 'ripple_info_PRE_unique.mat'), 'ripple_info_PRE', '-v7.3');

KDE_reactivation_ripples_PSTH_PRE = select_sws_rows(KDE_reactivation_ripples_PSTH_PRE, event_ids_first);
save(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_PRE_unique.mat'), 'KDE_reactivation_ripples_PSTH_PRE', '-v7.3');
fprintf('PRE: %d of %d events kept (%d bilateral pairs merged)\n', numel(event_ids_first), numel(info.session_id), ...
    sum(event_ids_first ~= event_ids_second));
end

function S = select_sws_rows(S, rows)
% Keep SWS ripple fields (and nan_mask) at the given rows; drop awake fields
fields = fieldnames(S);
for i = 1:numel(fields)
    f = fields{i};
    if contains(f, 'awake')
        S = rmfield(S, f);
    else
        S.(f) = S.(f)(rows, :);
    end
end
S.event_id = rows(:)';
end
