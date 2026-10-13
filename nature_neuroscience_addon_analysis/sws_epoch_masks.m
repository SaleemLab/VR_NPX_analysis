function E = sws_epoch_masks(info_PRE, info_POST, analysis_folder, cap_min)
% Ripple event selection by cumulative SWS epoch for the non-duplicated PRE and
% POST events (ripple_info_PRE_unique / ripple_info).
%   all     : every SWS ripple
%   first30 : ripples before cap_min of cumulative SWS (all SWS if less than cap_min)
%   matched : sessions with both sleeps only, ripples before min(cap_min, PRE SWS,
%             POST SWS) of cumulative SWS in both sleeps
% Cumulative SWS from behavioural_state_merged_all_PRE/POST.mat (SWS{k} = [onset offset]
% of the k-th session of that sleep; k-th PRE session = k-th sorted PRE session_id,
% mapped_session_id 11-22). Event peak times from ripples_all_PRE/POST.mat (first
% ripple of a merged bilateral pair).
% E(1) = PRE, E(2) = POST, fields: name, sws_session_id, sws_bouts, sws_duration (s),
% peaktime (s, per event), sel.(epoch) (logical per event), used_min.(epoch) (per session).
rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
names = ["PRE", "POST"];
infos = {info_PRE, info_POST};

for iS = 1:2
    E(iS).name = names(iS); %#ok<AGROW>
    B = load(fullfile(analysis_folder, sprintf('behavioural_state_merged_all_%s.mat', names(iS))), ...
        'behavioural_state_merged_all');
    ids = unique(infos{iS}.session_id(:));
    assert(numel(B.behavioural_state_merged_all.SWS) == numel(ids), '%s SWS sessions do not match ripple_info', names(iS));
    E(iS).sws_session_id = ids;
    E(iS).sws_bouts = B.behavioural_state_merged_all.SWS(:);
    E(iS).sws_duration = cellfun(@(iv) sum(iv(:, 2) - iv(:, 1)), E(iS).sws_bouts);
end

% PRE: ripples_all_index = row within probe hemisphere_id of ripples_all_PRE
load(fullfile(rev_dir, 'ripples_all_PRE.mat'), 'ripples_all');
info = info_PRE;
E(1).peaktime = NaN(numel(info.session_id), 1);
for h = 1:2
    in_h = info.hemisphere_id == h;
    idx = info.ripples_all_index(in_h);
    assert(isequal(ripples_all(h).mapped_session_id(idx), info.session_id(in_h)), 'PRE ripples_all rows do not match');
    E(1).peaktime(in_h) = ripples_all(h).peaktimes(idx);
end
% POST: event_id = row of [probe 1 SWS ripples; probe 2 SWS ripples] of ripples_all_POST
load(fullfile(rev_dir, 'ripples_all_POST.mat'), 'ripples_all');
sws1 = ripples_all(1).SWS_index == 1; sws2 = ripples_all(2).SWS_index == 1;
pk = [ripples_all(1).peaktimes(sws1); ripples_all(2).peaktimes(sws2)];
sid = [ripples_all(1).session_count(sws1); ripples_all(2).session_count(sws2)];
assert(isequal(sid(info_POST.event_id(:)), info_POST.session_id(:)), 'POST ripples_all rows do not match');
E(2).peaktime = pk(info_POST.event_id(:));

paired_ids = intersect(E(1).sws_session_id, E(2).sws_session_id);
for iS = 1:2
    sess = infos{iS}.session_id(:);
    nS = numel(E(iS).sws_session_id);
    E(iS).sel = struct('all', true(size(sess)), 'first30', false(size(sess)), 'matched', false(size(sess)));
    E(iS).used_min = struct('all', E(iS).sws_duration / 60, 'first30', NaN(nS, 1), 'matched', NaN(nS, 1));
    for k = 1:nS
        s = E(iS).sws_session_id(k);
        in_s = sess == s;
        [t_cut, E(iS).used_min.first30(k)] = sws_cutoff(E(iS).sws_bouts{k}, cap_min * 60);
        E(iS).sel.first30(in_s) = E(iS).peaktime(in_s) <= t_cut;
        if ismember(s, paired_ids)
            cap = min([cap_min, E(1).sws_duration(E(1).sws_session_id == s) / 60, ...
                E(2).sws_duration(E(2).sws_session_id == s) / 60]);
            [t_cut, E(iS).used_min.matched(k)] = sws_cutoff(E(iS).sws_bouts{k}, cap * 60);
            E(iS).sel.matched(in_s) = E(iS).peaktime(in_s) <= t_cut;
        end
    end
end
end

function [t_cut, used_min] = sws_cutoff(bouts, cap_s)
% Time at which cumulative SWS (bouts [onset offset], s) reaches cap_s; Inf if never
bouts = sortrows(bouts, 1);
cum = cumsum(bouts(:, 2) - bouts(:, 1));
if cum(end) <= cap_s + 1e-6
    t_cut = Inf; used_min = cum(end) / 60;
else
    k = find(cum >= cap_s, 1);
    t_cut = bouts(k, 2) - (cum(k) - cap_s);
    used_min = cap_s / 60;
end
end
