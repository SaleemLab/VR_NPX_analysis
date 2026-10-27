function Cp = ar_cell_properties(Ds, T, P)
% Per-cell properties of one RUN1 session (cells of Ds.cell_id), from in-lap
% running samples (speed > P.bayes_speed, all laps; ar_bayes_lap_stats):
%   region, r_T1, r_T2 (mean in-lap rate per track, Hz)
%   pref        (r_T1 - r_T2) / (r_T1 + r_T2)
%   pref_start  the same over 0-20 cm
%   peak_pos    peak of the rate map on the preferred track (cm, bin centre)
%   field_group "start" (< 40 cm), "middle" (40-100), "end" (>= 100)
%   w           log(max(r_T1, floor) / max(r_T2, floor)): linear weight of the
%               cell in the track log odds (floor = P.rate_floor)
%   session_rate mean rate over the whole RUN1 session (Hz); high_rate = > 10 Hz
nPos = numel(P.pos_edges) - 1;
[occ, cnt] = ar_bayes_lap_stats(Ds, T, P);
nC = numel(Ds.cell_id);
r = zeros(nC, 2, nPos); rm = zeros(nC, 2); rs = zeros(nC, 2);
for track = 1:2
    laps = T.track == track;
    o = sum(occ(laps, :), 1);
    c = sum(cnt(:, :, laps), 3);
    rr = c ./ o; rr(:, o == 0) = NaN;
    r(:, track, :) = reshape(rr, nC, 1, nPos);
    rm(:, track) = sum(c, 2) / sum(o);
    rs(:, track) = sum(c(:, 1:2), 2) / sum(o(1:2));
end
Cp = table();
Cp.cell_id = Ds.cell_id(:);
Cp.region = string(Ds.region(:));
Cp.r_T1 = rm(:, 1); Cp.r_T2 = rm(:, 2);
Cp.pref = (rm(:, 1) - rm(:, 2)) ./ (rm(:, 1) + rm(:, 2));
Cp.pref_start = (rs(:, 1) - rs(:, 2)) ./ (rs(:, 1) + rs(:, 2));
pt = 1 + (Cp.pref < 0);
Cp.peak_pos = NaN(nC, 1);
for i = 1:nC
    [~, k] = max(squeeze(r(i, pt(i), :)));
    Cp.peak_pos(i) = (k - 0.5) * (P.pos_edges(2) - P.pos_edges(1));
end
Cp.field_group = repmat("middle", nC, 1);
Cp.field_group(Cp.peak_pos < 40) = "start";
Cp.field_group(Cp.peak_pos >= 100) = "end";
Cp.w = log(max(rm(:, 1), P.rate_floor) ./ max(rm(:, 2), P.rate_floor));
[~, ci] = ismember(Ds.spike_id, Ds.cell_id);
Cp.session_rate = accumarray(ci(ci > 0), 1, [nC 1]) / (Ds.tvec(end) - Ds.tvec(1));
Cp.high_rate = Cp.session_rate > 10;
end
