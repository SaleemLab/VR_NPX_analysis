function maps = ar_rate_maps(occ, cnt, train, lap_track, nPos)
% Rate maps (cells x [T1 bins, T2 bins], Hz) from the training laps' running
% occupancy and spike counts (ar_bayes_lap_stats); unvisited bins are 0
maps = zeros(size(cnt, 1), 2 * nPos);
for track = 1:2
    laps = train & lap_track == track;
    o = sum(occ(laps, :), 1);
    c = sum(cnt(:, :, laps), 3);
    rate = c ./ o; rate(:, o == 0) = 0;
    maps(:, (track - 1) * nPos + (1:nPos)) = rate;
end
end
