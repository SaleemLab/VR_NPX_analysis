function [beh, spk, folds] = RUN1_session_behaviour(Ds, nsession, P)
% Behaviour, spikes and cross-validation folds of one RUN1 session (one element of
% RUN1_decoding_input.mat), shared by run_Bayesian_RUN_decoding_CV.m and
% run_Lseq_RUN_decoding_CV.m.
%   P.pos_edges       : position bin edges (cm), e.g. 0:10:140
%   P.speed_threshold : running threshold (cm/s)
%   P.nFolds          : number of CV folds (laps split per track)
% beh  : per behaviour sample (tvec): lap, track, pos_bin, running
% spk  : per spike: t, cell (index into Ds.cell_id), lap, track, pos_bin, running
% folds: fold id per lap; laps of each track are randomly split into nFolds groups of
%        (near) equal size (seed 1000*nsession + track)

tvec = Ds.tvec;
nLaps = size(Ds.lap_times, 1);
beh.tvec = tvec;
beh.lap = zeros(size(tvec));
for k = 1:nLaps
    beh.lap(tvec >= Ds.lap_times(k, 1) & tvec < Ds.lap_times(k, 2)) = k;
end
beh.track = zeros(size(tvec));
beh.track(beh.lap > 0) = Ds.lap_track(beh.lap(beh.lap > 0));
beh.pos_bin = discretize(Ds.position, P.pos_edges);
beh.speed = Ds.speed;
beh.running = beh.lap > 0 & ~isnan(beh.pos_bin) & Ds.speed > P.speed_threshold;
beh.dt = [diff(tvec); median(diff(tvec))];

[~, cell_idx] = ismember(Ds.spike_id, Ds.cell_id);
keep = cell_idx > 0 & Ds.spike_times >= tvec(1) & Ds.spike_times <= tvec(end);
spk.t = Ds.spike_times(keep);
spk.cell = cell_idx(keep);
sample = interp1(tvec, (1:numel(tvec))', spk.t, 'nearest');
spk.lap = beh.lap(sample);
spk.track = beh.track(sample);
spk.pos_bin = beh.pos_bin(sample);
spk.running = beh.running(sample);
% spikes in a lap but on a sample just outside its window keep lap 0
in_lap = spk.lap > 0;
spk.running(in_lap) = spk.running(in_lap) & spk.t(in_lap) >= Ds.lap_times(spk.lap(in_lap), 1) & spk.t(in_lap) < Ds.lap_times(spk.lap(in_lap), 2);

folds = zeros(nLaps, 1);
for track = 1:2
    laps = find(Ds.lap_track == track);
    s = RandStream('mrg32k3a', 'Seed', 1000 * nsession + track);
    order = laps(randperm(s, numel(laps)));
    folds(order) = mod(0:numel(order) - 1, P.nFolds) + 1;
end
end
