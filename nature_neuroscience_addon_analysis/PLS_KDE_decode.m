function R = PLS_KDE_decode(dec, spikes, event_times)
% Decode track probabilities in sleep events with a decoder from PLS_KDE_train
% (minimal version of the decoding part of PLS_KDE_Reactivation_decoding.m).
% Each event window is split into 20 ms bins stepped by 10 ms; spike counts are
% z-scored per neuron across all bins of all events, the per-bin mean across
% neurons is subtracted, then projected on dec.W and evaluated under the Track 1
% and Track 2 KDEs (bandwidth dec.bandwidth).
%   spikes      : [time id], ids 1..dec.nNeurons (same neurons and numbering as
%                 used for training)
%   event_times : [start end] per event (e.g. [onset-1 onset+1])
% Output R: event_id (per bin), event_bins (per bin [start end]),
%   T1_probability, T2_probability (per bin), bias = T1 ./ (T1 + T2).
bin_size = 0.02;
step = 0.01;

edges = cell(size(event_times, 1), 1); ids = edges;
for i = 1:size(event_times, 1)
    num_bins = floor((event_times(i, 2) - event_times(i, 1)) / step);
    edges{i} = linspace(event_times(i, 1), event_times(i, 1) + num_bins * step, num_bins + 1)';
    ids{i} = i * ones(num_bins + 1, 1);
end
bins = vertcat(edges{:});
bins(:, 2) = bins + bin_size;
event_id = vertcat(ids{:});

% Count per event: windows of neighbouring events overlap, so the concatenated
% bin list is not sorted, which CountInIntervals requires
unit_spikes = arrayfun(@(u) spikes(spikes(:, 2) == u, 1), 1:dec.nNeurons, 'UniformOutput', false);
counts = cell(size(event_times, 1), 1);
parfor i = 1:size(event_times, 1)
    b = bins(event_id == i, :);
    c = zeros(size(b, 1), numel(unit_spikes));
    for unit = 1:numel(unit_spikes)
        c(:, unit) = CountInIntervals(unit_spikes{unit}, b);
    end
    counts{i} = c;
end
counts = vertcat(counts{:});
n = zscore(counts, 0, 1);
n(~isfinite(n)) = 0;          % neurons silent in all event bins contribute nothing
n_residual = n - mean(n, 2);

proj = n_residual * dec.W;
R.event_id = event_id';
R.event_bins = bins;
R.T1_probability = mvksdensity(dec.X_T1_proj, proj, 'Bandwidth', dec.bandwidth);
R.T2_probability = mvksdensity(dec.X_T2_proj, proj, 'Bandwidth', dec.bandwidth);
R.bias = R.T1_probability ./ (R.T1_probability + R.T2_probability);
end
