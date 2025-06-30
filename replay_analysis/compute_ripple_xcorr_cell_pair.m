function [pair_corr, pair_pval,lags] = compute_ripple_xcorr_cell_pair(spike_times, spike_ids, V1_id, HPC_id, event_times, varargin)
% Faster version of ripple cross-correlation with parfor support

% --- Parse inputs ---
p = inputParser;
addParameter(p, 'window', 0.1);
addParameter(p, 'step', 0.02);
addParameter(p, 'lag_range', [-1 1]);
addParameter(p, 'shuffle_option', 1);

parse(p, varargin{:});
window = p.Results.window;
step = p.Results.step;
lag_range = p.Results.lag_range;
shuffle_option = p.Results.shuffle_option;
% --- Define lags ---


if length(lag_range) ==1
    lags = lag_range;
else
    lags = lag_range(1):step:lag_range(2);
    
end
numLags = length(lags);


nHPC = length(HPC_id);
nV1  = length(V1_id);
nRipples = length(event_times);

% --- Preallocate output ---
pair_corr = nan(nHPC, nV1, numLags);

% --- Pre-extract spike times per unit ---
HPC_spikes = cell(nHPC, 1);
V1_spikes  = cell(nV1, 1);
for i = 1:nHPC
    HPC_spikes{i} = spike_times(spike_ids == HPC_id(i));
end
for j = 1:nV1
    V1_spikes{j} = spike_times(spike_ids == V1_id(j));
end

% --- Parallel loop across pairs ---
% Preallocate output as cell for parfor safety
pair_corr_cell = cell(nHPC * nV1, 1);
pair_pval = nan(nHPC * nV1, 1);

parfor pair_idx = 1:(nHPC * nV1)
    % tic
    [iHPC, iV1] = ind2sub([nHPC, nV1], pair_idx);
    spikes_H = HPC_spikes{iHPC};
    spikes_V = V1_spikes{iV1};

    corr_vec = nan(1, numLags);
    corr_vec_shuffled = [];
    tenp=[];
    for iLag = 1:numLags
        lag = lags(iLag);

        hpc_counts = zeros(nRipples, 1);
        v1_counts  = zeros(nRipples, 1);

        for iEvt = 1:nRipples
            t_evt = event_times(iEvt);

            hpc_win = [t_evt, t_evt + window];
            v1_win = [t_evt + lag, t_evt + lag + window];

            hpc_counts(iEvt) = sum(spikes_H >= hpc_win(1) & spikes_H < hpc_win(2));
            v1_counts(iEvt)  = sum(spikes_V >= v1_win(1) & spikes_V < v1_win(2));
        end

        if any(hpc_counts) && any(v1_counts)
            [corr_vec(iLag) tenp]= corr(hpc_counts, v1_counts, 'type', 'Spearman');
            if lag == 0
                pair_pval(pair_idx) = tenp;
            end
        end

        if shuffle_option == 1
            if lag == 0
                for nshuffle = 1:1000
                    s = RandStream('philox4x32_10', 'Seed', nshuffle);
                    bins_to_shift = randsample(s,1:nRipples,1);
                    hpc_counts_shuffled = circshift(hpc_counts,bins_to_shift);
                    [corr_vec_shuffled(nshuffle)]= corr(hpc_counts_shuffled, v1_counts, 'type', 'Spearman','Rows', 'complete');
                end
                pair_pval(pair_idx) = mean(corr_vec(iLag) > corr_vec_shuffled)

            end
        end
    end

    pair_corr_cell{pair_idx} = corr_vec;
    % toc
end

% Convert back to 3D array
pair_corr = nan(nHPC, nV1, numLags);
for pair_idx = 1:(nHPC * nV1)
    [iHPC, iV1] = ind2sub([nHPC, nV1], pair_idx);
    pair_corr(iHPC, iV1, :) = pair_corr_cell{pair_idx};

end


pair_pval = reshape(pair_pval,nHPC,nV1);
end
