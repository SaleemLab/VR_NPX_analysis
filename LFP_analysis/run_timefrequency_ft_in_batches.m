function [FFT,timebin,freqs] = run_timefrequency_ft_in_batches(cfg, data, batch_size)
% Batch process trials in FieldTrip frequency analysis
%
% cfg        - your FieldTrip config (including all mtmconvol options)
% data       - FieldTrip data structure (with .trial, .time, etc.)
% batch_size - number of trials per batch (e.g., 500)

nTrials = numel(data.trial);
batch_starts = 1:batch_size:nTrials;

TFR_batch = [];
FFT = [];  % initialize output
global ft_default
ft_default.warning = 'no';   % turns off ft_warning
ft_default.feedback = 'no';  % suppresses ft_progress and others
warning('off', 'all')

for i = 1:length(batch_starts)
    idx_start = batch_starts(i);
    idx_end   = min(idx_start + batch_size - 1, nTrials);

    fprintf('Processing trials %d to %d...\n', idx_start, idx_end);

    % Extract batch
    batch_data = data;
    batch_data.trial = data.trial(idx_start:idx_end);
    batch_data.time  = data.time(idx_start:idx_end);

    clear TFR_batch
    % Run frequency analysis
    % warning('off', 'all');  % turn off all warnings
    TFR_batch = ft_freqanalysis(cfg, batch_data);
    % warning('on', 'all');  % turn off all warnings

    FFT_batch = single(permute(TFR_batch.fourierspctrm(:, :,:,TFR_batch.time <= 2 &   TFR_batch.time >= -2), [2 3 4 1])); % grab and save only 2 second windows
    
    if isempty(FFT)
        FFT = FFT_batch;
    else
        FFT = cat(4, FFT, FFT_batch);
    end
end
% FFT = single(FFT);
timebin = TFR_batch.time;
freqs = TFR_batch.freq;
warning('on', 'all')
fprintf('Done. All batches processed and concatenated.\n');
end



