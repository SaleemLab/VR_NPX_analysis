function [probability, lags] = calculate_time_lagged_probability(event_times1, event_times2, bin_size, time_window)
    lags = time_window(1):bin_size:time_window(end);
    probability = zeros(size(lags));
    
    for i = 1:length(event_times2)
        diffs = event_times1 - event_times2(i); % Compute time differences
        hist_counts = histcounts(diffs, [lags - bin_size/2, lags(end) + bin_size/2]); % Bin time differences
        probability = probability + hist_counts;
    end
    
    % Normalize by total number of reference events
    if ~isempty(event_times2)
        probability = probability / length(event_times2);
    end
end