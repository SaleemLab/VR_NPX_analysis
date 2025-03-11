function [probabilities,binnedArray] = calculate_event_probability(event_A, event_B, time_windows,shuffle_options)
% Calculate the probability distribution of observing event A relative to event B for a range of time windows.
%
% Parameters:
% event_A (array): Timestamps of event A
% event_B (array): Timestamps of event B
% time_windows (array): Array of time windows to evaluate
%
% Returns:
% probabilities (array): Array of probabilities for each time window


probabilities = zeros(1, length(time_windows) - 1);
binnedArray = zeros(size(event_B,1), length(time_windows)-1);
for i = 1:length(event_B)
    b = event_B(i);

    % Count occurrences of event_B within each bin
    if shuffle_options ==1
        s = RandStream('mrg32k3a','Seed',i); % Set random seed for resampling
        bins_to_shift = datasample(s,1:length(time_windows),1);
        bins= circshift(time_windows,bins_to_shift);

    else
        bins = time_windows;
    end

    % Check if there are any event A timestamps within the current time window of event B
    for j = 1:length(bins)-1
        probabilities(j) = probabilities(j) + sum((event_A >= b + bins(j)) & (event_A <= b + bins(j + 1)));
        binnedArray(i,j) = sum((event_A >= b + bins(j)) & (event_A <= b + bins(j + 1)));
    end
end

probabilities = probabilities ./ length(event_B);
end