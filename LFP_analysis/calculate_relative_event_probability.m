function probabilities = calculate_relative_event_probability(event_A, event_B,num_bins)
% Calculate the probability distribution of observing event B relative to event A for a range of normalized time bins.
%
% Parameters:
% event_A (array): Timestamps of event A (nevent x 2, where nevent,1 is onset and nevent,2 is offset)
% event_B (array): Timestamps of event B (onset times)
% num_bins (int): Number of bins to divide the relative duration of event_A
%
% Returns:
% probabilities (array): Array of probabilities for each bin
% outside_count (int): Count of event_B occurrences outside of event_A

% Initialize the probabilities array
probabilities = zeros(1, num_bins);
outside_count = 0;

% Iterate through each event_A
for i = 1:size(event_A, 1)
    onset_A = event_A(i, 1);
    offset_A = event_A(i, 2);
    duration_A = offset_A - onset_A;

    % Normalize event_B times relative to the duration of event_A
    relative_times = (event_B - onset_A) / duration_A;

    % Count occurrences of event_B within each bin
    for j = 1:num_bins
        bin_start = (j - 1) / num_bins;
        bin_end = j / num_bins;
        probabilities(j) = probabilities(j) + sum(relative_times >= bin_start & relative_times < bin_end);
    end

    % Count occurrences of event_B outside of event_A
    % outside_count = outside_count + sum(relative_times < 0 | relative_times > 1);
end

% Normalize probabilities by the number of event_B
probabilities = probabilities / size(event_B, 1);
end