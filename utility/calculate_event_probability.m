function probabilities = calculate_event_probability(event_A, event_B, time_windows)
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

    for j = 1:length(time_windows) - 1
        count_A_within_window = 0;
        for i = 1:length(event_B)
            b = event_B(i);
            % Check if there are any event A timestamps within the current time window of event B
            count_A_within_window = count_A_within_window + sum((event_A >= b + time_windows(j)) & (event_A <= b + time_windows(j + 1)));
        end
        % Calculate probability for the current time window
        probabilities(j) = count_A_within_window / length(event_B);
    end
end