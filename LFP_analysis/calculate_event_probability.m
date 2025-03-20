function [probabilities,binnedArray,event_index] = calculate_event_probability(event_A, event_B, time_windows,shuffle_options)
% Calculate the probability distribution of observing event A relative to event B for a range of time windows.
%
% Parameters:
% event_A (array): Timestamps of event A
% event_B (array): Timestamps of event B
% time_windows (array): Array of time windows to evaluate
%
% Returns:
% probabilities (array): Array of probabilities for each time window
% event_index (:,1) : Event A index (e.g. ripple)
% event_index (:,2) : Event B index (e.g.UP or DOWN)
% event_index (:,3) : Time of Event A retive to event B

bins_centre = linspace( time_windows(1)+mean(diff(time_windows))/2,time_windows(end)-mean(diff(time_windows))/2,length(time_windows)-1);
count = 1;
probabilities = zeros(1, length(time_windows) - 1);
binnedArray = zeros(size(event_B,1), length(time_windows)-1);
event_index = [];

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
    if size(event_A,2)==1
        for j = 1:length(bins)-1
            probabilities(j) = probabilities(j) + sum((event_A >= b + bins(j)) & (event_A <= b + bins(j + 1)));
            binnedArray(i,j) = sum((event_A >= b + bins(j)) & (event_A <= b + bins(j + 1)));

            if sum((event_A >= b + bins(j)) & (event_A <= b + bins(j + 1)))>0

                event_index(count,1) = find((event_A >= b + bins(j)) & (event_A <= b + bins(j + 1)));
                event_index(count,2) = i;
                event_index(count,3) = bins_centre(j);
                count = count +1;
            end
        end
    else
        for j = 1:length(bins)-1
            probabilities(j) = probabilities(j) + sum((event_A(:,end) >= b + bins(j)) & (event_A(:,1) <= b + bins(j + 1)));
            binnedArray(i,j) = sum((event_A(:,end) >= b + bins(j)) & (event_A(:,1) <= b + bins(j + 1)));

            if sum((event_A(:,end) >= b + bins(j)) & (event_A(:,1) <= b + bins(j + 1)))>0
                index = find((event_A(:,end) >= b + bins(j)) & (event_A(:,1) <= b + bins(j + 1)));

                for k = 1:length(index)
                    event_index(count,1) = index(k);
                    event_index(count,2) = i;
                    event_index(count,3) = bins_centre(j);
                    count = count +1;
                end
            end
        end
    end
end

probabilities = probabilities ./ length(event_B);
end