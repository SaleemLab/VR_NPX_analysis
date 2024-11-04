function [intervals] = MergeIntervals(intervals1,intervals2)
% Merge two intervals 
% . For example, if list 1 has [0 10; 45 55] and list two got [7 20; 32 36; 40 65].
% Then the output will be [0 20; 32 36; 40 65];


% Combine the two lists into one
combined_list = [intervals1; intervals2];

% Sort the combined list by the start of the intervals
combined_list = sortrows(combined_list);

% Initialize an empty array to store the merged intervals
merged_intervals = [];

% Initialize the first interval
current_interval = combined_list(1, :);

% Iterate through the combined list and merge intervals
for i = 2:size(combined_list, 1)
    if combined_list(i, 1) <= current_interval(2)
        % If the current interval overlaps or is adjacent to the next interval, merge them
        current_interval(2) = max(current_interval(2), combined_list(i, 2));
    else
        % If the current interval does not overlap, add it to the merged intervals
        merged_intervals = [merged_intervals; current_interval];
        current_interval = combined_list(i, :);
    end
end

% Add the last interval
intervals = [merged_intervals; current_interval];