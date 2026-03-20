function [intervals] = OverlapIntervals(intervals1,intervals2)
% Find overlaps between two intervals 
% For example, if list 1 has [0 20; 25 45] and list two got [10 30; 31 34; 35 55]. 
% Then the output is [10 20; 31 34; 35 45];

intervals = [];

% Iterate through each interval in list1
for i = 1:size(intervals1, 1)
    list1 = intervals1(i, :);
    
    % Iterate through each interval in list2
    for j = 1:size(intervals2, 1)
        list2 = intervals2(j, :);
        
        % Check if there is an overlap
        start_overlap = max(list1(1), list2(1));
        end_overlap = min(list1(2), list2(2));
        
        if start_overlap < end_overlap
            % Store the overlapping interval
            intervals = [intervals; start_overlap, end_overlap];
        end
    end
end