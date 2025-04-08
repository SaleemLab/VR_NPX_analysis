function [probabilities,event_index,normalized_duration,binnedArray] = calculate_relative_event_probability(event_A, event_B,num_bins,shuffle_options)
% Calculate the probability distribution of observing event B relative to event A for a range of normalized time bins.
%
% Parameters:
% event_A (array): Timestamps of event A (nevent x 2, where nevent,1 is onset and nevent,2 is offset)
% event_B (array): Timestamps of event B (onset times)
% num_bins (int): Number of bins to divide the relative duration of event_A
% shuffle_options: 0 or 1 for circularly shifting event count within each event.

% Returns:
% probabilities (array): Array of probabilities for each bin
% outside_count (int): Count of event_B occurrences outside of event_A
% event_index (array): Array indicating which event_A each event_B falls into
% normalized_duration (array): Array indicating the normalized time of each event_B within event_A

% Initialize the probabilities array
probabilities = zeros(1, num_bins);
binnedArray = zeros(size(event_A,1), num_bins);
outside_count = 0;
event_index = [];
normalized_duration = nan(1,length(event_B));
count = 1;
bin_edges = linspace(0, 1, num_bins+1);
bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;

% Iterate through each event_A
for i = 1:size(event_A, 1)
    onset_A = event_A(i, 1);
    offset_A = event_A(i, 2);
    duration_A = offset_A - onset_A;

    % Normalize event_B times relative to the duration of event_A
    relative_times = (event_B - onset_A) / duration_A;
    relative_times = relative_times(relative_times(:,end)>=0 & relative_times(:,1)<=1,:);
    
    if isempty(relative_times)
        continue
    end
    % Find event_Bs within this event_A
    % in_window_idx = find(relative_times >= 0 & relative_times < 1);
    % relative_times = relative_times(relative_times >= 0 & relative_times < 1);

    % Count occurrences of event_B within each bin
    if shuffle_options ~=0
        s = RandStream('mrg32k3a','Seed',i+10000*shuffle_options); % Set random seed for resampling
        bins_to_shift = datasample(s,1:num_bins,1);
        bins= circshift(1:num_bins,bins_to_shift);
         
    else
        bins = 1:num_bins;
    end

    for j = 1:num_bins
        n = bins(j);
        bin_start = (n - 1) / num_bins;
        bin_end = n / num_bins;
        probabilities(j) = probabilities(j) + sum(relative_times(:,end) >= bin_start & relative_times(:,1) < bin_end);

        binnedArray(i,j) =sum(relative_times(:,end) >= bin_start & relative_times(:,1) < bin_end);

        if sum(relative_times(:,end) >= bin_start & relative_times(:,1) < bin_end)>0
            index = find(relative_times(:,end) >= bin_start & relative_times(:,1) < bin_end);
            for k = 1:length(index)
                event_index(count,1) = index(k);
                event_index(count,2) = i;
                event_index(count,3) = bin_start+(bin_end-bin_start)/2;
                count = count +1;
            end
        end
    end

    % Count occurrences of event_B outside of event_A
    % outside_count = outside_count + sum(relative_times < 0 | relative_times > 1);




    % % Record event_B info
    % for k = 1:length(event_B)
    %     if relative_times(k) >= 0 && relative_times(k) <= 1
    %         event_index(k) = i; % Record event_A index
    %         normalized_duration(k) = relative_times(k); % Record normalized time
    %     end
    % end
end

% temp = event_index;
% [~,index]=find(~isnan(temp));
% 
% event_index=[];
% event_index(:,1) = index; % Record event_A index
% event_index(:,2) = temp(index); % Record event_A index

% Normalize probabilities by the number of event_B
probabilities = probabilities / size(event_A, 1);
end

