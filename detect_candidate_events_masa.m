function [replay,reactivations] = detect_candidate_events_masa(tvec,raw_LFP,mua_zscore,spike_times,peripherals,zscore_min,zscore_max,options)
%   based on Davidson, Kloosterman, and Wilson(2009,Neuron) criteria.  
%   - events detected by MUA only
%   - zscore of 3 minimum, and edges detected by zscore of 0.  
%     modificiation- If event is too long, then edge criteria shifted to
%     zscore of 0.25, then 0.5
%   - events within 50 ms are combined together
%   - replay events less than and equal 100ms are removed. 
% thus bayesian decoded candidate replay events have at least 5
% consecutive 20 ms bins, and 5 different neurons.

% things to decide:
%   - durations: shortest 40/50ms, above 100ms use for decoding. stored as
%   'reactivations' right now these are not merged or otherwise processed
%   - speed filter: 5cm/s
%   - theta filter (may have some bleeding from sharp wave band filtering) added 
%     for ranking: with csc high theta low ripple power
%   - adaptive thresholding: tag those that have been truncated and have a
%   look: start and end may have different min zscore as a result

parameters= list_of_parameters;
SR = 1/mean(diff(tvec));

% Ripple band
parameters = list_of_parameters;
filter_type  = 'bandpass';
filter_width = [125 300];                 % range of frequencies in Hz you want to filter between
filter_order = round(6*SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
norm_freq_range = filter_width/(SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_ripple = fir1(filter_order, norm_freq_range,filter_type);

ripple_LFP = filtfilt(b_ripple,1,raw_LFP);
ripple_zscore = zscore(abs(hilbert(ripple_LFP)));

% theta band
parameters = list_of_parameters;
filter_type  = 'bandpass';
filter_width = [4 12];                 % range of frequencies in Hz you want to filter between
filter_order = round(6*SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
norm_freq_range = filter_width/(SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_theta = fir1(filter_order, norm_freq_range,filter_type);

theta_LFP = filtfilt(b_theta,1,raw_LFP);
theta_zscore = zscore(abs(hilbert(theta_LFP)));
% ripple_zscore=interp1(LFP_time,ripple_zscore,MUA_time,'linear'); %adapt to MUA_time bins time steps
% theta_zscore= interp1(LFP_time,theta_zscore,MUA_time,'linear'); %adapt to MUA_time bins time steps

[replay,reactivations]=replay_search(tvec,spike_times,mua_zscore,ripple_zscore,theta_zscore,peripherals,parameters.min_zscore,parameters.max_zscore,options); % 0,3
% save extracted_candidate_events replay reactivations
end

function [replay, reactivations]=replay_search(time,spike_times,mua_zscore,ripple_zscore,theta_zscore,peripherals,zscore_min,zscore_max,options)
% spike_times = [];
% load extracted_clusters;
% load extracted_place_fields;
% load extracted_position;   

parameters= list_of_parameters;
mua_idx = find(mua_zscore>=zscore_max); %finds indices where MUA is above set threshold
max_search_length = parameters.max_search_length;  %300 ms max search
mua_bursts = find(diff([0 time(mua_idx)])>=0.01);  %find at least a 10 ms gap between bursts to define a burst

for k=1:length(mua_bursts)
    replay.error(k) = 0;
    replay.thresh_adapted_start(k) = 0;
    replay.thresh_adapted_end(k) = 0;
    
    % Finds start and end indices for the current MUA burst
    burst_index = mua_idx(mua_bursts(k)); %finds index of the current MUA peak
    burst_start = mua_idx(mua_bursts(k))-max_search_length;  %sets start of MUA burst
    if burst_start<1
        burst_start = 1;
    end
    burst_end = mua_idx(mua_bursts(k))+max_search_length;  %sets end of MUA burst
    if burst_end>length(time)
        burst_end = length(time);
    end
    
    % Set replay event START based on MUA zscore criteria
    MUAburst_1sthalf_idxs = burst_start:burst_index;    
    start_idx = max(find(mua_zscore(MUAburst_1sthalf_idxs)<zscore_min)); %find start of replay event (when zscore first crosses minimum point)    
    if isempty(start_idx)
        start_idx = max(find(mua_zscore(MUAburst_1sthalf_idxs)<=(zscore_min+0.25))); %if event is too wide, use a stricter criteria
        replay.thresh_adapted_start(k) = 0.25; % adapted to 0.25
    end
    if isempty(start_idx)
        start_idx = max(find(mua_zscore(MUAburst_1sthalf_idxs)<=(zscore_min+0.5))); %if event is too wide, use an even stricter criteria
        replay.thresh_adapted_start(k) = 0.5; % adapted to 0.5
    end
    if isempty(start_idx)
        replay.onset(k) = time(MUAburst_1sthalf_idxs(1));
        replay.error(k) = 1;
        start_idx = 1;
        %burst_start = burst_start-max_search_length; %negative?
    else
        replay.onset(k) = time(MUAburst_1sthalf_idxs(start_idx));
    end
    
    % Set replay event END based on MUA zscore criteria
    MUAburst_2ndhalf_idxs = burst_index:burst_end;    
    end_idx = min(find(mua_zscore(MUAburst_2ndhalf_idxs)<zscore_min));  %find end of replay event (when zscore crosses minimum point)
    if isempty(end_idx)
         end_idx = min(find(mua_zscore(MUAburst_2ndhalf_idxs)<=(zscore_min+0.25))); %if event is too wide, use a stricter criteria
         replay.thresh_adapted_end(k) = 0.25;
    end
      if isempty(end_idx)
         end_idx = min(find(mua_zscore(MUAburst_2ndhalf_idxs)<=(zscore_min+0.5))); %if event is too wide, use a stricter criteria
         replay.thresh_adapted_end(k) = 0.5;
     end
    if isempty(end_idx)
        replay.offset(k) = time(MUAburst_2ndhalf_idxs(end));
        replay.error(k) = 1;
        end_idx = length(MUAburst_2ndhalf_idxs);
        %burst_end = burst_end+max_search_length;
    else
        replay.offset(k) = time(MUAburst_2ndhalf_idxs(end_idx));
    end
end

[~,idx,~] = unique(replay.onset);
replay.error = replay.error(idx);  
replay.thresh_adapted_start = replay.thresh_adapted_start(idx);
replay.thresh_adapted_end = replay.thresh_adapted_end(idx);
replay.onset = replay.onset(idx);
replay.offset = replay.offset(idx);

[~,idx,~] = unique(replay.offset);
replay.error = replay.error(idx);  
replay.thresh_adapted_start = replay.thresh_adapted_start(idx);
replay.thresh_adapted_end = replay.thresh_adapted_end(idx);
replay.onset = replay.onset(idx);
replay.offset = replay.offset(idx);

%%% Eliminate short singleton events or very long events
replay.duration = replay.offset-replay.onset;
%     % have softer criteria for reactivations 50ms
% reactivations = replay;
% short_event_idx = find(reactivations.duration< parameters.min_react_duration);
% reactivations.onset(short_event_idx)   = [];
% reactivations.offset(short_event_idx)  = [];
% reactivations.duration(short_event_idx)= [];
% reactivations.error(short_event_idx)   = [];
% reactivations.thresh_adapted_start(short_event_idx) = [];
% reactivations.thresh_adapted_end(short_event_idx)   = [];

    % and stronger ones for candidate replay events 100ms
short_replay_idx = find(replay.duration< parameters.min_event_duration);  %eliminate <=100 ms bursts
replay.onset(short_replay_idx)    = [];
replay.offset(short_replay_idx)   = [];
replay.duration(short_replay_idx) = [];
replay.error(short_replay_idx)    = [];
replay.thresh_adapted_start(short_replay_idx) = [];
replay.thresh_adapted_end(short_replay_idx)   = [];

% merge next replay event if less than 50 ms away 
replay_raw = replay;
replay_raw.merged = 1:length(replay.onset);
merge_idx = find(replay.onset(2:end)-replay.offset(1:end-1)<0.05);
while ~isempty(merge_idx) %while next replay event is less than 50 ms away 
    replay.offset(merge_idx(1)) = replay.offset(merge_idx(1)+1); %change offset to offset of next event
    replay.duration(merge_idx(1)) = replay.offset(merge_idx(1))-replay.onset(merge_idx(1));  %update duration value
    replay.error(merge_idx(1)) = replay.error(merge_idx(1))+replay.error(merge_idx(1)+1); %if there is an error in either event, make sure this is copied/accumulated
    %remove indice of second replay event in merger, across all fields
    replay.duration(merge_idx(1)+1) =[];    
    replay.onset(merge_idx(1)+1)  = [];
    replay.offset(merge_idx(1)+1) = [];
    replay.error(merge_idx(1)+1)  = [];
    %these variables are for tracking where the merged events occured
    replay.thresh_adapted_start(merge_idx(1)+1) = [];
    replay.thresh_adapted_end(merge_idx(1)) = [];
    replay_raw.merged(merge_idx(1)+1) = merge_idx(1);
    merge_idx = find(replay.onset(2:end)-replay.offset(1:end-1)<0.05);  %see if any other replay events are separated by less than 50 ms
end

reactivations = replay;
% check for long events (more than 1 second) but leave them
long_events = find(replay.duration>parameters.max_event_duration);
number_of_long_events=length(long_events)

% Sorts spike times for all cells
spike_id_sorted = spike_times(:,1);
spike_times_sorted = spike_times(:,2);

%  for j = 1:length(place_fields.good_place_cells)
%       index = find(clusters.spike_id==place_fields.good_place_cells(j));
%       spike_id_sorted = [spike_id_sorted j*ones(1,length(index))];
%       spike_times_sorted = [spike_times_sorted clusters.spike_times(index)'];   
%  end

% Find replay event features
for i = 1:length(replay.onset)
    replay.zscore(i) = max(mua_zscore(find(time>replay.onset(i) & time<replay.offset(i))));
    replay.ripple_peak(i) = max(ripple_zscore(find(time>replay.onset(i) & time<replay.offset(i))));
    replay.mean_theta(i) = mean(theta_zscore(find(time>replay.onset(i) & time<replay.offset(i))));
    speed = interp1(peripherals.sglxTime,peripherals.speed,[replay.onset(i):.01:replay.offset(i)],'nearest');
    replay.speed(i) = median(speed);
    spikes = spike_id_sorted(find(spike_times_sorted>=replay.onset(i) & spike_times_sorted<=replay.offset(i)));
    replay.spike_count(i)  = length(spikes);
    replay.neuron_count(i) = length(unique(spikes));  %all good cells- good place field on at least one track
     %find midpoint (minimum z-score in middle third of event)
    event_time_indices = find(time>replay.onset(i) & time<replay.offset(i));
    middle_third_idx = ceil(event_time_indices+length(event_time_indices)/3):floor(event_time_indices+2*length(event_time_indices)/3);
    [~,minMUA_idx] = min(mua_zscore(middle_third_idx)); %find min MUA for the middle third
    replay.midpoint(i) = time(middle_third_idx(minMUA_idx)); %timestamp with min MUA in the middle third of the event
end
reactivations = replay;
% 
% % Find reactivation event features
% for i = 1:length(reactivations.onset)
%     reactivations.zscore(i) = max(mua_zscore(find(time>reactivations.onset(i) & time<reactivations.offset(i))));
%     reactivations.ripple_peak(i) = max(ripple_zscore(find(time>reactivations.onset(i) & time<reactivations.offset(i))));
%     reactivations.mean_theta(i) = mean(theta_zscore(find(time>reactivations.onset(i) & time<reactivations.offset(i))));
%     speed = interp1(peripherals.sglxTime,peripherals.speed,[reactivations.onset(i):.01:reactivations.offset(i)],'nearest');
%     reactivations.speed(i) = median(speed);
%     spikes = spike_id_sorted(find(spike_times_sorted>=reactivations.onset(i) & spike_times_sorted<=reactivations.offset(i)));
%     reactivations.spike_count(i)  = length(spikes);
%     reactivations.neuron_count(i) = length(unique(spikes));  %all good cells- good place field on at least one track
%      %find midpoint (minimum z-score in middle third of event)
%     event_time_indices = find(time>reactivations.onset(i) & time<reactivations.offset(i));
%     middle_third_idx = ceil(event_time_indices+length(event_time_indices)/3):floor(event_time_indices+2*length(event_time_indices)/3);
%     [~,minMUA_idx] = min(mua_zscore(middle_third_idx)); %find min MUA for the middle third
%     reactivations.midpoint(i) = time(middle_third_idx(minMUA_idx)); %timestamp with min MUA in the middle third of the event
% end

% index to ignore
high_v_idx = find(reactivations.speed>5);
% [~,idx,~] = unique(reactivations.onset);
reactivations.onset(high_v_idx)    = [];
reactivations.offset(high_v_idx)   = [];
reactivations.duration(high_v_idx) = [];
reactivations.error(high_v_idx)    = [];
reactivations.zscore(high_v_idx)   = [];
reactivations.ripple_peak(high_v_idx) = [];
reactivations.speed(high_v_idx) = [];
reactivations.spike_count(high_v_idx)  = [];
reactivations.neuron_count(high_v_idx) = [];
reactivations.midpoint(high_v_idx)=[];
reactivations.mean_theta(high_v_idx)   = [];
reactivations.thresh_adapted_start(high_v_idx) = [];
reactivations.thresh_adapted_end(high_v_idx)   = [];

reactivations.zscore_min = zscore_min;
reactivations.zscore_max = zscore_max;

% speed, neuron count and mua zscore threshold
high_v_idx = find(replay.zscore<zscore_max | replay.speed>5 | replay.neuron_count<5);
replay.onset(high_v_idx)    = [];
replay.offset(high_v_idx)   = [];
replay.duration(high_v_idx) = [];
replay.error(high_v_idx)    = [];
replay.zscore(high_v_idx)   = [];
replay.ripple_peak(high_v_idx) = [];
replay.speed(high_v_idx) = [];
replay.spike_count(high_v_idx)  = [];
replay.neuron_count(high_v_idx) = [];
replay.midpoint(high_v_idx)=[];
replay.mean_theta(high_v_idx)   = [];
replay.thresh_adapted_start(high_v_idx) = [];
replay.thresh_adapted_end(high_v_idx)   = [];

replay.zscore_min = zscore_min;
replay.zscore_max = zscore_max;
end

