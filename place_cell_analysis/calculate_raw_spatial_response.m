function spatial_response = calculate_raw_spatial_response(spike_id,cluster_id,spike_times,tvec,position,speed,track_ID_all,start_time_all,end_time_all,bin_size)
no_cluster = length(cluster_id);
spatial_response = cell(no_cluster,2);
t_bin = mean(diff(tvec));
position_edges = 0:bin_size:140;
position_centres= bin_size/2:bin_size:140-bin_size/2;

%   convert spikes time to corresponding postions
spike_position = interp1(tvec,position,spike_times,'nearest');
spike_speed = interp1(tvec,speed,spike_times,'nearest');

no_lap = size(start_time_all,1);
event_position = zeros(size(start_time_all));

position_bin_time = zeros(no_lap,length(position_centres));
for iLap = 1:no_lap
    spike_times_lap_index = spike_times <= end_time_all(iLap)...
        & spike_times >= start_time_all(iLap) & spike_speed > 5;

    spike_position(spike_times_lap_index) = spike_position(spike_times_lap_index)+1000*(iLap);
    event_position(iLap,1) = (iLap)*1000;
    position_bin_time(iLap,:) = t_bin.*histcounts(position(tvec>=start_time_all(iLap) ...
        & tvec <=end_time_all(iLap) & speed > 5 ),position_edges);
end
for track_id = 1:max(track_ID_all)
    temp_event_position = event_position(track_ID_all==track_id);

    cluster_spike_id = cell(size(cluster_id));
    no_cluster = length(cluster_id);

    parfor iCluster = 1:no_cluster
        cluster_spike_id = spike_id == cluster_id(iCluster);

        [~,~,binnedArray] = spatial_psth(spike_position(cluster_spike_id),temp_event_position, [position_edges(1),position_edges(end)], bin_size,position_bin_time(track_ID_all==track_id,:));
        binnedArray(isnan(binnedArray)) = 0;
        binnedArray(isinf(binnedArray)) = 0;
        spatial_response{iCluster,track_id} = binnedArray;
    end
end
