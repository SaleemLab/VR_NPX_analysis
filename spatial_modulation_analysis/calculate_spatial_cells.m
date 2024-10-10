function place_fields = calculate_spatial_cells(clusters,Task_info,Behaviour,x_window,x_bin_width)

t_bin = mean(diff(Behaviour.tvec));
position_edges = x_window(1):x_bin_width:x_window(2);
%   convert spikes time to corresponding postions
spike_position = interp1(Behaviour.tvec,Behaviour.position,clusters.spike_times,'nearest');
spike_speed = interp1(Behaviour.tvec,Behaviour.speed,clusters.spike_times,'nearest');

% spike_position_original = interp1(Behaviour.tvec,Behaviour.position,clusters.spike_times,'nearest');
% spike_position_original(spike_speed<=5)=nan;
% spike_track_original = interp1(Behaviour.tvec,Behaviour.track_ID,clusters.spike_times,'nearest');

no_lap = size(Task_info.start_time_all,1);
event_position = zeros(size(Task_info.start_time_all));

position_bin_time = zeros(no_lap,(x_window(2)-x_window(1))/x_bin_width);
for iLap = 1:no_lap
    spike_times_lap_index = clusters.spike_times <= Task_info.end_time_all(iLap)...
        & clusters.spike_times >= Task_info.start_time_all(iLap) & spike_speed > 5;

    spike_position(spike_times_lap_index) = spike_position(spike_times_lap_index)+1000*(iLap);
    event_position(iLap,1) = (iLap)*1000;
    position_bin_time(iLap,:) = t_bin.*histcounts(Behaviour.position(Behaviour.tvec>=Task_info.start_time_all(iLap) ...
        & Behaviour.tvec <=Task_info.end_time_all(iLap) & Behaviour.speed > 5 ),position_edges);
end

track1_event_position = event_position(Task_info.track_ID_all==1);
track2_event_position = event_position(Task_info.track_ID_all==2);


cluster_spike_id = cell(size(clusters.cluster_id));
no_cluster = length(clusters.cluster_id);

track1_ID = find(Task_info.track_ID_all == 1);
track2_ID = find(Task_info.track_ID_all == 2);

place_fields = [];
for track_id = 1:max(Behaviour.track_ID)
    place_fields(track_id).x_bin_edges = position_edges;
    place_fields(track_id).x_bin_centres = [(position_edges(2)-x_bin_width/2):x_bin_width:(position_edges(end-1)+x_bin_width/2)];

    place_fields(track_id).raw = cell(size(clusters.cluster_id));
    place_fields(track_id).skaggs_info = cell(size(clusters.cluster_id));
    place_fields(track_id).spatial_xcorr = cell(size(clusters.cluster_id));
    place_fields(track_id).cluster_id = clusters.cluster_id;
    place_fields(track_id).peak_depth = clusters.peak_depth;
    place_fields(track_id).peak_channel = clusters.peak_channel;
    place_fields(track_id).cell_type = clusters.cell_type;
    place_fields(track_id).peak_channel_waveforms = clusters.peak_channel_waveforms;
    place_fields(track_id).probe_id = clusters.probe_id;

    if isfield(clusters,'probe_hemisphere')
        place_fields(track_id).probe_hemisphere = clusters.probe_hemisphere;
    end
    
    place_fields(track_id).dwell_map = position_bin_time(Task_info.track_ID_all==1,:);

end


raw1 = cell(size(clusters.cluster_id));
raw2 = cell(size(clusters.cluster_id));
spatial_xcorr1 = cell(size(clusters.cluster_id));
spatial_xcorr2 = cell(size(clusters.cluster_id));
skaggs_info1 = nan(size(clusters.cluster_id));
skaggs_info2 = nan(size(clusters.cluster_id));

% list_of_parameters
% Create smoothing filter (gamma)
if x_bin_width > 2
    w= [1 1];  %moving average filter of 2 sample, will be become a filter of [0.25 0.5 0.25] with filtfilt
else
    w= gausswin(10);
end

w = w./sum(w); %make sure smoothing filter sums to 1



parfor iCluster = 1:no_cluster
    cluster_spike_id = clusters.spike_id == clusters.cluster_id(iCluster);

    [psth_track1,bins,binnedArray] = spatial_psth(spike_position(cluster_spike_id),track1_event_position, x_window, x_bin_width,position_bin_time(Task_info.track_ID_all==1,:));
    binnedArray(isnan(binnedArray)) = 0;
    raw1{iCluster} = binnedArray;
    %     raw1{iCluster} = binnedArray;

     % Track 1 ratemap
    ratemap= mean(raw1{iCluster});
    %need dwell map
    dwellmap= sum(place_fields(1).dwell_map);
    %remove locations where the firing rate is 0
    dwellmap(ratemap==0)=[];
    ratemap(ratemap==0)=[];

    % probability of being at location x
    prob_x= dwellmap/nansum(dwellmap);
    % overall firing rate on track
    meanFR= nansum(ratemap.*prob_x);
    norm_rate= ratemap/meanFR;
    log_norm= log2(norm_rate);
    skaggs_info1(iCluster)= nansum(prob_x.*norm_rate.*log_norm);

    [psth_track2,bins,binnedArray] = spatial_psth(spike_position(cluster_spike_id),track2_event_position, x_window, x_bin_width,position_bin_time(Task_info.track_ID_all==2,:));
    binnedArray(isnan(binnedArray)) = 0;
    raw2{iCluster} = binnedArray;
    %     raw2{iCluster} = binnedArray;

    % Track 2 ratemap
    ratemap= mean(raw2{iCluster});
    %need dwell map
    dwellmap= sum(place_fields(2).dwell_map);
    %remove locations where the firing rate is 0
    dwellmap(ratemap==0)=[];
    ratemap(ratemap==0)=[];

    % probability of being at location x
    prob_x= dwellmap/nansum(dwellmap);
    % overall firing rate on track
    meanFR= nansum(ratemap.*prob_x);
    norm_rate= ratemap/meanFR;
    log_norm= log2(norm_rate);
    skaggs_info2(iCluster)= nansum(prob_x.*norm_rate.*log_norm);

    % Compute autocorrelation of the firing rate time series
    [track1_CCG, lags] = xcorr(psth_track1,psth_track1,'coeff');
    spatial_xcorr1{iCluster} = [track1_CCG; lags];
    [track2_CCG, lags] = xcorr(psth_track2,psth_track2,'coeff');
    spatial_xcorr2{iCluster} = [track2_CCG; lags];
end

place_fields(1).raw = raw1;
place_fields(2).raw = raw2;
place_fields(1).spatial_xcorr = spatial_xcorr1;
place_fields(2).spatial_xcorr = spatial_xcorr2;
place_fields(1).skaggs_info = skaggs_info1;
place_fields(2).skaggs_info = skaggs_info2;
