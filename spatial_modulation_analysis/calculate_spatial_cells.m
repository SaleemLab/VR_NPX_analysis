function place_fields = calculate_spatial_cells(clusters,tvec,position,speed,track_ID_all,start_time_all,end_time_all,x_window,x_bin_width)
% Mainly for Masa 2 track place fields structures for Bayesian deocidng

t_bin = mean(diff(tvec));
position_edges = x_window(1):x_bin_width:x_window(2);
%   convert spikes time to corresponding postions
spike_position = interp1(tvec,position,clusters.spike_times,'nearest');
spike_speed = interp1(tvec,speed,clusters.spike_times,'nearest');

% spike_position_original = interp1(Behaviour.tvec,Behaviour.position,clusters.spike_times,'nearest');
% spike_position_original(spike_speed<=5)=nan;
% spike_track_original = interp1(Behaviour.tvec,Behaviour.track_ID,clusters.spike_times,'nearest');

no_lap = size(start_time_all,1);
event_position = zeros(size(start_time_all));

position_bin_time = zeros(no_lap,(x_window(2)-x_window(1))/x_bin_width);
for iLap = 1:no_lap
    spike_times_lap_index = clusters.spike_times <= end_time_all(iLap)...
        & clusters.spike_times >= start_time_all(iLap) & spike_speed > 1;

    spike_position(spike_times_lap_index) = spike_position(spike_times_lap_index)+1000*(iLap);
    event_position(iLap,1) = (iLap)*1000;
    position_bin_time(iLap,:) = t_bin.*histcounts(position(tvec>=start_time_all(iLap) ...
        & tvec <=end_time_all(iLap) & speed > 1 ),position_edges);
end

track1_event_position = event_position(track_ID_all==1);
track2_event_position = event_position(track_ID_all==2);


cluster_spike_id = cell(size(clusters.cluster_id));
no_cluster = length(clusters.cluster_id);

track1_ID = find(track_ID_all == 1);
track2_ID = find(track_ID_all == 2);

place_fields = [];
for track_id = 1:max(track_ID_all)
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
    place_fields(track_id).odd_even_stability = clusters.odd_even_stability(:,track_id);
    place_fields(track_id).first_second_stability = clusters.odd_even_stability(:,track_id);
    place_fields(track_id).peak_percentile = clusters.peak_percentile(:,track_id);
    place_fields(track_id).within_track_corr = clusters.within_track_corr(:,track_id);
    place_fields(track_id).across_track_corr = clusters.across_track_corr;
    
%     place_fields(track_id).probe_id = clusters.probe_id;

    if isfield(clusters,'probe_hemisphere')
        place_fields(track_id).probe_hemisphere = clusters.probe_hemisphere;
    end
    
    place_fields(track_id).dwell_map = position_bin_time(track_ID_all==1,:);

end

place_fields(track_id).across_track_corr = clusters.across_track_corr;

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

    [psth_track1,bins,binnedArray] = spatial_psth(spike_position(cluster_spike_id),track1_event_position, x_window, x_bin_width,position_bin_time(track_ID_all==1,:));
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

    [psth_track2,bins,binnedArray] = spatial_psth(spike_position(cluster_spike_id),track2_event_position, x_window, x_bin_width,position_bin_time(track_ID_all==2,:));
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



for track_id = 1:max(track_ID_all)
    %             place_fields_BAYESIAN(track_id).good_place_cells_LIBERAL= ...
    %                 find(HPC_clusters_RUN.peak_percentile(:,track_id)>0.95&HPC_clusters_RUN.odd_even_stability(:,track_id)>0.95);

    ratemap_matrix = [place_fields(track_id).raw{:}];
    ratemap_matrix = reshape(ratemap_matrix,size(place_fields(track_id).raw{1},1),[],length(place_fields(track_id).raw));%laps X position bins X cells
    average_maps= squeeze(mean(ratemap_matrix,1));
    place_fields(track_id).template=average_maps';% need ncell X nposition

    if isfield(clusters,'odd_even_stability')
        place_fields(track_id).good_place_cells_LIBERAL= ...
            find(clusters.odd_even_stability(:,track_id)>0.95);

        [~,peak_locations] = max(average_maps(:,place_fields(track_id).good_place_cells_LIBERAL));
        %     unsorted_cells(track_id,:) =
        [~,sort_id] = sort(peak_locations);
        place_fields(track_id).sorted_good_place_cells_LIBERAL=place_fields(track_id).good_place_cells_LIBERAL(sort_id);
    end

    %             for iCell = 1:length(place_fields_BAYESIAN(track_id).raw)
    %                 place_fields_BAYESIAN(track_id).template{iCell}=average_maps(:,iCell);% average map for decoding
    %             end


end

if isfield(clusters,'odd_even_stability')
    place_fields(1).all_good_place_cells_LIBERAL = unique([place_fields(1).good_place_cells_LIBERAL; place_fields(2).good_place_cells_LIBERAL]);
    place_fields(2).all_good_place_cells_LIBERAL = unique([place_fields(1).good_place_cells_LIBERAL; place_fields(2).good_place_cells_LIBERAL]);
end
