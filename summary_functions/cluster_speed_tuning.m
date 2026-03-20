function speed_response = cluster_speed_tuning(spike_times,tvec,speed,track_ID_all,start_time_all,end_time_all)

    % doing block-based distance calculation 
    t_bin = 1/60;  
    
    no_lap = length(start_time_all);
    % no_rewarded_lap = length(rewarded_lap_id);
    w = gausswin(15);
    w = w / sum(w);
    speed(isnan(speed)) = 0;
    speed_smoothed = filtfilt(w,1,speed')';
    tvec_edges = [tvec(1)-1/120 (tvec(1:end-1) + tvec(2:end))/2 tvec(end)+1/120];

    
    track1_time = zeros(size(tvec));
    track2_time = zeros(size(tvec));
    
    %find track 1, track 2, non track times
    for iLap = 1:no_lap
        if iLap < no_lap
            if track_ID_all(iLap) == 1
            track1_time_tmp = tvec >= start_time_all(iLap) ... 
                & tvec <= end_time_all(iLap);
            track1_time = track1_time | track1_time_tmp;
            else
                track2_time_tmp = tvec >= start_time_all(iLap) ... 
                & tvec <= end_time_all(iLap);
                track2_time = track2_time | track2_time_tmp;
            end

        end
    end
    non_track_time = (~track1_time)& (~track2_time);



    firing_rate(1,:) = histcounts(spike_times,tvec_edges)*60;
    firing_rate = filtfilt(w,1,firing_rate')';
    speed_count_edges = 0:2:50;
    %find the time bin indices corresponding to different speed bins
    speed_count_indices = discretize(speed_smoothed,speed_count_edges); 
    track1_speed_response = nan([1 length(speed_count_edges)-1]);
    track2_speed_response = nan([1 length(speed_count_edges)-1]);
    non_track_speed_response = nan([1 length(speed_count_edges)-1]);
    for iBin = 1:length(speed_count_edges)-1
        speed_count_index = speed_count_indices ==iBin;
        % calculate average firing rate at different speeds
        track1_speed_response(:,iBin) = sum(firing_rate(:,track1_time & speed_count_index),2)/sum(track1_time & speed_count_index);
        track2_speed_response(:,iBin) = sum(firing_rate(:,track2_time & speed_count_index),2)/sum(track2_time & speed_count_index);
        non_track_speed_response(:,iBin) = sum(firing_rate(:,non_track_time & speed_count_index),2)/sum(non_track_time & speed_count_index);
    end
    
    
   speed_response{1} = track1_speed_response;
   speed_response{2} = track2_speed_response;
   speed_response{3} = non_track_speed_response;
    
    
    
    
  

end