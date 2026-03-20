function spatial_rsp = cluster_spatial_responses(spike_times,tvec,position,speed,track_ID_all,start_time_all,end_time_all,method,delay,normalization)
if nargin < 10
normalization = 0;
end
if nargin < 9
    delay = 0;
end
spike_times = spike_times-delay;

w = gausswin(11);
w = w / sum(w);
speed(isnan(speed)) = 0;
speed= filtfilt(w,1,speed')';
spike_speed = interp1(tvec,speed,spike_times,'nearest');
distance = speed*(1/60);
distance(isnan(distance)) = 0;
%bin spike times into bins of tvec from beahviour at 60HZ
cluster_responses = zeros(1,length(tvec));
tvec_edges = [tvec(1)-1/(60*2) tvec+1/(60*2)];

cluster_responses(1,:) = histcounts(spike_times,tvec_edges)*60;
%smoothing filter ~200ms Gaussian window 11 timepoints at 60Hz

smoothed_cluster_responses = filtfilt(w,1,cluster_responses')';
if normalization
    smoothed_cluster_responses = normalize(smoothed_cluster_responses,2,'range');
end
%bin them based on laps
t1_index = find(track_ID_all == 1);
bin_no = 140; no_t1_laps = length(t1_index);
t2_index = find(track_ID_all == 2);
no_t2_laps = length(t2_index);

switch method
    case 'within'

        bin_edges = 0:140; bin_centres = 0.5:1:139.5;
        t1_spatial_responses = zeros(no_t1_laps,bin_no);

        for iLap = 1: length(t1_index)
            lap_start_time =  start_time_all(t1_index(iLap));
            lap_end_time =  end_time_all(t1_index(iLap));
            lap_tvec_index =  tvec >= lap_start_time &  tvec <= lap_end_time;
            lap_position =  position(lap_tvec_index);
            lap_response = smoothed_cluster_responses(:,lap_tvec_index);
            position_index = discretize(lap_position,bin_edges);
            temp_spatial_responses = zeros(1, length(bin_centres));
            for iBin = 1:length(bin_centres)
                bin_index = position_index == iBin & speed(lap_tvec_index) > 1;
                no_bin_index = sum(bin_index);
                temp_spatial_responses(iBin) = sum(lap_response(:,bin_index),2)/no_bin_index;
            end
            t1_spatial_responses(iLap, :) = temp_spatial_responses;
        end


        spatial_rsp{1,1} = t1_spatial_responses;
        if no_t2_laps > 0

            t2_spatial_responses = zeros([no_t2_laps,bin_no]);
            for iLap = 1: length(t2_index)
                lap_start_time =  start_time_all(t2_index(iLap));
                lap_end_time =  end_time_all(t2_index(iLap));
                lap_tvec_index =  tvec >= lap_start_time &  tvec <= lap_end_time;
                lap_position =  position(lap_tvec_index);
                lap_response = smoothed_cluster_responses(:,lap_tvec_index);
                position_index = discretize(lap_position,bin_edges);
                for iBin = 1:length(bin_centres)
                    bin_index = position_index == iBin & speed(lap_tvec_index) > 1;
                    no_bin_index = sum(bin_index);
                    t2_spatial_responses(iLap, iBin) = sum(lap_response(:,bin_index),2)/no_bin_index;
                end

            end

            spatial_rsp{1,2} = t2_spatial_responses;
        end

    case 'extension'

        for iLap = 1: length(t1_index)

            if iLap < length(t1_index)
                lap_start_time =  start_time_all(t1_index(iLap));
                lap_end_time =  start_time_all(t1_index(iLap+1));


            else
                lap_start_time =  start_time_all(t1_index(iLap));
                lap_end_time =  end_time_all(t1_index(iLap));


            end

            lap_tvec_index = find(tvec >=lap_start_time & tvec < lap_end_time);
            lap_position = cumsum(distance(lap_tvec_index));

            max_distance = max(lap_position);
            if max_distance > 0
                bin_edges = 0:ceil(max_distance);
                lap_response = smoothed_cluster_responses(:,lap_tvec_index);
                position_index = discretize(lap_position,bin_edges);
                t1_spatial_responses_tmp = zeros(1,length(bin_edges)-1);
                for iBin = 1:length(bin_edges)-1
                    bin_index = position_index == iBin & speed(lap_tvec_index) > 1;
                    no_bin_index = sum(bin_index);
                    t1_spatial_responses_tmp(:,iBin) = sum(lap_response(:,bin_index),2)/no_bin_index;
                end
                t1_spatial_responses_lap{iLap,1} = t1_spatial_responses_tmp;
            else
                t1_spatial_responses_lap{iLap,1} = zeros(1,1);
            end

        end

        spatial_rsp{1,1} = t1_spatial_responses_lap;
        if no_t2_laps > 0

            for iLap = 1: length(t2_index)

                if iLap < length(t2_index)
                    lap_start_time =  start_time_all(t2_index(iLap));
                    lap_end_time =  start_time_all(t2_index(iLap+1));
                else
                    lap_start_time =  start_time_all(t2_index(iLap));
                    lap_end_time =  end_time_all(t2_index(iLap));
                end

                lap_tvec_index = find(tvec >=lap_start_time & tvec < lap_end_time);
                lap_position = cumsum(distance(lap_tvec_index));

                max_distance = max(lap_position);
                if max_distance > 0
                    bin_edges = 0:ceil(max_distance);
                    lap_response = smoothed_cluster_responses(:,lap_tvec_index);
                    position_index = discretize(lap_position,bin_edges);
                    t2_spatial_responses_tmp = zeros([1,length(bin_edges)-1]);
                    for iBin = 1:length(bin_edges)-1
                        bin_index = position_index == iBin & speed(lap_tvec_index) > 1;
                        no_bin_index = sum(bin_index);
                        t2_spatial_responses_tmp(:,iBin) = sum(lap_response(:,bin_index),2)/no_bin_index;
                    end
                    t2_spatial_responses_lap{iLap,1} = t2_spatial_responses_tmp;
                else
                    t2_spatial_responses_lap{iLap,1} = zeros(1,1);
                end

            end

            spatial_rsp{1,2} = t2_spatial_responses_lap;
        end

end