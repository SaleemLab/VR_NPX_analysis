function [theta_modulation] = phase_precession_absolute_location(LFP_tvec,CA1_LFP,place_fields,clusters,Task_info,Behaviour,options)

if isfield(clusters,'merged_spike_id')
    clusters.cluster_id = unique(clusters.merged_cluster_id);
    for ncell = 1:length(unique(clusters.merged_cluster_id))
        tempt_peak_channel = clusters.peak_channel(clusters.merged_cluster_id == clusters.cluster_id(ncell));
        tempt_peak_depth = clusters.peak_depth(clusters.merged_cluster_id == clusters.cluster_id(ncell));
        tempt_peak_waveform = clusters.peak_channel_waveforms(clusters.merged_cluster_id == clusters.cluster_id(ncell),:);
        tempt_cell_type = clusters.cell_type(clusters.merged_cluster_id == clusters.cluster_id(ncell));

        if length(tempt_peak_depth)== 2
            tempt_peak_channel = tempt_peak_channel(1);
            tempt_peak_depth = tempt_peak_depth(1);
            tempt_peak_waveform = tempt_peak_waveform(1,:);
            tempt_cell_type = tempt_cell_type(1);
        else % find median peak depth assign that value to the unit
            [~,index]= min(tempt_peak_depth - median(tempt_peak_depth));
            tempt_peak_channel = tempt_peak_channel(index);
            tempt_peak_depth = tempt_peak_depth(index);
            tempt_peak_waveform = tempt_peak_waveform(index,:);
            tempt_cell_type = tempt_cell_type(index);
        end

        merged_peak_channel(ncell) = tempt_peak_channel;
        merged_peak_depth(ncell) = tempt_peak_depth;
        merged_peak_waveform(ncell,:) = tempt_peak_waveform;
        merged_cell_type(ncell,:) = tempt_cell_type;
    end

    clusters.peak_channel = merged_peak_channel;
    clusters.peak_depth = merged_peak_depth;
    clusters.peak_channel_waveforms = merged_peak_waveform;
    clusters.cell_type = merged_cell_type;

    clusters.spike_id = clusters.merged_spike_id;
    clusters.cluster_id = unique(clusters.merged_cluster_id);
end

% Find hilbert transform to extract phase info

[b a] = butter(3,[6/625 9/625],'bandpass');
theta_phase = (angle(hilbert(filtfilt(b,a,CA1_LFP))));
theta_phase_unwrap = unwrap(theta_phase); % unwrap for interpolation
LFP_speed = interp1(Behaviour.tvec,Behaviour.speed,LFP_tvec,'linear');
theta_phase_unwrap(LFP_speed<5) = nan; % speed filtered Theta LFP


spatial_cell_index = unique([find(place_fields(1).peak_percentile>0.95 & place_fields(1).odd_even_stability>0.95)...
    find(place_fields(2).peak_percentile>0.95 & place_fields(2).odd_even_stability>0.95)]);
good_cell_index = intersect(spatial_cell_index,find(ismember(place_fields(1).cluster_id,clusters.cluster_id)));
good_cell_index = good_cell_index(1:5);
t_bin = mean(diff(Behaviour.tvec));
sprintf('Theta modulation analysis for %i spatial cells',length(good_cell_index))

count =1;
for track = 1 : length(place_fields)

    track_times = Behaviour.tvec;
    track_linear = Behaviour.position;
    track_linear(Behaviour.speed<5)=nan;
    track_linear(Behaviour.track_ID~=track) = nan;

    track_times_temp = track_times;
    track_times_temp(isnan(track_linear))=nan;

    theta_phase_temp = wrapTo2Pi(interp1(LFP_tvec,theta_phase_unwrap,track_times_temp,'linear'));
    %     theta_phase_temp = wrapTo2Pi(theta_phase_unwrap);

    position_edges = 0:2:140;
    phase_edges = 0:20:360;
    phase_edges = pi/180.*phase_edges;

    %2D position X phase occupancy map
    [N,Xedges,Yedges,binX,binY] = histcounts2(track_linear,theta_phase_temp,position_edges,phase_edges);
    occupancy_map = t_bin.*N;

    start_time_this_track = Task_info.start_time_all(Task_info.track_ID_all==track);
    end_time_this_track = Task_info.end_time_all(Task_info.track_ID_all==track);

    no_lap = size(start_time_this_track,1);
    %     lap_phase_occupancy = zeros(no_lap,length(phase_edges)-1);
    %     theta_phase1= []; %odd lap theta phase
    %     theta_phase2 = [];% even lap theta phase
    %      theta_phase3= [];
    %     for iLap = 1:no_lap
    %         if  mod(iLap, 2) == 0
    %             theta_phase2 = [theta_phase2 theta_phase_unwrap(LFP_tvec>=start_time_this_track(iLap) ...
    %                 & LFP_tvec <end_time_this_track(iLap))];
    %         else
    %              theta_phase1 = [theta_phase1 theta_phase_unwrap(LFP_tvec>=start_time_this_track(iLap) ...
    %                 & LFP_tvec<end_time_this_track(iLap))];
    %         end
    %
    %         theta_phase3 = [theta_phase3 theta_phase_unwrap(LFP_tvec>=start_time_this_track(iLap) ...
    %             & LFP_tvec<end_time_this_track(iLap))];
    %
    %         lap_phase_occupancy(iLap,:) = t_bin.*histcounts(theta_phase_temp(Behaviour.tvec>=start_time_this_track(iLap) ...
    %             & Behaviour.tvec <=end_time_this_track(iLap) & Behaviour.speed > 5 ),phase_edges);
    %     end

    phase_occupancy_map = t_bin.*histcounts(theta_phase_temp,phase_edges);
    % phase_occupancy_map = mean(diff(LFP_tvec)).*histcounts(theta_phase_temp,phase_edges);
    %     phase_occupancy_map1 = mean(diff(LFP_tvec)).*histcounts(wrapTo2Pi(theta_phase1),phase_edges);
    %     phase_occupancy_map2 = mean(diff(LFP_tvec)).*histcounts(wrapTo2Pi(theta_phase2),phase_edges);
    %     phase_occupancy_map3 = mean(diff(LFP_tvec)).*histcounts(wrapTo2Pi(theta_phase3),phase_edges);

    % 2D occupancy map for spatial bin and phase bin



    % For each good cell in the track
    for ncell = 1 : length(good_cell_index)
        tic
        %% Find location for each spike independently for either direction

        theta_modulation(track).cluster_id(ncell) = place_fields(1).cluster_id(good_cell_index(ncell));
        theta_modulation(track).peak_percentile(ncell) = place_fields(track).peak_percentile(good_cell_index(ncell));
        theta_modulation(track).odd_even_stability(ncell) = place_fields(track).odd_even_stability(good_cell_index(ncell));
        theta_modulation(track).t1_t2_remapping(ncell) = place_fields(track).t1_t2_remapping(good_cell_index(ncell));
        
        place_cell_idx   = clusters.spike_id == place_fields(1).cluster_id(good_cell_index(ncell));  % find indices for this unit
        spike_times = clusters.spike_times(place_cell_idx);  % extract spike times for this place cell

        spike_position = interp1(track_times,track_linear,spike_times,'nearest');
        spike_phases = interp1(LFP_tvec,theta_phase_unwrap,spike_times,'linear');
        spike_phases = wrapTo2Pi(spike_phases);

        spike_speed = interp1(Behaviour.tvec,Behaviour.speed,spike_times,'nearest');

        spike_position = spike_position(spike_speed>5);
        spike_times = spike_times(spike_speed>5);
        spike_phases = spike_phases(spike_speed>5);

        [spike_count_map,Xedges,Yedges,binX,binY] = histcounts2(spike_position,spike_phases,position_edges,phase_edges);


        position_phase_map = spike_count_map./occupancy_map;
        position_phase_map(isinf(position_phase_map))=0;
        position_phase_map(isnan(position_phase_map))=0;


        gaussianWindow = gausswin(8/mean(diff(position_edges))*2.5*2+1);% 4 datapoint or 8cm position bin for one SD (alpha is fixed at 2.5)
        gaussianWindow = gaussianWindow/sum(gaussianWindow);

        theta_modulation(track).position_phase_map{ncell} = imgaussfilt(position_phase_map,[4,2]);
        theta_modulation(track).theta_phase_map{ncell} = histcounts(spike_phases,phase_edges)./phase_occupancy_map;
        theta_modulation(track).theta_modulation_index(ncell) = (max(theta_modulation(track).theta_phase_map{ncell})-min(theta_modulation(track).theta_phase_map{ncell}))/mean(theta_modulation(track).theta_phase_map{ncell});

        theta_modulation(track).position_phase_xcorr_map{ncell} = [];
        theta_modulation(track).position_phase_map_1D_smoothed{ncell}=[];

        average_map = filtfilt(gaussianWindow,1,mean(position_phase_map,2));

        for nphase = 1:length(phase_edges)-1
            % xxx(nphase,:) = filtfilt(gaussianWindow,1,position_phase_map(:,nphase));
            temp_map = filtfilt(gaussianWindow,1,position_phase_map(:,nphase));
            %             clip_size = (length(temp_map)-size(position_phase_map,1))/2;
            theta_modulation(track).position_phase_map_1D_smoothed{ncell}(:,nphase) = temp_map;
            [r,lags] = xcorr(temp_map-mean(temp_map)',average_map-mean(average_map)'); % zero mean for cross correlation


            %         r(lags<=80 & lags>=-80)
            theta_modulation(track).position_phase_xcorr_map{ncell} = [theta_modulation(track).position_phase_xcorr_map{ncell} r/max(r)];
        end

        [~,peak_index] = max(theta_modulation(track).position_phase_xcorr_map{ncell});
        % Define a sinusoidal model
        sinusoid = fittype('a*sin(b*x + c) + d', 'independent', 'x');

        % Set initial guess for parameters a (amplitude), b (frequency), c (phase offset), and d (vertical offset)
        bin_centres = phase_edges(1:end-1) + diff(phase_edges)/2;
        bin_centres = 180*bin_centres'./pi;

        startPoints = [max(peak_index)-min(peak_index), 180*mean(diff(bin_centres))./pi, 0, mean(peak_index)];

        % Fit the sinusoid to your data
        coeffs = fit(bin_centres, peak_index', sinusoid, 'Start', startPoints);

        % Extract the amplitude and phase offset
        amplitude = coeffs.a;
        phase_offset = coeffs.c;

%         % Plot the original data and the fit
%         figure;
%         plot(bin_centres, peak_index, 'bo');
%         hold on;
%         fitted_values = coeffs.a*sin(coeffs.b*bin_centres + coeffs.c) + coeffs.d;
%         plot(bin_centres, fitted_values, 'r-');
%         xlabel('Theta Phase Bin');
%         ylabel('Maximum Index');
%         legend('Data', 'Sinusoidal Fit');

        % Save variables
        theta_modulation(track).spike_times_phase_position{ncell} = [spike_times(spike_position<=139.5),spike_phases(spike_position<=139.5),spike_position(spike_position<=139.5)];
        theta_modulation(track).sin_fit_phase_offset(ncell) = coeffs.c;
        theta_modulation(track).sin_fit_amplitude(ncell) = coeffs.a;

        parfor nshuffle = 1:500
            s = RandStream('mrg32k3a','Seed',1000+nshuffle+ncell*100); % Set random seed for resampling

            if rand(s,1)>0.5
                s = RandStream('mrg32k3a','Seed',2000+nshuffle+ncell*100)
                shift_dir = 6+10*rand(s,1);
            else
                s = RandStream('mrg32k3a','Seed',2000+nshuffle+ncell*100)
                shift_dir = -6-10*rand(s,1);
            end
            spike_times_shuffled = spike_times + shift_dir;
            spike_times_shuffled(spike_times_shuffled>max(spike_times)) = spike_times_shuffled(spike_times_shuffled>max(spike_times))-max(spike_times)+min(spike_times);
            spike_phases_shuffled = interp1(LFP_tvec,theta_phase_unwrap,spike_times_shuffled,'linear');
            spike_phases_shuffled = wrapTo2Pi(spike_phases_shuffled);

            [spike_count_map_shuffled,Xedges,Yedges,binX,binY] = histcounts2(spike_position,spike_phases_shuffled,position_edges,phase_edges);
            phase_map_shuffled = histcounts(spike_phases_shuffled,phase_edges)./phase_occupancy_map;
            theta_modulation_shuffled(nshuffle) = (max(phase_map_shuffled)-min(phase_map_shuffled))/mean(phase_map_shuffled);


            position_phase_map_shuffled = spike_count_map_shuffled./occupancy_map;
            position_phase_map_shuffled(isinf(position_phase_map_shuffled))=0;
            position_phase_map_shuffled(isnan(position_phase_map_shuffled))=0;

            position_phase_xcorr_map_shuffled = [];
            for nphase = 1:length(phase_edges)-1

                temp_map = filtfilt(gaussianWindow,1,position_phase_map_shuffled(:,nphase));
                %             clip_size = (length(temp_map)-size(position_phase_map,1))/2;
                [r,lags] = xcorr(position_phase_map_shuffled(:,nphase),mean(position_phase_map_shuffled));
                %         r(lags<=80 & lags>=-80)
                position_phase_xcorr_map_shuffled = [ position_phase_xcorr_map_shuffled r/max(r)];
            end

            [~,peak_index] = max(position_phase_xcorr_map_shuffled);
            % Define a sinusoidal model
            sinusoid = fittype('a*sin(b*x + c) + d', 'independent', 'x');

            % Set initial guess for parameters a (amplitude), b (frequency), c (phase offset), and d (vertical offset)
            startPoints = [max(peak_index)-min(peak_index), 2*pi/mean(diff(bin_centres)), 0, mean(peak_index)];

            % Fit the sinusoid to your data
            coeffs = fit(bin_centres, peak_index', sinusoid, 'Start', startPoints);

            % Extract the amplitude and phase offset
            amplitude_shuffled(nshuffle) = coeffs.a;
            %             phase_offset_shuffled(nshuffle) = coeffs.c;
        end

        theta_modulation(track).phase_offset_percentile(ncell) = coeffs.c;
        theta_modulation(track).phase_amplitude_percentile(ncell) = sum(theta_modulation(track).sin_fit_amplitude(ncell)>amplitude_shuffled)/500;
        theta_modulation(track).theta_modulation_percentile(ncell) = sum(theta_modulation(track).theta_modulation_index(ncell)>theta_modulation_shuffled)/500;

        
                if count == 1
                    fig = figure;
                    fig.Position = [273 234 1040 720];
                    %             fig.Name = [273 234 1040 720];
        
                elseif count == 3
                    count = 1;
                    fig = figure;
                    fig.Position = [273 234 1040 720];
                    %             fig.Name = [273 234 1040 720];
        
                end
        
                subplot(3,3,count)
                imagesc(position_edges(1:end-1)+1,180/pi.*phase_edges(2:end),theta_modulation(track).position_phase_map{ncell}')
                colorbar
                colormap(flip(gray))
                xticks(position_edges(1:10:end-1)+1)
                yticks(180/pi.*phase_edges(2:3:end))
                title(sprintf('Track %i cell %i %s',track,theta_modulation(1).cluster_id(good_cell_index(ncell)),clusters.region(clusters.cluster_id==place_fields(1).cluster_id(good_cell_index(ncell)))))
                xlabel('Spatial location')
                ylabel('Theta phase')
        
                subplot(3,3,count+3)
                polarhistogram(theta_modulation(track).spike_times_phase_position{ncell}(:,2),0:2*pi/18:2*pi)
        
                subplot(3,3,count+6)
                lags= -69:69;
                imagesc(lags,180/pi.*phase_edges(2:end), theta_modulation(track).position_phase_xcorr_map{ncell}')
                hold on;xline(0,'r')
                colorbar
                colormap(flip(gray))
                xticks(lags(1:10:end-1)+1)
                yticks(180/pi.*phase_edges(2:3:end))
                xlabel('Spatial shift')
                ylabel('Theta phase')
                count = count + 1;

        %% phase precession vs location on track

        % circular-linear correlation coefficent from CircStat MATLABToolbox (Berens, 2009). Method first described in Kempter et al 2012
        % output of circ_corrcl is correlation coefficient and pval
        % needs input in radians
        %         if ~isempty(place_cell_times) % only useful for hippocampal cell
        %             [TPP(track).circ_lin_corr(ncell),TPP(track).circ_lin_PVAL(ncell)] = ...
        %                 circ_corrcl(TPP(track).spike_phases{ncell},TPP(track).spike_positions{ncell});
        %         end


        toc
    end

end

% save('extracted_phase_precession_absolute_location','TPP','half_laps_times')
%

end




function raster_plot(spike_times,y,col,height)

x2(1:3:length(spike_times)*3)=spike_times;
x2(2:3:length(spike_times)*3)=spike_times;
x2(3:3:length(spike_times)*3)=NaN;
y2(1:3:length(spike_times)*3)=y;
y2(2:3:length(spike_times)*3)=y+height;
y2(3:3:length(spike_times)*3)=NaN;
if isempty(col)
    plot(x2,y2,'linewidth',2);
else
    plot(x2,y2,'color',col,'linewidth',2);
end
end

function half_laps_times = extract_running_laps(position,lap_times)

parameters = list_of_parameters;

for track = 1 : length(lap_times)

    % Get half lap start and end time
    half_laps_timestamps = [lap_times(track).halfLaps_start' lap_times(track).halfLaps_stop'];

    % Split into 2 directions
    direction1 = half_laps_timestamps([1:2:size(half_laps_timestamps,1)],:);
    direction2 = half_laps_timestamps([2:2:size(half_laps_timestamps,1)],:);

    % Find these times in the position.data and get the indices
    direction1_idx = interp1(position.linear(track).timestamps,1:length(position.linear(track).timestamps),direction1,'nearest');
    direction2_idx = interp1(position.linear(track).timestamps,1:length(position.linear(track).timestamps),direction2,'nearest');

    % Turn this into logical index
    dir1_idx = zeros(length(position.linear(track).timestamps),1);
    for n = 1:size(direction1,1)
        dir1_idx(direction1_idx(n,1):direction1_idx(n,2)) = 1;
    end
    half_laps_times(track).direction_idx_1 = logical(dir1_idx);

    dir2_idx = zeros(length(position.linear(track).timestamps),1);
    for n = 1:size(direction2,1)
        dir2_idx(direction2_idx(n,1):direction2_idx(n,2)) = 1;
    end
    half_laps_times(track).direction_idx_2 = logical(dir2_idx);

    % Filter by speed
    speed_thresh = parameters.speed_threshold;   % arbitrarily chosen
    running_idx  =  position.v_cm(position.linear(track).clean_track_Indices) > speed_thresh;

    % remove portions along laps where animal not running
    half_laps_times(track).direction_idx_1(not(running_idx)) = 0;
    half_laps_times(track).direction_idx_2(not(running_idx)) = 0;

    % Save other variables
    half_laps_times(track).running_idx = running_idx;
    half_laps_times(track).direction1 = direction1;
    half_laps_times(track).direction2 = direction2;

end
end