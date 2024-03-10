function plot_theta_modulation(clusters,place_fields,theta_modulation_L,theta_modulation_R,options)

if ~isempty(theta_modulation_L) & ~isempty(theta_modulation_R)
    no_of_probes = 2;
elseif ~isempty(theta_modulation_L) | ~isempty(theta_modulation_R)
    no_of_probes = 1;
end

spatial_cell_index = unique([find(place_fields(1).peak_percentile>0.95 & place_fields(1).odd_even_stability>0.95)...
    find(place_fields(2).peak_percentile>0.95 & place_fields(2).odd_even_stability>0.95)]);

good_cell_index = intersect(spatial_cell_index,find(ismember(place_fields(1).cluster_id,clusters.cluster_id)));

% good_cell_index = good_cell_index(1:2);
sprintf('Theta modulation analysis for %i spatial cells',length(good_cell_index))

position_edges = 0:2:140;
phase_edges = 0:20:360;
phase_edges = pi/180.*phase_edges;

for ncell = 1 : length(good_cell_index)
    count = 1;
    %     for nprobe = 1:no_of_probes
    for track = 1 : length(place_fields)
        % For each good cell in the track


        if count == 1
            fig = figure;
            fig.Position = [273 234 1040 720];
            fig.Name = [273 234 1040 720];

        elseif count == 3
            count = 1;
            fig = figure;
            fig.Position = [273 234 1040 720];
            fig.Name = [273 234 1040 720];
        end

        if ~isempty(theta_modulation_L)
            cell_index = find(theta_modulation_L(1).cluster_id==place_fields(1).cluster_id(good_cell_index(ncell)));

            subplot(3,4,count)
            imagesc(position_edges(1:end-1)+1,180/pi.*phase_edges(2:end),theta_modulation_L(track).position_phase_map{cell_index}')
            colorbar
            colormap(flip(gray))
            xticks(position_edges(1:10:end-1)+1)
            yticks(180/pi.*phase_edges(2:3:end))
            title(sprintf('Track %i cell %i %s',track,theta_modulation_L(1).cluster_id(cell_index),clusters.region(clusters.cluster_id==place_fields(1).cluster_id(cell_index))))
            xlabel('Spatial location')
            ylabel('Theta phase')

            subplot(3,3,count+3)
            polarhistogram(theta_modulation_L(track).spike_times_phase_position{ncell}(:,2),0:2*pi/18:2*pi)

            subplot(3,3,count+6)
            lags= -69:69;
            imagesc(lags,180/pi.*phase_edges(2:end), theta_modulation_L(track).position_phase_xcorr_map{cell_index}')
            hold on;xline(0,'r')
            colorbar
            colormap(flip(gray))
            xticks(lags(1:10:end-1)+1)
            yticks(180/pi.*phase_edges(2:3:end))
            xlabel('Spatial shift')
            ylabel('Theta phase')
            count = count + 1;
        end

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
% end
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