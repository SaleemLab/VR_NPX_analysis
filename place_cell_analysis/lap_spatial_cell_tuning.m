function spatial_tunning_lap = lap_spatial_cell_tuning(clusters,place_fields_all,place_fields_even,place_fields_odd,position,lap_times,x_bins_width,options)
% Function to plot spatial cell activity
% Plot 1 is each cell's spatial ratemap for each lap
% Plot 2 is each cell's spatial tuning stability by plotting place field
% calculated by using all laps, odd no laps and even no laps.

% x_bins_width = 5;

%             clusters = HPC_clusters;
% clusters = CA1_clusters;
%             clusters = V1_clusters;

%     place_fields_even = CA1_place_fields_even.probe(nprobe);
%     place_fields_odd = CA1_place_fields_odd.probe(nprobe);
%     place_fields_all = CA1_place_fields.probe(nprobe);
if isfield(options,'probe_combined') && options.probe_combined == 1
    nprobe = 'combined';
else
    nprobe = options.probe_no;
end

% Grabbing ratemap for each lap for each cell
tic
spatial_tunning_lap = [];
if ~ischar(nprobe)
    for track_id = 1:2
        for nlap = 1:length(lap_times(track_id).start)
            field_this_lap = get_lap_place_fields_masa(x_bins_width,position,place_fields_all,...
                clusters.probe(nprobe),track_id,lap_times(track_id).start(nlap),lap_times(track_id).end(nlap));

            spatial_tunning_lap{track_id}{nlap} = field_this_lap.raw;
        end
    end
else
    for track_id = 1:2
        for nlap = 1:length(lap_times(track_id).start)
            field_this_lap = get_lap_place_fields_masa(x_bins_width,position,place_fields_all,...
                clusters,track_id,lap_times(track_id).start(nlap),lap_times(track_id).end(nlap));

            spatial_tunning_lap{track_id}{nlap} = field_this_lap.raw;
        end
    end
end