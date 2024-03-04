function [TPP, half_laps_times] = phase_precession_absolute_location(CA1_LFP,)

load('extracted_CSC.mat')
load('extracted_place_fields.mat')
load('extracted_clusters.mat')
load('extracted_position.mat')
load('extracted_laps.mat')

%Get timestamps for each half lap
half_laps_times = extract_running_laps(position,lap_times);
        
% Find hilbert transform to extract phase info
hilb = hilbert(CSC(4).theta);
theta_phase = angle(hilb);
theta_phase_unwrap = unwrap(theta_phase); % unwrap for interpolation


for track = 1 : length(place_fields.track)
    
    track_linear = position.linear(track).linear(~isnan(position.linear(track).linear));
    track_times = position.linear(track).timestamps;
    track_x = position.x(position.linear(track).clean_track_Indices);
    track_y = position.y(position.linear(track).clean_track_Indices);

    % For each good cell in the track 
    for pc = 1 : length(place_fields.track(track).good_cells)
        
        %% Find location for each spike independently for either direction
        
        TPP(track).cell_id(pc) = place_fields.track(track).good_cells(pc); 
        place_cell_idx   = clusters.spike_id == place_fields.track(track).good_cells(pc);  % find indices for this unit
        place_cell_times = clusters.spike_times(place_cell_idx);  % extract spike times for this place cell
        place_cell_times_track = place_cell_times(place_cell_times > track_times(1) & place_cell_times < track_times(end));
        
        % find the indices of the spike times in the position time frame; not the most acurate but use to filter spikes
        position_t_idx = interp1(track_times,1:length(track_times),place_cell_times_track,'nearest');
                
        % remove any spike times where position is not known
        position_t_idx(isnan(position_t_idx)) = [];

        % separate spikes for each direction
        % below gives 1 if spike occurs when direction_idx is 1, gives 0 if direction_idx is 0
        direction_test_1 = logical(half_laps_times(track).direction_idx_1(position_t_idx));
        direction_test_2 = logical(half_laps_times(track).direction_idx_2(position_t_idx));
        
        place_cell_times_1 = place_cell_times_track;     % spike times
        place_cell_times_2 = place_cell_times_track;
        place_cell_times_1(not(direction_test_1)) = [];  % removes the spike that don't occur in required direction
        place_cell_times_2(not(direction_test_2)) = [];
        
        % Save variables
        TPP(track).place_cell_times_track{pc} = place_cell_times_track;
        TPP(track).place_cell_times_1{pc} = place_cell_times_1;
        TPP(track).place_cell_times_2{pc} = place_cell_times_2;
        
        
        %% Find spike phase and position for each direction

        %  find phase for each spike      
        spike_phases_1 = interp1(CSC(1).CSCtime,theta_phase_unwrap,place_cell_times_1,'linear');
        TPP(track).spike_positions_1{pc} = interp1(track_times,track_linear,place_cell_times_1,'linear');
        TPP(track).spike_x_1{pc}         = interp1(track_times,track_x,place_cell_times_1,'linear');
        TPP(track).spike_y_1{pc}         = interp1(track_times,track_y,place_cell_times_1,'linear');
        
        spike_phases_2 = interp1(CSC(1).CSCtime,theta_phase_unwrap,place_cell_times_2,'linear');
        TPP(track).spike_positions_2{pc} = interp1(track_times,track_linear,place_cell_times_2,'linear');
        TPP(track).spike_x_2{pc}         = interp1(track_times,track_x,place_cell_times_2,'linear');
        TPP(track).spike_y_2{pc}         = interp1(track_times,track_y,place_cell_times_2,'linear');
        
        TPP(track).spike_phases_1{pc} = wrapToPi(spike_phases_1);
        TPP(track).spike_phases_2{pc} = wrapToPi(spike_phases_2);

        %% population analysis for phase precession vs location on track
        
        % circular-linear correlation coefficent from CircStat MATLABToolbox (Berens, 2009). Method first described in Kempter et al 2012
        % output of circ_corrcl is correlation coefficient and pval
        % needs input in radians
        if ~isempty(place_cell_times_1)
        [TPP(track).circ_lin_corr_dir1(pc),TPP(track).circ_lin_PVAL_dir1(pc)] = ...
            circ_corrcl(TPP(track).spike_phases_1{pc},TPP(track).spike_positions_1{pc});
        end
        if ~isempty(place_cell_times_2)
        [TPP(track).circ_lin_corr_dir2(pc),TPP(track).circ_lin_PVAL_dir2(pc)]= ...
            circ_corrcl(TPP(track).spike_phases_2{pc},TPP(track).spike_positions_2{pc});
        end
        
    end
    
end

save('extracted_phase_precession_absolute_location','TPP','half_laps_times')



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