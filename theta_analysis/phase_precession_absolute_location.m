function [TPP, half_laps_times] = phase_precession_absolute_location(LFP_tvec,CA1_LFP,place_fields,clusters,Task_info,Behaviour,options)

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

%Get timestamps for each half lap
half_laps_times = [];
        
% Find hilbert transform to extract phase info
hilb = hilbert(CA1_LFP);
theta_phase = angle(hilb);
theta_phase_unwrap = unwrap(theta_phase); % unwrap for interpolation

spatial_cell_index = unique([find(place_fields(1).peak_percentile>0.95 & place_fields(1).odd_even_stability>0.95)...
    find(place_fields(2).peak_percentile>0.95 & place_fields(2).odd_even_stability>0.95)]);
good_cell_index = intersect(spatial_cell_index,find(ismember(place_fields(1).cluster_id,clusters.cluster_id)));

for track = 1 : length(place_fields)
    
    track_linear = Behaviour.position;
    track_linear(Behaviour.track_ID~=track) = nan;

    track_times = Behaviour.tvec;

    % For each good cell in the track 
    for ncell = 1 : length(good_cell_index)
        
        %% Find location for each spike independently for either direction
        
        TPP(track).cell_id(ncell) = place_fields(1).cluster_id(good_cell_index(ncell)); 
        place_cell_idx   = clusters.spike_id == place_fields(1).cluster_id(good_cell_index(ncell));  % find indices for this unit
        place_cell_times = clusters.spike_times(place_cell_idx);  % extract spike times for this place cell
        place_cell_times_track = place_cell_times(place_cell_times > track_times(1) & place_cell_times < track_times(end));
        
        % find the indices of the spike times in the position time frame; not the most acurate but use to filter spikes
        position_t_idx = interp1(track_times,1:length(track_times),place_cell_times_track,'nearest');
                
        % remove any spike times where position is not known
        position_t_idx(isnan(position_t_idx)) = [];

        place_cell_times = place_cell_times_track;     % spike times

        % Save variables
        TPP(track).place_cell_times_track{ncell} = place_cell_times_track;
        TPP(track).place_cell_times{ncell} = place_cell_times;

        
        
        %% Find spike phase and position for each direction

        %  find phase for each spike      
        spike_phases = interp1(LFP_tvec,theta_phase_unwrap,place_cell_times,'linear');
%         polarhistogram(spike_phases,30)
        TPP(track).spike_positions{ncell} = interp1(track_times,track_linear,place_cell_times,'linear');
        TPP(track).spike_phases{ncell} = wrapToPi(spike_phases);

        %% population analysis for phase precession vs location on track
        
        % circular-linear correlation coefficent from CircStat MATLABToolbox (Berens, 2009). Method first described in Kempter et al 2012
        % output of circ_corrcl is correlation coefficient and pval
        % needs input in radians
        if ~isempty(place_cell_times)
            [TPP(track).circ_lin_corr(ncell),TPP(track).circ_lin_PVAL(ncell)] = ...
                circ_corrcl(TPP(track).spike_phases{ncell},TPP(track).spike_positions{ncell});
        end
        
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