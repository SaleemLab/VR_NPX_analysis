% GET PLACE FIELDS WITHIN LAPS
% loads list_of_parameters.m, extracted_clusters.mat,extracted_position.mat, extracted_waveform.mat
% uses function extract_laps.m, skaggs_information.m
% INPUTS: 
    % track_id: track being analysed
    % lap_start & lap_end: set laps to analyise by inputing start and end lap IDs 
    % bayesian_option: enter 1 or 0 for calculating place fields with wider x_bin required for bayesian decoding 

function lap_directional_place_fields = get_directional_lap_place_fields(track_id,lap_start,lap_end,bayesian_option,lap_option)
    
    load('extracted_directional_clusters.mat');
    load('extracted_clusters.mat');
    load('extracted_laps.mat');
    load('extracted_position.mat');
    if exist('extracted_waveforms.mat','file')
        load('extracted_waveforms.mat');
    else
        disp('no extracted_waveforms.mat file');
        allclusters_waveform=[];
    end

    % Set parameters
    parameters=list_of_parameters;
    if bayesian_option == 0
        x_bins_width = parameters.x_bins_width;
        w= gausswin(parameters.place_field_smoothing);       
    else                                                    % If calculating place fields to use in bayesian decoding
        x_bins_width = parameters.x_bins_width_bayesian;  
        w= [1 1]; %moving average filter of 2 sample, will be become a filter of [0.25 0.5 0.25] with filtfilt
    end
    w = w./sum(w); %make sure smoothing filter sums to 1 
    
    % Run threshold on pyramidal cells: half-width amplitude
    if ~isempty(allclusters_waveform)       
        PC_indices = [allclusters_waveform.half_width] > parameters.half_width_threshold; % cells that pass treshold of pyramidal cell half width
        pyramidal_cells = [allclusters_waveform(PC_indices).converted_ID];
    end
     
    
    for dir = 1 : 2 %for each direction
        %Get mean firing rate from original data set to be able to find interneurons
        if strcmp(bayesian_option,'Y')
            load extracted_directional_place_fields_BAYESIAN.mat
            place_fields.mean_rate = directional_place_fields_BAYESIAN(dir).place_fields.mean_rate;
            clear directional_place_fields_BAYESIAN
        else
            load extracted_directional_place_fields.mat
            place_fields.mean_rate = directional_place_fields(dir).place_fields.mean_rate;
            clear directional_place_fields
        end
        
    % Place field calculation
    place_fields.time_window = [lap_times(track_id).(sprintf('%s',lap_option,'Laps_start'))(lap_start) lap_times(track_id).(sprintf('%s',lap_option,'Laps_stop'))(lap_end)];
    time_bin_width = position.t(2)-position.t(1); %sets time_bin 
   
    position_index = isnan(position.linear(track_id).linear);
    position_speed = abs(position.v_cm);
    position_speed(position_index) = NaN;  %make sure speed is NaN if position is NaN
    
    position_during_spike = interp1(position.t,position.linear(track_id).linear,directional_clusters(track_id).(strcat('spikes_dir',num2str(dir)))(:,2),'nearest');  %interpolates position into spike time
    speed_during_spike = interp1(position.t,position_speed,directional_clusters(track_id).(strcat('spikes_dir',num2str(dir)))(:,2),'nearest');
    
    x_bin_edges = 0:x_bins_width:100*position.linear(track_id).length; % forces x_bins to be from 0 to 200cm
    x_bin_centres = (x_bin_edges(2)-x_bins_width/2):x_bins_width:(x_bin_edges(end-1)+x_bins_width/2);   
    x_bins = 0:x_bins_width:(100*position.linear(track_id).length); %bin position
    
    % Time spent at each x_bin (speed filtered)
    x_hist = time_bin_width.*histcounts(position.linear(track_id).linear(find(position.t>place_fields.time_window(1) &...
        position.t<place_fields.time_window(2) & position_speed>parameters.speed_threshold_laps...
        & position_speed<parameters.speed_threshold_max)),x_bin_edges); 
    
    place_fields.x_bin_centres = x_bin_centres;
    place_fields.x_bin_edges = x_bin_edges;
    place_fields.x_bins = x_bins;
    place_fields.x_bins_width = x_bins_width;
    place_fields.dwell_map = x_hist;

    for j=1:max(clusters.id_conversion(:,1))
        % Number of spikes per bin within time window (speed filtered)
        place_fields.spike_hist{j} = histcounts(position_during_spike(find(directional_clusters(track_id).(strcat('spikes_dir',num2str(dir)))(:,1) == j & ...
            directional_clusters(track_id).(strcat('spikes_dir',num2str(dir)))(:,2) > place_fields.time_window(1) & ...
            directional_clusters(track_id).(strcat('spikes_dir',num2str(dir)))(:,2) < place_fields.time_window(2) & ...
            speed_during_spike>parameters.speed_threshold_laps &...
            speed_during_spike<parameters.speed_threshold_max)),x_bin_edges);
        
        place_fields.raw{j} = place_fields.spike_hist{j}./x_hist; % place field calculation
        place_fields.raw{j}(find(isnan(place_fields.raw{j})==1)) = 0; %replace NaNs for zeros
        
        % zero bins with 0 dwell time, but make sure no spikes occurred
        non_visited_bins = find(x_hist==0);
        if sum(place_fields.spike_hist{j}(non_visited_bins))>0
            disp('ERROR: x_hist is zero, but spike histogram is not');
        else
            place_fields.raw{j}(non_visited_bins)= 0;
        end
        place_fields.non_visited_bins = non_visited_bins; %NaNs that have been replaced by O
        
        % Get place field information
        place_fields.smooth{j}         = filtfilt(w,1,place_fields.raw{j}); %smooth pl field
        place_fields.centre_of_mass(j) = sum(place_fields.smooth{j}.*x_bin_centres)/sum(place_fields.smooth{j}); %averaged center
        [place_fields.peak(j) , index] = max(place_fields.smooth{j}); %peak of smoothed place field and index of peak (center)
        place_fields.centre(j)   =  x_bin_centres(index);
        place_fields.raw_peak(j) = max(place_fields.raw{j}); % raw pl field peak
        place_fields.mean_rate_lap(j)  = sum(place_fields.spike_hist{j})/(place_fields.time_window(2)-place_fields.time_window(1));
        place_fields.half_max_width(j) = x_bins_width*half_max_width(place_fields.smooth{j}); %finds half width of smoothed place field
    end
    
    %calculate skagges information
    place_fields.skaggs_info= skaggs_information(place_fields);

    % Find cells that pass the 'Place cell' thresholds -
    % both peak of smoothed place field or peak of raw place field need to be above the respective thresholds
    putative_place_cells = find((place_fields.peak >= parameters.min_smooth_peak...
        & place_fields.raw_peak >= parameters.min_raw_peak)...
        & place_fields.mean_rate_lap <= parameters.max_mean_rate...
        & place_fields.mean_rate <= parameters.max_mean_rate...
        & place_fields.skaggs_info > 0);
    
    % Set a less conservative criteria for place cells, having to pass either peak firing rate thresholds (smoothed PF and raw PF)
    putative_place_cells_LIBERAL = find(place_fields.peak >= parameters.min_smooth_peak...
        | place_fields.raw_peak >= parameters.min_raw_peak...
        & place_fields.mean_rate_lap <= parameters.max_mean_rate...
        & place_fields.mean_rate <= parameters.max_mean_rate...
        & place_fields.skaggs_info > 0);
    
    if ~isempty(allclusters_waveform)
        place_fields.good_cells = intersect(putative_place_cells,pyramidal_cells); % Check that the cells that passed the threshold are pyramidal cells
        place_fields.good_cells_LIBERAL = intersect(putative_place_cells_LIBERAL,pyramidal_cells); 
    else
        place_fields.good_cells = putative_place_cells;
        place_fields.good_cells_LIBERAL = putative_place_cells_LIBERAL;
    end
    
    % Sort place fields to the location of their peak
    [~,index] = sort(place_fields.centre);
    place_fields.sorted = index;
    [~,index1] = sort(place_fields.centre(place_fields.good_cells));
    place_fields.sorted_good_cells =  place_fields.good_cells(index1);
    [~,index2] = sort(place_fields.centre(place_fields.good_cells_LIBERAL));
    place_fields.sorted_good_cells_LIBERAL = place_fields.good_cells_LIBERAL(index2);
    
    %%% Classify other cells as other cells, interneurons and pyramidal cells
    
    interneurons=[]; %interneurons classfication
    interneurons = [interneurons find(place_fields.mean_rate > parameters.max_mean_rate)];        
    place_fields.interneurons = unique(interneurons);
    
    putative_pyramidal_cells=[];  % putative pyramidal cells that pass the 'Pyramidal type' threshold (but not need to be place cells)
    putative_pyramidal_cells = [putative_pyramidal_cells find(place_fields.peak >= place_fields.mean_rate_lap <= parameters.max_mean_rate)];
    place_fields.pyramidal_cells = intersect(putative_pyramidal_cells,pyramidal_cells);
    
    other_cells = setdiff(1:max(clusters.id_conversion(:,1)),place_fields.good_cells,'stable'); %find the excluded putative pyramidal cells
    place_fields.other_cells = setdiff(other_cells,interneurons,'stable'); %remove also the interneurons
    
     lap_directional_place_fields(dir).place_fields = place_fields;

    end
end

function half_width = half_max_width(place_field)
 
    %interpolate place field to get better resolution
    new_step_size = 0.1;  %decrease value to get finer resolution interpolation of place field
    place_field_resampled = interp1(1:length(place_field),place_field,1:new_step_size:length(place_field),'linear');

    [peak,index] = max(place_field_resampled); %finds smoothed place field peak firing rate (FR)
    for i = index : length(place_field_resampled)
        if place_field_resampled(i)<peak/2 %finds the point after the peak where the FR is half the peak FR
            break;
        end
    end
    for j = index : -1 : 1 %finds the point before the peak where the FR is half the peak FR
        if place_field_resampled(j)<peak/2
            break;
        end
    end
    half_width = new_step_size*(i-j); %distance between half-peaks
    %(calculated in indicies of original place field, but converted to distance in cm in function above)

end