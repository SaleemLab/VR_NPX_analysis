    
function tracks_DecodingError = bayesian_decoding_error_stability
% Looks at place cell stability using bayesian decoding error. For each track, decodes each lap using the 4 laps (or to the last lap in
% short-exposures (T2)), and gets the mean and median decoding error (per lap and as a single value).
% MH, 06.02.20


    % Set parameters and load files
    parameters=list_of_parameters;
    load('extracted_clusters.mat');
    load('extracted_position.mat');
    load('extracted_laps.mat');
    
    for track = 1: length(lap_times)
        % Template laps
        lap_start(1) = 2;
        lap_end(1) = 10;
        
        % Decoded laps
        lap_start(2) = 1;
        lap_end(2) = 8;
        
        % Time window to be decoded
        timeWindow_decodedLaps = [lap_times(track).start(lap_start(2)), lap_times(track).end(lap_end(2))];
        
        % Calculates place field inside template lap(s) selected
        plFields_template.track = get_lap_place_fields_masa(track,lap_start(1),lap_end(1),1); %using Bayesian place fields
        
       % Bayesian decoding
        bayesian_spike_count = spike_count(plFields_template,timeWindow_decodedLaps(1),timeWindow_decodedLaps(2),0,[]);
        estimated_position = bayesian_decoding(plFields_template,bayesian_spike_count,0,[]); 
        bayesian_decodingError = decoding_error(estimated_position,plFields_template,track,track,timeWindow_decodedLaps(1),timeWindow_decodedLaps(2),'N');    % calculate decoding error
       
        out(track).tracks_compared = mat2str([track track]);
        out(track).laps_jump = 1; 
        out(track).decoded_laps = [lap_start(2) lap_end(2)];
        out(track).template_laps = [lap_start(1) lap_end(1)];
        out(track).decoding_error = bayesian_decodingError;

    end
    
    tracks_DecodingError = out;    
    
end
