function extract_place_field_lap(bayesian_option)
% Extract place fields for each complete lap and each half lap. 
% MH_ 2020

load('extracted_laps.mat')

for track = 1 : length(lap_times)
    
    for i = 1 : length(lap_times(track).completeLaps_start)
       
        disp([num2str(i) ' out of ' num2str(length(lap_times(track).completeLaps_start))])

        % Extract place field of each complete lap
        lap_start_time = lap_times(track).completeLap_id(i);
        lap_end_time = lap_times(track).completeLap_id(i);
        
        place_fields = get_lap_place_fields(track,lap_start_time,lap_end_time,bayesian_option,'complete');
        directional_place_fields = get_directional_lap_place_fields(track,lap_start_time,lap_end_time,bayesian_option,'complete');

        lap_place_fields(track).Complete_Lap{i} = place_fields;
        lap_directional_place_fields(track).dir1.Complete_Lap{i} =  directional_place_fields(1).place_fields;
        lap_directional_place_fields(track).dir2.Complete_Lap{i} =  directional_place_fields(2).place_fields;

        clear place_fields directional_place_fields

    end
    
    for i = 1 : length(lap_times(track).halfLaps_start)
        disp([num2str(i) ' out of ' num2str(length(lap_times(track).halfLaps_start))])
         % Extract place field of each half lap
        lap_start_time = lap_times(track).halfLap_id(i);
        lap_end_time = lap_times(track).halfLap_id(i);
        
        place_fields = get_lap_place_fields(track,lap_start_time,lap_end_time,bayesian_option,'half');
        directional_place_fields = get_directional_lap_place_fields(track,lap_start_time,lap_end_time,bayesian_option,'half');

        lap_place_fields(track).half_Lap{i} = place_fields;
        lap_directional_place_fields(track).dir1.half_Lap{i} =  directional_place_fields(1).place_fields;
        lap_directional_place_fields(track).dir2.half_Lap{i} =  directional_place_fields(2).place_fields;

        clear place_fields directional_place_fields
    end
    
end 

%SAVE
    if bayesian_option == 1
        lap_place_fields_BAYESIAN = lap_place_fields;
        lap_directional_place_fields_BAYESIAN = lap_directional_place_fields;
        save('extracted_lap_place_fields_BAYESIAN','lap_place_fields_BAYESIAN','-v7.3')
        save('extracted_directional_lap_place_fields_BAYESIAN.mat','lap_directional_place_fields_BAYESIAN','-v7.3')
    else
        save('extracted_lap_place_fields.mat','lap_directional_place_fields','-v7.3')
        save('extracted_directional_lap_place_fields.mat','lap_directional_place_fields','-v7.3')
    end



end