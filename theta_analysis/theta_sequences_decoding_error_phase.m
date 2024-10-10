% THETA SEQUENCE DECODING ERRROR
% MH, 2020
% Finds real position at the time of the theta sequence, and calculates the decoding error during the theta sequence
% INPUTS:
    % decoded_thetaSeq - uses decoded theta sequences that have passed velocity and number of active units thresholds. 
    %                   If empty, loads the file in the folder
    % save_option - 1 to save structure
     

function thetaSequences_decodingError = theta_sequences_decoding_error_phase(decoded_thetaSeq,save_option)

if isempty(decoded_thetaSeq)
    load Theta\decoded_theta_sequences.mat
end

load extracted_directional_place_fields.mat
directional_place_fields_BAYESIAN = directional_place_fields;
load extracted_position.mat
parameters = list_of_parameters;

fields = fieldnames(decoded_thetaSeq);
fields = fields(cellfun(@(s) isempty(strmatch ('phase', s)), fieldnames(decoded_thetaSeq)));
for d = 1 : length(fields) % for each direction
    
    place_fields =  directional_place_fields_BAYESIAN(d).place_fields;
    
    for t = 1 : length(place_fields.track) % for each track
        
        c = 1;
        for s = 1 : size(decoded_thetaSeq.(sprintf('%s',fields{d}))(t).theta_sequences,2) %for each theta sequence
            
            % decoded positions in theta sequence
            decoded_position = decoded_thetaSeq.(sprintf('%s',fields{d}))(t).theta_sequences(s).decoded_position;
            decoded_time = decoded_thetaSeq.(sprintf('%s',fields{d}))(t).theta_sequences(s).timebins_centre;
            
             % Find actual position centered on cycle trough (estimated peak of place field)
             linear_position =  position.linear(t).linear(~isnan(position.linear(t).linear));
             [~,time_idx] = min(abs(position.linear(t).timestamps - decoded_thetaSeq.(sprintf('%s',fields{d}))(t).theta_sequences(s).theta_cycle_trough_time));
             real_position = linear_position(time_idx); % position at the theta trough
            
            % conversion from actual position in each time bin to the closest binned position value
            [~,idx] = min(abs(directional_place_fields(d).place_fields.track(t).x_bin_centres - real_position));
            converted_position = place_fields.track(t).x_bin_centres(idx);
                  
            % Finds the maximum decoded probability for each time bin
            all_maxProbability = [];
            decoded_positions =[];
            for jj = 1 : size(decoded_position,2)
                [maxProbability_bin,indx_maxProbability_bin] =  max(decoded_position(:,jj)); %find the value and index of maximum probability per time bin(each column)
                all_maxProbability = [all_maxProbability, maxProbability_bin]; % save all max probabilities
                decoded_positions = [decoded_positions, place_fields.track(t).x_bin_centres(indx_maxProbability_bin)]; %find the maximum probability decoded position per time bin
            end
            
            % Normalize decoded position to current position
            position_window =  decoded_positions - converted_position;
             
            % Save in structure
            bayesian_decodingError(t).thetaSequence(c).index_from_theta_window = decoded_thetaSeq.(sprintf('%s',fields{d}))(t).theta_sequences(s).index_from_theta_windows;
            bayesian_decodingError(t).thetaSequence(c).theta_window_edges = [decoded_thetaSeq.(sprintf('%s',fields{d}))(t).theta_sequences(s).theta_cycle_peaks_times(1) decoded_thetaSeq.(sprintf('%s',fields{d}))(t).theta_sequences(s).theta_cycle_peaks_times(2)];
            bayesian_decodingError(t).thetaSequence(c).theta_window_trough = decoded_thetaSeq.(sprintf('%s',fields{d}))(t).theta_sequences(s).theta_cycle_trough_time;
            bayesian_decodingError(t).thetaSequence(c).decoded_positions= decoded_positions;
            bayesian_decodingError(t).thetaSequence(c).real_position = real_position;
            bayesian_decodingError(t).thetaSequence(c).real_position_window = position_window;
            bayesian_decodingError(t).thetaSequence(c).decoded_time = decoded_time;

            % Calculate decoding error with and without threshold for low probabilities
            bayesian_decodingError(t).thetaSequence(c).decodingErrors = abs(bayesian_decodingError(t).thetaSequence(c).decoded_positions - bayesian_decodingError(t).thetaSequence(c).real_position);   % difference between actual and decoded position
            bayesian_decodingError(t).thetaSequence(c).weighted_decodingErrors = bayesian_decodingError(t).thetaSequence(c).decodingErrors.*all_maxProbability; % multiply each diff value by the max probability on that time bin
            
            bayesian_decodingError(t).thetaSequence(c).mean_DecodingError = mean(bayesian_decodingError(t).thetaSequence(c).decodingErrors); %median decoding error
            bayesian_decodingError(t).thetaSequence(c).median_DecodingError = median(bayesian_decodingError(t).thetaSequence(c).decodingErrors); %mean decoding error
            bayesian_decodingError(t).thetaSequence(c).mean_weighted_decodingError = mean(bayesian_decodingError(t).thetaSequence(c).weighted_decodingErrors)/sum(all_maxProbability); %mean weighted decoding error
            bayesian_decodingError(t).thetaSequence(c).median_weighted_decodingError = median(bayesian_decodingError(t).thetaSequence(c).weighted_decodingErrors)/sum(all_maxProbability); %median weighted decoding error
            
            bayesian_decodingError(t).all_decoded_errors(c,:) =  bayesian_decodingError(t).thetaSequence(c).decodingErrors;
            bayesian_decodingError(t).all_decoded_weighted_errors(c,:) =  bayesian_decodingError(t).thetaSequence(c).weighted_decodingErrors;

            c = c+1;

        end
    end
    thetaSeq(d).bayesian_decodingError = bayesian_decodingError;    
end

thetaSequences_decodingError = thetaSeq;

if save_option == 1
    save Theta\thetaSequences_decodingError thetaSequences_decodingError
end
end