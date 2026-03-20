% THETA SEQUENCES DETECTION AND DECODING
% MH 2020
% Bins theta based on spike phases. Converts phase bins to time, and decodes entire tracks using the new bins. Then divides in theta windows
% (peak to peak), and applies thresholds to discard the noisy cycles. Finally, split the theta windows based on direction
% Loads: extracted_directional_place_fields_BAYESIAN 

function decoded_thetaSeq = theta_sequences_detection_decoding(phase_shift)

load extracted_directional_place_fields.mat

% Extracts spike phase to create bins. Then counts spike in each bin
disp('Spike count...');
[theta_spike_count] = theta_phase_spike_count(phase_shift); % using phase bins
save Theta\theta_seq_bayesian_spike_count theta_spike_count


% Run bayesian decoding with place fields from each direction
disp('Decoding position...');
estimated_position.dir1 = bayesian_decoding(directional_place_fields,theta_spike_count,1,'N');
estimated_position.dir2 = bayesian_decoding(directional_place_fields,theta_spike_count,2,'N');

save('Theta\estimated_position_theta_sequences','estimated_position','-v7.3')

% Apply thresholds to theta windows (peak to peak)
theta_windows = extract_theta_window(theta_spike_count,'Y');

% Split the decoded tracks in theta cycles, and classify according to direction
for track = 1 : length(theta_windows.track)
    c_1 = 1;
    c_2 = 1;
    for s = 1 : length(theta_windows.track(track).theta_windows) %for each theta window
        if theta_windows.track(track).theta_windows(s,3) ~= 0
            % Find indices of time bins edges and centres within each theta window (should match)
            idcs = find(estimated_position.dir1(track).replay_events(track).replay_time_edges >= theta_windows.track(track).theta_windows(s,1)  &...
                estimated_position.dir1(track).replay_events(track).replay_time_edges  <= theta_windows.track(track).theta_windows(s,2));
            idcs_cent = find(estimated_position.dir1(track).replay_events(track).replay_time_centered > theta_windows.track(track).theta_windows(s,1)  &...
                estimated_position.dir1(track).replay_events(track).replay_time_centered  < theta_windows.track(track).theta_windows(s,2));
            % If theta window is in direction 1 
            if theta_windows.track(track).theta_windows(s,3) == 1
                decoded_thetaSeq.direction1(track).theta_sequences(c_1).theta_cycle_peaks_times = [theta_windows.track(track).theta_windows(s,1) theta_windows.track(track).theta_windows(s,2)];
                decoded_thetaSeq.direction1(track).theta_sequences(c_1).theta_cycle_trough_time = theta_windows.track(track).theta_troughs_timestamps(s,1);
                decoded_thetaSeq.direction1(track).theta_sequences(c_1).timebins_edges = estimated_position.dir1(track).replay_events(track).replay_time_edges(idcs);
                decoded_thetaSeq.direction1(track).theta_sequences(c_1).timebins_centre = estimated_position.dir1(track).replay_events(track).replay_time_centered(idcs_cent);
                decoded_thetaSeq.direction1(track).theta_sequences(c_1).timebins_edge_index = idcs;
                decoded_thetaSeq.direction1(track).theta_sequences(c_1).timebins_centre_index = idcs_cent;
                decoded_thetaSeq.direction1(track).theta_sequences(c_1).decoded_position = estimated_position.dir1(track).replay_events(track).replay_OneTrack(:,idcs_cent); % normalized by all tracks
                decoded_thetaSeq.direction1(track).theta_sequences(c_1).raw_decoded_position = estimated_position.dir1(track).replay_events(track).replay_raw(:,idcs_cent); 
                decoded_thetaSeq.direction1(track).theta_sequences(c_1).index_from_theta_windows = s; %index of theta event in theta windows structure
                decoded_thetaSeq.direction1(track).track_active_units_ID = estimated_position.dir1(track).active_units_ID;
                c_1 = c_1+1;
            % If theta window is in direction 2 
            elseif theta_windows.track(track).theta_windows(s,3) == 2
                decoded_thetaSeq.direction2(track).theta_sequences(c_2).theta_cycle_peaks_times = [theta_windows.track(track).theta_windows(s,1) theta_windows.track(track).theta_windows(s,2)];
                decoded_thetaSeq.direction2(track).theta_sequences(c_2).theta_cycle_trough_time = theta_windows.track(track).theta_troughs_timestamps(s,1);
                decoded_thetaSeq.direction2(track).theta_sequences(c_2).timebins_edges = estimated_position.dir2(track).replay_events(track).replay_time_edges(idcs);
                decoded_thetaSeq.direction2(track).theta_sequences(c_2).timebins_centre = estimated_position.dir2(track).replay_events(track).replay_time_centered(idcs_cent);
                decoded_thetaSeq.direction2(track).theta_sequences(c_2).timebins_edge_index = idcs;
                decoded_thetaSeq.direction2(track).theta_sequences(c_2).timebins_centre_index = idcs_cent;
                decoded_thetaSeq.direction2(track).theta_sequences(c_2).decoded_position = estimated_position.dir2(track).replay_events(track).replay_OneTrack(:,idcs_cent); % normalized by all tracks
                decoded_thetaSeq.direction2(track).theta_sequences(c_2).raw_decoded_position = estimated_position.dir2(track).replay_events(track).replay_raw(:,idcs_cent); 
                decoded_thetaSeq.direction2(track).theta_sequences(c_2).index_from_theta_windows = s; %index of theta event in theta windows structure
                decoded_thetaSeq.direction2(track).track_active_units_ID = estimated_position.dir2(track).active_units_ID;
                c_2 = c_2 +1;
            end
        end
    end
end


% Apply last threshold: position threshold
[decoded_thetaSeq,~]= theta_sequences_position_threshold(decoded_thetaSeq);

% Saves structure
decoded_thetaSeq.phase_shift = phase_shift;
save('Theta\decoded_theta_sequences','decoded_thetaSeq','-v7.3')

%%%% Figure for Sanity check
% Plot real times of trough and theta window edges, and then plot the time bins used for decoding
% Blue circles: Theta window edges; Green circle: central trough; Magenta
% stars: time bins centers; Black crosses: time bin edges; 
s = 3;
load('Theta\theta_peak_trough.mat')
figure;
plot(decoded_thetaSeq.direction1(1).theta_sequences(s).theta_cycle_peaks_times,ones(1,2),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',6)
hold on
plot(decoded_thetaSeq.direction1(1).theta_sequences(s).theta_cycle_trough_time,1,'o','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',6)
plot(decoded_thetaSeq.direction1(1).theta_sequences(s).timebins_edges,ones(1,length(decoded_thetaSeq.direction1(1).theta_sequences(s).timebins_edges)),'xk','MarkerSize',8)
plot(decoded_thetaSeq.direction1(1).theta_sequences(s).timebins_centre,ones(1,length(decoded_thetaSeq.direction1(1).theta_sequences(s).timebins_centre)),'*m','MarkerSize',6)

[~,idx] = min(abs(theta_troughs(:,4) - decoded_thetaSeq.direction1(1).theta_sequences(s).theta_cycle_trough_time));
upper = theta_peaks(theta_peaks(:,4) > theta_troughs(idx,4),4); 
lower = theta_peaks(theta_peaks(:,4) < theta_troughs(idx,4),4); 
surround_peaks = [lower(end) upper(1)];
plot([surround_peaks(1) surround_peaks(2)],ones(1,2),'dr','MarkerFaceColor','r','MarkerSize',6)

legend('Theta Time Window Edges - CSC Time','Theta Trough Time - CSC Time','Decoded Time Bin Centres - Position time','Decoded Time Bin Edges - Position time',...
    'Surrounding peaks - CSC Time');

end




function  [decoded_thetaSeq,deleted_decoded_thetaSeq]= theta_sequences_position_threshold(decoded_thetaSeq)
% Set position threshold: theta sequences only analysed when the animal is in the middle of track (1/6 to 5/6) - Feng et al (2015, J Neuro)

if isempty(decoded_thetaSeq)
    load Theta\decoded_theta_sequences.mat
end
load extracted_position.mat

parameters = list_of_parameters;

% Exclude reward sites, both are 20 cm (total length 200cm)
lower_thresh = 20; 
high_thresh = 180;

for d = 1 : length(fieldnames(decoded_thetaSeq)) % for each direction
    
    for t = 1 : length(decoded_thetaSeq.(strcat('direction',num2str(d)))) % for each track
        indices_to_remove = [];
        
        for s = 1 : size(decoded_thetaSeq.(strcat('direction',num2str(d)))(t).theta_sequences,2) %for each theta sequence
            
            linear_position =  position.linear(t).linear(~isnan(position.linear(t).linear));
            [~,time_idx] = min(abs(position.linear(t).timestamps - decoded_thetaSeq.(strcat('direction',num2str(d)))(t).theta_sequences(s).theta_cycle_trough_time));
            real_position(s) = linear_position(time_idx);
            
            % If animal position is at the start or end of track, save index
            if real_position(s) < lower_thresh | real_position(s) > high_thresh
                indices_to_remove = [indices_to_remove s];
            end
        end
        
        %%% Uncomment for safety check figure
        %figure;
        %histogram(real_position)
        %xlabel('Position (cm)'); ylabel('number of theta windows')
        
        clear real_position
        % For checking, save only the deleted sequences
        deleted_decoded_thetaSeq.(strcat('direction',num2str(d)))(t).theta_sequences = decoded_thetaSeq.(strcat('direction',num2str(d)))(t).theta_sequences(indices_to_remove);
        % Delete indices
        decoded_thetaSeq.(strcat('direction',num2str(d)))(t).theta_sequences(indices_to_remove) = [];
        
    end
    
end

end