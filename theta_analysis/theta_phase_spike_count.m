% SPIKE COUNT FOR BAYESIAN DECODING FOR THETA SEQUENCES
% Method to extract spike count using time bins based on theta phase.
% First,extracts theta phase during running and bins it. Then finds spike count for each phase bin. 
% Finally, converts phase bin edges back to time and calculates spike count used in bayesian decoding.
% INPUTS:
    % phase_shift. Default is 1 (no shift). Factor by which the phase bin edge is shifted based on the function 'check_phase_range_window.mat'
% OUTPUTS:
    % bayesian_spike_count structure. Contains time vector with time bin edges, time bin centres and matrix of spike count per place field
    % theta_bins_width_concat. All time bins widths concatenated in an array
% Loads: 'extracted_clusters', 'extracted_directional_place_fields','extracted_position','extracted_CSC.mat'

function [bayesian_spike_count] = theta_phase_spike_count(phase_shift)

load extracted_CSC.mat
load extracted_clusters.mat
load extracted_position.mat
parameters = list_of_parameters;

phase_bin_width = parameters.phase_bin_width;
if isempty(phase_shift)
    phase_shift = 1;
end

theta_CSC_time = CSC(4).CSCtime;
theta_CSC = CSC(4).theta;
clear CSC
tic
for track = 1 : length(position.linear)
    % Get start and end times per track
    start_times(track) = min(position.linear(track).timestamps);
    end_times(track) = max(position.linear(track).timestamps);
    track_indices = find(theta_CSC_time >= start_times(track) & theta_CSC_time <= end_times(track)); %time track indices
    
    % Get theta phase
    hilb = hilbert(theta_CSC(track_indices));%hilbert of CSC theta track
    theta_phase = angle(hilb); % 0 is equivalent to peak of theta, while the trough is approx PI
    temp_theta_unwrap = unwrap(theta_phase); % unwrap for interpolation

    % Check that theta phase has been extracted properly. If there was a high frequency event during theta, this might cause an error during
    % the extraction which will be detected by having some values with descending order.
    if issorted(temp_theta_unwrap) == 0
        
        % Find error indices
        backwards_indices = find(diff(temp_theta_unwrap)<0);
        clear temp_theta_unwrap
        % Find start and end indices of all error segments
        jump_indices = find(diff(backwards_indices)>2);
        start_idx = []; end_idx = [];
        end_idx = [end_idx; backwards_indices(jump_indices); backwards_indices(end)];
        start_idx = [start_idx; backwards_indices(1); backwards_indices(jump_indices+1)];
        % Find indx of the closest trough and peak surrounding the error segment
        [~,locstr] = findpeaks(-theta_phase);
        [~,locspeak] = findpeaks(theta_phase);
        all_peaks_troughs = [];
       
        %%% SAFETY CHECK FIGURE
%         timeaxes = 1: length(theta_phase);
%         f1 = figure;
%         ax1 = subplot(2,1,1);
%         plot(timeaxes, theta_phase,'b');
%         hold on
%         plot(timeaxes(backwards_indices),theta_phase(backwards_indices),'go')
%         plot(timeaxes(start_idx),theta_phase(start_idx),'o','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','r')
%         plot(timeaxes(end_idx),theta_phase(end_idx),'o','MarkerSize',8,'MarkerFaceColor','y','MarkerEdgeColor','y')

        
        for s = 1 : length(start_idx)
            trough_locs = locstr(locstr < start_idx(s)); %get trough indices previous to the current start idx
            [~,trix] = min(abs(start_idx(s)-trough_locs)); %indx of closest trough
            closest_trough = trough_locs(trix);
            peaks_locs = locspeak(locspeak > end_idx(s)); %get peak indices after to the current end idx
            [~,pkix] = min(abs(end_idx(s)-peaks_locs)); %indx of closest peak
            closest_peak = peaks_locs(pkix);
            clear trough_locs peaks_locs
            
            % Delete points in between trough and peak and interpolate. Do it in two steps, from trough to 0 (or closest to 0), and from
            % 0 to peak. Done because 0 corresponds at the peak phase in the real theta signal, and we don't want to mess with spikes between theta windows (peak to peak)
            if ~isempty(closest_peak) && ~isempty(closest_trough)
                segment = theta_phase(closest_trough:closest_peak);
                [~, closest_to_zero] = min(abs(segment));
                new_points_trough = linspace(theta_phase(closest_trough),segment(closest_to_zero),length(1:closest_to_zero));
                new_points_peak = linspace(segment(closest_to_zero),theta_phase(closest_peak),length(closest_to_zero:length(segment)));
                if issorted([new_points_trough(1:end-1) new_points_peak])
                    theta_phase(closest_trough:closest_peak) = [new_points_trough(1:end-1) new_points_peak];
                else
                    new_points = linspace(theta_phase(closest_trough),theta_phase(closest_peak),length(closest_trough:closest_peak)); %do not break in zero
                    theta_phase(closest_trough:closest_peak) = new_points;
                end
            elseif isempty(closest_trough)%if it's the first cycle
                segment = theta_phase(1:closest_peak);
                [~, closest_to_zero] = min(abs(segment));
                temp_points_trough = linspace(theta_phase(1),segment(closest_to_zero),length(1:closest_to_zero));
                new_points_trough = theta_phase(1): abs(mean(diff(temp_points_trough))):theta_phase(1)+(abs(mean(diff(temp_points_trough)))*(length(temp_points_trough)-1));
                temp_points_peak = linspace(segment(closest_to_zero),theta_phase(closest_peak),length(closest_to_zero:length(segment)));
                new_points_peak = segment(closest_to_zero): abs(mean(diff(temp_points_peak))):segment(closest_to_zero)+(abs(mean(diff(temp_points_peak)))*(length(temp_points_peak)-1));
                if issorted([new_points_trough(1:end-1) new_points_peak])
                    theta_phase (1:closest_peak) = [new_points_trough(1:end-1) new_points_peak];
                else
                    new_points = linspace(theta_phase(1),theta_phase(closest_peak),length(1:closest_peak)); %do not break in zero
                    theta_phase(1:closest_peak) = new_points;
                end
            else %if it's the last cycle
                segment = theta_phase(start_idx(s):end);
                [~, closest_to_zero] = min(abs(segment));
                temp_points_trough = linspace(theta_phase(start_idx(s)),segment(closest_to_zero),length(1:closest_to_zero));
                new_points_trough = theta_phase(start_idx(s)): abs(mean(diff(temp_points_trough))):theta_phase(start_idx(s))+(abs(mean(diff(temp_points_trough)))*(length(temp_points_trough)-1));
                if segment(closest_to_zero) ~= theta_phase(end)
                    temp_points_peak = linspace(segment(closest_to_zero),theta_phase(end),length(closest_to_zero:length(segment)));
                    new_points_peak = segment(closest_to_zero): abs(mean(diff(temp_points_peak))):segment(closest_to_zero)+(abs(mean(diff(temp_points_peak)))*(length(temp_points_peak)-1));
                    if issorted([new_points_trough(1:end-1) new_points_peak])
                        theta_phase(start_idx(s):end) = [new_points_trough(1:end-1) new_points_peak];
                    else
                        new_points = linspace(theta_phase(start_idx(s)),theta_phase(end),length(theta_phase(start_idx(s):end)));
                        theta_phase(start_idx(s):end) = new_points;
                    end
                else
                    theta_phase(start_idx(s):end) = new_points_trough;
                end
            end
        end
        
        % figure(f1)
        % plot(ax1,timeaxes(all_peaks_troughs),theta_phase(all_peaks_troughs),'o','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k');
        % plot(timeaxes(start_idx(s):end), new_points,'ko')
    end
    
    %%% EXTRA SAFETY CHECK FIGURE  %%%%
    %figure(f1); ax2= subplot(2,1,2); plot(theta_phase); linkaxes([ax1 ax2])
    
    % GET TIME BINS FROM SPIKE PHASES
    theta_phase_unwrap{track} = unwrap(theta_phase); % Get new unwrap phase for interpolation
    clear hilb theta_phase start_idx % Clear space
    track_spike_times{track} = clusters.spike_times(clusters.spike_times >= start_times(track) & clusters.spike_times <= end_times(track));    % Get spikes times per track
    spike_phases = interp1(theta_CSC_time(track_indices),theta_phase_unwrap{track},track_spike_times{track},'linear','extrap');  % Interpolate phases for each spike time (get phase for each spike)
    phase_bin_edges{track} = -phase_bin_width*phase_shift : phase_bin_width : spike_phases(end);% Get phase bin centers and edges
    %phase_bin_edges{track} = 0 : phase_bin_width : spike_phases(end);
    time_bin_edges{track} = interp1(theta_phase_unwrap{track},theta_CSC_time(track_indices),phase_bin_edges{track},'linear','extrap'); % Interpolate phase bins back to time (output time bins width will be different)
    time_bins_width{track} = diff(time_bin_edges{track});  % Find time bins width

   
    %%% EXTRA SAFETY CHECK FIGURE  %%%%
%     figure
%     plot(theta_CSC_time(track_indices),theta_phase_unwrap{track},'k*','MarkerSize',1)
%     hold on
%     plot(track_spike_times{track},spike_phases,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4)
%     plot(time_bin_edges{track},phase_bin_edges{track},'d','MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',10)
%     legend({'Theta phase per time','Spike phase','Time bin edge based on phase'})
%     xlabel('Time(ms)')
%     ylabel('Theta phase')   
%     figure; plot(theta_CSC_time(track_indices), theta_phase*700,'b');
%     hold on
%     plot(theta_CSC_time(track_indices),theta_CSC(track_indices),'k')
%     plot(time_bin_edges{track},zeros(1,length(time_bin_edges{track})),'d','MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',10)
%     plot(track_spike_times{track},spike_phases,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4)
%      
     clear spike_phases track_indices
end

clear position theta_CSC CSC_track_times theta_CSC_time
toc

%%%% SPIKE COUNT  %%%%
load extracted_directional_place_fields.mat

%select place cells to use for theta sequence analysis
place_field_index = unique([directional_place_fields(1).place_fields.good_place_cells directional_place_fields(2).place_fields.good_place_cells]);
clear directional_place_fields

bayesian_spike_count.n.replay= [];    bayesian_spike_count.replay_events_indices =[];   replay_edges_concat=[];  
theta_bins_width_concat = [];
for i = 1 : length(time_bin_edges)
    % mantain same variable names than in spike_count
    
    t.replay_edges = cell2mat(time_bin_edges(i)); % takes each replay time vector separately
    % Takes time vectors and centres each bin
    bayesian_spike_count.replay_events(i).replay_time_edges = t.replay_edges;
    bayesian_spike_count.replay_events(i).replay_time_centered = t.replay_edges(1:end-1)+ time_bins_width{i}/2; %centres of bins
    replay_edges_concat = [replay_edges_concat t.replay_edges];
    theta_bins_width_concat = [theta_bins_width_concat time_bins_width{i},NaN]; 
    %concatenate indicies for each replay event, put NaN in between events so that spikes between edges of replay events are ignored
    bayesian_spike_count.replay_events_indices = [bayesian_spike_count.replay_events_indices, ones(1,length(bayesian_spike_count.replay_events(i).replay_time_centered))*i,NaN]; 
end

% Spike histogram per time bin for each place field
% performed on bins across all replay events (time bin between replay events gets ignored later because replay_event_indice is set to NaN
for k = 1 :length(place_field_index)
    bayesian_spike_count.n.replay(k,:) = histcounts(clusters.spike_times(find(clusters.spike_id==place_field_index(k))),replay_edges_concat);
end
bayesian_spike_count.theta_bins_width_concat = theta_bins_width_concat;

% FINAL SAFETY CHECK
 % figure
% plot(CSC(4).CSCtime,CSC(4).theta)
% hold on
% for track = 1 : 4
%     plot(theta_CSC_time(track_indices),theta_CSC(track_indices))
%     plot(bayesian_spike_count.replay_events(track).replay_time_edges,zeros(1,length(bayesian_spike_count.replay_events(track).replay_time_edges)),'o')
% end
end