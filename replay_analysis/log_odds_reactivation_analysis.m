function [place_fields_BAYESIAN decoded_replay_events,probability_ratio_original,decoded_replay_events_shuffled] = log_odds_reactivation_analysis(clusters,place_fields_BAYESIAN,replay,position,modify,stimulus_name,timebin)
% Input:
% Clusters
% place_fields_BAYESIAN
% replay (e.g. replay = reactivations.probe(nprobe))

%%% Decoded Track Replay Structure
% Extracts spikes for each replay event and decodes it. Saves in a
% structure called repay_track. Inside each field (each track) the replay events are saved.
% Loads: extracted_replay_events, extracted_clusters and extracted_place_fields_BAYESIAN

% REPLAY EVENTS STRUCTURE
%replay_events is an empty template for replay event analysis.
%Each track will create its own field
place_cell_index = place_fields_BAYESIAN(1).all_good_place_cells_LIBERAL;

replay_events = struct('replay_id',{},...%the id of the candidate replay events in chronological order
    'spikes',{});%column 1 is spike id, column 2 is spike time

% TAKE SPIKES FROM ONLY common good cells on both tracks

sorted_spikes = zeros(size(clusters.spike_id));
sorted_spikes(:,1) = clusters.spike_id;
sorted_spikes(:,2) = clusters.spike_times;
all_units = unique(clusters.spike_id);


non_good = setdiff(all_units,place_fields_BAYESIAN(1).cluster_id(place_cell_index));% Find cells that are not good place cells

for i = 1 : length(non_good)
    non_good_indices = find(sorted_spikes(:,1)== non_good(i));
    sorted_spikes(non_good_indices,:) = [];% remove spikes from unwanted cells
    non_good_indices =[];
end

num_spikes = length(sorted_spikes);
num_units = length(place_cell_index);

% EXTRACT SPIKES IN REPLAY EVENTS
num_replay = size(replay.onset, 2);
current_replay = 1;
current_replay_spikes = [];

for i = 1 : num_spikes
    % Collect spike data during replay
    if sorted_spikes(i,2) > replay.offset(current_replay)
        replay_events(current_replay).replay_id = current_replay;
        replay_events(current_replay).spikes = current_replay_spikes;
        current_replay = current_replay + 1;
        if current_replay > num_replay
            break
        end
        current_replay_spikes = [];
    end

    if sorted_spikes(i,2) >= replay.onset(current_replay)
        % If spike happens during replay, records it as replay spike
        current_replay_spikes = [current_replay_spikes; sorted_spikes(i,:)];
    end

    if i == num_spikes && current_replay == num_replay && isempty(current_replay_spikes)
        % for the last replay event, if it is empty, write empty
        replay_events(current_replay).replay_id = current_replay;
        replay_events(current_replay).spikes = current_replay_spikes;
    end

end

num_replay_events = length(replay_events);

if length(replay_events) > 1
    msg = [num2str(num_replay_events), ' candidate events.'];
    disp(msg);
end

% Save all replay events all tracks
for j = 1:length(place_fields_BAYESIAN)
    decoded_replay_events(j).replay_events = replay_events;
end

%%%%%% Place field modification if needed %%%%%%%%%%%
if strcmp(modify,'global_remapped')
    global_remapped_place_fields_id = [];
%     bayesian_spike_count = 'replayEvents_bayesian_spike_count';

    for event = 1:length(decoded_replay_events(1).replay_events)
        % global remapping by swapping the Cell ID for good cells
        global_remapped_place_fields = [];


        for track_id = 1:2
            s = RandStream('mrg32k3a','Seed',track_id+100*event); % Set random seed for resampling
            global_remapped_place_fields{track_id} = place_fields_BAYESIAN(track_id).raw;

            random_cell_index = randperm(s,length(place_fields_BAYESIAN(track_id).good_place_cells_LIBERAL));
            random_cell = place_fields_BAYESIAN(track_id).good_place_cells_LIBERAL(random_cell_index);
            %                         place_fields_BAYESIAN.track(track_id).random_cell = random_cell;
            original_cell = place_fields_BAYESIAN(track_id).good_place_cells_LIBERAL;

            for j=1:length(random_cell) %only swap good cells
                global_remapped_place_fields{track_id}{original_cell(j)}=place_fields_BAYESIAN(track_id).raw{random_cell(j)};
            end

            place_fields_BAYESIAN(track_id).global_remapped{event} = global_remapped_place_fields{track_id};

            % Also save the shuffled place cell id for subsequent spearsman
            % analysis
            global_remapped_place_fields_id{event}{track_id}(1,:) = original_cell;
            global_remapped_place_fields_id{event}{track_id}(2,:) = random_cell;
        end

    end

%     if contains(stimulus_name,'Masa2tracks')
%         save(sprintf('global_remapped_place_fields_id%s.mat',erase(stimulus_name,'Masa2tracks')),'global_remapped_place_fields_id');
%     else
%         save('global_remapped_place_fields_id.mat','global_remapped_place_fields_id');
%     end
end


%%%%%% BAYESIAN DECODING ON REPLAY EVENTS %%%%%%
replay_starts = replay.onset;
replay_ends = replay.offset;

% Get time vectors for bayesian decoding and matrix with spike count
% disp('Spike count...');
replayEvents_bayesian_spike_count = create_spike_count_masa(place_fields_BAYESIAN,clusters,...
    replay_starts,replay_ends,[],timebin);

%     save replayEvents_bayesian_spike_count replayEvents_bayesian_spike_count

%     replay_decoding_split_events

% Run bayesian decoding
fprintf('Decoding %s...\n',modify)

bayesian_spike_count = replayEvents_bayesian_spike_count;
if contains(stimulus_name,'RUN')|contains(stimulus_name,'Track')
    estimated_position = bayesian_decoding(place_fields_BAYESIAN,replayEvents_bayesian_spike_count,position,[],modify,[],timebin);
else
    estimated_position = bayesian_decoding_POST(place_fields_BAYESIAN,replayEvents_bayesian_spike_count,position,[],modify,[],timebin);
%      estimated_position = bayesian_decoding(place_fields_BAYESIAN,replayEvents_bayesian_spike_count,position,[],modify,[],timebin);
end
% Save in structure

for i = 1 : num_replay_events
    for j = 1:length(place_fields_BAYESIAN)
        decoded_replay_events(j).replay_events(i).replay = estimated_position(j).replay_events(i).replay;
%         decoded_replay_events(j).replay_bias = estimated_position(j).replay_bias;
        decoded_replay_events(j).replay_events(i).timebins_edges = estimated_position(j).replay_events(i).replay_time_edges;
        decoded_replay_events(j).replay_events(i).timebins_centre = estimated_position(j).replay_events(i).replay_time_centered;
        decoded_replay_events(j).replay_events(i).timebins_index = 1:length(estimated_position(j).replay_events(i).replay_time_centered);
        decoded_replay_events(j).replay_events(i).decoded_position = estimated_position(j).replay_events(i).replay; % normalized by all tracks
        if contains(stimulus_name,'RUN') % Only if it is awake replay during running
            decoded_replay_events(j).replay_events(i).replay_actual_position = estimated_position(j).replay_events(i).replay_actual_position;
        end
        decoded_replay_events(j).replay_events(i).probability_ratio  = estimated_position(j).replay_events(i).probability_ratio;
    end

end

for event = 1:length(estimated_position(1).replay_events)
    probability_ratio_original{1}(1,event) = estimated_position(1).replay_events(event).probability_ratio;
    probability_ratio_original{1}(2,event) = estimated_position(2).replay_events(event).probability_ratio;
end

fprintf('ratemap shuffle %s...\n',modify)

parfor nshuffle = 1:1000
    estimated_position_ratemap_shuffled = [];
    if contains(stimulus_name,'RUN')
        estimated_position_ratemap_shuffled = bayesian_decoding(place_fields_BAYESIAN,replayEvents_bayesian_spike_count,position,'ratemap shuffle',modify,nshuffle,timebin);
    else
        estimated_position_ratemap_shuffled = bayesian_decoding_POST(place_fields_BAYESIAN,replayEvents_bayesian_spike_count,position,'ratemap shuffle',modify,nshuffle,timebin);
    end
    %         estimated_position_ratemap_shuffled = log_odds_bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,place_cell_index,timebin,[],'ratemap shuffle','','N');
    for event = 1:length(estimated_position_ratemap_shuffled(1).replay_events)
        ratemap_shuffled_probability_ratio{nshuffle}(1,event,1) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio;
        ratemap_shuffled_probability_ratio{nshuffle}(2,event,1) = (estimated_position_ratemap_shuffled(2).replay_events(event).probability_ratio);
        %             ratemap_shuffled_probability_ratio{nshuffle}(1,event,2) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_first;
        %             ratemap_shuffled_probability_ratio{nshuffle}(2,event,2) = 1/(estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_first);
        %             ratemap_shuffled_probability_ratio{nshuffle}(1,event,3) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_second;
        %             ratemap_shuffled_probability_ratio{nshuffle}(2,event,3) = 1/(estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_second);
    end
    decoded_replay_events_shuffled{nshuffle} = estimated_position_ratemap_shuffled;
end

probability_ratio_original{2} = ratemap_shuffled_probability_ratio;

end