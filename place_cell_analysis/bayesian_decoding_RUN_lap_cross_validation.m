function [probability_ratio_RUN_lap estimated_position_lap_CV] = bayesian_decoding_RUN_lap_cross_validation(clusters,place_fields_all,position,lap_times)

% Bayesian decoding 10 fold cross validation
% output - probability_ratio_RUN_lap{1} is original and
% probability_ratio_RUN_lap{2} is shuffled
% probability_ratio_RUN_lap{1}{physical track}{decoded track}{lap id}

time_bin = 0.1;
cv_groups = [];
for track_id = 1:length(place_fields_all.track)

    % create 10 groups with random lap id for 10 fold cross
    % validation bayesian decoding
    s = RandStream('mrg32k3a','Seed',track_id); % Set random seed for resampling
    random_lap_index = randperm(s,length(lap_times(track_id).start));

    number_per_group = floor(length(random_lap_index)/10);

    startIndex = 1;
    for groupIndex = 1:10
        endIndex = startIndex + number_per_group - 1;
        cv_groups{track_id}{groupIndex} = random_lap_index(startIndex:endIndex);
        startIndex = endIndex + 1;
    end

    if startIndex <= length(random_lap_index)
        for groupIndex = 1:10
            if startIndex > length(random_lap_index)
                break
            end
            cv_groups{track_id}{groupIndex} = [cv_groups{track_id}{groupIndex} random_lap_index(startIndex)];
            startIndex = startIndex + 1;
        end
    end
end


clear estimated_position_lap_CV
for track_id = 1:length(place_fields_all.track)
    estimated_position_lap_CV(track_id).lap(length(lap_times(track_id).start)).track(1:2) = struct();
end
%         estimated_position_lap.track(1:2) = struct('track_template',struct([]));
%         estimated_position_lap.track(2) = struct('track_template',struct([]));
%         estimated_position_lap = struct('track',[]);
%         estimated_position_lap.track(2) = [];
probability_ratio_RUN_lap = [];

% bayesian decoding
for groupIndex = 1:10
    tic
    place_fields_option = [];
    % use 9/10 of data for place field construction
    place_fields_option.track(1).lap_id = setxor(cv_groups{1}{groupIndex},1:length(lap_times(1).start));
    place_fields_option.track(2).lap_id = setxor(cv_groups{2}{groupIndex},1:length(lap_times(2).start));
    place_fields_training = calculate_place_fields_masa_NPX(10,position,clusters,place_fields_option);
    place_fields_training.good_place_cells_LIBERAL = place_fields_all.good_place_cells_LIBERAL;
    place_fields_training.good_place_cells = place_fields_all.good_place_cells_LIBERAL;

    % remaining laps to test
    for track_id = 1:length(place_fields_all.track)
        for nlap = 1:length(cv_groups{track_id}{groupIndex})
            temp_estimated_position = [];

            if lap_times(track_id).end(cv_groups{track_id}{groupIndex}(nlap)) - lap_times(track_id).start(cv_groups{track_id}{groupIndex}(nlap)) < 2
                start_time = lap_times(track_id).start(cv_groups{track_id}{groupIndex}(nlap));
                end_time = lap_times(track_id).end(cv_groups{track_id}{groupIndex}(nlap)) + 2;
            else
                start_time = lap_times(track_id).start(cv_groups{track_id}{groupIndex}(nlap));
                end_time = lap_times(track_id).end(cv_groups{track_id}{groupIndex}(nlap));
            end
            bayesian_spike_count_RUN = create_spike_count_masa(place_fields_all,clusters,...
                start_time,end_time,[]);

            temp_estimated_position = bayesian_decoding(place_fields_training,bayesian_spike_count_RUN,position,[],[],[],time_bin);
            % first track is track to decode and second track is place field template track
            estimated_position_lap_CV(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track = temp_estimated_position;

            probability_ratio_RUN_lap{1}{track_id}{1}{cv_groups{track_id}{groupIndex}(nlap)} = temp_estimated_position(1).probability_ratio;
            probability_ratio_RUN_lap{1}{track_id}{2}{cv_groups{track_id}{groupIndex}(nlap)} = temp_estimated_position(2).probability_ratio;

            disp('Ratemap shuffle (original) shuffle')
            parfor nshuffle = 1:1000
                estimated_position_ratemap_shuffled = [];
                estimated_position_ratemap_shuffled = bayesian_decoding(place_fields_training,bayesian_spike_count_RUN,position,'ratemap shuffle',[],nshuffle,time_bin);
                %         estimated_position_ratemap_shuffled = log_odds_bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,place_cell_index,timebin,[],'ratemap shuffle','','N');

                ratemap_shuffled_probability_ratio{nshuffle}{track_id}{1}{cv_groups{track_id}{groupIndex}(nlap)} = estimated_position_ratemap_shuffled(1).probability_ratio;
                ratemap_shuffled_probability_ratio{nshuffle}{track_id}{2}{cv_groups{track_id}{groupIndex}(nlap)} = estimated_position_ratemap_shuffled(2).probability_ratio;
                %             ratemap_shuffled_probability_ratio{nshuffle}(1,event,2) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_first;
                %             ratemap_shuffled_probability_ratio{nshuffle}(2,event,2) = 1/(estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_first);
                %             ratemap_shuffled_probability_ratio{nshuffle}(1,event,3) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_second;
                %             ratemap_shuffled_probability_ratio{nshuffle}(2,event,3) = 1/(estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_second);
            end                
        end


    end
    toc
end

probability_ratio_RUN_lap{2} = ratemap_shuffled_probability_ratio;

temp_estimated_position = [];

