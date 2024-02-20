function [probability_ratio_RUN_lap estimated_position_lap_CV] = bayesian_decoding_RUN_lap_cross_validation(clusters,place_fields_all,Behaviour,Task_info)

% Bayesian decoding 10 fold cross validation
% output - probability_ratio_RUN_lap{1} is original and
% probability_ratio_RUN_lap{2} is shuffled
% probability_ratio_RUN_lap{1}{physical track}{decoded track}{lap id}
lap_times = [];
position = [];

time_bin = 0.1;
cv_groups = [];
for track_id = 1:length(place_fields_all)

    % create 10 groups with random lap id for 10 fold cross
    % validation bayesian decoding
    s = RandStream('mrg32k3a','Seed',track_id); % Set random seed for resampling
    random_lap_index = randperm(s,sum(Task_info.track_ID_all == track_id));

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

    for groupIndex = 1:10
        cv_groups{track_id}{groupIndex} = sort(cv_groups{track_id}{groupIndex});
    end
end


clear estimated_position_lap_CV
for track_id = 1:length(place_fields_all)
    % Actual track -> laps -> track template
    estimated_position_lap_CV(track_id).lap(sum(Task_info.track_ID_all == track_id)).track(1:2) = struct();
end
%         estimated_position_lap.track(1:2) = struct('track_template',struct([]));
%         estimated_position_lap.track(2) = struct('track_template',struct([]));
%         estimated_position_lap = struct('track',[]);
%         estimated_position_lap.track(2) = [];
probability_ratio_RUN_lap = [];
place_fields_BAYESIAN = calculate_spatial_cells(clusters,Task_info,Behaviour,[0 140],10); % use 10 cm bin for Bayesian decoding

good_cell_index = unique([find(place_fields_all(1).peak_percentile>0.99 & place_fields_all(1).odd_even_stability>0.99)...
    find(place_fields_all(2).peak_percentile>0.99 & place_fields_all(2).odd_even_stability>0.99)]);

% Laps
track_laps{1} = find(Task_info.track_ID_all == 1);
track_laps{2} = find(Task_info.track_ID_all == 2);

% bayesian decoding
for groupIndex = 1:10
    tic
    
    ratemap_matrix = [];
    training_laps = [];

    % use 9/10 of data for place field construction
    training_laps{1} = setxor(cv_groups{1}{groupIndex},1:sum(Task_info.track_ID_all == 1));
    training_laps{2} = setxor(cv_groups{2}{groupIndex},1:sum(Task_info.track_ID_all == 2));

    for track_id = 1:length(place_fields_all)
        ratemap_matrix = [place_fields_BAYESIAN(track_id).raw{good_cell_index}];

        ratemap_matrix = reshape(ratemap_matrix,size(place_fields_BAYESIAN(track_id).raw{1},1),[],length(good_cell_index));%laps X position bins X cells
        template_maps = normalize(squeeze(mean(ratemap_matrix(training_laps{track_id},:,:),1)),'range');
        template_maps(isnan(template_maps))=0;

        place_fields_BAYESIAN(track_id).template = template_maps';% ncell X position bin (for decoding)
        place_fields_BAYESIAN(track_id).good_cells = good_cell_index;
        %     place_fields_training = calculate_place_fields_masa_NPX(10,position,clusters,place_fields_option);
    end

   
    for track_id = 1:length(place_fields_all)
        % remaining laps to test
        start_time = [];
        end_time = [];

        for nlap = 1:length(cv_groups{track_id}{groupIndex})
            temp_estimated_position = [];

            if Task_info.end_time_all(track_laps{track_id}(cv_groups{track_id}{groupIndex}(nlap))) ...
                    - Task_info.start_time_all(track_laps{track_id}(cv_groups{track_id}{groupIndex}(nlap))) < 2

                start_time(nlap) = Task_info.start_time_all(track_laps{track_id}(cv_groups{track_id}{groupIndex}(nlap)));
                end_time(nlap) = Task_info.end_time_all(track_laps{track_id}(cv_groups{track_id}{groupIndex}(nlap))) + 2;
            else
                start_time(nlap) = Task_info.start_time_all(track_laps{track_id}(cv_groups{track_id}{groupIndex}(nlap)));
                end_time(nlap) = Task_info.end_time_all(track_laps{track_id}(cv_groups{track_id}{groupIndex}(nlap)));
            end
        end

        bayesian_spike_count_RUN = create_spike_count_masa(place_fields_BAYESIAN,clusters,...
            start_time,end_time,[]);

        for nlap = 1:length(cv_groups{track_id}{groupIndex})
            temp_estimated_position = bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count_RUN,Behaviour,[],[],[],time_bin);
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

