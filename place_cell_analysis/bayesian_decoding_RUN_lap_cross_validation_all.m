function [probability_ratio_RUN_lap estimated_position_lap_CV estimated_position_lap_CV_shuffled] = bayesian_decoding_RUN_lap_cross_validation_all(clusters,place_fields_all,Behaviour,Task_info,options)

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
place_fields_BAYESIAN = calculate_spatial_cells(clusters,Task_info,Behaviour,[0 140],10); % use 10 cm bin for Bayesian decoding


spatial_cell_index = unique([find(place_fields_all(1).peak_percentile>0.95 & place_fields_all(1).odd_even_stability>0.95)...
    find(place_fields_all(2).peak_percentile>0.95 & place_fields_all(2).odd_even_stability>0.95)]);
% spatial_cell_index = unique([find(place_fields_all(1).odd_even_stability>0.95)...
%     find(place_fields_all(2).odd_even_stability>0.95)]);
good_cell_index = intersect(spatial_cell_index,find(ismember(place_fields_all(1).cluster_id, clusters.cluster_id)));
good_cell_index = find(ismember(place_fields_BAYESIAN(1).cluster_id,place_fields_all(1).cluster_id(good_cell_index)));
speed_range = 5:5:40;

if isempty(good_cell_index)

    disp('No spatial cells')
else
    % good_cell_index = unique([find(place_fields_all(1).peak_percentile>=0 )...
    %     find(place_fields_all(2).peak_percentile>=0)]);
    % Laps
    track_laps{1} = find(Task_info.track_ID_all == 1);
    track_laps{2} = find(Task_info.track_ID_all == 2);

    % Position/Speed shuffle
    for track_id = 1:length(place_fields_all)
        ratemap_matrix = [place_fields_BAYESIAN(track_id).raw{good_cell_index}];

        ratemap_matrix = reshape(ratemap_matrix,size(place_fields_BAYESIAN(track_id).raw{1},1),[],length(good_cell_index));%laps X position bins X cells
        %         template_maps = normalize(squeeze(mean(ratemap_matrix,1)),'range');
        template_maps = squeeze(mean(ratemap_matrix(:,:,:),1));
        template_maps(isnan(template_maps))=0;

        place_fields_BAYESIAN(track_id).template = template_maps';% ncell X position bin (for decoding)
        place_fields_BAYESIAN(track_id).good_cells = good_cell_index;
        %     place_fields_training = calculate_place_fields_masa_NPX(10,position,clusters,place_fields_option);
    end

    place_fields_BAYESIAN(track_id).good_cells = good_cell_index;
    bayesian_spike_count_RUN_WHOLE{1} = create_spike_count_masa(place_fields_BAYESIAN,clusters,...
        Task_info.start_time_all(1 == Task_info.track_ID_all)',Task_info.end_time_all(1 == Task_info.track_ID_all),[]);
    estimated_position_WHOLE{1} = bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count_RUN_WHOLE{1},Behaviour,[],[],[],time_bin);

    bayesian_spike_count_RUN_WHOLE{2} = create_spike_count_masa(place_fields_BAYESIAN,clusters,...
        Task_info.start_time_all(2 == Task_info.track_ID_all)',Task_info.end_time_all(2 == Task_info.track_ID_all),[]);
    estimated_position_WHOLE{2} = bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count_RUN_WHOLE{2},Behaviour,[],[],[],time_bin);

    clear estimated_position_lap_CV
    % bayesian decoding
    for groupIndex = 1:10
        tic

        ratemap_matrix = [];
        training_laps = [];

        % use 9/10 of data for place field construction
        training_laps{1} = setxor(cv_groups{1}{groupIndex},1:sum(Task_info.track_ID_all == 1));
        training_laps{2} = setxor(cv_groups{2}{groupIndex},1:sum(Task_info.track_ID_all == 2));

        %     training_laps{1} = cv_groups{1}{groupIndex};
        %     training_laps{2} = cv_groups{2}{groupIndex};

        for track_id = 1:length(place_fields_all)
            ratemap_matrix = [place_fields_BAYESIAN(track_id).raw{good_cell_index}];

            ratemap_matrix = reshape(ratemap_matrix,size(place_fields_BAYESIAN(track_id).raw{1},1),[],length(good_cell_index));%laps X position bins X cells
            %         template_maps = normalize(squeeze(mean(ratemap_matrix,1)),'range');
            template_maps = squeeze(mean(ratemap_matrix(training_laps{track_id},:,:),1));
            template_maps(isnan(template_maps))=0;

            place_fields_BAYESIAN(track_id).template = template_maps';% ncell X position bin (for decoding)
            place_fields_BAYESIAN(track_id).good_cells = good_cell_index;
            %     place_fields_training = calculate_place_fields_masa_NPX(10,position,clusters,place_fields_option);
        end


        for track_id = 1:length(place_fields_all)
            % remaining laps to test

            bayesian_spike_count_RUN = create_spike_count_masa(place_fields_BAYESIAN,clusters,...
                Task_info.start_time_all(track_laps{track_id}(cv_groups{track_id}{groupIndex}))',Task_info.end_time_all(track_laps{track_id}(cv_groups{track_id}{groupIndex})),[]);

            temp_estimated_position = bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count_RUN,Behaviour,[],[],[],time_bin);
            % first track is track to decode and second track is place field template track
            %             estimated_position_lap_CV(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track = temp_estimated_position;

            disp('Ratemap shuffle (original) shuffle')
            parfor nshuffle = 1:1000
                estimated_position_ratemap_shuffled = [];

                estimated_position_ratemap_shuffled = bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count_RUN,Behaviour,'ratemap shuffle',[],nshuffle,time_bin);
                %         estimated_position_ratemap_shuffled = log_odds_bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,place_cell_index,timebin,[],'ratemap shuffle','','N');
                for nlap = 1:length(cv_groups{track_id}{groupIndex})
                    ratemap_shuffled_probability_ratio{nshuffle}{track_id}{1}{cv_groups{track_id}{groupIndex}(nlap)}  = estimated_position_ratemap_shuffled(1).laps(nlap).probability_ratio;
                    ratemap_shuffled_probability_ratio{nshuffle}{track_id}{1}{cv_groups{track_id}{groupIndex}(nlap)}  =  estimated_position_ratemap_shuffled(2).laps(nlap).probability_ratio;
                end
            end

            for nlap = 1:length(cv_groups{track_id}{groupIndex})
                if length(cv_groups{track_id}{groupIndex})==1
                    estimated_position_lap_CV(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track(1) = temp_estimated_position(1);
                    estimated_position_lap_CV(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track(2) = temp_estimated_position(2);

                    probability_ratio_RUN_lap{1}{track_id}{1}(cv_groups{track_id}{groupIndex}(nlap)) = temp_estimated_position(1).probability_ratio;
                    probability_ratio_RUN_lap{1}{track_id}{2}(cv_groups{track_id}{groupIndex}(nlap)) =  temp_estimated_position(2).probability_ratio;
                else
                    estimated_position_lap_CV(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track(1) = temp_estimated_position(1).laps(nlap);
                    estimated_position_lap_CV(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track(2) = temp_estimated_position(2).laps(nlap);

                    probability_ratio_RUN_lap{1}{track_id}{1}(cv_groups{track_id}{groupIndex}(nlap)) = temp_estimated_position(1).laps(nlap).probability_ratio;
                    probability_ratio_RUN_lap{1}{track_id}{2}(cv_groups{track_id}{groupIndex}(nlap)) =  temp_estimated_position(2).laps(nlap).probability_ratio;
                end
            end


            disp('Speed/position shuffle')
            [N,edges,xbin_WHOLE] = histcounts(estimated_position_WHOLE{track_id}(track_id).run_actual_position,estimated_position_WHOLE{track_id}(track_id).position_bin_centres);
            [N,edges,speed_bin_WHOLE] = histcounts(estimated_position_WHOLE{track_id}(track_id).actual_run_speed,speed_range);

            parfor nshuffle = 1:1000
                bayesian_spike_count_RUN_shuffled = bayesian_spike_count_RUN;
                temp_estimated_position = [];
                for nlap = 1:length(cv_groups{track_id}{groupIndex})
                    lap_indcies = find(bayesian_spike_count_RUN.lap_indices==nlap);
                    [N,edges,xbin] = histcounts(estimated_position_WHOLE{track_id}(track_id).laps(cv_groups{track_id}{groupIndex}(nlap)).run_actual_position,estimated_position_WHOLE{track_id}(1).position_bin_centres);
                    [N,edges,speed_bin] = histcounts(estimated_position_WHOLE{track_id}(track_id).laps(cv_groups{track_id}{groupIndex}(nlap)).actual_run_speed,speed_range);
                    for i = unique(xbin)
                        for j = unique(speed_bin(xbin==i))
                            s = RandStream('mrg32k3a','Seed',nlap*1000+1000*nshuffle+i*100+j*100); % Set random seed for resampling
                            this_time_bin = find(xbin_WHOLE==i&speed_bin_WHOLE~=j);
                            if isempty(this_time_bin)
                                continue
                            end

                            if sum(xbin==i&speed_bin==j) == 0
                                continue
                            end
                            y = randsample(s,this_time_bin,sum(xbin==i&speed_bin==j),'true');
                            bayesian_spike_count_RUN_shuffled.n.run(:,lap_indcies( xbin==i&speed_bin==j)) = bayesian_spike_count_RUN_WHOLE{track_id}.n.run(:,y);
                        end
                    end
                end

                temp_estimated_position = bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count_RUN_shuffled,Behaviour,[],[],[],time_bin);


                for nlap = 1:length(cv_groups{track_id}{groupIndex})
                    if length(cv_groups{track_id}{groupIndex})==1
                        estimated_position_lap_CV_shuffled{nshuffle}(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track(1).run = temp_estimated_position(1).run;
                        estimated_position_lap_CV_shuffled{nshuffle}(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track(2).run = temp_estimated_position(2).run;
                        estimated_position_lap_CV_shuffled{nshuffle}(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track(1).run_error = temp_estimated_position(1).run_error;
                        estimated_position_lap_CV_shuffled{nshuffle}(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track(2).run_error = temp_estimated_position(2).run_error;
                    else
                        estimated_position_lap_CV_shuffled{nshuffle}(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track(1).run = temp_estimated_position(1).laps(nlap).run;
                        estimated_position_lap_CV_shuffled{nshuffle}(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track(2).run = temp_estimated_position(2).laps(nlap).run;
                        estimated_position_lap_CV_shuffled{nshuffle}(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track(1).run_error = temp_estimated_position(1).laps(nlap).run_error;
                        estimated_position_lap_CV_shuffled{nshuffle}(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track(2).run_error = temp_estimated_position(2).laps(nlap).run_error;
                    end
                end
            end


            % first track is track to decode and second track is place field template track
            %             estimated_position_lap_CV(track_id).lap(cv_groups{track_id}{groupIndex}(nlap)).track = temp_estimated_position;

        end
        toc
    end

    probability_ratio_RUN_lap{2} = ratemap_shuffled_probability_ratio;

    % temp_estimated_position = [];
end
end