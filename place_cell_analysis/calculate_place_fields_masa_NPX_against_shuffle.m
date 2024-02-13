function place_fields = calculate_place_fields_masa_NPX_against_shuffle(clusters,Task_info,Behaviour,x_window,x_bin_width,lap_option)
% INPUTS:
%   x_bin_width: enter width value (5 for fine resolution, 14 for bayesian
%   decoding).
% loads list_of_parameters.m, extracted_clusters.mat,extracted_position.mat, extracted_waveform.mat
% uses function skaggs_information.m
% Main difference is to use bin centres rather than edges

%% Place field calculation

place_fields = calculate_spatial_cells(clusters,Task_info,Behaviour,x_window,x_bin_width);

disp('Circularly shifting position for each lap and re-calculate spatial firing rate')


tic
place_fields_shuffled = cell(1000,1);
position_shuffled = cell(1000,1);

parfor nshuffle = 1:1000
    %     tic
    disp(sprintf('Shuffle %i',nshuffle))
    position_shuffled{nshuffle} = Behaviour.position;

    for nlap = 1:length(Task_info.start_time_all)
        s = RandStream('mrg32k3a','Seed',nlap+nshuffle*1000); % Set random seed for resampling

        position_shuffled{nshuffle}(Behaviour.tvec  >= Task_info.start_time_all(nlap) & Behaviour.tvec <= Task_info.end_time_all(nlap)) =...
            circshift(Behaviour.position(Behaviour.tvec  >= Task_info.start_time_all(nlap) & Behaviour.tvec <= Task_info.end_time_all(nlap))...
            ,randi(s,length(Behaviour.tvec  >= Task_info.start_time_all(nlap) & Behaviour.tvec <= Task_info.end_time_all(nlap))));
    end

    Behaviour_temp = Behaviour;
    Behaviour_temp.position = position_shuffled{nshuffle};

    place_fields_shuffled{nshuffle} = calculate_spatial_cells(clusters,Task_info,Behaviour_temp,x_window,x_bin_width);
    %     toc
end
toc

% save('place_fields_shuffled.mat','place_fields_shuffled','-v7.3')
%% Quantify stability


%     if ~isempty(lap_option)
%         if contains(lap_option,'even')
%             selected_laps = 2:2:length(start_times);
%         elseif contains(lap_option,'odd')
%             selected_laps = 1:2:length(start_times);
%         elseif contains(lap_option,'first')
%             selected_laps = 1:2:length(start_times);
%         elseif contains(lap_option,'second')
%             selected_laps = 1:2:length(start_times);
%
%         elseif isnumeric(lap_option) % can input one lap or multiple laps (based on block or behaviour)
%             selected_laps = lap_option;
%         end
%     end

%         for nshuffle = 1:length(place_fields_shuffled)
%             shuffled_ratemap(nshuffle,:,:) = place_fields_shuffled{nshuffle}(track_id).raw{iCluster};
%         end
%
%         odd_map_shuffled = squeeze(mean(shuffled_ratemap(:,1:2:length(start_times),:),2));
%         even_map_shuffled = squeeze(mean(shuffled_ratemap(:,2:2:length(start_times),:),2));
%
%         odd_map = mean(place_fields_shuffled{nshuffle}(track_id).raw{iCluster}(1:2:length(start_times),:));
%         even_map = mean(place_fields_shuffled{nshuffle}(track_id).raw{iCluster}(2:2:length(start_times),:));
%
%         first_half_map_shuffled = squeeze(mean(shuffled_ratemap(:,1:2:length(start_times),:),2));
%         second_half_map_shuffled = squeeze(mean(shuffled_ratemap(:,2:2:length(start_times),:),2));
%
%         first_half_map = mean(place_fields(track_id).raw{iCluster}(1:length(start_times)/2,:));
%         second_half_map = mean(place_fields(track_id).raw{iCluster}(length(start_times)/2+1:end,:));
tic
% Spatial cell stability and peak
for iCluster = 1:length(place_fields(track_id).raw)
tic
    % Stability
    for track_id = 1:max(Behaviour.track_ID)
        start_times = Task_info.start_time_all(Task_info.track_ID_all == track_id);

        lap_correlation = corr(normalize(place_fields(track_id).raw{iCluster}','range'),...
            normalize(place_fields(track_id).raw{iCluster}','range')); % lap by lap correlation

        first_second_corr = mean(mean(lap_correlation(1:2:end,2:2:end),'omitnan'),'omitnan'); % First half vs Second half
        odd_even_corr = mean(mean(lap_correlation(1:round(length(start_times)/2),round(length(start_times)/2)+1:end),'omitnan'),'omitnan'); % Odd laps vs even laps

        place_fields(track_id).within_track_correlation(iCluster,:,:) = lap_correlation;
        place_fields(track_id).first_second_corr(iCluster) = first_second_corr;
        place_fields(track_id).odd_even_corr(iCluster) = odd_even_corr;

        for nshuffle = 1:1000
            peak_shuffled(nshuffle) = max(mean(place_fields_shuffled{nshuffle}(track_id).raw{iCluster}));
            skaggs_info_shuffled(nshuffle) = place_fields_shuffled{nshuffle}(track_id).skaggs_info(iCluster);

            lap_correlation = corr(normalize(place_fields_shuffled{nshuffle}(track_id).raw{iCluster}','range'),...
                normalize(place_fields_shuffled{nshuffle}(track_id).raw{iCluster}','range'));

            first_second_corr_shuffled(nshuffle) = mean(mean(lap_correlation(1:2:end,2:2:end),'omitnan'),'omitnan');
            odd_even_corr_shuffled(nshuffle) = mean(mean(lap_correlation(1:round(length(start_times)/2),round(length(start_times)/2)+1:end),'omitnan'),'omitnan');
        end

        place_fields(track_id).first_second_corr_shuffled(iCluster,:) = first_second_corr_shuffled;
        place_fields(track_id).first_second_stability(iCluster) = sum(first_second_corr > first_second_corr_shuffled)/length(first_second_corr_shuffled);

        place_fields(track_id).odd_even_corr_shuffled(iCluster,:) = odd_even_corr_shuffled;
        place_fields(track_id).odd_even_stability(iCluster) = sum(odd_even_corr > odd_even_corr_shuffled)/length(odd_even_corr_shuffled);

        place_fields(track_id).raw_peak(iCluster) = max(mean(place_fields(track_id).raw{iCluster})); % peak FR
        place_fields(track_id).peak_percentile(iCluster) = sum(max(mean(place_fields(track_id).raw{iCluster}))>peak_shuffled)/length(peak_shuffled); % peak FR relative to peak of shuffled data

        place_fields(track_id).skaggs_percentile(iCluster) = sum(place_fields(track_id).skaggs_info(iCluster)>skaggs_info_shuffled)/length(skaggs_info_shuffled);
    end

    % Remapping (correlation between track 1 and track 2)
    lap_correlation = corr(normalize(place_fields(1).raw{iCluster}','range'),...
        normalize(place_fields(2).raw{iCluster}','range'));
    lap_correlation(lap_correlation==1) = nan;
    t1_t2_corr = mean(mean(lap_correlation,'omitnan'),'omitnan');

    place_fields(track_id).across_tracks_correlation(iCluster,:,:) = lap_correlation;
    place_fields(track_id).t1_t2_corr(iCluster) = t1_t2_corr;

    for nshuffle = 1:1000
        lap_correlation = corr(normalize(place_fields_shuffled{nshuffle}(1).raw{iCluster}','range'),...
            normalize(place_fields_shuffled{nshuffle}(2).raw{iCluster}','range'));
        lap_correlation(lap_correlation==1) = nan;
        t1_t2_corr_shuffled(nshuffle) = mean(mean(lap_correlation,'omitnan'),'omitnan');
    end

    place_fields(track_id).t1_t2_corr_shuffled(iCluster,:) = t1_t2_corr_shuffled;
    place_fields(track_id).t1_t2_remapping(iCluster) = sum(t1_t2_corr > t1_t2_corr_shuffled)/length(t1_t2_corr_shuffled);

    place_fields(1).t1_t2_remapping(iCluster) = sum(t1_t2_corr > t1_t2_corr_shuffled)/length(t1_t2_corr_shuffled);
    place_fields(1).t1_t2_remapping(iCluster) = sum(t1_t2_corr > t1_t2_corr_shuffled)/length(t1_t2_corr_shuffled);

    % Reliability/ explained variance (currently taking 80 seconds to run one cell)

    %     if  place_fields(1).first_second_stability(iCluster) > 0.99 | place_fields(2).first_second_stability(iCluster) > 0.99 % if stability better than shuffled data
    %         spikeTimes = clusters.spike_times(clusters.spike_id == place_fields(1).cluster_id(iCluster));
    %         [reliability,shuffled_reliability] = spatial_cell_reliability_analysis(spikeTimes,Behaviour,position_shuffled,Task_info,x_window,2,1);
    %     else
    spikeTimes = clusters.spike_times(clusters.spike_id == place_fields(1).cluster_id(iCluster));
    [reliability,~] = spatial_cell_reliability_analysis(spikeTimes,Behaviour,position_shuffled,Task_info,x_window,1,0);
    shuffled_reliability = zeros(1000,2,10); % 1000 shuffles x 2 tracks x 10 folds
    %     end
    reliability(isinf(reliability)) = nan;
    place_fields(1).explained_variance(iCluster,:) = reliability(1,:);
    place_fields(2).explained_variance(iCluster,:) = reliability(2,:);
%     place_fields(1).reliability(iCluster,:) = sum(mean(reliability(1,:),'omitnan')>mean(squeeze(shuffled_reliability(:,1,:)),2,'omitnan'))/1000;
%     place_fields(2).reliability(iCluster,:) = sum(mean(reliability(2,:),'omitnan')>mean(squeeze(shuffled_reliability(:,2,:)),2,'omitnan'))/1000;

%     place_fields(2).reliability(iCluster,:) = mean(reliability(2,:));
%     sum(t1_t2_corr > t1_t2_corr_shuffled)/length(t1_t2_corr_shuffled);
%     place_fields(1).explained_variance_shuffled(iCluster,:,:) = shuffled_reliability(:,1,:);
%     place_fields(2).explained_variance_shuffled(iCluster,:,:) = shuffled_reliability(:,2,:);
toc
end
toc

for track_id = 1:max(Behaviour.track_ID)
    % Find cells that pass the 'Place cell' thresholds -
    % both peak of smoothed place field or peak of raw place field need to be above the respective thresholds
    place_fields(track_id).explained_variance(isinf(place_fields(track_id).explained_variance)) = nan;

    putative_place_cells = find(place_fields(track_id).raw_peak >= 1 ...
        & mean(place_fields(track_id).explained_variance,2,'omitnan')' > 0 ...
        & place_fields(track_id).first_second_stability >= 0.99 ...
        & place_fields(track_id).peak_percentile >= 0.99); % rather than based on skaggs_info. Now based on if reliability and stability is significantly better than shuffle distribution.

    %         % Set a less conservative criteria for place cells, having to pass either peak firing rate thresholds (smoothed PF and raw PF)
    %         putative_place_cells_LIBERAL = find(place_fields.track(track_id).peak >= parameters.min_smooth_peak...
    %             | place_fields.track(track_id).raw_peak >= parameters.min_raw_peak...
    %             & place_fields.mean_rate <= parameters.max_mean_rate...
    %             & place_fields.track(track_id).skaggs_info > 0);
    %


    if ~isempty(clusters.cell_type)
        pyramidal_cells = clusters.cluster_id(clusters.cell_type == 1);
        place_fields.track(track_id).good_cells = intersect(putative_place_cells,pyramidal_cells); % Check that the cells that passed the threshold are pyramidal cells
        place_fields.track(track_id).good_cells_LIBERAL = putative_place_cells; % No cell type restrictions
    else
        place_fields.track(track_id).good_cells_LIBERAL = putative_place_cells;
    end
    %
    % Sort place fields according to the location of their peak
    all_maps = reshape([place_fields(track_id).raw{:}],[size(place_fields(track_id).raw{1}) length(place_fields(track_id).raw)]);
    average_maps = squeeze(mean(all_maps,1));
    [~,peak_id]=max(average_maps);
    [~,index] = sort(peak_id);
    place_fields(track_id).sorted = index;
    [~,index1] = sort(place_fields.track(track_id).centre(place_fields(track_id).good_cells));
    place_fields.track(track_id).sorted_good_cells = place_fields.track(track_id).good_cells(index1);
    [~,index2] = sort(place_fields.track(track_id).centre(place_fields.track(track_id).good_cells_LIBERAL));
    place_fields.track(track_id).sorted_good_cells_LIBERAL = place_fields.track(track_id).good_cells_LIBERAL(index2);
end
%% Classify cells as good place cells, interneuron, pyramidal cells & other cells

%interneurons classfication
interneurons = find(place_fields.mean_rate > parameters.max_mean_rate);
place_fields.interneurons=interneurons;

good_place_cells=[]; track=[];
for track_id=1:length(position.linear) %good cells classfication
    all_cells = 1:1:length(place_fields.track(track_id).raw);
    good_place_cells = [good_place_cells place_fields.track(track_id).sorted_good_cells];
    track =[track track_id*ones(size(place_fields.track(track_id).sorted_good_cells))];
end
place_fields.good_place_cells = unique(good_place_cells);
place_fields.all_cells = unique(all_cells);

good_place_cells_LIBERAL=[];
for track_id=1:length(position.linear) %good cells (liberal threshold) classfication
    good_place_cells_LIBERAL = [good_place_cells_LIBERAL place_fields.track(track_id).sorted_good_cells_LIBERAL];
end
place_fields.good_place_cells_LIBERAL = unique(good_place_cells_LIBERAL);

% cells that are unique for each track
unique_cells=[];
for track_id = 1:length(position.linear)
    place_fields.track(track_id).unique_cells = setdiff(good_place_cells(track==track_id),good_place_cells(track~=track_id),'stable');
    unique_cells = [unique_cells, place_fields.track(track_id).unique_cells];
end
place_fields.unique_cells = unique_cells;  % all cells that have good place fields only on a single track

% putative pyramidal cells classification:  pyramidal cells that pass the 'Pyramidal type' threshold (but not need to be place cells)
putative_pyramidal_cells = find(place_fields.mean_rate <= parameters.max_mean_rate);

if ~isempty(clusters.cell_type)
    place_fields.pyramidal_cells = intersect(putative_pyramidal_cells,pyramidal_cells);
else
    place_fields.pyramidal_cells = putative_pyramidal_cells;
end
place_fields.pyramidal_cells=unique(place_fields.pyramidal_cells);

other_cells = setdiff(1:length(unique(clusters.spike_id)),good_place_cells,'stable'); %find the excluded putative pyramidal cells
place_fields.other_cells = setdiff(other_cells,interneurons,'stable'); %remove also the interneurons

%save place fields (different filenames used based on x_bins_width chosen)
%     if x_bins_width== parameters.x_bins_width_bayesian
%         place_fields_BAYESIAN=place_fields;
%         save extracted_place_fields_BAYESIAN place_fields_BAYESIAN;
%     elseif x_bins_width== parameters.x_bins_width
%         save extracted_place_fields place_fields;
%     else disp('error: x_bin_width does not match expected value')
%     end


close(H)
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
