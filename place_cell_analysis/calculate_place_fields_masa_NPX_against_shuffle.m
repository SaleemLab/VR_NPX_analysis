function place_fields = calculate_place_fields_masa_NPX_against_shuffle(clusters,Task_info,Behaviour,x_window,x_bin_width,lap_option)
% INPUTS:
%   x_bin_width: enter width value (5 for fine resolution, 14 for bayesian
%   decoding).
% loads list_of_parameters.m, extracted_clusters.mat,extracted_position.mat, extracted_waveform.mat
% uses function skaggs_information.m
% Main difference is to use bin centres rather than edges

%% Place field calculation

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
for iCluster = 1:length(place_fields(1).raw)
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

    place_fields(1).across_tracks_correlation(iCluster,:,:) = lap_correlation;
    place_fields(1).t1_t2_corr(iCluster) = t1_t2_corr;
    place_fields(2).across_tracks_correlation(iCluster,:,:) = lap_correlation;
    place_fields(2).t1_t2_corr(iCluster) = t1_t2_corr;

    for nshuffle = 1:1000
        lap_correlation = corr(normalize(place_fields_shuffled{nshuffle}(1).raw{iCluster}','range'),...
            normalize(place_fields_shuffled{nshuffle}(2).raw{iCluster}','range'));
        lap_correlation(lap_correlation==1) = nan;
        t1_t2_corr_shuffled(nshuffle) = mean(mean(lap_correlation,'omitnan'),'omitnan');
    end

    place_fields(1).t1_t2_corr_shuffled(iCluster,:) = t1_t2_corr_shuffled;
    place_fields(1).t1_t2_remapping(iCluster) = sum(t1_t2_corr > t1_t2_corr_shuffled)/length(t1_t2_corr_shuffled);
    place_fields(2).t1_t2_corr_shuffled(iCluster,:) = t1_t2_corr_shuffled;
    place_fields(2).t1_t2_remapping(iCluster) = sum(t1_t2_corr > t1_t2_corr_shuffled)/length(t1_t2_corr_shuffled);

    % Reliability/ explained variance (currently taking 80 seconds to run one cell)

    %     if  place_fields(1).first_second_stability(iCluster) > 0.99 | place_fields(2).first_second_stability(iCluster) > 0.99 % if stability better than shuffled data
    %         spikeTimes = clusters.spike_times(clusters.spike_id == place_fields(1).cluster_id(iCluster));
    %         [reliability,shuffled_reliability] = spatial_cell_reliability_analysis(spikeTimes,Behaviour,position_shuffled,Task_info,x_window,2,1);
    %     else
    spikeTimes = clusters.spike_times(clusters.spike_id == place_fields(1).cluster_id(iCluster));
    [reliability,~] = spatial_cell_reliability_analysis(spikeTimes,Behaviour,position_shuffled,Task_info,x_window,2,0);
%     shuffled_reliability = zeros(1000,2,10); % 1000 shuffles x 2 tracks x 10 folds
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

    putative_place_cells = find(place_fields(track_id).raw_peak >= 1 ...
        & place_fields(track_id).first_second_stability >= 0.99 ...
        & place_fields(track_id).peak_percentile >= 0.99 ...
        & place_fields(track_id).skaggs_percentile >= 0.99); % If peak FR, skaggs_info and 1st vs 2nd half stability are significantly better than shuffle distribution.

%         & mean(place_fields(track_id).explained_variance,2,'omitnan')' > 0 ...


    %         % Set a less conservative criteria for place cells, having to pass either peak firing rate thresholds (smoothed PF and raw PF)
    %         putative_place_cells_LIBERAL = find(place_fields.track(track_id).peak >= parameters.min_smooth_peak...
    %             | place_fields.track(track_id).raw_peak >= parameters.min_raw_peak...
    %             & place_fields.mean_rate <= parameters.max_mean_rate...
    %             & place_fields.track(track_id).skaggs_info > 0);
    %


    if ~isempty(clusters.cell_type)
        pyramidal_cells = find(clusters.cell_type == 1)';
        place_fields(track_id).good_cells = intersect(putative_place_cells,pyramidal_cells); % Check that the cells that passed the threshold are pyramidal cells
        place_fields(track_id).good_cells_LIBERAL = putative_place_cells; % No cell type restrictions
    else
        place_fields(track_id).good_cells_LIBERAL = putative_place_cells; % No cell type restrictions
    end
    %
    % Sort place fields according to the location of their peak
    all_maps = reshape([place_fields(track_id).raw{:}],[size(place_fields(track_id).raw{1}) length(place_fields(track_id).raw)]);
    average_maps = squeeze(mean(all_maps,1));
    [~,peak_id]=max(average_maps);
    [~,index] = sort(peak_id);
    place_fields(track_id).sorted = index;

    all_maps = reshape([place_fields(track_id).raw{place_fields(track_id).good_cells}],[size(place_fields(track_id).raw{1}) length(place_fields(track_id).good_cells)]);
    average_maps = squeeze(mean(all_maps,1));
    [~,peak_id]=max(average_maps);
    [~,index1] = sort(peak_id);
    place_fields(track_id).sorted_good_cells = place_fields(track_id).good_cells(index1);

    all_maps = reshape([place_fields(track_id).raw{place_fields(track_id).good_cells_LIBERAL}],[size(place_fields(track_id).raw{1}) length(place_fields(track_id).good_cells_LIBERAL)]);
    average_maps = squeeze(mean(all_maps,1));
    [~,peak_id]=max(average_maps);
    [~,index2] = sort(peak_id);
    place_fields(track_id).sorted_good_cells_LIBERAL = place_fields(track_id).good_cells_LIBERAL(index2);
end

for track_id = 1:2
    place_fields(track_id).all_good_cells_LIBERAL = unique([place_fields(:).good_cells_LIBERAL]);
    place_fields(track_id).common_good_cells_LIBERAL = intersect(place_fields(1).good_cells_LIBERAL,place_fields(2).good_cells_LIBERAL);
    place_fields(track_id).all_good_cells = unique([place_fields(:).good_cells]);
    place_fields(track_id).common_good_cells = intersect(place_fields(1).good_cells,place_fields(2).good_cells);
end

[C,i1,i2] = setxor(place_fields(1).good_cells_LIBERAL,place_fields(2).good_cells_LIBERAL);
place_fields(1).unique_good_cells_LIBERAL = place_fields(1).good_cells_LIBERAL(i1);
place_fields(2).unique_good_cells_LIBERAL = place_fields(2).good_cells_LIBERAL(i2);

[C,i1,i2] = setxor(place_fields(1).good_cells,place_fields(2).good_cells);
place_fields(1).unique_good_cells = place_fields(1).good_cells(i1);
place_fields(2).unique_good_cells = place_fields(2).good_cells(i2);

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
