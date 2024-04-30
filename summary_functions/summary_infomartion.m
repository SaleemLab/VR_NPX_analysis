%%merge all sessions
addpath(genpath('/Users/atom/Documents/GitHub/V1_MEC_acute_MAT'))
% Load clusters_all_all from base_folder
base_folder = fullfile('/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data');
load(fullfile(base_folder, 'clusters_all.mat'));
load(fullfile(base_folder,'spatial_responses_V1_HVA.mat'))

% select clusters you want to plot for population analysis e.g. only plot neurons from V1
V1_index = (clusters_all.region == "V1");
HVA_index =  clusters_all.region =='HVA';
spatially_tuned_neurons = clusters_all.odd_even_stability >=0.95 ...
    & clusters_all.peak_percentile >=0.95;
overall_cluster_index = V1_index & spatially_tuned_neurons(:,1) | HVA_index;
cluster_id = clusters_all.cluster_id(overall_cluster_index);
% Call population_spatial_summary function with cluster_id
% population_spatial_summary(clusters_all, cluster_id, spatial_response,spatial_response_extended);
clusters_all.raw_peak(isnan(clusters_all.raw_peak)) = 0;
ids_index = zeros(size(cluster_id));
ids_to_check = 100:631;
ids_index(ids_to_check) = 1;
mean_spatial_response = cellfun(@(x) mean(x,'omitnan'), spatial_response, 'UniformOutput', false);
peak_spatial_response = cell2mat(cellfun(@(x) max(x,[],2,'omitnan'), mean_spatial_response,'UniformOutput', false));
low_firing_neuron_filter = peak_spatial_response(:,1) < 1 & peak_spatial_response(:,2) < 1;
cluster_summary(clusters_all,cluster_id(~low_firing_neuron_filter & ids_index),...
    'within',spatial_response(~low_firing_neuron_filter & ids_index,:), ...
    spatial_response_extended(~low_firing_neuron_filter & ids_index,:));

%2303205200367,2303202200449

%% plot speed profiles of each session

for iS = 1:10
    figure; % create a new figure
    speed = clusters_all.speed{iS};
    tvec = clusters_all.tvec{iS};
    start_time_all = clusters_all.start_time_all{iS};
    end_time_all = clusters_all.end_time_all{iS};
    w = gausswin(15);
    w = w / sum(w);
    speed(isnan(speed)) = 0;
    speed_smoothed = filtfilt(w,1,speed')';

    plot(tvec, speed); % plot tvec against speed
    hold on; % keep the plot for adding patches
    plot(tvec,speed_smoothed)

    % loop over start and end times
    for i = 1:length(start_time_all)
        % find the indices in tvec corresponding to the start and end times
        start_idx = find(tvec >= start_time_all(i), 1, 'first');
        end_idx = find(tvec <= end_time_all(i), 1, 'last');

        % create a patch
        patch('XData', [tvec(start_idx) tvec(end_idx) tvec(end_idx) tvec(start_idx)], ...
            'YData', [min(speed) min(speed) max(speed) max(speed)], ...
            'FaceColor', 'b', 'FaceAlpha', 0.1);
    end

    hold off; % release the plot

end

%% test if smoothing changes the extended positions
speed = clusters_all.speed{5};
tvec = clusters_all.tvec{5};
start_time_all = clusters_all.start_time_all{5};
end_time_all = clusters_all.end_time_all{5};
position = clusters_all.position{5};
track_ID_all = clusters_all.track_ID_all{5};
t1_index = find(track_ID_all == 1);
bin_no = 140; no_t1_laps = length(t1_index);
t2_index = find(track_ID_all == 2);
no_t2_laps = length(t2_index);
w = gausswin(15);
w = w / sum(w);
speed(isnan(speed)) = 0;
speed_smoothed = filtfilt(w,1,speed')';
distance = speed_smoothed*(1/60);
distance(isnan(distance)) = 0;
figure;
for iLap = 1: 16
    subplot(4,4,iLap)
    if iLap < length(t1_index)
        lap_start_time =  start_time_all(t1_index(iLap));
        lap_end_time =  start_time_all(t1_index(iLap+1));


    else
        lap_start_time =  start_time_all(t1_index(iLap));
        lap_end_time =  end_time_all(t1_index(iLap));


    end

    lap_tvec_index = find(tvec >=lap_start_time & tvec < lap_end_time);
    lap_position = cumsum(distance(lap_tvec_index));
    hold on;
    plot(position(lap_tvec_index));
    plot(lap_position);

    max_distance = max(lap_position);

    bin_edges = 0:ceil(max_distance);

    position_index = discretize(lap_position,bin_edges);

end