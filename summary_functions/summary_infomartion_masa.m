%%merge all sessions
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis\'))
% Load clusters_all_all from base_folder
base_folder = fullfile('Z:\ibn-vision\USERS\Masa');
load(fullfile(base_folder, 'clusters_all.mat'));
load(fullfile(base_folder,'spatial_responses_V1.mat'))

% select clusters you want to plot for population analysis e.g. only plot neurons from V1
V1_index = contains(clusters_all.region,"V1");
spatially_tuned_neurons = clusters_all.odd_even_stability >=0.95 ...
    & clusters_all.peak_percentile >=0.95;
overall_cluster_index = V1_index & spatially_tuned_neurons(:,1);
cluster_id = clusters_all.cluster_id(overall_cluster_index);
% Call population_spatial_summary function with cluster_id
% population_spatial_summary(clusters_all, cluster_id, spatial_response,spatial_response_extended);
clusters_all.raw_peak(isnan(clusters_all.raw_peak)) = 0;
ids_index = zeros(size(cluster_id));
ids_to_check = 1:97;
ids_index(ids_to_check) = 1;
mean_spatial_response = cellfun(@(x) mean(x,'omitnan'), spatial_response, 'UniformOutput', false);
peak_spatial_response = cell2mat(cellfun(@(x) max(x,[],2,'omitnan'), mean_spatial_response,'UniformOutput', false));
low_firing_neuron_filter = peak_spatial_response(:,1) < 1 & peak_spatial_response(:,2) < 1;
cluster_summary(clusters_all,cluster_id(~low_firing_neuron_filter & ids_index),...
    'within',spatial_response(~low_firing_neuron_filter & ids_index,:), ...
    spatial_response_extended(~low_firing_neuron_filter & ids_index,:));

%% calculate and plot V1 population SMIs
SMI = nan(length(cluster_id),5);
for iC = 1:length(cluster_id)
    current_cluster = cluster_spatial_info(clusters_all, cluster_id(iC),...
     spatial_response(iC,:), spatial_response_extended(iC,:));
    SMI_tmp = current_cluster.SMI;
    SMI_t1_tmp = SMI_tmp{1,1};
    SMI_t2_tmp = SMI_tmp{1,2};
    SMI(iC,:) = [SMI_t1_tmp,SMI_t2_tmp];
end
SMI_t1 = SMI(:,1:2);
max_SMI_t1 = max(SMI_t1,[],2);
[SMI_count,SMI_edges] = histcounts(max_SMI_t1,100);
figure;
% no_neuron = length(SMI_tmp);
% cum_prob = cumsum(ones([no_neuron,1]).*(1/no_neuron));
plot(SMI_edges(2:end),cumsum(SMI_count)./sum(SMI_count),'LineWidth',3);
xlim([-1,1])
ylim([0,1])
xline(0,'r')
yline(0.5,'r')
SMI_t2 = SMI(:,3:4);
max_SMI_t2 = max(SMI_t2,[],2);
[SMI_count,SMI_edges] = histcounts(max_SMI_t2,100);
hold on;
% no_neuron = length(SMI_tmp);
% cum_prob = cumsum(ones([no_neuron,1]).*(1/no_neuron));
plot(SMI_edges(2:end),cumsum(SMI_count)./sum(SMI_count),'LineWidth',3);

SMI_t1 = SMI(:,1:2);
mean_SMI_t1 = mean(SMI_t1,2);
[SMI_count,SMI_edges] = histcounts(mean_SMI_t1,100);
figure;
% no_neuron = length(SMI_tmp);
% cum_prob = cumsum(ones([no_neuron,1]).*(1/no_neuron));
plot(SMI_edges(2:end),cumsum(SMI_count)./sum(SMI_count),'LineWidth',3);
xlim([-1,1])
ylim([0,1])
xline(0,'r')
yline(0.5,'r')
SMI_t2 = SMI(:,3:4);
mean_SMI_t2 = mean(SMI_t2,2);
[SMI_count,SMI_edges] = histcounts(mean_SMI_t2,100);
hold on;
% no_neuron = length(SMI_tmp);
% cum_prob = cumsum(ones([no_neuron,1]).*(1/no_neuron));
plot(SMI_edges(2:end),cumsum(SMI_count)./sum(SMI_count),'LineWidth',3);


figure;
for iT = 1:5
subplot(3,2,iT)
SMI_tmp = SMI(:,iT);
SMI_tmp = SMI_tmp(~isnan(SMI_tmp));
[SMI_count,SMI_edges] = histcounts(SMI_tmp,50);
% no_neuron = length(SMI_tmp);
% cum_prob = cumsum(ones([no_neuron,1]).*(1/no_neuron));
% plot(SMI_edges(2:end),cumsum(SMI_count)./sum(SMI_count),'LineWidth',3);
plot(SMI_edges(2:end),SMI_count)
hold on;
xline(0,'r')
yline(0.5,'r')
end
figure;
for iT = 1:5
SMI_tmp = SMI(:,iT);
SMI_tmp = SMI_tmp(~isnan(SMI_tmp));
[SMI_count,SMI_edges] = histcounts(SMI_tmp,30);
% no_neuron = length(SMI_tmp);
% cum_prob = cumsum(ones([no_neuron,1]).*(1/no_neuron));
% plot(SMI_edges(2:end),cumsum(SMI_count)./sum(SMI_count),'LineWidth',3);
% histogram(SMI_tmp,50)
plot(SMI_edges(2:end),SMI_count./sum(SMI_count),'LineWidth',3)
hold on;
end
xlim([-1,1])
ylim([0,0.22])
xline(0,'r')
yline(0.5,'r')
legend({'t1 A','t1 B','t2 C','t2 B','t2 Y'})
title('Spatial Modulation Index Cumulative Probability')
%2303205200367,2303202200449
ids_to_check = ismember(cluster_id,interesting_ids);
interesting_ids = [
2303202100449;
2303202100462;
2303202100515;
2303202100516;
2303202100527;
2303202100547;
2303202100557;
2303202100632;
2303202100637;
2303202100651;
2303202100669;
2303203100477;
2303203100488;
2303203100489;
2303203100497;
2303203100525;
2303203100526;
2303203100571;
2303203100574;
2303203100579;
2303203100599;
2303203100604;
2303203100605;
2303203100608;
2303203100610;
2303203100615;
2303203100631;
2303203100636;
2303204100531;
2303204100657;
2303204100675;
2303204100700;
2303204100708;
2303204100723;
2303205100364;
2303205100399;
2303205100404;
2303205100417;
2303205100429;
2303205100439;
2303205100440;
2303205100486;
2303205100499;
2303205100527;
2303205100532;
2303401100388;
2303401100391;
2303401100406;
2303401100432;
2303401100434;
2303401100440;
2303401100456;
2303401100469;
2303401100499;
2303401100503;
2303401100528;
2303401100580;
2303401100582;
2303402100176;
2303402100186;
2303402100187;
2303402100189;
2303402100200;
2303402100229;
2303402100268;
2303402100274;
2303402100296;
2303402100319;
2303402100320;
2303402100362;
2303701100688;
2303701100709;
2303701100722;
2303702100577;
2303702100585;
2303702100595;
2303702100600;
2303702100665;
2303702100668;
2303702100670;
2303703100506;
2303703100515;
2303703100552;
2303703100581;
2303703100595;
2303703100643;
2303703100662;
2303704100394;
2303801100308;
2303801100322;
2303801100332;
2303801100395;
2303801100429;
2303801100467;
2303802100257;
2303802100263;
2303802100354;
2303802100359];


%% test the effect of delay
cluster_index = find(clusters_all.cluster_id == 2303202100557);
iS = clusters_all.session_count(cluster_index);
spike_times = clusters_all.spike_times{cluster_index};
tvec = clusters_all.tvec{iS};
start_time_all = clusters_all.start_time_all{iS};
end_time_all = clusters_all.end_time_all{iS};
track_ID_all = clusters_all.track_ID_all{iS};
position = clusters_all.position{iS};
speed = clusters_all.speed{iS};
delay = 0:0.04:0.56;
for iD = 1:15
cluster_sr = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'within',delay(iD));
cluster_sr_ex = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'extension',delay(iD));
cluster_summary(clusters_all,2303202100637,'within',cluster_sr,cluster_sr_ex,...
    'plot_behaviour',0,'break_loop',1);
end

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