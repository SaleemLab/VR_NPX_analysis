function fig = population_spatial_summary(clusters_all,cluster_ids,method)
    if nargin < 3
        method = 'default';
    end
% this function aims to summarize the spatial aspects of selected population from cluster_ids
cluster_ids_index = ismember(clusters_all.cluster_id,cluster_ids);

all_spatial_response_t1_psth = zeros(length(cluster_ids),140);
all_spatial_response_t2_psth = zeros(length(cluster_ids),140);
odd_spatial_response_t1_psth = zeros(length(cluster_ids),140);
even_spatial_response_t1_psth = zeros(length(cluster_ids),140);
odd_spatial_response_t2_psth = zeros(length(cluster_ids),140);
even_spatial_response_t2_psth = zeros(length(cluster_ids),140);
sessions = unique(clusters_all.session_count(cluster_ids_t1_index));
for iS = 1:length(sessions)
    this_session_t1 = clusters_all.session_count == sessions(iS);
    this_session_t2 = clusters_all.session_count == sessions(iS);
    rewarded_laps = clusters_all.rewarded_lap_id{sessions(iS)};
    track_ID_all = clusters_all.track_ID_all{sessions(iS)};
    rewarded_laps_index = zeros(length(track_ID_all),1);
    rewarded_laps_index(rewarded_laps) = 1;
    rewarded_laps_index = logical(rewarded_laps_index);
    t1_index = ones(sum(track_ID_all == 1),1);
    t2_index = ones(sum(track_ID_all == 2),1);
    rewarded_laps_t1 = t1_index & rewarded_laps_index(track_ID_all == 1);
    rewarded_laps_t2 = t2_index & rewarded_laps_index(track_ID_all == 2);
    switch method
        case 'default'
            laps_index_t1 = rewarded_laps_t1;
            laps_index_t2 = rewarded_laps_t2;
        case 'active'

    end

    laps_t1_tmp = 1:length(laps_index_t1);
    laps_t2_tmp = 1:length(laps_index_t2);
    odd_laps_t1 = mod(laps_t1_tmp,2) == 1;
    even_laps_t1 = mod(laps_t1_tmp,2) == 0;
    odd_laps_t2 = mod(laps_t2_tmp,2) == 1;
    even_laps_t2 = mod(laps_t2_tmp,2) == 0;
    clusters_t1_in_this_session_index = this_session_t1 & cluster_ids_t1_index;
    clusters_t2_in_this_session_index = this_session_t2 & cluster_ids_t2_index;
    spatial_response_t1_this_session = cell2mat(clusters_all.spatial_response(clusters_t1_in_this_session_index));
    spatial_response_t2_this_session = cell2mat(clusters_all.spatial_response(clusters_t2_in_this_session_index));
    clusters_ids_index_t1 = ismember(cluster_ids,clusters_all.cluster_id(clusters_t1_in_this_session_index));
    clusters_ids_index_t2 = ismember(cluster_ids,clusters_all.cluster_id(clusters_t2_in_this_session_index));
    all_spatial_response_t1_psth(clusters_ids_index_t1,:) = squeeze(mean(spatial_response_t1_this_session(:,laps_index_t1,:),2,'omitnan')); 
    all_spatial_response_t2_psth(clusters_ids_index_t2,:) = squeeze(mean(spatial_response_t2_this_session(:,laps_index_t2,:),2,'omitnan'));
    odd_spatial_response_t1_psth(clusters_ids_index_t1,:) = squeeze(mean(spatial_response_t1_this_session(:,laps_index_t1(odd_laps_t1),:),2,'omitnan'));
    even_spatial_response_t1_psth(clusters_ids_index_t1,:) = squeeze(mean(spatial_response_t1_this_session(:,laps_index_t1(even_laps_t1),:),2,'omitnan'));
    odd_spatial_response_t2_psth(clusters_ids_index_t2,:) = squeeze(mean(spatial_response_t2_this_session(:,laps_index_t2(odd_laps_t2),:),2,'omitnan'));
    even_spatial_response_t2_psth(clusters_ids_index_t2,:) = squeeze(mean(spatial_response_t2_this_session(:,laps_index_t2(even_laps_t2),:),2,'omitnan'));
end
% sort the spatial response by the peak location (or anything you want in the future)

[~,peak_location_t1] = max(all_spatial_response_t1_psth,[],2);
[~,peak_location_t2] = max(all_spatial_response_t2_psth,[],2);
[~,sort_index_t1] = sort(peak_location_t1);
[~,sort_index_t2] = sort(peak_location_t2);
[~,odd_peak_location_t1] = max(odd_spatial_response_t1_psth,[],2);
[~,odd_peak_location_t2] = max(odd_spatial_response_t2_psth,[],2);
[~,even_peak_location_t1] = max(even_spatial_response_t1_psth,[],2);
[~,even_peak_location_t2] = max(even_spatial_response_t2_psth,[],2);
[~,sort_index_odd_t1] = sort(odd_peak_location_t1);
[~,sort_index_even_t1] = sort(even_peak_location_t1);
[~,sort_index_odd_t2] = sort(odd_peak_location_t2);
[~,sort_index_even_t2] = sort(even_peak_location_t2);
all_spatial_response_t1_psth = all_spatial_response_t1_psth(sort_index_t1,:);
all_spatial_response_t2_psth = all_spatial_response_t2_psth(sort_index_t2,:);
odd_spatial_response_t1_psth = odd_spatial_response_t1_psth(sort_index_even_t1,:);
even_spatial_response_t1_psth = even_spatial_response_t1_psth(sort_index_odd_t1,:);
odd_spatial_response_t2_psth = odd_spatial_response_t2_psth(sort_index_even_t2,:);
even_spatial_response_t2_psth = even_spatial_response_t2_psth(sort_index_odd_t2,:);
fig = figure;
subplot(4,5,[1 2 6 7])
% plot t1 and t2 snake plots of the same population
h1 = imagesc(normalize(flip(all_spatial_response_t1_psth),2,'scale'));
set(h1,'ButtonDownFcn',{@mycallback,cluster_ids(sort_index_t1)});
clim([0 1]);colormap(flip(gray(50)));colorbar
xlabel('position(cm)')
ylabel('units')
title('t1')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

subplot(4,5,[11 12 16 17])
h2 = imagesc(normalize(flip(all_spatial_response_t2_psth),2,'scale'));
set(h2,'ButtonDownFcn',{@mycallback,cluster_ids(sort_index_t2)});
clim([0 1]);colormap(flip(gray(50)));colorbar
xlabel('position(cm)')
ylabel('units')
title('t2')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

subplot(4,5,3)
h3 = imagesc(normalize(flip(odd_spatial_response_t1_psth),2,'scale'));
set(h3,'ButtonDownFcn',{@mycallback,cluster_ids(sort_index_even_t1)});
clim([0 1]);colormap(flip(gray(50)));
ylabel('units')
title('odd sorted by even')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

subplot(4,5,8)
h4 =imagesc(normalize(flip(even_spatial_response_t1_psth),2,'scale'));
set(h4,'ButtonDownFcn',{@mycallback,cluster_ids(sort_index_odd_t1)});
clim([0 1]);colormap(flip(gray(50)));
ylabel('units')
title('even sorted by odd')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

subplot(4,5,13)
h5 = imagesc(normalize(flip(odd_spatial_response_t2_psth),2,'scale'));
set(h5,'ButtonDownFcn',{@mycallback,cluster_ids(sort_index_even_t2)});
clim([0 1]);colormap(flip(gray(50)));
ylabel('units')
title('odd sorted by even')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

subplot(4,5,18)
h6 = imagesc(normalize(flip(even_spatial_response_t2_psth),2,'scale'));
set(h6,'ButtonDownFcn',{@mycallback,cluster_ids(sort_index_odd_t2)});
clim([0 1]);colormap(flip(gray(50)));
ylabel('units')
title('even sorted by odd')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

%% plot the peak in the middle and -40 to 40cm from the peak
peak_response_t1 = zeros(length(cluster_ids),81);
peak_response_t2 = zeros(length(cluster_ids),81);
for iC = 1:length(cluster_ids)
peak_location_t1 = find(all_spatial_response_t1_psth(iC,:) == max(all_spatial_response_t1_psth(iC,:)));
peak_location_t2 = find(all_spatial_response_t2_psth(iC,:) == max(all_spatial_response_t2_psth(iC,:)));

start_t1 = max(1, peak_location_t1 - 40);
end_t1 = min(size(all_spatial_response_t1_psth, 2), peak_location_t1 + 40);
start_t2 = max(1, peak_location_t2 - 40);
end_t2 = min(size(all_spatial_response_t2_psth, 2), peak_location_t2 + 40);

peak_response_t1(iC,(41-(peak_location_t1-start_t1)):(41+end_t1-peak_location_t1)) = all_spatial_response_t1_psth(iC,start_t1:end_t1);
peak_response_t2(iC,(41-(peak_location_t2-start_t2)):(41+end_t2-peak_location_t2)) = all_spatial_response_t2_psth(iC,start_t2:end_t2);
end
subplot(4,5,[4 5 9 10])
h7 = imagesc(normalize(flip(peak_response_t1),2,'scale'));
set(h7,'ButtonDownFcn',{@mycallback,cluster_ids(sort_index_t1)});
clim([0 1]);colormap(flip(gray(50)));colorbar
xlabel('distance from peak (cm)')
ylabel('units')
title('t1 distance from peak')

subplot(4,5,[14 15 19 20])
h8 = imagesc(normalize(flip(peak_response_t2),2,'scale'));
set(h8,'ButtonDownFcn',{@mycallback,cluster_ids(sort_index_t2)});
clim([0 1]);colormap(flip(gray(50)));colorbar
xlabel('distance from peak (cm)')
ylabel('units')
title('t2 distance from peak')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

function mycallback(~, ~, sort_index)
    ax = gca;
    point = get(ax, 'CurrentPoint');
    row = round(point(1, 2));
    if row >= 1 && row <= length(sort_index)
        disp(['You clicked on row ', num2str(row), ', cluster ID: ', num2str(sort_index(row))]);
        cluster_summary(place_fields,sort_index(row),'within')
    end
end
end 