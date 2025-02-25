function plot_spatial_heatmap_both_track(spike_times,spike_id,Task_info,Behaviour,window,psthBinSize,varargin)

% Default values
p = inputParser;
addParameter(p,'place_fields',[],@isstruct) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'unit_depth',[],@isnumeric) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'unit_region',[],@isstring) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'unit_id',unique(spike_id),@isnumeric) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)

% assign parameters (either defaults or given)
parse(p,varargin{:});
place_fields = p.Results.place_fields;
unit_depth = p.Results.unit_depth;
unit_region = p.Results.unit_region;
unit_id = p.Results.unit_id;

if ~isempty(place_fields) % if place fields, only select spatially tuned cells
    spatial_cell_index = unique([find(place_fields(1).peak_percentile>0.95 & place_fields(1).odd_even_stability>0.95)...
        find(place_fields(2).peak_percentile>0.95 & place_fields(2).odd_even_stability>0.95)]);
    good_cell_index = intersect(spatial_cell_index,find(ismember(place_fields(1).cluster_id,unit_id)));

    [Lia,Locb] = ismember(unit_id,place_fields(1).cluster_id(good_cell_index));
    unit_depth = unit_depth(Lia);
    unit_region = unit_region(Lia);
    unit_id = unit_id(Lia);
end

t_bin = mean(diff(Behaviour.tvec));
no_lap = size(Task_info.start_time_all,1);
position_edges = window(1):psthBinSize:window(2);
%   convert spikes time to corresponding postions
spike_position = interp1(Behaviour.tvec,Behaviour.position,spike_times,'nearest');
spike_speed = interp1(Behaviour.tvec,Behaviour.speed,spike_times,'nearest');

event_position = zeros(size(Task_info.start_time_all));

position_bin_time = zeros(no_lap,(window(2)-window(1))/psthBinSize);
for iLap = 1:no_lap
%     if iLap < no_lap
%         spike_times_lap_index = spike_times < Task_info.start_time_all(iLap+1)...
%             & spike_times >= Task_info.start_time_all(iLap) & spike_speed > 5;
%         position_bin_time(iLap,:) = t_bin.*histcounts(Behaviour.position(Behaviour.tvec>=Task_info.start_time_all(iLap) ...
%             & Behaviour.tvec <Task_info.start_time_all(iLap+1) & Behaviour.speed > 5 ),position_edges);
%     else
        spike_times_lap_index = spike_times <= Task_info.end_time_all(iLap)...
            & spike_times >= Task_info.start_time_all(iLap) & spike_speed > 5;
        position_bin_time(iLap,:) = t_bin.*histcounts(Behaviour.position(Behaviour.tvec>=Task_info.start_time_all(iLap) ...
            & Behaviour.tvec <=Task_info.end_time_all(iLap) & Behaviour.speed > 5 ),position_edges);
%     end
    spike_position(spike_times_lap_index) = spike_position(spike_times_lap_index)+1000*(iLap);
    event_position(iLap,1) = (iLap)*1000;

end

track1_event_position = event_position(Task_info.track_ID_all==1);
track2_event_position = event_position(Task_info.track_ID_all==2);

cluster_spike_id = cell(size(unit_id));
no_cluster = length(unit_id);

track1_ID = find(Task_info.track_ID_all == 1);
track2_ID = find(Task_info.track_ID_all == 2);


% Define Gaussian window for smoothing
gaussianWindow = gausswin(5);

% Normalize to have an area of 1 (i.e., to be a probability distribution)
gaussianWindow = gaussianWindow / sum(gaussianWindow);

colour_lines = {[145,191,219]/255,[69,117,180]/255,[215,48,39]/255,[252,141,89]/255};



for iCluster = 1:no_cluster
    fig = figure;
    fig.Position = [109 135 1000 600];
    fig.Name = sprintf('%s unit %i at %i um',unit_region(iCluster),unit_id(iCluster),unit_depth(iCluster));

    subplot(1,3,1)
    cluster_spike_id{iCluster} = spike_id == unit_id(iCluster);

    [~,~,rasterX,rasterY,~,~] = psthAndBA(spike_position(cluster_spike_id{iCluster}),event_position, window, psthBinSize/20);
    rasterX_reshaped = reshape(rasterX,[3 length(rasterX)/3]);
    rasterY_reshaped = reshape(rasterY,[3 length(rasterY)/3]);
    track1_raster_id = ismember(rasterY_reshaped(1,:),track1_ID);
    track2_raster_id = ismember(rasterY_reshaped(1,:),track2_ID);
    track1_rasterX = reshape(rasterX_reshaped(:,track1_raster_id),[1 sum(track1_raster_id)*3]);
    track1_rasterY = reshape(rasterY_reshaped(:,track1_raster_id),[1 sum(track1_raster_id)*3]);
    track2_rasterX = reshape(rasterX_reshaped(:,track2_raster_id),[1 sum(track2_raster_id)*3]);
    track2_rasterY = reshape(rasterY_reshaped(:,track2_raster_id),[1 sum(track2_raster_id)*3]);
    plot(track1_rasterX,track1_rasterY,'LineWidth',0.5,'Color',colour_lines{2})


    hold on;
    plot(track2_rasterX,track2_rasterY,'LineWidth', 0.5,'Color',colour_lines{3})

    yline((find(diff(Task_info.track_ID_all)==-1))+1,'LineWidth',1,'Color',[0.5 0.5 0.5])
    yline((find(diff(Task_info.track_ID_all)==1))+1,'LineWidth',1,'Color',[0.5 0.5 0.5])
    xline(140, 'LineWidth',1,'Color',[0.5 0.5 0.5])
    ylim([0 length(Task_info.start_time_all)])
    xlim(window)
    ylabel('lap')
    xlabel('position')
    %     title(sprintf('%s unit %i at %i um',unit_region(iCluster),unit_id(iCluster),unit_depth(iCluster)))
    xticks([0 30 50 70 90 110])
%     xticklabels([0 20 40 60 80 100 120 140])
    set(gca,'TickDir','out','box','off','Color','none','FontSize',12)



    [psth_track1,bins,binnedArray1] = spatial_psth(spike_position(cluster_spike_id{iCluster}),track1_event_position, window, psthBinSize,position_bin_time(Task_info.track_ID_all==1,:));
    [psth_track2,bins,binnedArray2] = spatial_psth(spike_position(cluster_spike_id{iCluster}),track2_event_position, window, psthBinSize,position_bin_time(Task_info.track_ID_all==2,:));
    ratemaps_track1 = binnedArray1/psthBinSize;
    ratemaps_track2 = binnedArray2/psthBinSize;

    for nevent = 1:size(ratemaps_track1,1)
        ratemaps_track1(nevent,:) = conv(ratemaps_track1(nevent,:),gaussianWindow,'same');
    end

    for nevent = 1:size(ratemaps_track2,1)
        ratemaps_track2(nevent,:) = conv(ratemaps_track2(nevent,:),gaussianWindow,'same');
    end

    average_map_track1 = mean(ratemaps_track1,'omitnan');
    %     average_map_track1_odd = mean(ratemaps_track1(1:2:end,:),'omitnan');
    %     average_map_track1_even = mean(ratemaps_track1(2:2:end,:),'omitnan');

    average_map_track2 = mean(ratemaps_track2,'omitnan');
    %     average_map_track2_odd = mean(ratemaps_track2(1:2:end,:),'omitnan');
    %     average_map_track2_even = mean(ratemaps_track2(2:2:end,:),'omitnan');


    subplot(3,3,[2])
    h(1)=plot(bins,average_map_track1,'LineWidth',2,'Color',colour_lines{2});
    map_error = std(ratemaps_track1(:,:),'omitnan')./sqrt(size(ratemaps_track1(:,:),1));
    hold on
    % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
    patch([bins fliplr(bins)],[average_map_track1+map_error fliplr(average_map_track1-map_error)],colour_lines{2},'FaceAlpha','0.3','LineStyle','none');

    h(2)=plot(bins,average_map_track2,'LineWidth',2,'Color',colour_lines{3});
    map_error = std(ratemaps_track1(:,:),'omitnan')./sqrt(size(ratemaps_track1(:,:),1));
    hold on
    % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
    patch([bins fliplr(bins)],[average_map_track2+map_error fliplr(average_map_track2-map_error)],colour_lines{3},'FaceAlpha','0.3','LineStyle','none');
    xticks([0 30 50 70 90 110])
    set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

%     xline([30 50 70 90 110],'LineWidth',1,'Color',[0.5 0.5 0.5])
    legend([h(1:2)],{'Track L','Track R'},'Color','none','box','off','Units', 'normalized','Position', [0.60 0.85 0.1 0.1])

    xlim(window)
    xlabel('position (cm)')
    ylabel('firing rate (Hz)')

    reshape([ratemaps_track1; ratemaps_track2])))

    subplot(6,3,[11 14 17])
    imagesc(ratemaps_track1)
    clim([0 )
    set(gca, 'YDir', 'normal')
    colormap(flipud(gray))
    colorbar
    set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
    ylabel('lap')
    xlabel('position')
    xticks([0 30 50 70 90 110]/psthBinSize)
    xticklabels([0 30 50 70 90 110])

    title('Track L')
    subplot(6,3,[12 15 18])
    imagesc(ratemaps_track2)
    set(gca, 'YDir', 'normal')
    clim([0 (max(max([ratemaps_track1; ratemaps_track2])))])
    colormap(flipud(gray))
    colorbar
    ylabel('lap')
    xlabel('position')
    xticks([0 30 50 70 90 110]/psthBinSize)
    xticklabels([0 30 50 70 90 110])
    set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

    title('Track R')
    sgt = sgtitle(sprintf('%s unit %i at %i um',unit_region(iCluster),unit_id(iCluster),unit_depth(iCluster)))
    sgt.FontSize = 15;

end


