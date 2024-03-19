function plot_perievent_spiketimes(spike_times,spike_id,Task_info,Behaviour,subplot_xy,window,psthBinSize,varargin)

% Default values
p = inputParser;
addParameter(p,'event_times',Task_info.end_time_all,@isnumeric) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'event_label',[],@isstr) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'event_id',[],@isnumeric) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)


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

event_label = p.Results.event_label;
event_times = p.Results.event_times;
event_id =  p.Results.event_id;

if ~isempty(place_fields) % if place fields, only select spatially tuned cells
    spatial_cell_index = unique([find(place_fields(1).peak_percentile>0.95 & place_fields(1).odd_even_stability>0.95)...
        find(place_fields(2).peak_percentile>0.95 & place_fields(2).odd_even_stability>0.95)]);
    good_cell_index = intersect(spatial_cell_index,find(ismember(place_fields(1).cluster_id,unit_id)));

    [Lia,Locb] = ismember(unit_id,place_fields(1).cluster_id(good_cell_index));
    unit_depth = unit_depth(Lia);
    unit_region = unit_region(Lia);
    unit_id = unit_id(Lia);
end


no_subplot_x = subplot_xy(1); %no of subplot in one figure columns
no_subplot_y = subplot_xy(2); %no of subplot in one figure rows

t_bin = mean(diff(Behaviour.tvec));
no_events = size(event_times,1);
time_edges = window(1):psthBinSize:window(2);
% spike_speed = interp1(Behaviour.tvec,Behaviour.speed,spike_times,'nearest');
spike_times_events = spike_times;

for nevent = 1:no_events
    if nevent < no_events
        spike_times_event_index = spike_times < event_times(nevent+1)+window(1) ...
            & spike_times >=  (event_times(nevent)+window(1));
    else
        spike_times_event_index = spike_times >=  (event_times(nevent)+window(1));
    end
    spike_times_events(spike_times_event_index) = spike_times_events(spike_times_event_index)+100000*(nevent);
    event_times(nevent,1) = event_times(nevent,1)+(nevent)*100000;

end

if isempty(event_label)
    event_label = 'event';
    event_id = ones(length(event_times),1);
end

cluster_spike_id = cell(size(unit_id));
no_cluster = length(unit_id);

no_subplot = no_subplot_x*no_subplot_y;
track1_ID = find(event_id == 1);
track2_ID = find(event_id == 2);


% Define Gaussian window for smoothing
gaussianWindow = gausswin(0.2*1/psthBinSize);

% Normalize to have an area of 1 (i.e., to be a probability distribution)
gaussianWindow = gaussianWindow / sum(gaussianWindow);


colour_lines = {[145,191,219]/255,[69,117,180]/255,[215,48,39]/255,[252,141,89]/255};
for iPlot = 1: ceil(no_cluster/(no_subplot))
    fig = figure;
    fig.Position = [109 77 1426 878];
    fig.Name = sprintf('%s events raster plot %s %i',event_label,unit_region(1+(iPlot-1)*no_subplot),iPlot);

    for iCluster = 1:no_subplot
        if iCluster+(iPlot-1)*no_subplot > no_cluster
            break
        end
        subplot_scale = 0:2;
        subplot(no_subplot_y*5,no_subplot_x,iCluster+(floor((iCluster-1)/no_subplot_x)+subplot_scale)*no_subplot_x+floor((iCluster-1)/no_subplot_x)*no_subplot)
        cluster_spike_id{iCluster+(iPlot-1)*no_subplot} = spike_id == unit_id(iCluster+(iPlot-1)*no_subplot);

        [~,~,rasterX,rasterY,~,~] = psthAndBA(spike_times_events(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),event_times, window, 0.001);
        rasterX_reshaped = reshape(rasterX,[3 length(rasterX)/3]);
        rasterY_reshaped = reshape(rasterY,[3 length(rasterY)/3]);
        event1_raster_id = ismember(rasterY_reshaped(1,:), find(event_id == 1));
        event2_raster_id = ismember(rasterY_reshaped(1,:), find(event_id == 2));
        event1_rasterX = reshape(rasterX_reshaped(:,event1_raster_id),[1 sum(event1_raster_id)*3]);
        event1_rasterY = reshape(rasterY_reshaped(:,event1_raster_id),[1 sum(event1_raster_id)*3]);
        event2_rasterX = reshape(rasterX_reshaped(:,event2_raster_id),[1 sum(event2_raster_id)*3]);
        event2_rasterY = reshape(rasterY_reshaped(:,event2_raster_id),[1 sum(event2_raster_id)*3]);
        plot(event1_rasterX,event1_rasterY,'LineWidth',0.5)
        hold on;
        plot(event2_rasterX,event2_rasterY,'LineWidth', 0.5)

        xline(0, 'LineWidth',1,'Color',[0.5 0.5 0.5])
        ylim([0 length(event_id)])
        xlim(window)
        ylabel('Event')
        xlabel('Time')
        title(sprintf('%s unit %i at %i um',unit_region(iCluster+(iPlot-1)*no_subplot),unit_id(iCluster+(iPlot-1)*no_subplot),unit_depth(iCluster+(iPlot-1)*no_subplot)))

        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

        subplot_scale = 3:4;
        subplot(no_subplot_y*5,no_subplot_x,iCluster+(floor((iCluster-1)/no_subplot_x)+subplot_scale)*no_subplot_x+floor((iCluster-1)/no_subplot_x)*no_subplot)


        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray1] = psthAndBA(spike_times_events(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),event_times(event_id==1), window, psthBinSize);
        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray2] = psthAndBA(spike_times_events(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),event_times(event_id==2), window, psthBinSize);
        for nevent = 1:size(binnedArray1,1)
            psth_track1(nevent,:) = conv(binnedArray1(nevent,:)/psthBinSize,gaussianWindow,'same');
        end
        
        for nevent = 1:size(binnedArray2,1)
            psth_track2(nevent,:) = conv(binnedArray2(nevent,:)/psthBinSize,gaussianWindow,'same');
        end


        average_map_track1 = mean(psth_track1,'omitnan');
        average_map_track1_odd = mean(psth_track1(1:2:end,:),'omitnan');
        average_map_track1_even = mean(psth_track1(2:2:end,:),'omitnan');

        average_map_track2 = mean(psth_track2,'omitnan');
        average_map_track2_odd = mean(psth_track2(1:2:end,:),'omitnan');
        average_map_track2_even = mean(psth_track2(2:2:end,:),'omitnan');

        h(1)=plot(bins,average_map_track1_odd,'LineWidth',2,'Color',colour_lines{1});
        map_error = std(psth_track1(1:2:end,:))./sqrt(size(psth_track1(1:2:end,:),1));
        hold on
        % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
        patch([bins fliplr(bins)],[average_map_track1_odd+map_error fliplr(average_map_track1_odd-map_error)],colour_lines{1},'FaceAlpha','0.3','LineStyle','none');

        h(2)=plot(bins,average_map_track1_even,'LineWidth',2,'Color',colour_lines{2});
        map_error = std(psth_track1(2:2:end,:))./sqrt(size(psth_track1(2:2:end,:),1));
        hold on
        % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
        patch([bins fliplr(bins)],[average_map_track1_even+map_error fliplr(average_map_track1_even-map_error)],colour_lines{2},'FaceAlpha','0.3','LineStyle','none');

        h(3)=plot(bins,average_map_track2_odd,'LineWidth',2,'Color',colour_lines{3});
        map_error = std(psth_track2(1:2:end,:))./sqrt(size(psth_track2(1:2:end,:),1));
        hold on
        % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
        patch([bins fliplr(bins)],[average_map_track2_odd+map_error fliplr(average_map_track2_odd-map_error)],colour_lines{3},'FaceAlpha','0.3','LineStyle','none');

        h(4)=plot(bins,average_map_track2_even,'LineWidth',2,'Color',colour_lines{4});
        map_error = std(psth_track2(2:2:end,:))./sqrt(size(psth_track2(2:2:end,:),1));
        hold on
        % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
        patch([bins fliplr(bins)],[average_map_track2_even+map_error fliplr(average_map_track2_even-map_error)],colour_lines{4},'FaceAlpha','0.3','LineStyle','none');


        xline([0],'LineWidth',1,'Color',[0.5 0.5 0.5])
        if iCluster == no_subplot
            legend([h(1:4)],{'Track 1 odd','Track 1 even','Track 2 odd','Track 2 even'},'Color','none')
        end
        xlim(window)
        xlabel('Time')

        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

    end
end