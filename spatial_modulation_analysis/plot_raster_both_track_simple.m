function plot_raster_both_track_simple(spike_times,spike_id,Task_info,Behaviour,subplot_xy,window,psthBinSize,varargin)

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
no_subplot_x = subplot_xy(1); %no of subplot in one figure columns
no_subplot_y = subplot_xy(2); %no of subplot in one figure rows

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

no_subplot = no_subplot_x*no_subplot_y;
track1_ID = find(Task_info.track_ID_all == 1);
track2_ID = find(Task_info.track_ID_all == 2);


% Define Gaussian window for smoothing
gaussianWindow = gausswin(5);

% Normalize to have an area of 1 (i.e., to be a probability distribution)
gaussianWindow = gaussianWindow / sum(gaussianWindow);


colour_lines = {[145,191,219]/255,[69,117,180]/255,[215,48,39]/255,[252,141,89]/255};
for iPlot = 1: ceil(no_cluster/(no_subplot))
    fig = figure;
    fig.Position = [109 77 1426 878];
    fig.Name = sprintf('Spatial modulation raster plot %s %i',unit_region(1+(iPlot-1)*no_subplot),iPlot);

    for iCluster = 1:no_subplot
        if iCluster+(iPlot-1)*no_subplot > no_cluster
            break
        end
        subplot_scale = 0:2;
        subplot(no_subplot_y*5,no_subplot_x,iCluster+(floor((iCluster-1)/no_subplot_x)+subplot_scale)*no_subplot_x+floor((iCluster-1)/no_subplot_x)*no_subplot)
        cluster_spike_id{iCluster+(iPlot-1)*no_subplot} = spike_id == unit_id(iCluster+(iPlot-1)*no_subplot);

        [~,~,rasterX,rasterY,~,~] = psthAndBA(spike_position(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),event_position, window, psthBinSize/20);
        rasterX_reshaped = reshape(rasterX,[3 length(rasterX)/3]);
        rasterY_reshaped = reshape(rasterY,[3 length(rasterY)/3]);
        track1_raster_id = ismember(rasterY_reshaped(1,:),track1_ID);
        track2_raster_id = ismember(rasterY_reshaped(1,:),track2_ID);
        track1_rasterX = reshape(rasterX_reshaped(:,track1_raster_id),[1 sum(track1_raster_id)*3]);
        track1_rasterY = reshape(rasterY_reshaped(:,track1_raster_id),[1 sum(track1_raster_id)*3]);
        track2_rasterX = reshape(rasterX_reshaped(:,track2_raster_id),[1 sum(track2_raster_id)*3]);
        track2_rasterY = reshape(rasterY_reshaped(:,track2_raster_id),[1 sum(track2_raster_id)*3]);
        plot(track1_rasterX,track1_rasterY,'LineWidth',0.5)


        hold on;
        plot(track2_rasterX,track2_rasterY,'LineWidth', 0.5)

        yline((find(diff(Task_info.track_ID_all)==-1))+1,'LineWidth',1,'Color',[0.5 0.5 0.5])
        yline((find(diff(Task_info.track_ID_all)==1))+1,'LineWidth',1,'Color',[0.5 0.5 0.5])
        xline(140, 'LineWidth',1,'Color',[0.5 0.5 0.5])
        ylim([0 length(Task_info.start_time_all)])
        xlim(window)
        ylabel('lap')
        xlabel('position')
        title(sprintf('%s unit %i at %i um',unit_region(iCluster+(iPlot-1)*no_subplot),unit_id(iCluster+(iPlot-1)*no_subplot),unit_depth(iCluster+(iPlot-1)*no_subplot)))

        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

        subplot_scale = 3:4;
        subplot(no_subplot_y*5,no_subplot_x,iCluster+(floor((iCluster-1)/no_subplot_x)+subplot_scale)*no_subplot_x+floor((iCluster-1)/no_subplot_x)*no_subplot)

        if isempty(place_fields)
            [psth_track1,bins,binnedArray1] = spatial_psth(spike_position(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),track1_event_position, window, psthBinSize,position_bin_time(Task_info.track_ID_all==1,:));
            [psth_track2,bins,binnedArray2] = spatial_psth(spike_position(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),track2_event_position, window, psthBinSize,position_bin_time(Task_info.track_ID_all==2,:));
            ratemaps_track1 = binnedArray1/psthBinSize;
            ratemaps_track2 = binnedArray2/psthBinSize;
        else
            ratemaps_track1 = place_fields(1).raw{place_fields(1).cluster_id==unit_id(iCluster+(iPlot-1)*no_subplot)};
            ratemaps_track2 = place_fields(2).raw{place_fields(2).cluster_id==unit_id(iCluster+(iPlot-1)*no_subplot)};
            bins = place_fields(1).x_bin_centres;
        end

        for nevent = 1:size(ratemaps_track1,1)
            ratemaps_track1(nevent,:) = conv(ratemaps_track1(nevent,:),gaussianWindow,'same');
        end

        for nevent = 1:size(ratemaps_track2,1)
            ratemaps_track2(nevent,:) = conv(ratemaps_track2(nevent,:),gaussianWindow,'same');
        end

        average_map_track1 = mean(ratemaps_track1,'omitnan');
%         average_map_track1_odd = mean(ratemaps_track1(1:2:end,:),'omitnan');
%         average_map_track1_even = mean(ratemaps_track1(2:2:end,:),'omitnan');

        average_map_track2 = mean(ratemaps_track2,'omitnan');
%         average_map_track2_odd = mean(ratemaps_track2(1:2:end,:),'omitnan');
%         average_map_track2_even = mean(ratemaps_track2(2:2:end,:),'omitnan');


        h(1)=plot(bins,average_map_track1,'LineWidth',2,'Color',colour_lines{1});
        map_error = std(ratemaps_track1)./sqrt(size(ratemaps_track1,1));
        hold on
        % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
        patch([bins fliplr(bins)],[average_map_track1+map_error fliplr(average_map_track1-map_error)],colour_lines{1},'FaceAlpha','0.3','LineStyle','none');

        [peak1,peak1_index] = max(average_map_track1(bins >= 46 & bins <= 86));
        [peak2,peak2_index] = max(average_map_track1(bins > 86 & bins <= 126));

        if  peak1 > peak2
            peak1_index = peak1_index+46/psthBinSize;
            peak2_index = peak1_index+40/psthBinSize;
            peak2 = average_map_track1(peak1_index+40/psthBinSize);

            SMI(1,iCluster+(iPlot-1)*no_subplot)=(peak1-peak2)/(peak1+peak2);
        elseif  peak1 < peak2
            peak2_index = peak2_index+86/psthBinSize;
            peak1_index = peak2_index-40/psthBinSize;
            peak1 = average_map_track1(peak2_index-40/psthBinSize);
            SMI(1,iCluster+(iPlot-1)*no_subplot)=(peak2-peak1)/(peak1+peak2);
        elseif  peak1 == peak2
            SMI(1,iCluster+(iPlot-1)*no_subplot)=0;
        end

        xline(bins(peak1_index),'LineWidth',2,'Color',colour_lines{1},'LineStyle','--')
        xline(bins(peak2_index),'LineWidth',2,'Color',colour_lines{1},'LineStyle','--')

        %         if ismember(unit_id(iCluster+(iPlot-1)*no_subplot),good_cell{1})
        %             xline(place_fields(1).x_bin_centres(peak1_index),'LineWidth',2,'Color',colour_lines{1},'LineStyle','--')
        %             xline(place_fields(1).x_bin_centres(peak2_index),'LineWidth',2,'Color',colour_lines{1},'LineStyle','--')
        %
        %         else
        %             SMI(1,iCluster+(iPlot-1)*no_subplot)=nan;
        %         end

        h(2)=plot(bins,average_map_track2,'LineWidth',2,'Color',colour_lines{3});
        map_error = std(ratemaps_track2)./sqrt(size(ratemaps_track2,1));
        hold on
        % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
        patch([bins fliplr(bins)],[average_map_track2_odd+map_error fliplr(average_map_track2_odd-map_error)],colour_lines{3},'FaceAlpha','0.3','LineStyle','none');

        [peak1,peak1_index] = max(average_map_track2(bins >= 46 & bins <= 86));
        [peak2,peak2_index] = max(average_map_track2(bins > 86 & bins <= 126));

        if  peak1 > peak2
            peak1_index = peak1_index+46/psthBinSize;
            peak2_index = peak1_index+40/psthBinSize;
            peak2 = average_map_track2(peak2_index);
            SMI(2,iCluster+(iPlot-1)*no_subplot)=(peak1-peak2)/(peak1+peak2);
        elseif  peak1 < peak2
            peak2_index = peak2_index+86/psthBinSize;
            peak1_index = peak2_index-40/psthBinSize;
            peak1 = average_map_track2(peak1_index);
            SMI(2,iCluster+(iPlot-1)*no_subplot)=(peak2-peak1)/(peak1+peak2);
        elseif  peak1 == peak2
            SMI(2,iCluster+(iPlot-1)*no_subplot)=0;
        end

        xline(bins(peak1_index),'LineWidth',2,'Color',colour_lines{3},'LineStyle','--')
        xline(bins(peak2_index),'LineWidth',2,'Color',colour_lines{3},'LineStyle','--')

% 
%         if ismember(unit_id(iCluster+(iPlot-1)*no_subplot),place_fields(2).cluster_id(place_fields(2).good_cells_LIBERAL))
%             xline(bins(peak1_index),'LineWidth',2,'Color',colour_lines{3},'LineStyle','--')
%             xline(bins(peak2_index),'LineWidth',2,'Color',colour_lines{3},'LineStyle','--')
%         else
%             SMI(2,iCluster+(iPlot-1)*no_subplot)=nan;
%         end


        xline([30 50 70 90 110],'LineWidth',1,'Color',[0.5 0.5 0.5])
        if iCluster == no_subplot
            legend([h(1:4)],{'Track 1 odd','Track 1 even','Track 2 odd','Track 2 even'},'Color','none')
        end
        xlim(window)
        xlabel('position')

        title(sprintf('SMI T1=%.3f T2=%.3f',SMI(1,iCluster+(iPlot-1)*no_subplot),SMI(2,iCluster+(iPlot-1)*no_subplot)))
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

    end
end