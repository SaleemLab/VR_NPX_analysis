function plot_raster_both_track(cluster,Task_info,Behaviour,subplot_xy)
no_subplot_x = subplot_xy(1); %no of subplot in one figure columns
no_subplot_y = subplot_xy(2); %no of subplot in one figure rows


    %   convert spikes time to corresponding postions
    spike_position = interp1(Behaviour.tvec,Behaviour.position,cluster.spike_times,'nearest');
    spike_speed = interp1(Behaviour.tvec,Behaviour.speed,cluster.spike_times,'nearest');
    event_position = zeros(size(Task_info.start_time_all));
    for iLap = 1:size(Task_info.start_time_all,1)
        spike_times_lap_index = cluster.spike_times <= Task_info.end_time_all(iLap)...
            & cluster.spike_times >= Task_info.start_time_all(iLap) & spike_speed > 5;
        
        spike_position(spike_times_lap_index) = spike_position(spike_times_lap_index)+1000*(iLap);
        event_position(iLap,1) = (iLap)*1000;


        
    end
    track1_event_position = event_position(Task_info.track_ID_all==1);
    track2_event_position = event_position(Task_info.track_ID_all==2);
    window = [0,140];
    psthBinSize = 1;
    
    cluster_spike_id = cell(size(cluster.cluster_id));
    no_cluster = length(cluster.cluster_id);
    
    no_subplot = no_subplot_x*no_subplot_y;
    track1_ID = find(Task_info.track_ID_all == 1);
    track2_ID = find(Task_info.track_ID_all == 2);
    for iPlot = 1: ceil(no_cluster/(no_subplot))
        figure;
    for iCluster = 1:no_subplot
        if iCluster+(iPlot-1)*no_subplot > no_cluster
            break
        end
        subplot_scale = 0:3; 
        subplot(no_subplot_y*5,no_subplot_x,iCluster+(floor((iCluster-1)/no_subplot_x)+subplot_scale)*no_subplot_x+floor((iCluster-1)/no_subplot_x)*no_subplot)
        cluster_spike_id{iCluster+(iPlot-1)*no_subplot} = cluster.spike_id == cluster.cluster_id(iCluster+(iPlot-1)*no_subplot);
    
        [~,~,rasterX,rasterY,~,~] = psthAndBA(spike_position(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),event_position, window, psthBinSize);
        rasterX_reshaped = reshape(rasterX,[3 length(rasterX)/3]);
        rasterY_reshaped = reshape(rasterY,[3 length(rasterY)/3]);
        track1_raster_id = ismember(rasterY_reshaped(1,:),track1_ID);
        track2_raster_id = ismember(rasterY_reshaped(1,:),track2_ID);
        track1_rasterX = reshape(rasterX_reshaped(:,track1_raster_id),[1 sum(track1_raster_id)*3]);
        track1_rasterY = reshape(rasterY_reshaped(:,track1_raster_id),[1 sum(track1_raster_id)*3]);
        track2_rasterX = reshape(rasterX_reshaped(:,track2_raster_id),[1 sum(track2_raster_id)*3]);
        track2_rasterY = reshape(rasterY_reshaped(:,track2_raster_id),[1 sum(track2_raster_id)*3]);
        plot(track1_rasterX,track1_rasterY,'LineWidth',1)
        hold on; 
        plot(track2_rasterX,track2_rasterY,'LineWidth', 1)

            ylim([0 length(Task_info.start_time_all)])
            xlim(window)
            ylabel('lap')
            xlabel('position')
            title(['unit: ',num2str(iCluster+(iPlot-1)*no_subplot),'at depth ',num2str(cluster.peak_depth(iCluster+(iPlot-1)*no_subplot))])
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)


        subplot(no_subplot_y*5,no_subplot_x,iCluster+(floor((iCluster-1)/no_subplot_x)+4)*no_subplot_x+floor((iCluster-1)/no_subplot_x)*no_subplot)
        [psth_track1,bins,~,~,~,binnedarray_track1] = psthAndBA(spike_position(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),track1_event_position, window, psthBinSize);
        [psth_track2,bins,~,~,~,binnedarray_track2] = psthAndBA(spike_position(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),track2_event_position, window, psthBinSize);
        plot(bins,mean(binnedarray_track1),'LineWidth',2)
        hold on;
        plot(bins,mean(binnedarray_track2),'LineWidth',2)
        legend(['track 1'],['track 2'])
        xlim(window)
        xlabel('position')
        
        title('PSTH')
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
    end
end