function plot_raster_single_track(cluster,Task_info,Behaviour,subplot_xy)
no_subplot_x = subplot_xy(1); %no of subplot in one figure columns
no_subplot_y = subplot_xy(2); %no of subplot in one figure rows

for track_id = 1:2
    %   convert spikes time to corresponding postions
    spike_position = interp1(Behaviour.tvec,Behaviour.position,cluster.spike_times,'nearest');

    event_position = zeros(size(Task_info.start_time_all));
    track_start_time{track_id} = Task_info.start_time_all(Task_info.track_ID_all == track_id);
    track_end_time{track_id} = Task_info.end_time_all(Task_info.track_ID_all == track_id)
    for iLap = 1:size(Task_info.start_time_all(Task_info.track_ID_all == track_id),1)
        spike_times_lap_index = cluster.spike_times <= track_end_time{track_id}(iLap)...
            & cluster.spike_times >= track_start_time{track_id}(iLap);
        spike_position(spike_times_lap_index) = spike_position(spike_times_lap_index)+1000*(iLap-1);
        event_position(iLap,1) = (iLap-1)*1000;
    end
    window = [0,140];
    psthBinSize = 1;

    cluster_spike_id = cell(size(cluster.cluster_id));
    no_cluster = length(cluster.cluster_id);

    no_subplot = no_subplot_x*no_subplot_y;
    for iPlot = 1: ceil(no_cluster/(no_subplot))
        figure(iPlot);
        for iCluster = 1:no_subplot
            if iCluster+(iPlot-1)*no_subplot > no_cluster
                break
            end
            subplot(no_subplot_y*2,no_subplot_x,iCluster+(floor((iCluster-1)/no_subplot_x)-1+track_id)*no_subplot_x);
            cluster_spike_id{iCluster+(iPlot-1)*no_subplot} = cluster.spike_id == cluster.cluster_id(iCluster+(iPlot-1)*no_subplot);

            [~,~,rasterX,rasterY,~,~] = psthAndBA(spike_position(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),event_position, window, psthBinSize);
            plot(rasterX,rasterY,'LineWidth',1.5)
            ylim([0 length(track_start_time{track_id})])
            xlim(window)
            ylabel('lap')
            xlabel('position')
            title(['unit: ',num2str(iCluster+(iPlot-1)*no_subplot),' in Track: ',num2str(track_id)])
            set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        end
    end
end
end