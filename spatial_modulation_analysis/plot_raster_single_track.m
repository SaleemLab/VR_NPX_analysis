function plot_raster_single_track(cluster,Task_info,Behaviour,subplot_xy,window,psthBinSize)
no_subplot_x = subplot_xy(1); %no of subplot in one figure columns
no_subplot_y = subplot_xy(2); %no of subplot in one figure rows

    t_bin = 1/60;
    position_edges = window(1):psthBinSize:window(2);
    color = ['b';'r'];
    
for track_id = 1:2
    spike_position = interp1(Behaviour.tvec,Behaviour.position,cluster.spike_times,'nearest');
    spike_speed = interp1(Behaviour.tvec,Behaviour.speed,cluster.spike_times,'nearest');
    %   convert spikes time to corresponding postions
    no_lap = size(Task_info.start_time_all(Task_info.track_ID_all == track_id),1);
    event_position = zeros([no_lap,1]);
    track_start_time{track_id} = Task_info.start_time_all(Task_info.track_ID_all == track_id);
    track_end_time{track_id} = Task_info.end_time_all(Task_info.track_ID_all == track_id);
    position_bin_time = zeros(no_lap,(window(2)-window(1))/psthBinSize);
    for iLap = 1:size(Task_info.start_time_all(Task_info.track_ID_all == track_id),1)
        spike_times_lap_index = cluster.spike_times <= track_end_time{track_id}(iLap)...
            & cluster.spike_times >= track_start_time{track_id}(iLap) & spike_speed > 5;
        spike_position(spike_times_lap_index) = spike_position(spike_times_lap_index)+1000*iLap;
        event_position(iLap,1) = iLap*1000;
        position_bin_time(iLap,:) = t_bin.*histcounts(Behaviour.position(Behaviour.tvec>=track_start_time{track_id}(iLap) ...
            & Behaviour.tvec <=track_end_time{track_id}(iLap) & Behaviour.speed > 5 ),position_edges);
    end


    cluster_spike_id = cell(size(cluster.cluster_id));
    no_cluster = length(cluster.cluster_id);

    no_subplot = no_subplot_x*no_subplot_y;
    for iPlot = 1: ceil(no_cluster/(no_subplot))
        fig = figure(iPlot);
        fig.Position = [109 77 1426 878];
        
        for iCluster = 1:no_subplot
            if iCluster+(iPlot-1)*no_subplot > no_cluster
                break
            end
            subplot(no_subplot_y*5,no_subplot_x,[iCluster+(floor((iCluster-1)/no_subplot_x)-2+track_id*2)*no_subplot_x,iCluster+(floor((iCluster-1)/no_subplot_x)-1+track_id*2)*no_subplot_x]);
            cluster_spike_id{iCluster+(iPlot-1)*no_subplot} = cluster.spike_id == cluster.cluster_id(iCluster+(iPlot-1)*no_subplot);
            [psth,bin,rasterX,rasterY,~,~] = psthAndBA(spike_position(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),event_position, window, psthBinSize/20);
            plot(rasterX,rasterY,'LineWidth',1,'Color',color(track_id,:))
            yline(find(diff(Task_info.track_ID_all)==(-1)^track_id)+1,'LineWidth',1,'Color',[0.5 0.5 0.5])

            ylim([0 length(track_start_time{track_id})])
            xlim(window)
            ylabel('lap')
            xlabel('position')
            title(['unit: ',num2str(iCluster+(iPlot-1)*no_subplot),' in Track: ',num2str(track_id)])
            set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

            subplot(no_subplot_y*5,no_subplot_x,iCluster+(floor((iCluster-1)/no_subplot_x)-1+5)*no_subplot_x)
            [psth,bins] = spatial_psth(spike_position(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),event_position, window, psthBinSize,position_bin_time);
            hold on;
            
            plot(bins,psth,'LineWidth',2,'Color',color(track_id,:))
            xline([ 30 50 70 90 110],'LineWidth',1,'Color',[0.5 0.5 0.5])
            if iCluster == no_subplot
            legend(['track 1'],['track 2'],'bkgd','boxoff','Color','none')
            end
            xlim(window)
            xlabel('position')
            
            title('PSTH')
            set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        end
    end
end
end