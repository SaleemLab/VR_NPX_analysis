function SMI = calculate_spatial_modulation_index(cluster,Task_info,Behaviour,window,psthBinSize,varargin)

% Default values
p = inputParser;
addParameter(p,'place_fields',[],@isstruct) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn))
addParameter(p,'plot_option',0,@isnumeric) % Powers Selected frequency for plotting
addParameter(p,'subplot_xy',[],@isnumeric) % Powers Selected frequency for plotting

% assign parameters (either defaults or given)
parse(p,varargin{:});
place_fields = p.Results.place_fields;
plot_option = p.Results.plot_option;
subplot_xy = p.Results.subplot_xy;

% calculate_spatial_modulation_index(V1_clusters_L,place_fields_V1_L,Task_info,Behaviour,[5 1],[0 140],5,1)
no_subplot_x = subplot_xy(1); %no of subplot in one figure columns
no_subplot_y = subplot_xy(2); %no of subplot in one figure rows

t_bin = 1/60;
position_edges = window(1):psthBinSize:window(2);
%   convert spikes time to corresponding postions
spike_position = interp1(Behaviour.tvec,Behaviour.position,cluster.spike_times,'nearest');
spike_speed = interp1(Behaviour.tvec,Behaviour.speed,cluster.spike_times,'nearest');
no_lap = size(Task_info.start_time_all,1);
event_position = zeros(size(Task_info.start_time_all));

position_bin_time = zeros(no_lap,(window(2)-window(1))/psthBinSize);
for iLap = 1:no_lap
    spike_times_lap_index = cluster.spike_times <= Task_info.end_time_all(iLap)...
        & cluster.spike_times >= Task_info.start_time_all(iLap) & spike_speed > 5;

    spike_position(spike_times_lap_index) = spike_position(spike_times_lap_index)+1000*(iLap);
    event_position(iLap,1) = (iLap)*1000;
    position_bin_time(iLap,:) = t_bin.*histcounts(Behaviour.position(Behaviour.tvec>=Task_info.start_time_all(iLap) ...
        & Behaviour.tvec <=Task_info.end_time_all(iLap) & Behaviour.speed > 5 ),position_edges);
end


track1_event_position = event_position(Task_info.track_ID_all==1);
track2_event_position = event_position(Task_info.track_ID_all==2);


cluster_spike_id = cell(size(cluster.cluster_id));
no_cluster = length(cluster.cluster_id);

if plot_option == 1
    no_subplot = no_subplot_x*no_subplot_y;
    track1_ID = find(Task_info.track_ID_all == 1);
    track2_ID = find(Task_info.track_ID_all == 2);
    for iPlot = 1: ceil(no_cluster/(no_subplot))
        fig = figure;
        fig.Position = [109 77 1426 878];
        for iCluster = 1:no_subplot
            if iCluster+(iPlot-1)*no_subplot > no_cluster
                break
            end
            subplot_scale = 0:3;
            subplot(no_subplot_y*5,no_subplot_x,iCluster+(floor((iCluster-1)/no_subplot_x)+subplot_scale)*no_subplot_x+floor((iCluster-1)/no_subplot_x)*no_subplot)
            cluster_spike_id{iCluster+(iPlot-1)*no_subplot} = cluster.spike_id == cluster.cluster_id(iCluster+(iPlot-1)*no_subplot);

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

            yline(find(diff(Task_info.track_ID_all)==-1)+1,'LineWidth',1,'Color',[0.5 0.5 0.5])
            yline(find(diff(Task_info.track_ID_all)==1)+1,'LineWidth',1,'Color',[0.5 0.5 0.5])

            ylim([0 length(Task_info.start_time_all)])
            xlim(window)
            ylabel('lap')
            xlabel('position')
            title(['unit: ',num2str(iCluster+(iPlot-1)*no_subplot),'at depth ',num2str(cluster.peak_depth(iCluster+(iPlot-1)*no_subplot))])
            set(gca,'TickDir','out','box','off','Color','none','FontSize',12)


            subplot(no_subplot_y*5,no_subplot_x,iCluster+(floor((iCluster-1)/no_subplot_x)+4)*no_subplot_x+floor((iCluster-1)/no_subplot_x)*no_subplot)
            
            %         [psth_track1,bins] = spatial_psth(spike_position(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),track1_event_position, window, psthBinSize,position_bin_time(Task_info.track_ID_all==1,:));
            %         [psth_track2,bins] = spatial_psth(spike_position(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),track2_event_position, window, psthBinSize,position_bin_time(Task_info.track_ID_all==2,:));

            plot(bins,psth_track1,'LineWidth',2)

            hold on;
            plot(bins,psth_track2,'LineWidth',2)
            xline([30 50 70 90 110],'LineWidth',1,'Color',[0.5 0.5 0.5])
            if iCluster == no_subplot
                legend(['track 1'],['track 2'],'bkgd','boxoff','Color','none')
            end
            xlim(window)
            xlabel('position')

            title('PSTH')
            set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        end
    end

else

end

