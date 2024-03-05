function control_plots_phase_precession(LFP_tvec,CA1_LFP,place_fields,clusters,TPP,Task_info,Behaviour,options)

% Sanity check plots at different time points of theta phase precession extraction code
% Uses output from main code 'phase_precession_absolute_location', which is a structure called TPP
if isfield(clusters,'merged_spike_id')
    clusters.cluster_id = unique(clusters.merged_cluster_id);
    for ncell = 1:length(unique(clusters.merged_cluster_id))
        tempt_peak_channel = clusters.peak_channel(clusters.merged_cluster_id == clusters.cluster_id(ncell));
        tempt_peak_depth = clusters.peak_depth(clusters.merged_cluster_id == clusters.cluster_id(ncell));
        tempt_peak_waveform = clusters.peak_channel_waveforms(clusters.merged_cluster_id == clusters.cluster_id(ncell),:);
        tempt_cell_type = clusters.cell_type(clusters.merged_cluster_id == clusters.cluster_id(ncell));

        if length(tempt_peak_depth)== 2
            tempt_peak_channel = tempt_peak_channel(1);
            tempt_peak_depth = tempt_peak_depth(1);
            tempt_peak_waveform = tempt_peak_waveform(1,:);
            tempt_cell_type = tempt_cell_type(1);
        else % find median peak depth assign that value to the unit
            [~,index]= min(tempt_peak_depth - median(tempt_peak_depth));
            tempt_peak_channel = tempt_peak_channel(index);
            tempt_peak_depth = tempt_peak_depth(index);
            tempt_peak_waveform = tempt_peak_waveform(index,:);
            tempt_cell_type = tempt_cell_type(index);
        end

        merged_peak_channel(ncell) = tempt_peak_channel;
        merged_peak_depth(ncell) = tempt_peak_depth;
        merged_peak_waveform(ncell,:) = tempt_peak_waveform;
        merged_cell_type(ncell,:) = tempt_cell_type;
    end

    clusters.peak_channel = merged_peak_channel;
    clusters.peak_depth = merged_peak_depth;
    clusters.peak_channel_waveforms = merged_peak_waveform;
    clusters.cell_type = merged_cell_type;

    clusters.spike_id = clusters.merged_spike_id;
    clusters.cluster_id = unique(clusters.merged_cluster_id);
end

spatial_cell_index = unique([find(place_fields(1).peak_percentile>0.95 & place_fields(1).odd_even_stability>0.95)...
    find(place_fields(2).peak_percentile>0.95 & place_fields(2).odd_even_stability>0.95)]);
good_cell_index = intersect(spatial_cell_index,find(ismember(place_fields(1).cluster_id,clusters.cluster_id)));

for track = 1:length(place_fields)
    track_linear = Behaviour.position;
    track_linear(Behaviour.track_ID~=track) = nan;

    track_times = Behaviour.tvec;

    % For each good cell in the track
    for ncell = 1 : length(good_cell_index) % change here if you want to check just some cells - instead of looping through all the cells, add
        % a variable before the loop with some random cell IDs from the good cells

        % sanity check that the correct segments and spikes are being chosen
        % also compare with LFP
        figure

        subplot(211,'next','add')
        plot(track_times,track_linear)
        ylabel('Linearized position (cm)')
        raster_plot(TPP(track).place_cell_times_track{ncell},100,'r',25) % all spikes
        text(max(track_times)+5,max(track_linear)+20,'all spikes','Color','r','FontSize',12,'FontWeight','bold')
        title('Lap extraction with speed filter ( >5 cms^{-1}) and spike times')

        subplot(212,'next','add')
        plot(LFP_tvec,CA1_LFP,'color',[0.2 0.2 0.2 0.4],'LineWidth',1.5)
        raster_plot(TPP(track).place_cell_times{ncell},500,'k',300)
        xlabel('Time (s)')
        ylabel('Amplitude (V)')
        set(findall(gcf,'-property','FontSize'),'FontSize',16)


        %% plot the phase precession over position - plot 2 cycles for clarity
        figure
        subplot(121)
        plot([TPP(track).spike_positions_1{pc},TPP(track).spike_positions_1{pc}],[rad2deg(TPP(track).spike_phases_1{pc})+180,rad2deg(TPP(track).spike_phases_1{pc})+180+360],...
            'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3,'LineStyle','none')
        xlabel('Linearized position')
        ylabel('Phase (degrees)')
        title('Direction 1')
        subplot(122)
        plot([TPP(track).spike_positions_2{pc},TPP(track).spike_positions_2{pc}],[rad2deg(TPP(track).spike_phases_2{pc})+180,rad2deg(TPP(track).spike_phases_2{pc})+180+360],...
            'Marker','o','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',3,'LineStyle','none')
        xlabel('Linearized position')
        ylabel('Phase (degrees)')
        title('Direction 2')
        set(findall(gcf,'-property','FontSize'),'FontSize',16)

        % create colour gradient for 360 degrees
        n=360;
        red = (1/ceil(n/2)):(1/ceil(n/2)):1;
        nSones = ones(ceil(n/2),1);
        red = [red' ; nSones];
        blue = flipud(red);
        green = zeros(1,size(red,1));
        cols = [red green' blue];

        %% plot spikes over animal's trajectory, with colour representing phase
        figure
        subplot(121,'next','add')
        plot(track_x,track_y,'o','MarkerSize',2,...
            'MarkerEdgeColor','none','MarkerFaceColor',[0.3 0.3 0.3])
        spikes_x = TPP(track).spike_x_1{pc};
        spikes_y = TPP(track).spike_y_1{pc};
        spikes_pase = TPP(track).spike_phases_1{pc};
        for n=1:length(spikes_pase)
            plot(spikes_x(n),spikes_y(n),'o','MarkerSize',5,...
                'MarkerEdgeColor','none','MarkerFaceColor',cols(ceil(rad2deg(spikes_pase(n)))+180,:))
        end
        title('Direction 1')

        subplot(122,'next','add')
        plot(track_x,track_y,'o','MarkerSize',2,...
            'MarkerEdgeColor','none','MarkerFaceColor',[0.3 0.3 0.3])
        spikes_x = TPP(track).spike_x_2{pc};
        spikes_y = TPP(track).spike_y_2{pc};
        spike_phases_2 = TPP(track).spike_phases_2{pc};
        for n=1:length(spike_phases_2)
            plot(spikes_x(n),spikes_y(n),'o','MarkerSize',5,...
                'MarkerEdgeColor','none','MarkerFaceColor',cols(ceil(rad2deg(spike_phases_2(n)))+180,:))
        end
        title('Direction 2')
        c = colorbar;
        cmap = colormap(cols);
        c.Limits = [0 1];
        c.Label.String  = 'Phase';
        c.Ticks = [0,1];
        c.TickLabels = {'0','360'};
        set(findall(gcf,'-property','FontSize'),'FontSize',16)




        %% plot correlation values and phase precession for a specific cell

        spike_id = 13; %index = 2 of a good cell for DIR1
        %index = 1 of a bad cell for DIR1

        figure

        subplot(521,'next','add')
        histogram([TPP(track).circ_lin_PVAL_dir1(:)],'FaceColor','k','BinWidth',0.01)
        plot([0.05 0.05],ylim,'r--')
        title('direction 1')
        xlim([0 1])

        subplot(522,'next','add')
        histogram([TPP(track).circ_lin_PVAL_dir2(:)],'FaceColor','g','BinWidth',0.01)
        plot([0.05 0.05],ylim,'r--')
        title('direction 2')
        xlim([0 1])

        subplot(5,2,[3,5],'next','add')
        plot([TPP(track).circ_lin_PVAL_dir1(:)],[TPP(track).circ_lin_corr_dir1(:)],'k.','MarkerSize',10)
        plot(TPP(track).circ_lin_PVAL_dir1(spike_id),...
            TPP(track).circ_lin_corr_dir1(spike_id),'ro','MarkerSize',10)
        plot([0.05 0.05],[0 1],'r--')
        xlabel('p value')
        ylabel('correlation coefficient')

        subplot(5,2,[4,6],'next','add')
        plot([TPP(track).circ_lin_PVAL_dir2(:)],[TPP(track).circ_lin_corr_dir2(:)],'g.','MarkerSize',10)
        plot(TPP(track).circ_lin_PVAL_dir2(spike_id),...
            TPP(track).circ_lin_corr_dir2(spike_id),'ro','MarkerSize',10)
        plot([0.05 0.05],[0 1],'r--')
        xlabel('p value')
        ylabel('correlation coefficient')

        subplot(5,2,[7,9],'next','add')
        if ~isnan(TPP(track).circ_lin_corr_dir1(spike_id))
            plot([TPP(track).spike_positions_1{spike_id},TPP(track).spike_positions_1{spike_id}]...
                ,[rad2deg(TPP(track).spike_phases_1{spike_id})+180,rad2deg(TPP(track).spike_phases_1{spike_id})+180+360],...
                'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3,'LineStyle','none')
        end
        xlabel('Linearized position')
        ylabel('Phase (degrees)')

        subplot(5,2,[8,10],'next','add')
        if ~isnan(TPP(track).circ_lin_corr_dir2(spike_id))
            plot([TPP(track).spike_positions_2{spike_id},TPP(track).spike_positions_2{spike_id}]...
                ,[rad2deg(TPP(track).spike_phases_2{spike_id})+180,rad2deg(TPP(track).spike_phases_2{spike_id})+180+360],...
                'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3,'LineStyle','none')
        end
        xlabel('Linearized position')
        ylabel('Phase (degrees)')
        set(findall(gcf,'-property','FontSize'),'FontSize',16)





    end

end

end



function raster_plot(spike_times,y,col,height)

x2(1:3:length(spike_times)*3)=spike_times;
x2(2:3:length(spike_times)*3)=spike_times;
x2(3:3:length(spike_times)*3)=NaN;
y2(1:3:length(spike_times)*3)=y;
y2(2:3:length(spike_times)*3)=y+height;
y2(3:3:length(spike_times)*3)=NaN;
if isempty(col)
    plot(x2,y2,'linewidth',2);
else
    plot(x2,y2,'color',col,'linewidth',2);
end
end