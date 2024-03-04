function control_plots_phase_precession(track)

% Sanity check plots at different time points of theta phase precession extraction code
% Uses output from main code 'phase_precession_absolute_location', which is a structure called TPP

load('extracted_phase_precession_absolute_location.mat')
load('extracted_position.mat')
load('extracted_CSC.mat')
load('extracted_place_fields.mat')


track_linear = position.linear(track).linear(~isnan(position.linear(track).linear));
track_times = position.linear(track).timestamps;
track_x = position.x(position.linear(track).clean_track_Indices);
track_y = position.y(position.linear(track).clean_track_Indices);


% For each good cell in the track
for pc = 1 : length(place_fields.track(track).good_cells) % change here if you want to check just some cells - instead of looping through all the cells, add
                                                            % a variable before the loop with some random cell IDs from the good cells
        
    % sanity check that the correct segments and spikes are being chosen 
    % also compare with LFP
    figure
    
    subplot(211,'next','add')
    plot(track_times,track_linear)
    area(track_times,half_laps_times(track).direction_idx_1*230,'facecolor','k','facealpha',0.2,'LineStyle','none')
    area(track_times,half_laps_times(track).direction_idx_2*230,'facecolor','g','facealpha',0.2,'LineStyle','none')
    ylabel('Linearized position (cm)')
    raster_plot(TPP(track).place_cell_times_track{pc},100,'m',25) % all spikes
    raster_plot(TPP(track).place_cell_times_1{pc},130,'k',25) % spikes in direction 1 (only during running)
    raster_plot(TPP(track).place_cell_times_2{pc},60,'g',25) % spikes in direction 2 (only during running)
    text(max(track_times)+5,max(track_linear)+40,'spikes in direction 1','Color','k','FontSize',12,'FontWeight','bold')
    text(max(track_times)+5,max(track_linear)+20,'all spikes','Color','m','FontSize',12,'FontWeight','bold')
    text(max(track_times)+5,max(track_linear),'spikes in direction 2','Color','g','FontSize',12,'FontWeight','bold')
    title('Lap extraction with speed filter ( >5 cms^{-1}) and spike times')
    
    subplot(212,'next','add')
    plot(CSC(1).CSCtime,CSC(1).theta,'color',[0.2 0.2 0.2 0.4],'LineWidth',1.5)
    raster_plot(TPP(track).place_cell_times_1{pc},500,'k',300)
    raster_plot(TPP(track).place_cell_times_2{pc},500,'g',300)
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