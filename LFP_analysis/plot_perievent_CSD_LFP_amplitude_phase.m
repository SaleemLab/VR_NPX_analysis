function plot_perievent_CSD_LFP_amplitude_phase(lfpAvg,csd,power,chan_config,sorted_config,best_channels)
filter_type = lfpAvg.filter_type;
for type = 1:length(lfpAvg.filter_type)
    for event = 1:length(lfpAvg.event_group)
        fig = figure
        fig.Position = [334.7143 102.7143 1600 830]
        fig.Name = sprintf('%s event (%s filtered)',lfpAvg.event_group{event},filter_type{type});

        subplot(1,5,1);
        colour_line= {[177,0,38]/256,[213,62,79]/256,[252,141,89]/256,...
            [153,213,148]/256,[26,152,80]/256,[66,146,198]/256,[8,69,148]/256};
        freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','300-600 Hz'};
        selected_frequency = [1 2 6 7];
        %     subplot(1,4,2)
        hold on
        plot([0 1],[chan_config.Ks_ycoord(best_channels.first_in_brain_channel) chan_config.Ks_ycoord(best_channels.first_in_brain_channel)],'--k','LineWidth',2)
        plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel))],'--b','LineWidth',2)
        plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L5_channel))],'--c','LineWidth',2)
        plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel))],'--r','LineWidth',2)

        for n = selected_frequency
            pl(n) = plot(power(:,n)./max(power(:,n)),sorted_config.Ks_ycoord','Color',colour_line{n})

            %     p(n) = plot(power_differnece,1:96,colour_line{n})
            hold on
        end
        % legend('1-3','4-12','30-60','60-100','125-300')
        legend([pl(selected_frequency)],{freq_legends{selected_frequency}});
        ylim([0 4000])

        subplot(1,5,2)
%         fig = figure
%         fig.Position = [334.7143 102.7143 340 420]
        taxis = lfpAvg.(filter_type{type})(event).timestamps;
        for ch=1:size(lfpAvg.(filter_type{type})(event).filtered,2)

            %         sh_tmp = 1e0*(500000*lfpAvg.SO(event).raw(:,ch)) + sorted_config.Ks_ycoord(ch);
            sh_tmp = 1e0*(50000*lfpAvg.(filter_type{type})(event).filtered(:,ch)) + sorted_config.Ks_ycoord(ch);
            plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
            clear sh_tmp
        end
        % ylim([-1000 offset+1000]);
        xlim([taxis(1) taxis(end)]);
        xlabel('time (ms)');ylabel('depth(um)');title('LFP');
        plot([0 0],ylim,'--r');hold on;

        title(sprintf('LFP (%s filtered)',filter_type{type}))
        hold on
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(best_channels.first_in_brain_channel) chan_config.Ks_ycoord(best_channels.first_in_brain_channel)],'--k','LineWidth',2)
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel))],'--b','LineWidth',2)
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L5_channel))],'--c','LineWidth',2)
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel))],'--r','LineWidth',2)
        ylim([0 4000])
%         ylim([10*(best_channels.CA1_channel-15) 10*(best_channels.CA1_channel+15)])
%         xlim([-100 200])
%         set(gca,"TickDir","out",'box', 'off','Color','none','Fontsize',14)
% save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project\ripples',[])

        subplot(1,5,3);
        taxis = csd.(filter_type{type})(event).timestamps;
        cmax = max(max(csd.(filter_type{type})(event).data));
        contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.(filter_type{type})(event).data,1)); csd.(filter_type{type})(event).data'; zeros(1,size(csd.(filter_type{type})(event).data,1))],40,'LineColor','none');hold on;
        colormap jet; caxis([-cmax/2 cmax/2]);
        xlabel('time (s)');ylabel('channel');title('CSD');
        plot([0 0],[1 size(csd.(filter_type{type})(event).data,2)],'--k');hold on;
        % set(gca,'YDir','reverse')
        % yticks(1:length(sorted_config.Ks_ycoord))
        % yticklabels(sorted_config.Ks_ycoord)
        plot([0 0],ylim,'--r');hold on;
        title(sprintf('CSD (%s filtered)',filter_type{type}))
        colorbar
        hold on
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(best_channels.first_in_brain_channel) chan_config.Ks_ycoord(best_channels.first_in_brain_channel)],'--k','LineWidth',2)
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel))],'--b','LineWidth',2)
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L5_channel))],'--c','LineWidth',2)
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel))],'--r','LineWidth',2)
        ylim([0 4000])

        subplot(1,5,4);
        taxis = lfpAvg.(filter_type{type})(event).timestamps;
        cmax = max(max(lfpAvg.(filter_type{type})(event).power));
        % imagesc(taxis,sorted_config.Ks_ycoord,[lfpAvg_SO(1).power']);hold on;
        contourf(taxis,sorted_config.Ks_ycoord,[lfpAvg.(filter_type{type})(event).power'],40,'LineColor','none');hold on;
        colormap jet; caxis([-cmax/2 cmax/2]);
        xlabel('time (s)');ylabel('channel');title('CSD');
        plot([0 0],[1 size(lfpAvg.(filter_type{type})(event).power,2)],'--k');hold on;
        % set(gca,'YDir','reverse')
        % yticks(1:length(sorted_config.Ks_ycoord))
        % yticklabels(sorted_config.Ks_ycoord)
        plot([0 0],ylim,'--r');hold on;
        title(sprintf('LFP Power (%s filtered)',filter_type{type}))
        colorbar
        hold on
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(best_channels.first_in_brain_channel) chan_config.Ks_ycoord(best_channels.first_in_brain_channel)],'--k','LineWidth',2)
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel))],'--b','LineWidth',2)
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L5_channel))],'--c','LineWidth',2)
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel))],'--r','LineWidth',2)
        ylim([0 4000])

        subplot(1,5,5);
        taxis = lfpAvg.(filter_type{type})(event).timestamps;
        cmax = max(max(lfpAvg.(filter_type{type})(event).phase));
        %     imagesc(taxis,sorted_config.Ks_ycoord,[lfpAvg.SO(event).power']);hold on;
        contourf(taxis,sorted_config.Ks_ycoord,[lfpAvg.(filter_type{type})(event).phase'],40,'LineColor','none');hold on;
        colormap jet; caxis([-cmax cmax]);
        xlabel('time (s)');ylabel('channel');title('CSD');
        plot([0 0],[1 size(lfpAvg.(filter_type{type})(event).phase,2)],'--k');hold on;
        % set(gca,'YDir','reverse')
        % yticks(1:length(sorted_config.Ks_ycoord))
        % yticklabels(sorted_config.Ks_ycoord)
        plot([0 0],ylim,'--r');hold on;
        title(sprintf('LFP phase (%s filtered)',filter_type{type}))
        colorbar
        hold on
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(best_channels.first_in_brain_channel) chan_config.Ks_ycoord(best_channels.first_in_brain_channel)],'--k','LineWidth',2)
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel))],'--b','LineWidth',2)
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L5_channel))],'--c','LineWidth',2)
        plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel))],'--r','LineWidth',2)
        ylim([0 4000])
        
        set(gca,"TickDir","out",'box', 'off','Color','none','Fontsize',14)

        sgtitle(sprintf('%s event (%s filtered)',lfpAvg.event_group{event},filter_type{type}))
%         filename = sprintf('%s event (%s filtered).fig',lfpAvg.event_group{event},filter_type{type})
%         saveas(gcf,filename)
%         close
        
    end
end