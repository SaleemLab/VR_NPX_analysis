function plot_perievent_CSD_LFP_amplitude_phase(lfpAvg,csd,power,chan_config,sorted_config,best_channels,options)

if isfield(best_channels,'xcoord') %if different channels or depths for different columns/shanks
    col = find(best_channels.xcoord == lfpAvg.xcoord);
    fieldnames = fields(best_channels);  % Get the names of all fields in the structure
%     best_channels1 = best_channels;
    % Loop over all fields
    for i = 1:length(fieldnames)
        % Get the ith value from the current field
        if ~isnan(best_channels.(fieldnames{i})((col)))
            best_channels.(fieldnames{i}) = best_channels.(fieldnames{i})((col));
        else
            shank_id = ceil(best_channels.xcoord/250);
            best_channels.(fieldnames{i}) = median(best_channels.(fieldnames{i})(shank_id == shank_id(col)),'omitnan');
        end
    end
end

filter_type = lfpAvg.filter_type;

nprobe = options.probe_id+1;

for type = 1:length(lfpAvg.filter_type)
    for event = 1:length(lfpAvg.event_group)

        fig = figure
        fig.Position = [161 116 1700 830]
        if isfield(options,'CSD_V1_CA1_normalisation')
            fig.Name = sprintf('%s %s %s event (%s filtered CSD normalised) probe %i X coord %i',options.SUBJECT,options.SESSION,lfpAvg.event_group{event},filter_type{type},nprobe,unique(sorted_config.Ks_xcoord));
        else
            fig.Name = sprintf('%s %s %s event (%s filtered) probe %i X coord %i',options.SUBJECT,options.SESSION,lfpAvg.event_group{event},filter_type{type},nprobe,unique(sorted_config.Ks_xcoord));
        end
        subplot(1,5,1);
        colour_line= {[177,0,38]/256,[213,62,79]/256,[252,141,89]/256,...
            [153,213,148]/256,[26,152,80]/256,[66,146,198]/256,[8,69,148]/256};
        freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','300-600 Hz'};
        selected_frequency = [2,6,7];
        %     subplot(1,4,2)
        hold on
        plot([0 1],[best_channels.surface_depth best_channels.surface_depth],'--k','LineWidth',2)
        
        if isfield(best_channels,'L4_channel')
            if ~isempty(best_channels.L4_channel) % L4 is assumed to roughly 100 micron thick
                plot([0 1],[best_channels.L4_depth-50 best_channels.L4_depth-50],'--b','LineWidth',0.5)
                plot([0 1],[best_channels.L4_depth+50 best_channels.L4_depth+50],'--b','LineWidth',0.5)
                plot([0 1],[best_channels.L4_depth best_channels.L4_depth],'--b','LineWidth',2)
            end
        end

        if isfield(best_channels,'L5_channel')
            if ~isempty(best_channels.L5_channel)
                plot([0 1],[best_channels.L5_depth-100 best_channels.L5_depth-100],'--c','LineWidth',0.5)
                
                if ~isfield(best_channels,'L4_channel') | isempty(best_channels.L4_channel) | best_channels.L5_depth+100 < best_channels.L4_depth-60
                    % If L4 empty or upper bound of L5 does not overlap
                    % with lower bound of L4, L5 is assumed to be roughly
                    % 200 micron thick
                    plot([0 1],[best_channels.L5_depth+100 best_channels.L5_depth+100],'--c','LineWidth',0.5)
                else % otherwise L5 thickness limited by lower bound of L4
                    plot([0 1],[best_channels.L4_depth-60 best_channels.L4_depth-60],'--c','LineWidth',0.5)
                end
                
                plot([0 1],[best_channels.L5_depth best_channels.L5_depth],'--c','LineWidth',2)

            end
        end

        if isfield(best_channels,'CA1_channel')
            if ~isempty(best_channels.CA1_channel)
                plot([0 1],[best_channels.CA1_depth-100 best_channels.CA1_depth-100],'--r','LineWidth',0.5)
                plot([0 1],[best_channels.CA1_depth+100 best_channels.CA1_depth+100],'--r','LineWidth',0.5)
                plot([0 1],[best_channels.CA1_depth best_channels.CA1_depth],'--r','LineWidth',2)
            end
        end

        for n = selected_frequency
            pl(n) = plot(power(:,n)./max(power(:,n)),sorted_config.Ks_ycoord','Color',colour_line{n})

            %     p(n) = plot(power_differnece,1:96,colour_line{n})
            hold on
        end
        % legend('1-3','4-12','30-60','60-100','125-300')
        ylabel('probe depth (um)')
        legend([pl(selected_frequency)],{freq_legends{selected_frequency}},'Location','southeast','Color','none');
        ylim([0 max(chan_config.Ks_ycoord)*1.2])
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


        subplot(1,5,2)
        colororder({'k','r'})

        taxis = lfpAvg.(filter_type{type})(event).timestamps;
        for ch=1:size(lfpAvg.(filter_type{type})(event).filtered,2)

            %         sh_tmp = 1e0*(500000*lfpAvg.SO(event).raw(:,ch)) + sorted_config.Ks_ycoord(ch);
            sh_tmp = 1e0*(500000*lfpAvg.(filter_type{type})(event).filtered(:,ch)) + sorted_config.Ks_ycoord(ch);
            plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
            clear sh_tmp
        end
        % ylim([-1000 offset+1000]);
        xlim([taxis(1) taxis(end)]);
        xlabel('time (ms)');title('LFP');
        plot([0 0],ylim,'--r');hold on;

        title(sprintf('LFP (%s filtered)',filter_type{type}))

        hold on

        plot([min(taxis) max(taxis)],[best_channels.surface_depth best_channels.surface_depth],'--k','LineWidth',2)

         if isfield(best_channels,'L4_channel')
            if ~isempty(best_channels.L4_channel) % L4 is assumed to roughly 100 micron thick
                plot([min(taxis) max(taxis)],[best_channels.L4_depth-50 best_channels.L4_depth-50],'--b','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.L4_depth+50 best_channels.L4_depth+50],'--b','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.L4_depth best_channels.L4_depth],'--b','LineWidth',2)
            end
        end

        if isfield(best_channels,'L5_channel')
            if ~isempty(best_channels.L5_channel)
                plot([min(taxis) max(taxis)],[best_channels.L5_depth-100 best_channels.L5_depth-100],'--c','LineWidth',0.5)
                
                if ~isfield(best_channels,'L4_channel') | isempty(best_channels.L4_channel) | best_channels.L5_depth+100 < best_channels.L4_depth-60
                    % If L4 empty or upper bound of L5 does not overlap
                    % with lower bound of L4, L5 is assumed to be roughly
                    % 200 micron thick
                    plot([min(taxis) max(taxis)],[best_channels.L5_depth+100 best_channels.L5_depth+100],'--c','LineWidth',0.5)
                else % otherwise L5 thickness limited by lower bound of L4
                    plot([min(taxis) max(taxis)],[best_channels.L4_depth-60 best_channels.L4_depth-60],'--c','LineWidth',0.5)
                end
                
                plot([min(taxis) max(taxis)],[best_channels.L5_depth best_channels.L5_depth],'--c','LineWidth',2)

            end
        end

        if isfield(best_channels,'CA1_channel')
            if ~isempty(best_channels.CA1_channel)
                plot([min(taxis) max(taxis)],[best_channels.CA1_depth-100 best_channels.CA1_depth-100],'--r','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.CA1_depth+100 best_channels.CA1_depth+100],'--r','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.CA1_depth best_channels.CA1_depth],'--r','LineWidth',2)
            end
        end

        ylim([0 max(sorted_config.Ks_ycoord)*1.2])
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        yyaxis right
        ax = gca;
        ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'r';
        ylim([-best_channels.surface_depth max(chan_config.Ks_ycoord)*1.2-best_channels.surface_depth])


        subplot(1,5,3);
        taxis = csd.(filter_type{type})(event).timestamps;
        if isfield(options,'CSD_V1_CA1_normalisation')
            % Normalise CSD by V1 and CA1(and below)
%             cmax = max(max(csd.(filter_type{type})(event).data(:,:)));
%             csd.(filter_type{type})(event).data(:,:) = csd.(filter_type{type})(event).data(:,:)./cmax;

            cortex_channels = find(sorted_config.Ks_ycoord>best_channels.CA1_depth-200 & sorted_config.Ks_ycoord<best_channels.surface_depth);
            cmax = max(max(csd.(filter_type{type})(event).data(:,cortex_channels)));
            cmin = min(min(csd.(filter_type{type})(event).data(:,cortex_channels)));
            temp = csd.(filter_type{type})(event).data(:,cortex_channels);
      
            csd.(filter_type{type})(event).data(:,cortex_channels) = 2*(temp-cmin)./(cmax-cmin)-1;

            CA1_channels = find(sorted_config.Ks_ycoord>best_channels.CA1_depth-700 & sorted_config.Ks_ycoord<best_channels.CA1_depth-200);
            cmax = max(max(csd.(filter_type{type})(event).data(:,CA1_channels)));
            cmin = min(min(csd.(filter_type{type})(event).data(:,CA1_channels)));
            temp = csd.(filter_type{type})(event).data(:,CA1_channels);
            csd.(filter_type{type})(event).data(:,CA1_channels) = 2*(temp-cmin)./(cmax-cmin)-1;

            other_channels = find(sorted_config.Ks_ycoord>=0 & sorted_config.Ks_ycoord<best_channels.CA1_depth-700);
            cmax = max(max(csd.(filter_type{type})(event).data(:,other_channels)));
            cmin = min(min(csd.(filter_type{type})(event).data(:,other_channels)));
            temp = csd.(filter_type{type})(event).data(:,other_channels);
            csd.(filter_type{type})(event).data(:,other_channels) = 2*(temp-cmin)./(cmax-cmin)-1;

            cmax = max(max(csd.(filter_type{type})(event).data(:,1:end)));
            contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.(filter_type{type})(event).data,1)); csd.(filter_type{type})(event).data'; zeros(1,size(csd.(filter_type{type})(event).data,1))],40,'LineColor','none');hold on;
            colormap jet; caxis([-1 1]);
            xlabel('time (s)');title('CSD');
            plot([0 0],[1 size(csd.(filter_type{type})(event).data,2)],'--k');hold on;
            % set(gca,'YDir','reverse')
            % yticks(1:length(sorted_config.Ks_ycoord))
            % yticklabels(sorted_config.Ks_ycoord)
            plot([0 0],ylim,'--r');hold on;
            title(sprintf('Norm CSD (%s filtered)',filter_type{type}))
        else
            cmax = max(max(csd.(filter_type{type})(event).data(:,1:end)));
            contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.(filter_type{type})(event).data,1)); csd.(filter_type{type})(event).data'; zeros(1,size(csd.(filter_type{type})(event).data,1))],40,'LineColor','none');hold on;
            colormap jet; caxis([-cmax/2 cmax/2]);
            xlabel('time (s)');title('CSD');
            plot([0 0],[1 size(csd.(filter_type{type})(event).data,2)],'--k');hold on;
            % set(gca,'YDir','reverse')
            % yticks(1:length(sorted_config.Ks_ycoord))
            % yticklabels(sorted_config.Ks_ycoord)
            plot([0 0],ylim,'--r');hold on;
            title(sprintf('CSD (%s filtered)',filter_type{type}))
        end
        %

        colorbar
        hold on
        plot([min(taxis) max(taxis)],[best_channels.surface_depth best_channels.surface_depth],'--k','LineWidth',2)

         if isfield(best_channels,'L4_channel')
            if ~isempty(best_channels.L4_channel) % L4 is assumed to roughly 100 micron thick
                plot([min(taxis) max(taxis)],[best_channels.L4_depth-50 best_channels.L4_depth-50],'--b','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.L4_depth+50 best_channels.L4_depth+50],'--b','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.L4_depth best_channels.L4_depth],'--b','LineWidth',2)
            end
        end

        if isfield(best_channels,'L5_channel')
            if ~isempty(best_channels.L5_channel)
                plot([min(taxis) max(taxis)],[best_channels.L5_depth-100 best_channels.L5_depth-100],'--c','LineWidth',0.5)
                
                if ~isfield(best_channels,'L4_channel') | isempty(best_channels.L4_channel) | best_channels.L5_depth+100 < best_channels.L4_depth-60
                    % If L4 empty or upper bound of L5 does not overlap
                    % with lower bound of L4, L5 is assumed to be roughly
                    % 200 micron thick
                    plot([min(taxis) max(taxis)],[best_channels.L5_depth+100 best_channels.L5_depth+100],'--c','LineWidth',0.5)
                else % otherwise L5 thickness limited by lower bound of L4
                    plot([min(taxis) max(taxis)],[best_channels.L4_depth-60 best_channels.L4_depth-60],'--c','LineWidth',0.5)
                end
                
                plot([min(taxis) max(taxis)],[best_channels.L5_depth best_channels.L5_depth],'--c','LineWidth',2)

            end
        end

        if isfield(best_channels,'CA1_channel')
            if ~isempty(best_channels.CA1_channel)
                plot([min(taxis) max(taxis)],[best_channels.CA1_depth-100 best_channels.CA1_depth-100],'--r','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.CA1_depth+100 best_channels.CA1_depth+100],'--r','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.CA1_depth best_channels.CA1_depth],'--r','LineWidth',2)
            end
        end
        ylim([0 max(chan_config.Ks_ycoord)*1.2])

        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


        subplot(1,5,4);
        taxis = lfpAvg.(filter_type{type})(event).timestamps;
%         if isfield(options,'CSD_V1_CA1_normalisation')
%             % Normalise CSD by V1 and CA1(and below)
% 
%             cortex_channels = find(sorted_config.Ks_ycoord>best_channels.CA1_depth-200 & sorted_config.Ks_ycoord<best_channels.surface_depth);
%             cmax = max(max(lfpAvg.(filter_type{type})(event).power(:,cortex_channels)));
%             cmin = min(min(lfpAvg.(filter_type{type})(event).power(:,cortex_channels)));
%             temp = lfpAvg.(filter_type{type})(event).power(:,cortex_channels);
%             lfpAvg.(filter_type{type})(event).power(:,cortex_channels) = 2*(temp-cmin)./(cmax-cmin)-1;
% 
%             CA1_channels = find(sorted_config.Ks_ycoord>best_channels.CA1_depth-700 & sorted_config.Ks_ycoord<best_channels.CA1_depth-200);
%             cmax = max(max(lfpAvg.(filter_type{type})(event).power(:,CA1_channels)));
%             cmin = min(min(lfpAvg.(filter_type{type})(event).power(:,CA1_channels)));
%             temp = lfpAvg.(filter_type{type})(event).power(:,CA1_channels);
%             lfpAvg.(filter_type{type})(event).power(:,CA1_channels) = 2*(temp-cmin)./(cmax-cmin)-1;
% 
%             other_channels = find(sorted_config.Ks_ycoord>=0 & sorted_config.Ks_ycoord<best_channels.CA1_depth-700);
%             cmax = max(max(lfpAvg.(filter_type{type})(event).power(:,other_channels)));
%             cmin = min(min(lfpAvg.(filter_type{type})(event).power(:,other_channels)));
%             temp = lfpAvg.(filter_type{type})(event).power(:,other_channels);
%             lfpAvg.(filter_type{type})(event).power(:,other_channels) = 2*(temp-cmin)./(cmax-cmin)-1;
% 
% %             cmax = max(max(lfpAvg.(filter_type{type})(event).power(:,1:end)));
%             contourf(taxis,sorted_config.Ks_ycoord,[lfpAvg.(filter_type{type})(event).power'],40,'LineColor','none');hold on;
%             colormap jet; caxis([-1/2 1/2]);
%             xlabel('time (s)');title('CSD');
%             plot([0 0],[1 size(lfpAvg.(filter_type{type})(event).power,2)],'--k');hold on;
%             % set(gca,'YDir','reverse')
%             % yticks(1:length(sorted_config.Ks_ycoord))
%             % yticklabels(sorted_config.Ks_ycoord)
%             plot([0 0],ylim,'--r');hold on;
%             title(sprintf('Norm power (%s filtered)',filter_type{type}))
%         else
            cmax = max(max(lfpAvg.(filter_type{type})(event).power(:,1:end)));
            contourf(taxis,sorted_config.Ks_ycoord,[lfpAvg.(filter_type{type})(event).power'],40,'LineColor','none');hold on;
            colormap jet; caxis([-cmax/2 cmax/2]);
            xlabel('time (s)');title('CSD');
            plot([0 0],[1 size(lfpAvg.(filter_type{type})(event).power,2)],'--k');hold on;
            % set(gca,'YDir','reverse')
            % yticks(1:length(sorted_config.Ks_ycoord))
            % yticklabels(sorted_config.Ks_ycoord)
            plot([0 0],ylim,'--r');hold on;
            title(sprintf('LFP power (%s filtered)',filter_type{type}))
%         end
        %

        colorbar
        hold on
        plot([min(taxis) max(taxis)],[best_channels.surface_depth best_channels.surface_depth],'--k','LineWidth',2)

         if isfield(best_channels,'L4_channel')
            if ~isempty(best_channels.L4_channel) % L4 is assumed to roughly 100 micron thick
                plot([min(taxis) max(taxis)],[best_channels.L4_depth-50 best_channels.L4_depth-50],'--b','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.L4_depth+50 best_channels.L4_depth+50],'--b','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.L4_depth best_channels.L4_depth],'--b','LineWidth',2)
            end
        end

        if isfield(best_channels,'L5_channel')
            if ~isempty(best_channels.L5_channel)
                plot([min(taxis) max(taxis)],[best_channels.L5_depth-100 best_channels.L5_depth-100],'--c','LineWidth',0.5)
                
                if ~isfield(best_channels,'L4_channel') | isempty(best_channels.L4_channel) | best_channels.L5_depth+100 < best_channels.L4_depth-60
                    % If L4 empty or upper bound of L5 does not overlap
                    % with lower bound of L4, L5 is assumed to be roughly
                    % 200 micron thick
                    plot([min(taxis) max(taxis)],[best_channels.L5_depth+100 best_channels.L5_depth+100],'--c','LineWidth',0.5)
                else % otherwise L5 thickness limited by lower bound of L4
                    plot([min(taxis) max(taxis)],[best_channels.L4_depth-60 best_channels.L4_depth-60],'--c','LineWidth',0.5)
                end
                
                plot([min(taxis) max(taxis)],[best_channels.L5_depth best_channels.L5_depth],'--c','LineWidth',2)

            end
        end

        if isfield(best_channels,'CA1_channel')
            if ~isempty(best_channels.CA1_channel)
                plot([min(taxis) max(taxis)],[best_channels.CA1_depth-100 best_channels.CA1_depth-100],'--r','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.CA1_depth+100 best_channels.CA1_depth+100],'--r','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.CA1_depth best_channels.CA1_depth],'--r','LineWidth',2)
            end
        end
        ylim([0 max(chan_config.Ks_ycoord)*1.2])

        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


        subplot(1,5,5);
        taxis = lfpAvg.(filter_type{type})(event).timestamps;
%         if isfield(options,'CSD_V1_CA1_normalisation')
%             % Normalise CSD by V1 and CA1(and below)
%             cortex_channels = find(sorted_config.Ks_ycoord>best_channels.CA1_depth-200 & sorted_config.Ks_ycoord<best_channels.surface_depth);
%             cmax = max(max(lfpAvg.(filter_type{type})(event).phase(:,cortex_channels)));
%             cmin = min(min(lfpAvg.(filter_type{type})(event).phase(:,cortex_channels)));
%             temp = lfpAvg.(filter_type{type})(event).power(:,cortex_channels);
%             lfpAvg.(filter_type{type})(event).phase(:,cortex_channels) = 2*(temp-cmin)./(cmax-cmin)-1;
% 
%             CA1_channels = find(sorted_config.Ks_ycoord>best_channels.CA1_depth-700 & sorted_config.Ks_ycoord<best_channels.CA1_depth-200);
%             cmax = max(max(lfpAvg.(filter_type{type})(event).phase(:,CA1_channels)));
%             cmin = min(min(lfpAvg.(filter_type{type})(event).phase(:,CA1_channels)));
%             temp = lfpAvg.(filter_type{type})(event).phase(:,CA1_channels);
%             lfpAvg.(filter_type{type})(event).phase(:,CA1_channels) = 2*(temp-cmin)./(cmax-cmin)-1;
% 
%             other_channels = find(sorted_config.Ks_ycoord>=0 & sorted_config.Ks_ycoord<best_channels.CA1_depth-700);
%             cmax = max(max(lfpAvg.(filter_type{type})(event).phase(:,other_channels)));
%             cmin = min(min(lfpAvg.(filter_type{type})(event).phase(:,other_channels)));
%             temp = lfpAvg.(filter_type{type})(event).phase(:,other_channels);
%             lfpAvg.(filter_type{type})(event).phase(:,other_channels) = 2*(temp-cmin)./(cmax-cmin)-1;
% 
% %             cmax = max(max(lfpAvg.(filter_type{type})(event).phase(:,1:end)));
%             contourf(taxis,sorted_config.Ks_ycoord,[lfpAvg.(filter_type{type})(event).phase'],40,'LineColor','none');hold on;
%             colormap jet; caxis([-1/2 1/2]);
%             xlabel('time (s)');title('CSD');
%             plot([0 0],[1 size(lfpAvg.(filter_type{type})(event).phase,2)],'--k');hold on;
%             % set(gca,'YDir','reverse')
%             % yticks(1:length(sorted_config.Ks_ycoord))
%             % yticklabels(sorted_config.Ks_ycoord)
%             plot([0 0],ylim,'--r');hold on;
%             title(sprintf('Norm phase (%s filtered)',filter_type{type}))
%         else
            cmax = max(max(lfpAvg.(filter_type{type})(event).phase(:,1:end)));
            contourf(taxis,sorted_config.Ks_ycoord,[lfpAvg.(filter_type{type})(event).phase'],40,'LineColor','none');hold on;
            colormap jet; caxis([-cmax/2 cmax/2]);
            xlabel('time (s)');title('CSD');
            plot([0 0],[1 size(lfpAvg.(filter_type{type})(event).phase,2)],'--k');hold on;
            % set(gca,'YDir','reverse')
            % yticks(1:length(sorted_config.Ks_ycoord))
            % yticklabels(sorted_config.Ks_ycoord)
            plot([0 0],ylim,'--r');hold on;
            title(sprintf('LFP phase (%s filtered)',filter_type{type}))
%         end
        %

        colorbar
        hold on
        plot([min(taxis) max(taxis)],[best_channels.surface_depth best_channels.surface_depth],'--k','LineWidth',2)

         if isfield(best_channels,'L4_channel')
            if ~isempty(best_channels.L4_channel) % L4 is assumed to roughly 100 micron thick
                plot([min(taxis) max(taxis)],[best_channels.L4_depth-50 best_channels.L4_depth-50],'--b','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.L4_depth+50 best_channels.L4_depth+50],'--b','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.L4_depth best_channels.L4_depth],'--b','LineWidth',2)
            end
        end

        if isfield(best_channels,'L5_channel')
            if ~isempty(best_channels.L5_channel)
                plot([min(taxis) max(taxis)],[best_channels.L5_depth-100 best_channels.L5_depth-100],'--c','LineWidth',0.5)
                
                if ~isfield(best_channels,'L4_channel') | isempty(best_channels.L4_channel) | best_channels.L5_depth+100 < best_channels.L4_depth-60
                    % If L4 empty or upper bound of L5 does not overlap
                    % with lower bound of L4, L5 is assumed to be roughly
                    % 200 micron thick
                    plot([min(taxis) max(taxis)],[best_channels.L5_depth+100 best_channels.L5_depth+100],'--c','LineWidth',0.5)
                else % otherwise L5 thickness limited by lower bound of L4
                    plot([min(taxis) max(taxis)],[best_channels.L4_depth-60 best_channels.L4_depth-60],'--c','LineWidth',0.5)
                end
                
                plot([min(taxis) max(taxis)],[best_channels.L5_depth best_channels.L5_depth],'--c','LineWidth',2)

            end
        end

        if isfield(best_channels,'CA1_channel')
            if ~isempty(best_channels.CA1_channel)
                plot([min(taxis) max(taxis)],[best_channels.CA1_depth-100 best_channels.CA1_depth-100],'--r','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.CA1_depth+100 best_channels.CA1_depth+100],'--r','LineWidth',0.5)
                plot([min(taxis) max(taxis)],[best_channels.CA1_depth best_channels.CA1_depth],'--r','LineWidth',2)
            end
        end
        ylim([0 max(chan_config.Ks_ycoord)*1.2])

        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    end
end
