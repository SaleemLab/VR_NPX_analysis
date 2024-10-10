function plot_cluster_density_profile(power,chan_config,sorted_config,best_channels,clusters,options)

if isfield(best_channels,'xcoord') %if different channels or depths for different columns/shanks
    col = find(best_channels.xcoord == unique(sorted_config.Ks_xcoord));
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



if isfield(options,'probe_no')
    nprobe= options.probe_no;
else
    nprobe = 1;
end


fig = figure
fig.Position = [334.7143 102.7143 556 830]

if isfield(options,'spatial_cell_id')
    fig.Name = sprintf('%s %s spatial cluster density probe %i X coord %i',options.SUBJECT,options.SESSION,nprobe,unique(sorted_config.Ks_xcoord));
    sgtitle(sprintf('%s %s spatial cluster density probe %i X coord %i',options.SUBJECT,options.SESSION,nprobe,unique(sorted_config.Ks_xcoord)))
else
    fig.Name = sprintf('%s %s cluster density probe %i X coord %i',options.SUBJECT,options.SESSION,nprobe,unique(sorted_config.Ks_xcoord));
    sgtitle(sprintf('%s %s cluster density probe %i X coord %i',options.SUBJECT,options.SESSION,nprobe,unique(sorted_config.Ks_xcoord)))
end

subplot(1,2,1);
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


subplot(1,2,2);
colour_line= {[177,0,38]/256,[213,62,79]/256,[252,141,89]/256,...
    [153,213,148]/256,[26,152,80]/256,[66,146,198]/256,[8,69,148]/256};
freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','300-600 Hz'};
selected_frequency = [2,6,7];
%     subplot(1,4,2)

% select good units based on simple metrics (remove obvious bad clusters)
if isfield(clusters,'isi_viol')
    good_unit_index = find(clusters.isi_viol < 0.1 ...
        & clusters.amplitude_cutoff < 0.1 ...
        & clusters.amplitude > 50);

elseif isfield(clusters,'sliding_rp_violation')

    if isfield(clusters,'merged_spike_id')

        clusters.cluster_id = unique(clusters.merged_cluster_id);
        for ncell = 1:length(unique(clusters.merged_cluster_id))
            tempt_peak_channel = clusters.peak_channel(clusters.merged_cluster_id == clusters.cluster_id(ncell));
            tempt_peak_depth = clusters.peak_depth(clusters.merged_cluster_id == clusters.cluster_id(ncell));
            tempt_peak_waveform = clusters.peak_channel_waveforms(clusters.merged_cluster_id == clusters.cluster_id(ncell),:);
            tempt_cell_type = clusters.cell_type(clusters.merged_cluster_id == clusters.cluster_id(ncell));
            tempt_sliding_rp_violation = clusters.sliding_rp_violation(clusters.merged_cluster_id == clusters.cluster_id(ncell));
            tempt_amplitude = clusters.amplitude_median(clusters.merged_cluster_id == clusters.cluster_id(ncell));
            tempt_amplitude_cutoff = clusters.amplitude_cutoff(clusters.merged_cluster_id == clusters.cluster_id(ncell));

           
            if length(tempt_peak_depth)== 2
                tempt_peak_channel = tempt_peak_channel(1);
                tempt_peak_depth = tempt_peak_depth(1);
                tempt_peak_waveform = tempt_peak_waveform(1,:);
                tempt_cell_type = tempt_cell_type(1);
                tempt_sliding_rp_violation = tempt_sliding_rp_violation(1);
                tempt_amplitude = tempt_amplitude(1);
                tempt_amplitude_cutoff = tempt_amplitude_cutoff(1);

            else % find median peak depth assign that value to the unit
                [~,index]= min(tempt_peak_depth - median(tempt_peak_depth));
                tempt_peak_channel = tempt_peak_channel(index);
                tempt_peak_depth = tempt_peak_depth(index);
                tempt_peak_waveform = tempt_peak_waveform(index,:);
                tempt_cell_type = tempt_cell_type(index);
                tempt_sliding_rp_violation = tempt_sliding_rp_violation(index);
                tempt_amplitude = tempt_amplitude(index);
                tempt_amplitude_cutoff = tempt_amplitude_cutoff(index);
            end

            merged_peak_channel(ncell) = tempt_peak_channel;
            merged_peak_depth(ncell) = tempt_peak_depth;
            merged_peak_waveform(ncell,:) = tempt_peak_waveform;
            merged_cell_type(ncell) = tempt_cell_type;
            merged_sliding_rp_violation(ncell) = tempt_sliding_rp_violation;
            merged_amplitude(ncell) = tempt_amplitude;
            merged_amplitude_cutoff(ncell) = tempt_amplitude_cutoff;
        end

        clusters.peak_channel = merged_peak_channel;
        clusters.peak_depth = merged_peak_depth;
        clusters.peak_channel_waveforms = merged_peak_waveform;
        clusters.cell_type = merged_cell_type;
        clusters.sliding_rp_violation = merged_sliding_rp_violation;
        clusters.amplitude_median = merged_amplitude;
        clusters.amplitude_cutoff = merged_amplitude_cutoff;

        clusters.spike_id = clusters.merged_spike_id;
        clusters.cluster_id = unique(clusters.merged_cluster_id);
    end



    if isfield(options,'spatial_cell_id')
        good_unit_index = find(ismember(clusters(nprobe).cluster_id,options.spatial_cell_id));

    else
        good_unit_index = find(clusters.sliding_rp_violation < 0.1 ...
            & clusters.amplitude_cutoff < 0.1 ...
            & clusters.amplitude_median > 50);
    end

end

bin_edge = sorted_config.Ks_ycoord(1)-(sorted_config.Ks_ycoord(2)-sorted_config.Ks_ycoord(1))/2:...
    mean(diff(sorted_config.Ks_ycoord)):...
    max(sorted_config.Ks_ycoord)+(sorted_config.Ks_ycoord(end)-sorted_config.Ks_ycoord(end-1))/2;

% Find density of units for each cell type

% intersect(good_unit_index,find(clusters.cell_type==1)
pyramidal_density = histcounts(clusters.peak_depth(intersect(good_unit_index,find(clusters.cell_type==1))),bin_edge);
narrow_interneuron_density = histcounts(clusters.peak_depth(intersect(good_unit_index,find(clusters.cell_type==2))),bin_edge);
wide_interneuron_density = histcounts(clusters.peak_depth(intersect(good_unit_index,find(clusters.cell_type==3))),bin_edge);
bin_max = 1.2*max([pyramidal_density narrow_interneuron_density wide_interneuron_density]);

hold on
plot([0 bin_max],[best_channels.surface_depth best_channels.surface_depth],'--k','LineWidth',2)

if isfield(best_channels,'L4_channel')
    if ~isempty(best_channels.L4_channel) % L4 is assumed to roughly 100 micron thick
        plot([0 bin_max],[best_channels.L4_depth-50 best_channels.L4_depth-50],'--b','LineWidth',0.5)
        plot([0 bin_max],[best_channels.L4_depth+50 best_channels.L4_depth+50],'--b','LineWidth',0.5)
        plot([0 bin_max],[best_channels.L4_depth best_channels.L4_depth],'--b','LineWidth',2)
    end
end

if isfield(best_channels,'L5_channel')
    if ~isempty(best_channels.L5_channel)
        plot([0 bin_max],[best_channels.L5_depth-100 best_channels.L5_depth-100],'--c','LineWidth',0.5)

        if ~isfield(best_channels,'L4_channel') | isempty(best_channels.L4_channel) | best_channels.L5_depth+100 < best_channels.L4_depth-60
            % If L4 empty or upper bound of L5 does not overlap
            % with lower bound of L4, L5 is assumed to be roughly
            % 200 micron thick
            plot([0 bin_max],[best_channels.L5_depth+100 best_channels.L5_depth+100],'--c','LineWidth',0.5)
        else % otherwise L5 thickness limited by lower bound of L4
            plot([0 bin_max],[best_channels.L4_depth-60 best_channels.L4_depth-60],'--c','LineWidth',0.5)
        end

        plot([0 bin_max],[best_channels.L5_depth best_channels.L5_depth],'--c','LineWidth',2)

    end
end

if isfield(best_channels,'CA1_channel')
    if ~isempty(best_channels.CA1_channel)
        plot([0 bin_max],[best_channels.CA1_depth-100 best_channels.CA1_depth-100],'--r','LineWidth',0.5)
        plot([0 bin_max],[best_channels.CA1_depth+100 best_channels.CA1_depth+100],'--r','LineWidth',0.5)
        plot([0 bin_max],[best_channels.CA1_depth best_channels.CA1_depth],'--r','LineWidth',2)
    end
end



hold on
b(1) = barh(sorted_config.Ks_ycoord',pyramidal_density,'r','EdgeColor','none','FaceAlpha',0.3);
b(2) = barh(sorted_config.Ks_ycoord',narrow_interneuron_density,'b','EdgeColor','none','FaceAlpha',0.3);
b(3) = barh(sorted_config.Ks_ycoord',wide_interneuron_density,'k','EdgeColor','none','FaceAlpha',0.3);
legend(b(:),{'pyramidal','narrow','wide'},'Location','southeast','Color','none')
% legend('1-3','4-12','30-60','60-100','125-300')
ylabel('probe depth (um)')
legend([pl(selected_frequency)],{freq_legends{selected_frequency}},'Location','southeast','Color','none');
ylim([0 max(chan_config.Ks_ycoord)*1.2])
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

%         sgtitle(sprintf('%s %s %s event (%s filtered) probe %i',options.SUBJECT,options.SESSION,lfpAvg.event_group{event},filter_type{type},nprobe))
%         filename = sprintf('%s %s %s event (%s filtered) probe %i.pdf',options.SUBJECT,options.SESSION,lfpAvg.event_group{event},filter_type{type},nprobe)
%         saveas(gcf,filename)
%         filename = sprintf('%s %s %s event (%s filtered) probe %i.fig',options.SUBJECT,options.SESSION,lfpAvg.event_group{event},filter_type{type},nprobe)
%         saveas(gcf,filename)
%         close


%
%
% for type = 1:length(lfpAvg.filter_type)
%     for event = 1:length(lfpAvg.event_group)
%         fig = figure
%         fig.Position = [334.7143 102.7143 800 830]
%         fig.Name = sprintf('%s %s %s event (%s filtered) probe %i',options.SUBJECT,options.SESSION,lfpAvg.event_group{event},filter_type{type},nprobe);
%         subplot(1,3,1);
%         colour_line= {[177,0,38]/256,[213,62,79]/256,[252,141,89]/256,...
%             [153,213,148]/256,[26,152,80]/256,[66,146,198]/256,[8,69,148]/256};
%         freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','300-600 Hz'};
%         selected_frequency = [1 2 6 7];
%         %     subplot(1,4,2)
%         hold on
%         plot([0 1],[chan_config.Ks_ycoord(best_channels.first_in_brain_channel) chan_config.Ks_ycoord(best_channels.first_in_brain_channel)],'--k','LineWidth',2)
%
%         if ~isempty(best_channels.L4_channel)
%             plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)-5) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)-5)],'--b','LineWidth',0.5)
%             plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)+5) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)+5)],'--b','LineWidth',0.5)
%             plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel))],'--b','LineWidth',2)
%         end
%
%         if ~isempty(best_channels.L5_channel)
%             plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)-10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L5_channel)-10)],'--c','LineWidth',0.5)
%
%             if isempty(best_channels.L4_channel) | find(chan_config.Channel ==best_channels.L5_channel)+10 < find(chan_config.Channel ==best_channels.L4_channel)-6
%                 plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)+10) chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)+10)],'--c','LineWidth',0.5)
%             else
%                 plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L4_channel)-6) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)-6)],'--c','LineWidth',0.5)
%             end
%
%             plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L5_channel))],'--c','LineWidth',2)
%
%         end
%
%         if ~isempty(best_channels.CA1_channel)
%             plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)-10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)-10)],'--r','LineWidth',0.5)
%             plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)+10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)+10)],'--r','LineWidth',0.5)
%             plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel))],'--r','LineWidth',2)
%         end
%
%         for n = selected_frequency
%             pl(n) = plot(power(:,n)./max(power(:,n)),sorted_config.Ks_ycoord','Color',colour_line{n})
%
%             %     p(n) = plot(power_differnece,1:96,colour_line{n})
%             hold on
%         end
%         % legend('1-3','4-12','30-60','60-100','125-300')
%         ylabel('probe depth (um)')
%         legend([pl(selected_frequency)],{freq_legends{selected_frequency}},'Location','southeast','Color','none');
%         ylim([0 4000])
%         set(gca,"TickDir","out",'box', 'off','Color','none')
%
%
%         subplot(1,3,2)
%         colororder({'k','r'})
%
%         taxis = lfpAvg.(filter_type{type})(event).timestamps;
%         for ch=1:size(lfpAvg.(filter_type{type})(event).filtered,2)
%
%             %         sh_tmp = 1e0*(500000*lfpAvg.SO(event).raw(:,ch)) + sorted_config.Ks_ycoord(ch);
%             sh_tmp = 1e0*(500000*lfpAvg.(filter_type{type})(event).filtered(:,ch)) + sorted_config.Ks_ycoord(ch);
%             plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
%             clear sh_tmp
%         end
%         % ylim([-1000 offset+1000]);
%         xlim([taxis(1) taxis(end)]);
%         xlabel('time (ms)');title('LFP');
%         plot([0 0],ylim,'--r');hold on;
%
%         title(sprintf('LFP (%s filtered)',filter_type{type}))
%         hold on
%         plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(best_channels.first_in_brain_channel) chan_config.Ks_ycoord(best_channels.first_in_brain_channel)],'--k','LineWidth',2)
%
%         if ~isempty(best_channels.L4_channel)
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)-5) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)-5)],'--b','LineWidth',0.5)
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)+5) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)+5)],'--b','LineWidth',0.5)
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel))],'--b','LineWidth',2)
%         end
%
%         if ~isempty(best_channels.L5_channel)
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)-10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L5_channel)-10)],'--c','LineWidth',0.5)
%
%             if isempty(best_channels.L4_channel) | find(chan_config.Channel ==best_channels.L5_channel)+10 < find(chan_config.Channel ==best_channels.L4_channel)-6
%                 plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)+10) chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)+10)],'--c','LineWidth',0.5)
%             else
%                 plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L4_channel)-6) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)-6)],'--c','LineWidth',0.5)
%             end
%
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L5_channel))],'--c','LineWidth',2)
%         end
%         if ~isempty(best_channels.CA1_channel)
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)-10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)-10)],'--r','LineWidth',0.5)
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)+10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)+10)],'--r','LineWidth',0.5)
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel))],'--r','LineWidth',2)
%         end
%         ylim([0 4000])
%         set(gca,"TickDir","out",'box', 'off','Color','none')
%
%         yyaxis right
%         ylim([-chan_config.Ks_ycoord(best_channels.first_in_brain_channel) 4000-chan_config.Ks_ycoord(best_channels.first_in_brain_channel)])
% %         ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'r';
%
%         subplot(1,3,3);
%         taxis = csd.(filter_type{type})(event).timestamps;
%         cmax = max(max(csd.(filter_type{type})(event).data(:,1:50)));
%         contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.(filter_type{type})(event).data,1)); csd.(filter_type{type})(event).data'; zeros(1,size(csd.(filter_type{type})(event).data,1))],40,'LineColor','none');hold on;
%         colormap jet; caxis([-cmax/2 cmax/2]);
%         xlabel('time (s)');title('CSD');
%         plot([0 0],[1 size(csd.(filter_type{type})(event).data,2)],'--k');hold on;
%         % set(gca,'YDir','reverse')
%         % yticks(1:length(sorted_config.Ks_ycoord))
%         % yticklabels(sorted_config.Ks_ycoord)
%         plot([0 0],ylim,'--r');hold on;
%         title(sprintf('CSD (%s filtered)',filter_type{type}))
%         colorbar
%         hold on
%         plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(best_channels.first_in_brain_channel) chan_config.Ks_ycoord(best_channels.first_in_brain_channel)],'--k','LineWidth',2)
%
%         if ~isempty(best_channels.L4_channel)
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)-5) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)-5)],'--b','LineWidth',0.5)
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)+5) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)+5)],'--b','LineWidth',0.5)
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel))],'--b','LineWidth',2)
%         end
%         if ~isempty(best_channels.L5_channel)
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)-10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L5_channel)-10)],'--c','LineWidth',0.5)
%             if isempty(best_channels.L4_channel) | find(chan_config.Channel ==best_channels.L5_channel)+10 < find(chan_config.Channel ==best_channels.L4_channel)-6
%                 plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)+10) chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)+10)],'--c','LineWidth',0.5)
%             else
%                 plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L4_channel)-6) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L4_channel)-6)],'--c','LineWidth',0.5)
%             end
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels.L5_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.L5_channel))],'--c','LineWidth',2)
%         end
%         if ~isempty(best_channels.CA1_channel)
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)-10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)-10)],'--r','LineWidth',0.5)
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)+10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)+10)],'--r','LineWidth',0.5)
%             plot([min(taxis) max(taxis)],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels.CA1_channel))],'--r','LineWidth',2)
%         end
%         ylim([0 4000])
%
%         set(gca,"TickDir","out",'box', 'off','Color','none')
%
%
% %         sgtitle(sprintf('%s %s %s event (%s filtered) probe %i',options.SUBJECT,options.SESSION,lfpAvg.event_group{event},filter_type{type},nprobe))
% %         filename = sprintf('%s %s %s event (%s filtered) probe %i.pdf',options.SUBJECT,options.SESSION,lfpAvg.event_group{event},filter_type{type},nprobe)
% %         saveas(gcf,filename)
% %         filename = sprintf('%s %s %s event (%s filtered) probe %i.fig',options.SUBJECT,options.SESSION,lfpAvg.event_group{event},filter_type{type},nprobe)
% %         saveas(gcf,filename)
%         %         close
%
%     end
% end

