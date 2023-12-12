function place_fields_lap = plot_spatial_cell_tuning(clusters,place_fields_all,place_fields_even,place_fields_odd,position,lap_times,options)
% Function to plot spatial cell activity
% Plot 1 is each cell's spatial ratemap for each lap
% Plot 2 is each cell's spatial tuning stability by plotting place field
% calculated by using all laps, odd no laps and even no laps.
ROOTPATH = options.ROOTPATH;

if isfield(options,'region')
    region = options.region;
else
    region = [];
end

if isfield(options,'probe_combined') && options.probe_combined == 1
    nprobe = 'combined';
else
    nprobe = options.probe_no;
    probe_hemisphere = options.probe_hemisphere;
    probe_hemisphere_text = {'left','right'}
end
x_bins_width = 10;

%             clusters = HPC_clusters;
% clusters = CA1_clusters;
%             clusters = V1_clusters;

%     place_fields_even = CA1_place_fields_even.probe(nprobe);
%     place_fields_odd = CA1_place_fields_odd.probe(nprobe);
%     place_fields_all = CA1_place_fields.probe(nprobe);

% Grabbing ratemap for each lap for each cell
tic
place_fields_lap = [];
if ~ischar(nprobe)
    for track_id = 1:2
        for nlap = 1:length(lap_times(track_id).start)
            field_this_lap = get_lap_place_fields_masa(x_bins_width,position,place_fields_all,...
                clusters.probe(nprobe),track_id,lap_times(track_id).start(nlap),lap_times(track_id).end(nlap));

            place_fields_lap{track_id}{nlap} = field_this_lap.raw;
        end
    end
else
    for track_id = 1:2
        for nlap = 1:length(lap_times(track_id).start)
            field_this_lap = get_lap_place_fields_masa(x_bins_width,position,place_fields_all,...
                clusters,track_id,lap_times(track_id).start(nlap),lap_times(track_id).end(nlap));

            place_fields_lap{track_id}{nlap} = field_this_lap.raw;
        end
    end
end
toc
% Basic visualisation of each cell's spatial activity
for track_id = 1:2
    fig = figure(track_id)
    fig.Position = [300 150 1400 920]
    if ischar(nprobe)
        fig.Name = sprintf('%s %s %s spatial cell maps lap heatmap probes %s %i',options.SUBJECT,options.SESSION,region,nprobe,track_id);
    else
        fig.Name = sprintf('%s %s %s spatial cell maps lap heatmap probes %s %i',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere},track_id);
    end
    for ncell = 1:length(place_fields_all.all_cells)
        this_cell_place_field = [];
        for nlap = 1:length(lap_times(track_id).start)
            this_cell_place_field(nlap,:) = place_fields_lap{track_id}{nlap}{ncell};
        end

        subplot(ceil(sqrt(length(place_fields_all.all_cells))),ceil(sqrt(length(place_fields_all.all_cells))),ncell)
        imagesc(this_cell_place_field./max(this_cell_place_field')')
        %                     if sum(ncell ==
        %                         title(ncell,'Color','r')
        %                     else
        title(ncell,'Color','k')
        %                     end
        set(gca,"TickDir","out",'box', 'off','Color','none')
    end
    sgtitle(sprintf('Track %i',track_id))
end
% 
% 
ppf = 18;
pcount = 1;
nfigure= 3;
% Visualise each good place cell's spatial tuning stability
for ncell = place_fields_all.good_place_cells_LIBERAL  
    if pcount == 19
        nfigure = nfigure + 1
        pcount = 1
    end
    fig = figure(nfigure)
    fig.Position = [300 150 810 800];
    if ischar(nprobe)
        fig.Name = sprintf('%s %s %s spatial cell tuning probes %s (%i)',options.SUBJECT,options.SESSION,region,nprobe,nfigure);
    else
        fig.Name = sprintf('%s %s %s spatial cell tuning probes %s (%i)',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere},nfigure);
    end

    subplot(6,3,pcount)
    max_FR = max([place_fields_all.track(1).raw{ncell} place_fields_all.track(2).raw{ncell}]);
    plot(1:x_bins_width:140,place_fields_all.track(1).raw{ncell}/max_FR,'r'); hold on; plot(1:x_bins_width:140,place_fields_all.track(2).raw{ncell}/max_FR,'b')
    ylabel('Normalised FR')
    xlabel('Position bins')
    xticks([30 50 70 90 110 140]) % 1-5 landmarks and end
    xline(100,'r','reward')
    set(gca,"TickDir","out",'box', 'off','Color','none')
    title(sprintf('all laps cell %i',ncell))

    pcount = pcount + 1;
    subplot(6,3,pcount)
    max_FR = max([place_fields_odd.track(1).raw{ncell} place_fields_odd.track(2).raw{ncell}]);
    plot(1:x_bins_width:140,place_fields_odd.track(1).raw{ncell}/max_FR,'r'); hold on; plot(1:x_bins_width:140,place_fields_odd.track(2).raw{ncell}/max_FR,'b')
    title('odd laps')
    set(gca,"TickDir","out",'box', 'off','Color','none')
    xticks([30 50 70 90 110 140]) % 1-5 landmarks and end
    xline(100,'r')

    pcount = pcount + 1;
    subplot(6,3,pcount)
    max_FR = max([place_fields_even.track(1).raw{ncell} place_fields_even.track(2).raw{ncell}]);
    plot(1:x_bins_width:140,place_fields_even.track(1).raw{ncell}/max_FR,'r'); hold on; plot(1:x_bins_width:140,place_fields_even.track(2).raw{ncell}/max_FR,'b')
    title('even laps')
    pcount = pcount + 1;
    set(gca,"TickDir","out",'box', 'off','Color','none')
    xticks([30 50 70 90 110 140]) % 1-5 landmarks and end
    xline(100,'r')
end

%% plotting relative to shuffles 
% fig = figure;
% fig.Name = sprintf('%s %s %s spatial cell tuning curves vs shuffle probes %s',options.SUBJECT,options.SESSION,region,nprobe);
% subplot(2,2,1)
% hold on
% x = 1:x_bins_width:140
% 
% CI_shuffle = prctile(place_fields_all.track(1).raw_shuffled{150},[1 99]);
% plot(x, CI_shuffle(2,:), 'k--', 'LineWidth', 1);
% plot(x, CI_shuffle(1,:), 'k--', 'LineWidth', 1);
% x2 = [x, fliplr(x)];
% inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
% h(2) = fill(x2, inBetween, 'k','FaceAlpha',0.2);
% 
% h(1) = plot(1:x_bins_width:140,place_fields_all.track(1).raw{ncell},'r')
% 
% ylabel('FR')
% xlabel('Position bins')
% xticks([30 50 70 90 110 140]) % 1-5 landmarks and end
% xline(100,'r','reward')
% set(gca,"TickDir","out",'box', 'off','Color','none')
% legend([h(1) h(2)],{'Original','Shuffled'},'Color','none')
% title(sprintf('Track 1 Tuning curve cell %i',ncell))
% 
% 
% subplot(2,2,2)
% hold on
% x = 1:x_bins_width:140
% CI_shuffle = prctile(place_fields_all.track(2).raw_shuffled{150},[1 99]);
% plot(x, CI_shuffle(2,:), 'k--', 'LineWidth', 1);
% plot(x, CI_shuffle(1,:), 'k--', 'LineWidth', 1);
% x2 = [x, fliplr(x)];
% inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
% h(2) = fill(x2, inBetween, 'k','FaceAlpha',0.2);
% 
% h(1) = plot(1:x_bins_width:140,place_fields_all.track(2).raw{ncell},'b')
% 
% ylabel('FR')
% xlabel('Position bins')
% xticks([30 50 70 90 110 140]) % 1-5 landmarks and end
% xline(100,'r','reward')
% legend([h(1) h(2)],{'Original','Shuffled'},'Color','none')
% set(gca,"TickDir","out",'box', 'off','Color','none')
% title(sprintf('Track 2 Tuning curve cell %i',ncell))
fig = figure;
fig.Position = [500 70 1300 930];
if ischar(nprobe)
    sgtitle(sprintf('%s %s %s place cell maps sorted probes %s',options.SUBJECT,options.SESSION,region,nprobe))
    fig.Name = sprintf('%s %s %s place cell maps sorted probes %s',options.SUBJECT,options.SESSION,region,nprobe);
else
    sgtitle(sprintf('%s %s %s cell maps sorted probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere}))
    fig.Name = sprintf('%s %s %s cell maps sorted probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere})
end

c = 1;

for test = 1:3
    if test == 3
        place_fields = place_fields_all;
        title_text = 'All laps';
    elseif test == 1
        place_fields = place_fields_odd;
        %         sgtitle('Odd laps')
        title_text = 'Odd laps'
    elseif test == 2
        place_fields = place_fields_even;
        %         sgtitle('Even laps')
        title_text = 'Even laps'
    end

    for kk=1:length(place_fields_all.track)
        for j=1:length(place_fields_all.track)
            y_vector=[];
            matrix=[];
            normalized_matrix=[];
            if isempty(place_fields_all.track(j).sorted_good_cells_LIBERAL)
                continue
            end
            for ii=1:length(place_fields_all.track(j).sorted_good_cells_LIBERAL)
                %plot sorted
                %             matrix=[];
                %             normalized_matrix=[];
                matrix(ii,:)=place_fields.track(kk).raw{place_fields_all.track(j).sorted_good_cells_LIBERAL(ii)};
                normalized_matrix(ii,:)=(matrix(ii,:)-min(matrix(ii,:)))/(max(matrix(ii,:))-min(matrix(ii,:)));
                %                 subplot(length(place_fields.track),length(place_fields.track),c)
                subplot(3,4,c)
                plfield_row= normalized_matrix(ii,:)+(1.5*ii-1);
                plot(1:length(plfield_row),plfield_row,'k'); hold on;
                xx = [1:length(plfield_row), fliplr(1:length(plfield_row))];
                inBetween = [(1.5*ii-1)*ones(size(plfield_row)), fliplr(plfield_row)];
                fill(xx, inBetween,[139,0,0]/255);
                y_vector= [y_vector, 1.5*ii-1];
            end
            xlim([0 size(normalized_matrix,2)+2]);
            ylim([0 max(y_vector)+1.2]);
            yt=place_fields_all.track(j).sorted_good_cells_LIBERAL;
            set(gca,'ytick',y_vector);
            set(gca,'yticklabel',yt);
            xticks([30 50 70 90 110 140]/10)
            xticklabels([30 50 70 90 110 140])
            xline(10,'r')
            ylabel('Unit ID');
            xlabel('sorted linearized position (bins)');
            c=c+1;
            title([{[title_text,' place cells on track ' num2str(j)]} ; {['sorted by track ' num2str(kk)]}]);
            set(gca,"TickDir","out",'box', 'off','Color','none')
        end
    end

end




% if ischar(nprobe)
%     sgtitle(sprintf('%s %s place cell heatmaps sorted probes %s',options.SUBJECT,options.SESSION,nprobe))
% 
%     filename = sprintf('%s %s place cell heatmaps sorted probes %s.pdf',options.SUBJECT,options.SESSION,nprobe)
%     saveas(gcf,filename)
%     filename = sprintf('%s %s place cell heatmaps sorted probes %s.fig',options.SUBJECT,options.SESSION,nprobe)
%     saveas(gcf,filename)
% else
%     sgtitle(sprintf('%s %s place cell heatmaps sorted probe %i',options.SUBJECT,options.SESSION,nprobe))
% 
%     filename = sprintf('%s %s place cell heatmaps sorted probe %i.pdf',options.SUBJECT,options.SESSION,nprobe)
%     saveas(gcf,filename)
%     filename = sprintf('%s %s place cell heatmaps sorted probe %i.fig',options.SUBJECT,options.SESSION,nprobe)
%     saveas(gcf,filename)
% end
