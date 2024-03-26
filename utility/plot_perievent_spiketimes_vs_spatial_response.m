function plot_perievent_spiketimes_vs_spatial_response(spike_times,spike_id,Task_info,Behaviour,subplot_xy,window,psthBinSize,varargin)

% Default values
p = inputParser;
addParameter(p,'event_times',Task_info.end_time_all,@isnumeric) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'event_label',[],@isstr) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'event_id',[],@isnumeric) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)


addParameter(p,'place_fields',[],@isstruct) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'unit_depth',[],@isnumeric) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'unit_region',[],@isstring) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'unit_id',unique(spike_id),@isnumeric) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)

% assign parameters (either defaults or given)
parse(p,varargin{:});
place_fields = p.Results.place_fields;
unit_depth = p.Results.unit_depth;
unit_region = p.Results.unit_region;
unit_id = p.Results.unit_id;

event_label = p.Results.event_label;
event_times = p.Results.event_times;
event_id =  p.Results.event_id;

if ~isempty(place_fields) % if place fields, only select spatially tuned cells
    spatial_cell_index = unique([find(place_fields(1).peak_percentile>0.95 & place_fields(1).odd_even_stability>0.95)...
        find(place_fields(2).peak_percentile>0.95 & place_fields(2).odd_even_stability>0.95)]);
    good_cell_index = intersect(spatial_cell_index,find(ismember(place_fields(1).cluster_id,unit_id)));

    [Lia,Locb] = ismember(unit_id,place_fields(1).cluster_id(good_cell_index));
    unit_depth = unit_depth(Lia);
    unit_region = unit_region(Lia);
    unit_id = unit_id(Lia);
    x_bin_width = mean(diff(place_fields(1).x_bin_centres));
end


no_subplot_x = subplot_xy(1); %no of subplot in one figure columns
no_subplot_y = subplot_xy(2); %no of subplot in one figure rows

t_bin = mean(diff(Behaviour.tvec));
no_events = size(event_times,1);
time_edges = window(1):psthBinSize:window(2);
% spike_speed = interp1(Behaviour.tvec,Behaviour.speed,spike_times,'nearest');
spike_times_events = spike_times;
event_times_unchanged = event_times;

for nevent = 1:no_events
    if nevent < no_events
        spike_times_event_index = spike_times < event_times(nevent+1)+window(1) ...
            & spike_times >=  (event_times(nevent)+window(1));
    else
        spike_times_event_index = spike_times >=  (event_times(nevent)+window(1));
    end
    spike_times_events(spike_times_event_index) = spike_times_events(spike_times_event_index)+100000*(nevent);
    event_times(nevent,1) = event_times(nevent,1)+(nevent)*100000;

end

if isempty(event_label)
    event_label = 'event';
    event_id = ones(length(event_times),1);
end

cluster_spike_id = cell(size(unit_id));
no_cluster = length(unit_id);

no_subplot = no_subplot_x*no_subplot_y;
track1_ID = find(event_id == 1);
track2_ID = find(event_id == 2);


% Define Gaussian window for smoothing
gaussianWindow = gausswin(0.2*1/psthBinSize);
% Normalize to have an area of 1 (i.e., to be a probability distribution)
gaussianWindow = gaussianWindow / sum(gaussianWindow);


% Define Gaussian window for spatial smoothing
spatial_w = gausswin(11);
% Normalize to have an area of 1 (i.e., to be a probability distribution)
spatial_w = spatial_w / sum(spatial_w);


% interpolated tvec and track position
tvec = interp1(Behaviour.tvec,Behaviour.tvec,Behaviour.tvec(1):psthBinSize:Behaviour.tvec(end));
trackPosition = interp1(Behaviour.tvec,discretize(Behaviour.position,place_fields(1).x_bin_edges(1):x_bin_width:place_fields(1).x_bin_edges(end)),...
    Behaviour.tvec(1):psthBinSize:Behaviour.tvec(end),'previous');


colour_lines = {[145,191,219]/255,[69,117,180]/255,[215,48,39]/255,[252,141,89]/255};
for iPlot = 1: ceil(no_cluster/(no_subplot))
    fig = figure;
    fig.Position = [109 77 1700 878];
    fig.Name = sprintf('%s events raster plot %s %i',event_label,unit_region(1+(iPlot-1)*no_subplot),iPlot);

    for iCluster = 1:no_subplot
        if iCluster+(iPlot-1)*no_subplot > no_cluster
            break
        end
        subplot_scale = 0:2;
        subplot(no_subplot_y*5,no_subplot_x,iCluster+(floor((iCluster-1)/no_subplot_x)+subplot_scale)*no_subplot_x+floor((iCluster-1)/no_subplot_x)*no_subplot)
        cluster_spike_id{iCluster+(iPlot-1)*no_subplot} = spike_id == unit_id(iCluster+(iPlot-1)*no_subplot);

        [~,~,rasterX,rasterY,~,~] = psthAndBA(spike_times_events(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),event_times, window, 0.001);
        rasterX_reshaped = reshape(rasterX,[3 length(rasterX)/3]);
        rasterY_reshaped = reshape(rasterY,[3 length(rasterY)/3]);
        event1_raster_id = ismember(rasterY_reshaped(1,:), find(event_id == 1));
        event2_raster_id = ismember(rasterY_reshaped(1,:), find(event_id == 2));
        event1_rasterX = reshape(rasterX_reshaped(:,event1_raster_id),[1 sum(event1_raster_id)*3]);
        event1_rasterY = reshape(rasterY_reshaped(:,event1_raster_id),[1 sum(event1_raster_id)*3]);
        event2_rasterX = reshape(rasterX_reshaped(:,event2_raster_id),[1 sum(event2_raster_id)*3]);
        event2_rasterY = reshape(rasterY_reshaped(:,event2_raster_id),[1 sum(event2_raster_id)*3]);
        plot(event1_rasterX,event1_rasterY,'LineWidth',0.5,'Color',colour_lines{2})
        hold on;
        plot(event2_rasterX,event2_rasterY,'LineWidth', 0.5,'Color',colour_lines{3})

        xline(0, 'LineWidth',1,'Color',[0.5 0.5 0.5])
        ylim([0 length(event_id)])
        xlim(window)
        ylabel('Event')
%         xlabel('Time')
        title(sprintf('%s unit %i at %i um',unit_region(iCluster+(iPlot-1)*no_subplot),unit_id(iCluster+(iPlot-1)*no_subplot),unit_depth(iCluster+(iPlot-1)*no_subplot)))

        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

        subplot_scale = 3;
        subplot(no_subplot_y*5,no_subplot_x,iCluster+(floor((iCluster-1)/no_subplot_x)+subplot_scale)*no_subplot_x+floor((iCluster-1)/no_subplot_x)*no_subplot)

        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray1] = psthAndBA(spike_times_events(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),event_times(event_id==1), window, psthBinSize);
        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray2] = psthAndBA(spike_times_events(cluster_spike_id{iCluster+(iPlot-1)*no_subplot}),event_times(event_id==2), window, psthBinSize);


        if size(binnedArray1,1)==1
            psth_track1 = conv(binnedArray1/psthBinSize,gaussianWindow,'same');
        else
            for nevent = 1:size(binnedArray1,1)
                psth_track1(nevent,:) = conv(binnedArray1(nevent,:)/psthBinSize,gaussianWindow,'same');
            end
        end

        if size(binnedArray2,1)==1
            psth_track2 = conv(binnedArray2/psthBinSize,gaussianWindow,'same');
        else
            for nevent = 1:size(binnedArray2,1)
                psth_track2(nevent,:) = conv(binnedArray2(nevent,:)/psthBinSize,gaussianWindow,'same');
            end
        end

        average_map_track1 = mean(psth_track1,'omitnan');
%         average_map_track1_odd = conv(mean(psth_track1(1:2:end,:),'omitnan'), gaussianWindow, 'same');
%         average_map_track1_even = conv(mean(psth_track1(2:2:end,:),'omitnan'), gaussianWindow, 'same');

        average_map_track2 = mean(psth_track2,'omitnan');
%         average_map_track2_odd = conv(mean(psth_track2(1:2:end,:),'omitnan'), gaussianWindow, 'same');
%         average_map_track2_even = conv(mean(psth_track2(2:2:end,:),'omitnan'), gaussianWindow, 'same');


        h(1)=plot(bins,average_map_track1,'LineWidth',2,'Color',colour_lines{2});
        map_error = std(psth_track1)./sqrt(size(psth_track1,1));
        hold on
        % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
        patch([bins fliplr(bins)],[average_map_track1+map_error fliplr(average_map_track1-map_error)],colour_lines{2},'FaceAlpha','0.3','LineStyle','none');

%         h(2)=plot(bins,average_map_track1_even,'LineWidth',2,'Color',colour_lines{2});
%         map_error = std(average_map_track1_even)./sqrt(length(average_map_track1_even));
%         hold on
%         % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
%         patch([bins fliplr(bins)],[average_map_track1_even+map_error fliplr(average_map_track1_even-map_error)],colour_lines{2},'FaceAlpha','0.3','LineStyle','none');

        h(2)=plot(bins,average_map_track2,'LineWidth',2,'Color',colour_lines{3});
        map_error = std(psth_track2)./sqrt(length(psth_track2));
        hold on
        % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
        patch([bins fliplr(bins)],[average_map_track2+map_error fliplr(average_map_track2-map_error)],colour_lines{3},'FaceAlpha','0.3','LineStyle','none');

        %         h(4)=plot(bins,average_map_track2_even,'LineWidth',2,'Color',colour_lines{4});
        %         map_error = std(average_map_track2_even)./sqrt(length(average_map_track2_even));
        %         hold on
        %         % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
        %         patch([bins fliplr(bins)],[average_map_track2_even+map_error fliplr(average_map_track2_even-map_error)],colour_lines{4},'FaceAlpha','0.3','LineStyle','none');

        yHat = [];
        for track_id = 1:2
            position_this_event = [];
            responseProfile = conv(mean(place_fields(track_id).raw{place_fields(1).cluster_id==unit_id(iCluster+(iPlot-1)*no_subplot)}),spatial_w,'same');

            track_event_time = event_times_unchanged(event_id==track_id);

            for event = 1:length(track_event_time)
                if length(trackPosition(tvec>track_event_time(event)+window(1) & tvec<=track_event_time(event)+window(2))) == length(average_map_track1)
                    position_this_event(event,:) = trackPosition(tvec>track_event_time(event)+window(1) & tvec<=track_event_time(event)+window(2));
                elseif length(trackPosition(tvec>=track_event_time(event)+window(1) & tvec<=track_event_time(event)+window(2))) == length(average_map_track1)
                    position_this_event(event,:) = trackPosition(tvec>=track_event_time(event)+window(1) & tvec<=track_event_time(event)+window(2));
                elseif length(trackPosition(tvec>=track_event_time(event)+window(1) & tvec<track_event_time(event)+window(2))) == length(average_map_track1)
                    position_this_event(event,:) = trackPosition(tvec>=track_event_time(event)+window(1) & tvec<track_event_time(event)+window(2));
                elseif track_event_time(event)+window(2)>max(tvec)
                    temp = trackPosition(end);
                    position_this_event(event,:) = [trackPosition(tvec>track_event_time(event)+window(1) & tvec<track_event_time(event)+window(2))...
                        temp*ones(1,size(position_this_event,2)-sum(tvec>track_event_time(event)+window(1) & tvec<track_event_time(event)+window(2)))];
                else
                    position_this_event(event,:) = trackPosition(tvec>track_event_time(event)+window(1) & tvec<track_event_time(event)+window(2));
                end

                if contains(event_label,'Lap start')
                    position_this_event(event,bins<=0)=1;% Before lap start fix at location 1
                    position_this_event(event,bins<=0.1&position_this_event(event,:)>10)=1; % before lap start
                    position_this_event(isnan(position_this_event))=0;
                    position_jump = find([0 abs(diff(position_this_event(event,:)))]>5);
                    position_this_event(event,position_jump) =  position_this_event(event,position_jump+2);
                elseif contains(event_label,'Lap end')
                    position_this_event(event,bins>=0)=place_fields(1).x_bin_edges(end)/x_bin_width;
                end
                
                position_this_event(isnan(position_this_event))=0;
                position_jump = find([0 abs(diff(position_this_event(event,:)))]>3);
                position_this_event(event,position_jump) =  position_this_event(event,position_jump+2);

%                  position_this_event(position_this_event==70)=0;
                yHat{track_id}(event,:) = interp1(1:length(responseProfile), responseProfile, position_this_event(event,:)); % predicted firing rate
            end

            yHat{track_id}(isnan(yHat{track_id}))=0;
            average_resp_track{track_id} = mean( yHat{track_id},'omitnan');
        end

        h(3)=plot(bins,average_resp_track{1},'LineWidth',2,'Color',colour_lines{1});
        map_error = std(yHat{1})./sqrt(size(yHat{1},1));
        hold on
        % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
        patch([bins fliplr(bins)],[average_resp_track{1}+map_error fliplr(average_resp_track{1}-map_error)],colour_lines{1},'FaceAlpha','0.2','LineStyle','none');

        %         h(2)=plot(bins,average_map_track1_even,'LineWidth',2,'Color',colour_lines{2});
        %         map_error = std(average_map_track1_even)./sqrt(length(average_map_track1_even));
        %         hold on
        %         % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
        %         patch([bins fliplr(bins)],[average_map_track1_even+map_error fliplr(average_map_track1_even-map_error)],colour_lines{2},'FaceAlpha','0.3','LineStyle','none');

        h(4)=plot(bins,average_resp_track{2},'LineWidth',2,'Color',colour_lines{4});
        map_error = std(yHat{2})./sqrt(size(yHat{2},1));
        hold on
        % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
        patch([bins fliplr(bins)],[average_resp_track{2}+map_error fliplr(average_resp_track{2}-map_error)],colour_lines{4},'FaceAlpha','0.2','LineStyle','none');


        xline([0],'LineWidth',1,'Color',[0.5 0.5 0.5])
        if iCluster == no_subplot
            %             legend([h(1:4)],{'Track 1 odd','Track 1 even','Track 2 odd','Track 2 even'},'Color','none')
            legend([h(1:4)],{'Track Left event','Track Right event','Track Left spatial response','Track Right spatial response'},'Color','none')
        end
        xlim(window)
        xlabel('Time')
        ylabel('FR')
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

        % spatial map
        subplot_scale = 4;
        subplot(no_subplot_y*5,no_subplot_x,iCluster+(floor((iCluster-1)/no_subplot_x)+subplot_scale)*no_subplot_x+floor((iCluster-1)/no_subplot_x)*no_subplot)
        ratemaps_track1 = place_fields(1).raw{place_fields(1).cluster_id==unit_id(iCluster+(iPlot-1)*no_subplot)};
        ratemaps_track2 = place_fields(2).raw{place_fields(2).cluster_id==unit_id(iCluster+(iPlot-1)*no_subplot)};
        bins = place_fields(1).x_bin_centres;


        average_map_track1 = conv(mean(ratemaps_track1,'omitnan'), spatial_w, 'same');
%         average_map_track1_odd = conv(mean(psth_track1(1:2:end,:),'omitnan'), gaussianWindow, 'same');
%         average_map_track1_even = conv(mean(psth_track1(2:2:end,:),'omitnan'), gaussianWindow, 'same');

        average_map_track2 = conv(mean(ratemaps_track2,'omitnan'), spatial_w, 'same');
%         average_map_track2_odd = conv(mean(psth_track2(1:2:end,:),'omitnan'), gaussianWindow, 'same');
%         average_map_track2_even = conv(mean(psth_track2(2:2:end,:),'omitnan'), gaussianWindow, 'same');


        h(1)=plot(bins,average_map_track1,'LineWidth',2,'Color',colour_lines{2});
        map_error = std(ratemaps_track1)./sqrt(size(ratemaps_track1,1));
        hold on
        % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
        patch([bins fliplr(bins)],[average_map_track1+map_error fliplr(average_map_track1-map_error)],colour_lines{2},'FaceAlpha','0.3','LineStyle','none');

%         h(2)=plot(bins,average_map_track1_even,'LineWidth',2,'Color',colour_lines{2});
%         map_error = std(average_map_track1_even)./sqrt(length(average_map_track1_even));
%         hold on
%         % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
%         patch([bins fliplr(bins)],[average_map_track1_even+map_error fliplr(average_map_track1_even-map_error)],colour_lines{2},'FaceAlpha','0.3','LineStyle','none');

        h(2)=plot(bins,average_map_track2,'LineWidth',2,'Color',colour_lines{3});
        map_error = std(ratemaps_track2)./sqrt(size(ratemaps_track2,1));
        hold on
        % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
        patch([bins fliplr(bins)],[average_map_track2+map_error fliplr(average_map_track2-map_error)],colour_lines{3},'FaceAlpha','0.3','LineStyle','none');
        xlabel('Position')
        ylabel('FR')
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
    end
end