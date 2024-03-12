function plot_theta_modulation(clusters,place_fields,theta_modulation_L,theta_modulation_R,options)

if ~isempty(theta_modulation_L) & ~isempty(theta_modulation_R)
    no_of_probes = 2;
elseif ~isempty(theta_modulation_L) | ~isempty(theta_modulation_R)
    no_of_probes = 1;
end

spatial_cell_index = unique([find(place_fields(1).peak_percentile>0.95 & place_fields(1).odd_even_stability>0.95)...
    find(place_fields(2).peak_percentile>0.95 & place_fields(2).odd_even_stability>0.95)]);

good_cell_index = intersect(spatial_cell_index,find(ismember(place_fields(1).cluster_id,clusters.cluster_id)));

% good_cell_index = good_cell_index(1:2);
sprintf('Theta modulation analysis for %i spatial cells',length(good_cell_index))

position_edges = 0:2:140;
phase_edges = 0:20:360;
phase_edges = pi/180.*phase_edges;

position_bins = position_edges(1:end-1) + diff(position_edges)/2;
bin_centres = phase_edges(1:end-1) + diff(phase_edges)/2;
bin_centres = 180*bin_centres'./pi;

for ncell = 1 : length(good_cell_index)
    count = 1;
    %     for nprobe = 1:no_of_probes
    modulated = 0;
    
    if ~isempty(theta_modulation_L)
        cell_index = find(theta_modulation_L(1).cluster_id==place_fields(1).cluster_id(good_cell_index(ncell)));
        modulated = modulated + theta_modulation_L(1).theta_modulation_percentile(cell_index)>0.95;
        modulated = modulated + theta_modulation_L(2).theta_modulation_percentile(cell_index)>0.95;
    end

    if ~isempty(theta_modulation_R)
        cell_index = find(theta_modulation_R(1).cluster_id==place_fields(1).cluster_id(good_cell_index(ncell)));
        modulated = modulated + theta_modulation_R(1).theta_modulation_percentile(cell_index)>0.95;
        modulated = modulated + theta_modulation_R(2).theta_modulation_percentile(cell_index)>0.95;
    end

    if modulated == 0 % only plot theta modulated cells
        continue
    end

    for track = 1 : length(place_fields)
        % For each good cell in the track


        if count == 1
            fig = figure;
            fig.Position = [400 100 1350 880];
            fig.Name = sprintf('%s %s theta modulation cell %i %s depth %i um',options(1).SUBJECT,options(1).SESSION,...
                place_fields(1).cluster_id(good_cell_index(ncell)),clusters.region(clusters.cluster_id==place_fields(1).cluster_id(good_cell_index(ncell))),...
                place_fields(1).peak_depth(good_cell_index(ncell)));
            sgtitle(sprintf('theta modulation cell %i %s depth %i um',place_fields(1).cluster_id(good_cell_index(ncell)),...
                clusters.region(clusters.cluster_id==place_fields(1).cluster_id(good_cell_index(ncell))),place_fields(1).peak_depth(good_cell_index(ncell))))
        end

        if ~isempty(theta_modulation_L)
            cell_index = find(theta_modulation_L(1).cluster_id==place_fields(1).cluster_id(good_cell_index(ncell)));

            ax1 = subplot(4,4,count);
            imagesc(position_edges(1:end-1)+1,180/pi.*phase_edges(2:end),theta_modulation_L(track).position_phase_map{cell_index}')
            set(gca,'YDir','normal')
            colorbar
            colormap(ax1,flip(gray))
%             colormap(summer)
            xticks(position_edges(1:10:end-1)+1)
            yticks(180/pi.*phase_edges(2:3:end))
            title(sprintf('Track %i cell %i %s (L theta)',track,theta_modulation_L(1).cluster_id(cell_index),clusters.region(clusters.cluster_id==place_fields(1).cluster_id(cell_index))))
            xlabel('Spatial location')
            ylabel('Theta phase')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            subplot(4,4,count+4)
            plot(bin_centres,theta_modulation_L(track).theta_phase_map{cell_index})
            xticks(0:45:360)
            xlabel('theta phases')
            ylabel('FR')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            if theta_modulation_L(track).theta_modulation_percentile(cell_index) > 0.95
                title(sprintf('L Theta modulation index %.2f',theta_modulation_L(track).theta_modulation_index(cell_index)),Color='r')
            else
                title(sprintf('L Theta modulation index %.2f',theta_modulation_L(track).theta_modulation_index(cell_index)))
            end

            ax2 = subplot(4,4,count+8);
            lags= -2*69:2:69*2;
            imagesc(lags,180/pi.*phase_edges(2:end), theta_modulation_L(track).position_phase_xcorr_map{cell_index}')
            xlim([-30 30])
            set(gca,'YDir','normal')
            hold on;xline(0,'r')
            colorbar
            colormap(ax2,parula)
            %             colormap(flip(gray))
            xticks(lags(  find(lags==0)-30:4:find(lags==0)+30))
            yticks(180/pi.*phase_edges(2:3:end))
            xlabel('Spatial shift')
            ylabel('Theta phase')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            title(sprintf('Sin fit phase %.2f (ampitude percentile %.2f)',theta_modulation_L(track).sin_fit_phase_offset(cell_index),theta_modulation_L(track).phase_amplitude_percentile(cell_index)));


            subplot(4,4,count+12)
            
            gaussianWindow = gausswin(8/mean(diff(position_edges))*2.5*2+1);% 4 datapoint or 8cm position bin for one SD (alpha is fixed at 2.5)
            gaussianWindow = gaussianWindow/sum(gaussianWindow);
            average_map = filtfilt(gaussianWindow,1,mean(place_fields(track).raw{good_cell_index(ncell)}));

            [pks,locs] = findpeaks(average_map);% Find local maxima and ignore peaks at start or end
            TF = islocalmin(average_map);
            if ~isempty(locs)
                [~,index]=max(pks);
                peak_index= locs(index);
                TF = find(TF);
                lower_min = TF(find(TF<peak_index));
                if ~isempty(lower_min)
                    lower_min = lower_min(end);
                else
                    lower_min = 1;
                end

                upper_min= TF(find(TF>peak_index));

                if ~isempty(upper_min)
                    upper_min = upper_min(1);
                else
                    upper_min = 70;
                end
            else
                [~,peak_index]=max(average_map);
            end

            [peak_FR]=max(average_map);
            average_map = (average_map-min(average_map))/(peak_FR+min(average_map));
            average_map = average_map/max(average_map);


            thresholded = find((average_map(peak_index+1:end)<0.3) == 1);
            if ~isempty(thresholded)
                thresholded = thresholded(1);
                upper_bound =  peak_index + thresholded;
            else
                upper_bound = length(average_map);
            end

            if upper_bound>upper_min
                upper_bound = upper_min;
            end

            thresholded = find((average_map(1:peak_index-1)<0.3) == 1);
            if ~isempty(thresholded)
                lower_bound =  thresholded(end);
            else
                lower_bound = 0;
            end

            if lower_bound<lower_min
                lower_bound = lower_min;
            end

            spike_id = theta_modulation_L(track).spike_times_phase_position{cell_index}(:,3) > lower_bound*2 & theta_modulation_L(track).spike_times_phase_position{cell_index}(:,3) < upper_bound*2;
            if sum(spike_id)>10
                [circ_corr,pval] = circ_corrcl(wrapTo2Pi(theta_modulation_L(track).spike_times_phase_position{cell_index}(spike_id,2)),theta_modulation_L(track).spike_times_phase_position{cell_index}(spike_id,3));
                theta_modulation_L(track).circ_corr(cell_index)=circ_corr;
                theta_modulation_L(track).circ_corr_pval(cell_index)=pval;

                scatter(theta_modulation_L(track).spike_times_phase_position{cell_index}(spike_id,3),180/pi*wrapTo2Pi(theta_modulation_L(track).spike_times_phase_position{cell_index}(spike_id,2)),4,'r','filled','MarkerFaceAlpha',0.1)
                xlabel('position')
                ylabel('theta phase')

                hold on
                scatter(theta_modulation_L(track).spike_times_phase_position{cell_index}(~spike_id,3),180/pi*wrapTo2Pi(theta_modulation_L(track).spike_times_phase_position{cell_index}(~spike_id,2)),4,'k','filled','MarkerFaceAlpha',0.05)
                xlabel('position')
                ylabel('theta phase')
                ylim([0 380])
                set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

                set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

                plot(position_bins,average_map*380,'b')
                yyaxis right
                ax = gca;
                ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'b';
                %             yylabel('FR')
                ylim([min(filtfilt(gaussianWindow,1,mean(place_fields(track).raw{good_cell_index(ncell)}))) max(filtfilt(gaussianWindow,1,mean(place_fields(track).raw{good_cell_index(ncell)})))])
                ylabel('FR')

                if pval <=0.05
                    title(sprintf('circ corr %.2f %.2f',circ_corr,pval),Color='r');
                else
                    title(sprintf('circ corr %.2f %.2f',circ_corr,pval));
                end
            end
            count = count + 1;
        end



        if ~isempty(theta_modulation_R)
            cell_index = find(theta_modulation_R(1).cluster_id==place_fields(1).cluster_id(good_cell_index(ncell)));

            ax1 = subplot(4,4,count);
            imagesc(position_edges(1:end-1)+1,180/pi.*phase_edges(2:end),theta_modulation_R(track).position_phase_map{cell_index}')
            set(gca,'YDir','normal')
            colorbar
            colormap(ax1,flip(gray))
            xticks(position_edges(1:10:end-1)+1)
            yticks(180/pi.*phase_edges(2:3:end))
            title(sprintf('Track %i cell %i %s (R theta)',track,theta_modulation_R(1).cluster_id(cell_index),clusters.region(clusters.cluster_id==place_fields(1).cluster_id(cell_index))))
            xlabel('Spatial location')
            ylabel('Theta phase')
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            subplot(4,4,count+4)
            plot(bin_centres,theta_modulation_R(track).theta_phase_map{cell_index})
            xticks(0:45:360)
            xlabel('theta phases')
            ylabel('FR')

            if theta_modulation_R(track).theta_modulation_percentile(cell_index) > 0.95
                title(sprintf('R Theta modulation index %.2f',theta_modulation_R(track).theta_modulation_index(cell_index)),Color='r')
            else
                title(sprintf('R Theta modulation index %.2f',theta_modulation_R(track).theta_modulation_index(cell_index)))
            end
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            ax2 = subplot(4,4,count+8);
            lags= -2*69:2:69*2;
            imagesc(lags,180/pi.*phase_edges(2:end), theta_modulation_R(track).position_phase_xcorr_map{cell_index}')
            set(gca,'YDir','normal')
            xlim([-30 30])
            colormap(ax2,parula)
            hold on;xline(0,'r')
            colorbar
%             colormap(flip(gray))
            xticks(lags(  find(lags==0)-30:4:find(lags==0)+30))
            yticks(180/pi.*phase_edges(2:3:end))
            xlabel('Spatial shift')
            ylabel('Theta phase')

            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
            title(sprintf('Sin fit phase %.2f (ampitude percentile %.2f)',theta_modulation_R(track).sin_fit_phase_offset(cell_index),theta_modulation_R(track).phase_amplitude_percentile(cell_index)));


            subplot(4,4,count+12)
            gaussianWindow = gausswin(8/mean(diff(position_edges))*2.5*2+1);% 4 datapoint or 8cm position bin for one SD (alpha is fixed at 2.5)
            gaussianWindow = gaussianWindow/sum(gaussianWindow);
            average_map = filtfilt(gaussianWindow,1,mean(place_fields(track).raw{good_cell_index(ncell)}));

            [pks,locs] = findpeaks(average_map);% Find local maxima and ignore peaks at start or end
            TF = islocalmin(average_map);
            if ~isempty(locs)
                [~,index]=max(pks);
                peak_index= locs(index);
                TF = find(TF);
                lower_min = TF(find(TF<peak_index));
                if ~isempty(lower_min)
                    lower_min = lower_min(end);
                else
                    lower_min = 1;
                end

                upper_min= TF(find(TF>peak_index));

                if ~isempty(upper_min)
                    upper_min = upper_min(1);
                else
                    upper_min = 70;
                end
            else
                [~,peak_index]=max(average_map);

                if peak_index < 15 % avoid first peak due to track onset
                    [~,peak_index2]=max(average_map(15:end));
                    
                    peak_index = peak_index2+15;
                end
            end

            [peak_FR]=max(average_map);
            average_map = (average_map-min(average_map))/(peak_FR+min(average_map));
            average_map = average_map/max(average_map);

            
            thresholded = find((average_map(peak_index+1:end)<0.3) == 1);
            if ~isempty(thresholded)
                thresholded = thresholded(1);
                upper_bound =  peak_index + thresholded;
            else
                upper_bound = length(average_map);
            end

            if upper_bound>upper_min
                upper_bound = upper_min;
            end

            thresholded = find((average_map(1:peak_index-1)<0.3) == 1);
            if ~isempty(thresholded)
                lower_bound =  thresholded(end);
            else
                lower_bound = 0;
            end

            if lower_bound<lower_min
                lower_bound = lower_min;
            end


            spike_id = theta_modulation_R(track).spike_times_phase_position{cell_index}(:,3) > lower_bound*2 & theta_modulation_R(track).spike_times_phase_position{cell_index}(:,3) < upper_bound*2;
            if sum(spike_id)>10
                [circ_corr,pval] = circ_corrcl(wrapTo2Pi(theta_modulation_R(track).spike_times_phase_position{cell_index}(spike_id,2)),theta_modulation_R(track).spike_times_phase_position{cell_index}(spike_id,3));
                theta_modulation_R(track).circ_corr(cell_index)=circ_corr;
                theta_modulation_R(track).circ_corr_pval(cell_index)=pval;

                scatter(theta_modulation_R(track).spike_times_phase_position{cell_index}(spike_id,3),180/pi*wrapTo2Pi(theta_modulation_R(track).spike_times_phase_position{cell_index}(spike_id,2)),4,'r','filled','MarkerFaceAlpha',0.05)
                xlabel('position')
                ylabel('theta phase')

                hold on
                scatter(theta_modulation_R(track).spike_times_phase_position{cell_index}(~spike_id,3),180/pi*wrapTo2Pi(theta_modulation_R(track).spike_times_phase_position{cell_index}(~spike_id,2)),4,'k','filled','MarkerFaceAlpha',0.01)
                xlabel('position')
                ylabel('theta phase')
                ylim([0 380])
                set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

                set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

                plot(position_bins,average_map*380,'b')
                yyaxis right
                ax = gca;
                ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'b';
                %             yylabel('FR')
                ylim([min(filtfilt(gaussianWindow,1,mean(place_fields(track).raw{good_cell_index(ncell)}))) max(filtfilt(gaussianWindow,1,mean(place_fields(track).raw{good_cell_index(ncell)})))])
                ylabel('FR')

                if pval <=0.05
                    title(sprintf('circ corr %.2f %.2f',circ_corr,pval),Color='r');
                else
                    title(sprintf('circ corr %.2f %.2f',circ_corr,pval));
                end
            end
                        count = count + 1;
        end

        %% phase precession vs location on track

        % circular-linear correlation coefficent from CircStat MATLABToolbox (Berens, 2009). Method first described in Kempter et al 2012
        % output of circ_corrcl is correlation coefficient and pval
        % needs input in radians
        %         if ~isempty(place_cell_times) % only useful for hippocampal cell
        %             [TPP(track).circ_lin_corr(ncell),TPP(track).circ_lin_PVAL(ncell)] = ...
        
        %         end


    end

end

% save('extracted_phase_precession_absolute_location','TPP','half_laps_times')
%
end
% end



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

function half_laps_times = extract_running_laps(position,lap_times)

parameters = list_of_parameters;

for track = 1 : length(lap_times)

    % Get half lap start and end time
    half_laps_timestamps = [lap_times(track).halfLaps_start' lap_times(track).halfLaps_stop'];

    % Split into 2 directions
    direction1 = half_laps_timestamps([1:2:size(half_laps_timestamps,1)],:);
    direction2 = half_laps_timestamps([2:2:size(half_laps_timestamps,1)],:);

    % Find these times in the position.data and get the indices
    direction1_idx = interp1(position.linear(track).timestamps,1:length(position.linear(track).timestamps),direction1,'nearest');
    direction2_idx = interp1(position.linear(track).timestamps,1:length(position.linear(track).timestamps),direction2,'nearest');

    % Turn this into logical index
    dir1_idx = zeros(length(position.linear(track).timestamps),1);
    for n = 1:size(direction1,1)
        dir1_idx(direction1_idx(n,1):direction1_idx(n,2)) = 1;
    end
    half_laps_times(track).direction_idx_1 = logical(dir1_idx);

    dir2_idx = zeros(length(position.linear(track).timestamps),1);
    for n = 1:size(direction2,1)
        dir2_idx(direction2_idx(n,1):direction2_idx(n,2)) = 1;
    end
    half_laps_times(track).direction_idx_2 = logical(dir2_idx);

    % Filter by speed
    speed_thresh = parameters.speed_threshold;   % arbitrarily chosen
    running_idx  =  position.v_cm(position.linear(track).clean_track_Indices) > speed_thresh;

    % remove portions along laps where animal not running
    half_laps_times(track).direction_idx_1(not(running_idx)) = 0;
    half_laps_times(track).direction_idx_2(not(running_idx)) = 0;

    % Save other variables
    half_laps_times(track).running_idx = running_idx;
    half_laps_times(track).direction1 = direction1;
    half_laps_times(track).direction2 = direction2;

end
end