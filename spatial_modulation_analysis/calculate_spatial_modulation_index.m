function spatial_modulation = calculate_spatial_modulation_index(place_fields)


% Define Gaussian window for smoothing
gaussianWindow = gausswin(5);
% 8/mean(diff(position_edges))*2.5*2+1
% Normalize to have an area of 1 (i.e., to be a probability distribution)
gaussianWindow = gaussianWindow / sum(gaussianWindow);
spatial_modulation = [];
ratemaps_track1 = place_fields(1).raw{place_fields(1).cluster_id==cell_index(iCluster+(iPlot-1)*no_subplot)};
ratemaps_track2 = place_fields(2).raw{place_fields(2).cluster_id==cell_index(iCluster+(iPlot-1)*no_subplot)};

average_map_track1 = conv(mean(ratemaps_track1), gaussianWindow, 'same');
average_map_track1_odd = conv(mean(ratemaps_track1(1:2:end,:)), gaussianWindow, 'same');
average_map_track1_even = conv(mean(ratemaps_track1(2:2:end,:)), gaussianWindow, 'same');

average_map_track2 = conv(mean(ratemaps_track2), gaussianWindow, 'same');
average_map_track2_odd = conv(mean(ratemaps_track2(1:2:end,:)), gaussianWindow, 'same');
average_map_track2_even = conv(mean(ratemaps_track2(2:2:end,:)), gaussianWindow, 'same');

% (place_fields(2).x_bin_centres >= 50 & place_fields(2).x_bin_centres <= 130)

h(1)=plot(place_fields(1).x_bin_centres,average_map_track1_odd,'LineWidth',2,'Color',colour_lines{1});
map_error = std(average_map_track1_odd)./sqrt(length(average_map_track1_odd));
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([place_fields(1).x_bin_centres fliplr(place_fields(1).x_bin_centres)],[average_map_track1_odd+map_error fliplr(average_map_track1_odd-map_error)],colour_lines{1},'FaceAlpha','0.3','LineStyle','none');

[peak1,peak1_index] = max(average_map_track1(place_fields(1).x_bin_centres >= 46 & place_fields(1).x_bin_centres <= 86));
[peak2,peak2_index] = max(average_map_track1(place_fields(1).x_bin_centres > 86 & place_fields(1).x_bin_centres <= 126));

if  peak1 > peak2
    peak1_index = peak1_index+46/psthBinSize;
    peak2_index = peak1_index+40/psthBinSize;
    peak2 = average_map_track1(peak1_index+40/psthBinSize);

    SMI(1,iCluster+(iPlot-1)*no_subplot)=(peak1-peak2)/(peak1+peak2);
elseif  peak1 < peak2
    peak2_index = peak2_index+86/psthBinSize;
    peak1_index = peak2_index-40/psthBinSize;
    peak1 = average_map_track1(peak2_index-40/psthBinSize);
    SMI(1,iCluster+(iPlot-1)*no_subplot)=(peak2-peak1)/(peak1+peak2);
elseif  peak1 == peak2
    SMI(1,iCluster+(iPlot-1)*no_subplot)=0;
end


if ismember(cell_index(iCluster+(iPlot-1)*no_subplot),place_fields(1).cluster_id(place_fields(1).good_cells_LIBERAL))
    xline(place_fields(1).x_bin_centres(peak1_index),'LineWidth',2,'Color',colour_lines{1},'LineStyle','--')
    xline(place_fields(1).x_bin_centres(peak2_index),'LineWidth',2,'Color',colour_lines{1},'LineStyle','--')

else
    SMI(1,iCluster+(iPlot-1)*no_subplot)=nan;
end

h(2)=plot(place_fields(1).x_bin_centres,average_map_track1_even,'LineWidth',2,'Color',colour_lines{2});
map_error = std(average_map_track1_even)./sqrt(length(average_map_track1_even));
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([place_fields(1).x_bin_centres fliplr(place_fields(1).x_bin_centres)],[average_map_track1_even+map_error fliplr(average_map_track1_even-map_error)],colour_lines{2},'FaceAlpha','0.3','LineStyle','none');

h(3)=plot(place_fields(1).x_bin_centres,average_map_track2_odd,'LineWidth',2,'Color',colour_lines{3});
map_error = std(average_map_track2_odd)./sqrt(length(average_map_track2_odd));
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([place_fields(1).x_bin_centres fliplr(place_fields(1).x_bin_centres)],[average_map_track2_odd+map_error fliplr(average_map_track2_odd-map_error)],colour_lines{3},'FaceAlpha','0.3','LineStyle','none');

[peak1,peak1_index] = max(average_map_track2(place_fields(2).x_bin_centres >= 46 & place_fields(2).x_bin_centres <= 86));
[peak2,peak2_index] = max(average_map_track2(place_fields(2).x_bin_centres > 86 & place_fields(2).x_bin_centres <= 126));

if  peak1 > peak2
    peak1_index = peak1_index+46/psthBinSize;
    peak2_index = peak1_index+40/psthBinSize;
    peak2 = average_map_track2(peak2_index);
    SMI(2,iCluster+(iPlot-1)*no_subplot)=(peak1-peak2)/(peak1+peak2);
elseif  peak1 < peak2
    peak2_index = peak2_index+86/psthBinSize;
    peak1_index = peak2_index-40/psthBinSize;
    peak1 = average_map_track2(peak1_index);
    SMI(2,iCluster+(iPlot-1)*no_subplot)=(peak2-peak1)/(peak1+peak2);
elseif  peak1 == peak2
    SMI(2,iCluster+(iPlot-1)*no_subplot)=0;
end


