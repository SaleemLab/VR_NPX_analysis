function [psth,bin,binnedArray] = spatial_psth(spike_position,event_position, window, psthBinSize,position_bin_time)

no_bin = (window(2)-window(1))/psthBinSize;
binnedArray = zeros([length(event_position),no_bin]);
no_lap = length(event_position);
for iLap = 1:no_lap
    lap_position_edges = event_position(iLap)+window(1):psthBinSize:event_position(iLap)+window(2);
    spike_count_position(iLap,:) = histcounts(spike_position((spike_position >= event_position(iLap)+window(1)) & (spike_position <= event_position(iLap)+window(2))),lap_position_edges);
    
end

binnedArray = spike_count_position./position_bin_time;

psth = mean(binnedArray,'omitnan');
bin = (window(1)+psthBinSize/2):psthBinSize:(window(2)-psthBinSize/2);