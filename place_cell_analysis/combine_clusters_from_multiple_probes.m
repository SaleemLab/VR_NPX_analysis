function clusters_combined = combine_clusters_from_multiple_probes(clusters)
% Functions for combining clusters from multiple probes into one structure
% 

% give new cell id to cell from probe 2
for ncell = 1:length(clusters.probe(2).cell_type)
    if clusters.probe(2).id_conversion(ncell,2) < 10000
        clusters.probe(2).spike_id(clusters.probe(2).spike_id == clusters.probe(2).id_conversion(ncell,2))...
            = clusters.probe(2).id_conversion(ncell,2) + 10000;
        clusters.probe(2).id_conversion(ncell,2) = clusters.probe(2).id_conversion(ncell,2) + 10000;
    end
end

clusters_combined = [];
combined_spike_times = [clusters.probe(1).spike_times; clusters.probe(2).spike_times];
[~,spike_index]  = sort(combined_spike_times) ;
clusters_combined.spike_id = [clusters.probe(1).spike_id; clusters.probe(2).spike_id];
clusters_combined.spike_id = clusters_combined.spike_id(spike_index);
clusters_combined.spike_times = combined_spike_times(spike_index);
clusters_combined.peak_channel = [clusters.probe(1).peak_channel clusters.probe(2).peak_channel];
clusters_combined.peak_channel_waveform = [clusters.probe(1).peak_channel_waveforms; clusters.probe(2).peak_channel_waveforms];
clusters_combined.cell_type = [clusters.probe(1).cell_type clusters.probe(2).cell_type];
clusters_combined.probe_no = [ones(1,length(clusters.probe(1).cell_type)) 2*ones(1,length(clusters.probe(2).cell_type))];

clusters_combined.id_conversion = [clusters.probe(1).id_conversion; clusters.probe(2).id_conversion];
clusters_combined.id_conversion(:,1) = 1:length(clusters_combined.cell_type);
clusters_combined.MUA_zscore = [clusters.probe(1).MUA_zscore; clusters.probe(1).MUA_zscore];
clusters_combined.MUA_tvec = clusters.probe(1).MUA_tvec;

end