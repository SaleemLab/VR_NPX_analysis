function clusters_combined = combine_clusters_from_multiple_probes(clusters1,clusters2)
% Functions for combining clusters from multiple probes into one structure
%
all_fields = fieldnames(clusters1);

if isfield(clusters2,'cluster_id')
    if clusters2.cluster_id(1)<10000 % if unique id is already given
        for ncell =1:length(clusters2.cluster_id)
            % give cluster 2 spikes new cluster id
            clusters2.spike_id(clusters2.spike_id == clusters2.cluster_id(ncell)) = clusters2.cluster_id(ncell) + 10000;
            clusters2.cluster_id(ncell) = clusters2.cluster_id(ncell) + 10000;

        end
    end
end

if isfield(clusters2,'merged_cluster_id')
    merged_id = unique(clusters2.merged_cluster_id);
    for ncell =1:length(merged_id)
        clusters2.merged_spike_id(clusters2.merged_spike_id == merged_id(ncell)) = merged_id(ncell) + 10000;
    end

    for ncell =1:length(clusters2.merged_cluster_id)
        clusters2.merged_cluster_id(ncell) = clusters2.merged_cluster_id(ncell) + 10000;
    end
end

clusters_combined = [];

for n = 1:length(all_fields)

    clusters_combined.(all_fields{n}) = [clusters1.(all_fields{n}); clusters2.(all_fields{n})];
end

if isfield(clusters1,'spike_id')
    combined_spike_times = [clusters1.spike_times; clusters2.spike_times];
    [~,spike_index]  = sort(combined_spike_times) ;
    % clusters_combined.spike_id = [clusters1.spike_id; clusters2.spike_id];
    clusters_combined.spike_id = clusters_combined.spike_id(spike_index);
    clusters_combined.spike_times = clusters_combined.spike_times(spike_index);
end

if isfield(clusters1,'merged_spike_id')
    clusters_combined.merged_spike_id = clusters_combined.merged_spike_id(spike_index);
%     clusters_combined.spike_times = clusters_combined.spike_times(spike_index);
end


end