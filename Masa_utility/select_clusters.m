function selected_clusters = select_clusters(clusters,params)

selected_clusters = struct();
all_metrics = fieldnames(params);
cluster_filter_index = ones(size(clusters.cluster_id));
for iMetric = 1:size(all_metrics,1)
    temp_cluster_index = params.(all_metrics{iMetric})(clusters.(all_metrics{iMetric}));
    cluster_filter_index = cluster_filter_index &temp_cluster_index;
end

all_cluster_fields = fieldnames(clusters);

for iField = 1:size(all_cluster_fields,1)
    if contains(all_cluster_fields{iField},'timevec')
        selected_clusters.timevec = clusters.timevec;
    elseif contains(all_cluster_fields{iField},'spike_id')
        spike_id_selected = ismember(clusters.spike_id,clusters.cluster_id(cluster_filter_index));
        selected_clusters.spike_id = clusters.spike_id(spike_id_selected);
    elseif contains(all_cluster_fields{iField},'spike_times')
        selected_clusters.spike_times = clusters.spike_times(spike_id_selected);
    elseif contains(all_cluster_fields{iField},'probe_id')
        selected_clusters.probe_id = clusters.probe_id;
    elseif contains(all_cluster_fields{iField},'sorter')
        selected_clusters.sorter = clusters.sorter;
    elseif contains(all_cluster_fields{iField},'probe_hemisphere')
        selected_clusters.sorter = clusters.probe_hemisphere;
    else
        temp_cluster_field = clusters.(all_cluster_fields{iField});
        selected_clusters.(all_cluster_fields{iField}) = temp_cluster_field(cluster_filter_index,:);
    end
end