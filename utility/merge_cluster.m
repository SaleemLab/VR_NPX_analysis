function merged_clusters  = merge_cluster(clusters,match_ids)

% load clusters from a probe
%merged_id is the unit match suggestions
% find ids that occur more than once in merged_id

merged_clusters = clusters;

merged_id = match_ids(:,2); original_id = match_ids(:,1); unstable_id = match_ids(:,3); % convert from 0 based counting to 1 based counting

merged_clusters.merged_cluster_id = clusters.cluster_id;

merged_clusters.merged_spike_id = clusters.spike_id;

unique_ids = sort(unique(merged_id));

[id_counts, edges] = histcounts(merged_id,[unique_ids; unique_ids(end)+1]-0.5); % find ids that occur more than once

ids_to_merge = unique_ids(id_counts > 1); % these ids have other original ids merged to this id so will only look at these clusters

for number_id = 1:length(ids_to_merge)

    id_temp = ids_to_merge(number_id);
    
    original_ids_merged = original_id(merged_id == id_temp); % find the original ids of the merged ones

    merged_clusters.merged_spike_id(ismember(clusters.spike_id,original_ids_merged)) = id_temp;%now merge the spikes in the clusters from original ids to the merged one 

    merged_clusters.merged_cluster_id(ismember(clusters.cluster_id,original_ids_merged)) = id_temp; %convert cluster ids to merged ids
   
end

merged_clusters.unstable_ids = original_id(unstable_id);
















