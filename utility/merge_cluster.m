function merged_clusters  = merge_cluster(clusters,match_ids)

% load clusters from a probe
%merged_id is the unit match suggestions
% find ids that occur more than once in merged_id
merged_clusters = clusters;

if size(match_ids,2) < 4 % Within session merge
    merged_id = match_ids(:,2); original_id = match_ids(:,1); unstable_id = logical(match_ids(:,3)); % convert from 0 based counting to 1 based counting
else % Across session merge
    merged_id = match_ids(:,4); original_id = match_ids(:,1); unstable_id = logical(match_ids(:,3)); % convert from 0 based counting to 1 based counting
    within_session_merge_id = match_ids(:,2);
end

merged_clusters.merged_cluster_id = clusters.cluster_id;

merged_clusters.merged_spike_id = clusters.spike_id;

unique_ids = sort(unique(merged_id));

[id_counts, edges] = histcounts(merged_id,[unique_ids; unique_ids(end)+1]-0.5); % find ids that occur more than once

ids_to_merge = unique_ids(id_counts > 1); % these ids have other original ids merged to this id so will only look at these clusters

for number_id = 1:length(ids_to_merge)

    id_temp = ids_to_merge(number_id);

    if size(match_ids,2) > 3
        original_ids_merged = original_id(merged_id == id_temp); % find the original ids of the merged ones
        within_session_id_this_cluster = within_session_merge_id(merged_id == id_temp); % find the within session ids of the merged ones
        original_ids_merged(within_session_id_this_cluster ~= mode(within_session_id_this_cluster)) = []; % Remove cluster that is not the same cluster according to within session merge 

%         original_ids_merged
    else
        original_ids_merged = original_id(merged_id == id_temp); % find the original ids of the merged ones

    end

    merged_clusters.merged_spike_id(ismember(clusters.spike_id,original_ids_merged)) = id_temp;%now merge the spikes in the clusters from original ids to the merged one

    merged_clusters.merged_cluster_id(ismember(clusters.cluster_id,original_ids_merged)) = id_temp; %convert cluster ids to merged ids
end

[C,ia,ib] = intersect(merged_clusters.cluster_id,original_id); % find clusters that exist during this session
merged_clusters.unstable_ids = unstable_id(ib); % actually from select_clusters point of view, logical is better
















