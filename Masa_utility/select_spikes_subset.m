function [spike_times,spike_id]= select_spikes_subset(clusters,cell_index)

spike_times = [];
spike_id = [];
for ncell = cell_index
    original_id = clusters.id_conversion(ncell,2);
    spike_times = [spike_times; clusters.spike_times(clusters.spike_id == original_id)];
    spike_id = [spike_id; original_id*ones(sum(clusters.spike_id == original_id),1)];

    if isempty(clusters.spike_times(clusters.spike_id == original_id)) % Just as a place holder for spike id even when cell does not fire at all
        spike_times = [spike_times; 0];
        spike_id = [spike_id; original_id];
    end
end

[spike_times,index]= sort(spike_times);
spike_id = spike_id(index);

end