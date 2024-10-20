
%% Ripple (DOWN -> UP vs UP -> DOWN)

all_spike_data{1} = [superficial_clusters.spike_id superficial_clusters.spike_times];
all_spike_data{2} = [L4_clusters.spike_id L4_clusters.spike_times];
all_spike_data{3} = [L5_clusters.spike_id L5_clusters.spike_times];
all_spike_data{4} = [CA1_clusters.spike_id CA1_clusters.spike_times];

group_name = {'superficial','L4','L5','CA1'};

% clear all_spike_data
% V1_T1_clusters.spike_id = [];
% V1_T1_clusters.spike_times = [];
% for cell = 1:length(V1_track_1_cells)
%     V1_T1_clusters.spike_id = [V1_T1_clusters.spike_id; V1_clusters.spike_id((V1_clusters.spike_id == V1_track_1_cells(cell)))];
%     V1_T1_clusters.spike_times =[V1_T1_clusters.spike_times; V1_clusters.spike_times((V1_clusters.spike_id == V1_track_1_cells(cell)))];
% end
% 
% [V1_T1_clusters.spike_times index] = sort(V1_T1_clusters.spike_times);
% V1_T1_clusters.spike_id(index) = V1_T1_clusters.spike_id;
% 
% all_spike_data{1} = [V1_T1_clusters.spike_id V1_T1_clusters.spike_times];
% all_spike_data{2} = [CA1_clusters.spike_id CA1_clusters.spike_times];
% 
% group_name = {'V1 Track 1 cells','CA1'};

[ripples.DOWN_UP_transition,index] = RestrictInts(ripples.peaktimes,[slow_waves.ints.UP(:,1) slow_waves.ints.UP(:,1)+0.2]);
plot_perievent_spiketime_histogram(all_spike_data,ripples.DOWN_UP_transition,'group','by cell zscore','group_name',group_name,'event_name','DOWN-UP transition ripples','twin',[-0.5 0.5])

[ripples.UP_DOWN_transition,index] = RestrictInts(ripples.peaktimes,[slow_waves.ints.DOWN(:,1)-0.2 slow_waves.ints.DOWN(:,1)]);
plot_perievent_spiketime_histogram(all_spike_data,ripples.UP_DOWN_transition,'group','by cell zscore','group_name',group_name,'event_name','UP-DOWN transition ripples','twin',[-0.5 0.5])
plot_perievent_spiketime_histogram(all_spike_data,ripples.UP_DOWN_transition,'group','by region','group_name',group_name,'event_name','UP-DOWN transition ripples','twin',[-0.5 0.5])

[ripples.UP_state,index] = RestrictInts(ripples.peaktimes,[slow_waves.ints.UP(:,1) + 0.2 slow_waves.ints.UP(:,2) - 0.2]);
plot_perievent_spiketime_histogram(all_spike_data,ripples.UP_state,'group','by cell zscore','group_name',group_name,'event_name','UP state','twin',[-0.5 0.5])

[ripples.UP_state,index] = RestrictInts(ripples.peaktimes,[slow_waves.ints.UP(:,1) + 0.2 slow_waves.ints.UP(:,2) - 0.2]);
plot_perievent_spiketime_histogram(all_spike_data,ripples.UP_state,'group','by region','group_name',group_name,'event_name','UP state','twin',[-0.5 0.5])

[ripples.DOWN_state,index] = RestrictInts(ripples.peaktimes,[slow_waves.ints.DOWN(:,1) slow_waves.ints.DOWN(:,2)]);
plot_perievent_spiketime_histogram(all_spike_data,ripples.DOWN_state,'group','by cell zscore','group_name',group_name,'event_name','DOWN state','twin',[-0.5 0.5])

%% UP state with ripple and without ripple

[slow_waves.DOWN_UP_transition_ripple,index] = RestrictInts(slow_waves.ints.UP(:,1)+0.2,[ripples.SWS_onset ripples.SWS_onset+0.2]);
plot_perievent_spiketime_histogram(all_spike_data,slow_waves.ints.UP(index,1),'group','by region','group_name',group_name,'event_name','DOWN-UP transition with ripples','twin',[-0.5 0.5])
plot_perievent_spiketime_histogram(all_spike_data,slow_waves.ints.UP(~index,1),'group','by region','group_name',group_name,'event_name','DOWN-UP transition without ripples','twin',[-0.5 0.5])
plot_perievent_spiketime_histogram(all_spike_data,slow_waves.ints.UP(index,1),'group','by cell zscore','group_name',group_name,'event_name','DOWN-UP transition with ripples','twin',[-0.5 0.5])
plot_perievent_spiketime_histogram(all_spike_data,slow_waves.ints.UP(~index,1),'group','by cell zscore','group_name',group_name,'event_name','DOWN-UP transition without ripples','twin',[-0.5 0.5])

% [slow_waves.UP_ripple,index] = RestrictInts(slow_waves.ints.UP(:,1)+0.300,[ripples.SWS_onset ripples.SWS_onset+0.1]);
% plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_ripple,'group','by cell zscore','group_name',group_name,'event_name','UP with ripples','twin',[-0.5 0.5])
% plot_perievent_spiketime_histogram(all_spike_data,slow_waves.ints.UP(~index,1),'group','by cell zscore','group_name',group_name,'event_name','UP without ripples','twin',[-0.5 0.5])

[slow_waves.UP_DOWN_ripple,index] = RestrictInts(slow_waves.ints.DOWN(:,1),[ripples.SWS_onset ripples.SWS_onset+0.1]);
slow_waves.UP_DOWN_no_ripple = slow_waves.ints.DOWN(~index,1);
plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_DOWN_ripple,'group','by region','group_name',group_name,'event_name','UP-DOWN with ripples','twin',[-0.5 0.5])
plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_DOWN_no_ripple,'group','by region','group_name',group_name,'event_name','UP-DOWN without ripples','twin',[-0.5 0.5])
plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_DOWN_ripple,'group','by cell zscore','group_name',group_name,'event_name','UP-DOWN with ripples','twin',[-0.5 0.5])
plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_DOWN_no_ripple,'group','by cell zscore','group_name',group_name,'event_name','UP-DOWN without ripples','twin',[-0.5 0.5])

%% UP state sequence
% all_spike_data{1} = [V1_clusters.spike_id V1_clusters.spike_times];
% all_spike_data{2} = [CA1_clusters.spike_id CA1_clusters.spike_times];
% group_name = {'V1','CA1'};
UP_index = [];
DOWN_index = [];
index = [];
Not_index = [];
for event = 1:length(slow_waves.ints.UP)
    this_event_window = [slow_waves.ints.UP(event,1) slow_waves.ints.UP(event,2)-0.2];
    SWR_number(event) = sum(ripples.SWS_onset >  slow_waves.ints.UP(event,1) & ripples.SWS_onset <  slow_waves.ints.UP(event,2));

    if sum(ripples.SWS_onset >  slow_waves.ints.UP(event,1)+0.2 & ripples.SWS_onset <  slow_waves.ints.UP(event,2)-0.2) > 0
        index = [index event];
    end

    if sum(ripples.SWS_onset >  slow_waves.ints.UP(event,1) & ripples.SWS_onset <  slow_waves.ints.UP(event,1) + 0.2) > 0
        UP_index = [UP_index event];
    end

    if sum(ripples.SWS_onset >  slow_waves.ints.UP(event,2)-0.2 & ripples.SWS_onset <  slow_waves.ints.UP(event,2)) > 0
        DOWN_index = [DOWN_index event];
    end

    if sum(ripples.SWS_onset >  slow_waves.ints.UP(event,1) & ripples.SWS_onset <  slow_waves.ints.UP(event,2)) == 0
        Not_index = [Not_index event];
    end
end

UP_duration =  slow_waves.ints.UP(:,2) - slow_waves.ints.UP(:,1);
DOWN_duration = slow_waves.ints.UP(2:end,1) - slow_waves.ints.UP(1:end-1,2);

slow_waves.UP_ripple = slow_waves.ints.UP(index,1);
slow_waves.DOWN_UP_ripple = slow_waves.ints.UP(UP_index,1);
slow_waves.UP_DOWN_ripple = slow_waves.ints.UP(DOWN_index+1,1);

slow_waves.both_ripple = slow_waves.ints.UP(intersect([UP_index index],DOWN_index+1),1);
slow_waves.DOWN_without_ripple_1 = slow_waves.ints.UP(Not_index(1:2:end),2);
slow_waves.DOWN_without_ripple_2 = slow_waves.ints.UP(Not_index(2:2:end),2);
% slow_waves.DOWN_without_ripple_3 = slow_waves.ints.UP(Not_index(3:3:end),2);
slow_waves.UP_without_ripple = slow_waves.ints.UP(Not_index,1);
slow_waves.UP_without_ripple_1 = slow_waves.ints.UP(Not_index(1:2:end),1);
slow_waves.UP_without_ripple_2 = slow_waves.ints.UP(Not_index(2:2:end),1);
% slow_waves.UP_without_ripple_3 = slow_waves.ints.UP(Not_index(3:3:end),1);
% [slow_waves.UP_ripple,index] = RestrictInts([slow_waves.ints.UP(:,1) slow_waves.ints.UP(:,2)],[ripples.SWS_onset ripples.SWS_onset+0.1]);
PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_without_ripple_1,'group','by cell zscore','sort_option','latency','group_name',group_name,...
    'bin_size',0.001,'event_name','UP without ripples','twin',[0 0.15]);
cell_id = [];

for n = 1:length(PSTH)
    cell_id{n} = PSTH(n).cell_id;
end

PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_without_ripple_2,'group','by cell zscore','sort_option','sorted','sorted_cell_id',cell_id,...
    'group_name',group_name,'bin_size',0.001,'event_name','UP without ripples 2','twin',[0 0.15]);
PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_ripple,'group','by cell zscore','sort_option','sorted','sorted_cell_id',cell_id,...
    'group_name',group_name,'bin_size',0.001,'event_name','UP with ripples','twin',[0 0.15]);
PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.DOWN_UP_ripple,'group','by cell zscore','sort_option','sorted','sorted_cell_id',cell_id,...
    'group_name',group_name,'bin_size',0.001,'event_name','DOWN-UP with ripples','twin',[0 0.15]);



PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_DOWN_ripple,'group','by cell zscore','sort_option','sorted','sorted_cell_id',cell_id,...
    'group_name',group_name,'bin_size',0.001,'event_name','UP_DOWN with ripples','twin',[0 0.40]);
PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_DOWN_ripple,'group','by region','sort_option','sorted','sorted_cell_id',cell_id,...
    'group_name',group_name,'bin_size',0.001,'event_name','UP_DOWN with ripples','twin',[0 0.40]);



%%
% FR during UP -> DOWN?
% FR during DOWN -> UP?
PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_without_ripple,'group','by cell zscore','sort_option','latency','group_name',group_name,...
    'bin_size',0.01,'event_name','UP without ripples','twin',[0 0.15],'smooth_option',0);
cell_id = [];

% PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_without_ripple_1,'group','by cell zscore','sort_option','latency','group_name',group_name,...
%     'bin_size',0.01,'event_name','UP without ripples','twin',[0 0.15],'smooth_option',0);


for n = 1:length(PSTH)
    cell_id{n} = PSTH(n).cell_id;
end
PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_without_ripple_1,'group','by cell zscore','sort_option','sorted','sorted_cell_id',cell_id,'group_name',group_name,...
    'bin_size',0.01,'event_name','UP without ripples 1','twin',[0 0.15],'smooth_option',0);
PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_without_ripple_2,'group','by cell zscore','sort_option','sorted','sorted_cell_id',cell_id,'group_name',group_name,...
    'bin_size',0.01,'event_name','UP without ripples 2','twin',[0 0.15],'smooth_option',0);

PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_without_ripple_1,'group','by cell zscore','sort_option','sorted','sorted_cell_id',cell_id,'group_name',group_name,...
    'bin_size',0.1,'event_name','UP without ripples','twin',[-0.2 0.2],'smooth_option',0);
PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_without_ripple_2,'group','by cell zscore','sort_option','sorted','sorted_cell_id',cell_id,'group_name',group_name,...
    'bin_size',0.1,'event_name','UP without ripples 2','twin',[-0.2 0.2],'smooth_option',0);


PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_without_ripple_1,'group','by cell zscore','sort_option','sorted','sorted_cell_id',cell_id,'group_name',group_name,...
    'bin_size',0.2,'event_name','UP without ripples 1','twin',[-0.2 0.2],'smooth_option',0,'plot_option',1);
PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_without_ripple_2,'group','by cell zscore','sort_option','sorted','sorted_cell_id',cell_id,'group_name',group_name,...
    'bin_size',0.2,'event_name','UP without ripples 2','twin',[-0.2 0.2],'smooth_option',0,'plot_option',1);
PSTH = plot_perievent_spiketime_histogram(all_spike_data,slow_waves.UP_without_ripple_2,'group','by cell zscore','sort_option','sorted','sorted_cell_id',cell_id,'group_name',group_name,...
    'bin_size',0.2,'event_name','UP without ripples','twin',[-0.2 0.2],'smooth_option',0,'plot_option',1);

