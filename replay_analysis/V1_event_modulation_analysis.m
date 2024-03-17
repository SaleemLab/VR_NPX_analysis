function V1_event_modulation  = V1_event_modulation_analysis(spike_times,spike_id,Task_info,Behaviour,window,psthBinSize,varargin)

% Default values
p = inputParser;
addParameter(p,'event_times',Task_info.end_time_all,@isnumeric) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'event_label',[],@iscell) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'event_id',[],@isnumeric) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)


addParameter(p,'place_fields',[],@isstruct) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'unit_depth',[],@isnumeric) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'unit_region',[],@isstring) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)
addParameter(p,'unit_id',unique(spike_id),@isnumeric) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn)

% assign parameters (either defaults or given)
parse(p,varargin{:});
place_fields = p.Results.place_fields;
unit_depth = p.Results.unit_depth;
unit_region = p.Results.unit_region;
unit_id = p.Results.unit_id;

event_label = p.Results.event_label;
event_times = p.Results.event_times;
event_id =  p.Results.event_id;

if ~isempty(place_fields) % if place fields, only select spatially tuned cells
    spatial_cell_index = unique([find(place_fields(1).peak_percentile>0.95 & place_fields(1).odd_even_stability>0.95)...
        find(place_fields(2).peak_percentile>0.95 & place_fields(2).odd_even_stability>0.95)]);
    good_cell_index = intersect(spatial_cell_index,find(ismember(place_fields(1).cluster_id,unit_id)));

    [Lia,Locb] = ismember(unit_id,place_fields(1).cluster_id(good_cell_index));
    unit_depth = unit_depth(Lia);
    unit_region = unit_region(Lia);
    unit_id = unit_id(Lia);
end

t_bin = mean(diff(Behaviour.tvec));
no_events = size(event_times,1);
time_edges = window(1):psthBinSize:window(2);
% spike_speed = interp1(Behaviour.tvec,Behaviour.speed,spike_times,'nearest');
spike_times_events = spike_times;

for nevent = 1:no_events
    if nevent < no_events
        spike_times_event_index = spike_times < event_times(nevent+1)+window(1) ...
            & spike_times >=  (event_times(nevent)+window(1));
    else
        spike_times_event_index = spike_times >=  (event_times(nevent)+window(1));
    end
    spike_times_events(spike_times_event_index) = spike_times_events(spike_times_event_index)+100000*(nevent);
    event_times(nevent,1) = event_times(nevent,1)+(nevent)*100000;

end

if isempty(event_label)
    event_label{1} = 'event';
    event_id = ones(length(event_times),1);
end

cluster_spike_id = cell(size(unit_id));
no_cluster = length(unit_id);

track1_ID = find(event_id == 1);
track2_ID = find(event_id == 2);


% Define Gaussian window for smoothing
gaussianWindow = gausswin(0.2*1/psthBinSize);

% Normalize to have an area of 1 (i.e., to be a probability distribution)
gaussianWindow = gaussianWindow / sum(gaussianWindow);
tic
for iCell = 1:no_cluster

    V1_event_modulation(1).cluster_id(iCell) = unit_id(iCell);
    V1_event_modulation(1).region(iCell) = unit_region(iCell);
    V1_event_modulation(1).peak_depth(iCell) = unit_depth(iCell);

    V1_event_modulation(2).cluster_id(iCell) = unit_id(iCell);
    V1_event_modulation(2).region(iCell) = unit_region(iCell);
    V1_event_modulation(2).peak_depth(iCell) = unit_depth(iCell);

    cluster_spike_id{iCell} = spike_id == unit_id(iCell);

    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray1] = psthAndBA(spike_times_events(cluster_spike_id{iCell}),event_times(event_id==1), window, psthBinSize);
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray2] = psthAndBA(spike_times_events(cluster_spike_id{iCell}),event_times(event_id==2), window, psthBinSize);

    psth_track1 = binnedArray1/psthBinSize;
    psth_track1 = conv(mean(psth_track1,'omitnan'), gaussianWindow, 'same');
    psth_track2 = binnedArray2/psthBinSize;
    psth_track2 = conv(mean(psth_track2,'omitnan'), gaussianWindow, 'same');

    shiftedArrays1 = zeros(1000,size(binnedArray1,1),size(binnedArray1,2));
    shiftedArrays2 = zeros(1000,size(binnedArray2,1),size(binnedArray2,2));
    PSTH_shuffled1 = zeros(1000,size(binnedArray2,2));
    PSTH_shuffled2 = zeros(1000,size(binnedArray2,2));

    parfor nshuffle = 1:1000
        %         shiftedArrays1(:,:,nshuffle) = zeros(size(binnedArray1,1),size(binnedArray1,2));
        %         shiftedArrays2(:,:,nshuffle) = zeros(size(binnedArray2,1),size(binnedArray2,2));

        s = RandStream('mrg32k3a','Seed',nshuffle+iCell*1000); % Set random seed for resampling
        shift_values = round(size(binnedArray1,2)*rand(s,1,size(binnedArray1,1)))'; % Generate all shift values at once
        temp = arrayfun(@(nevent) circshift(binnedArray1(nevent,:), shift_values(nevent), 2), 1:size(binnedArray1,1), 'UniformOutput', false);
        shiftedArrays1(nshuffle,:,:) = cell2mat(temp');

        PSTH_shuffled1(nshuffle,:) = squeeze(mean(shiftedArrays1(nshuffle,:,:),'omitnan')/psthBinSize);
        PSTH_shuffled1(nshuffle,:)  = conv(PSTH_shuffled1(nshuffle,:) , gaussianWindow, 'same');

%         FR_difference_shuffled1(nshuffle) = (max(PSTH_shuffled1(nshuffle,:))-min(PSTH_shuffled1(nshuffle,:)))/mean(PSTH_shuffled1(nshuffle,:));
        %
        s = RandStream('mrg32k3a','Seed',1+nshuffle+iCell*1000); % Set random seed for resampling
        shift_values = round(size(binnedArray2,2)*rand(s,1,size(binnedArray2,1)))'; % Generate all shift values at once
        temp = arrayfun(@(nevent) circshift(binnedArray2(nevent,:), shift_values(nevent), 2), 1:size(binnedArray2,1), 'UniformOutput', false);
        shiftedArrays2(nshuffle,:,:) = cell2mat(temp');

        PSTH_shuffled2(nshuffle,:) = squeeze(mean(shiftedArrays2(nshuffle,:,:),'omitnan')/psthBinSize);
        PSTH_shuffled2(nshuffle,:)  = conv(PSTH_shuffled2(nshuffle,:) , gaussianWindow, 'same');
%         FR_difference_shuffled2(nshuffle) = (max(PSTH_shuffled2(nshuffle,:))-min(PSTH_shuffled2(nshuffle,:)))/mean(PSTH_shuffled2(nshuffle,:));
        %
    end

    twin = bins>-1 & bins<1;
    parfor nshuffle = 1:1000
        FR_difference_shuffled1(nshuffle) = (max(PSTH_shuffled1(nshuffle,twin))-min(PSTH_shuffled1(nshuffle,twin)))/mean(PSTH_shuffled1(nshuffle,twin));
        FR_difference_shuffled2(nshuffle) = (max(PSTH_shuffled2(nshuffle,twin))-min(PSTH_shuffled2(nshuffle,twin)))/mean(PSTH_shuffled2(nshuffle,twin));
    end

    V1_event_modulation(1).spike_count{iCell} = binnedArray1;
    V1_event_modulation(2).spike_count{iCell} = binnedArray2;

    V1_event_modulation(1).PSTH{iCell} = psth_track1;
    V1_event_modulation(2).PSTH{iCell} = psth_track2;

    V1_event_modulation(1).PSTH_shuffled{iCell} = PSTH_shuffled1;
    V1_event_modulation(2).PSTH_shuffled{iCell} = PSTH_shuffled2;

    V1_event_modulation(1).max_min_difference(iCell) = max(psth_track1)-min(psth_track1);
    V1_event_modulation(2).max_min_difference(iCell) = max(psth_track2)-min(psth_track2);

    FR_modulation1= (max(psth_track1(bins>-1 & bins<1))-min(psth_track1(bins>-1 & bins<1)))/mean(psth_track1(bins>-1 & bins<1));
    FR_modulation2 = (max(psth_track2(bins>-1 & bins<1))-min(psth_track2(bins>-1 & bins<1)))/mean(psth_track2(bins>-1 & bins<1));

    V1_event_modulation(1).V1_event_modulation(iCell) = FR_modulation1;
    V1_event_modulation(2).V1_event_modulation(iCell) = FR_modulation2;

    V1_event_modulation(1).V1_event_modulation_percentile(iCell) = sum(abs(FR_modulation1)>abs(FR_difference_shuffled1))/length(FR_difference_shuffled1); % abs FR change as a neuron can be activated or inhibited
    V1_event_modulation(2).V1_event_modulation_percentile(iCell) = sum(abs(FR_modulation2)>abs(FR_difference_shuffled2))/length(FR_difference_shuffled2);

    % PRE-V1_event modulation
    temp = psth_track1(bins<0 & bins>-1);
    tmep_shuffled = PSTH_shuffled1(:,bins<0 & bins>-1);
    V1_event_modulation(1).pre_V1_event_activation(iCell) = sum(max(temp)>max(tmep_shuffled'))/length(FR_difference_shuffled1);
    V1_event_modulation(1).pre_V1_event_inhibition(iCell) = sum(min(temp)<min(tmep_shuffled'))/length(FR_difference_shuffled1);

    temp = psth_track2(bins<0 & bins>-1);
    tmep_shuffled = PSTH_shuffled2(:,bins<0 & bins>-1);
    V1_event_modulation(2).pre_V1_event_activation(iCell) = sum(max(temp)>max(tmep_shuffled'))/length(FR_difference_shuffled1);
    V1_event_modulation(2).pre_V1_event_inhibition(iCell) = sum(min(temp)<min(tmep_shuffled'))/length(FR_difference_shuffled1);

    % POST-V1_event modulation
    temp = psth_track1(bins>0 & bins<1);
    tmep_shuffled = PSTH_shuffled1(:,bins>0 & bins<1);
    V1_event_modulation(1).post_V1_event_activation(iCell) = sum(max(temp)>max(tmep_shuffled'))/length(FR_difference_shuffled1);
    V1_event_modulation(1).post_V1_event_inhibition(iCell) = sum(min(temp)<min(tmep_shuffled'))/length(FR_difference_shuffled1);

    temp = psth_track2(bins>0 & bins<1);
    tmep_shuffled = PSTH_shuffled2(:,bins>0 & bins<1);
    V1_event_modulation(2).post_V1_event_activation(iCell) = sum(max(temp)>max(tmep_shuffled'))/length(FR_difference_shuffled1);
    V1_event_modulation(2).post_V1_event_inhibition(iCell) = sum(min(temp)<min(tmep_shuffled'))/length(FR_difference_shuffled1);


end
toc
end