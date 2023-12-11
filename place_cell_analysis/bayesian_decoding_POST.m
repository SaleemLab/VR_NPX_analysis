function estimated_position = bayesian_decoding_POST(place_fields_BAYESIAN,bayesian_spike_count,position,option,modify,nseed,varargin)
% INPUTS:
   % place fields: matrix or []. Place fields that want to be used as a template in the decoding. If empty, it will load place fields for the whole session (extracted_place_fields_BAYESIAN.mat)
   % Save_option: enter 'Y' for saving. Else, won't save.
% OUTPUTS:
    % estimated_position structure. 
% Loads: 'extracted_position','extracted_place_fields_BAYESIAN','bayesian_spike_count'
% This code is modified for NPX 2 probe analysis. Things will change
% 26/09/2023 Masa place_field_index changed to good place cells LIBERAL


parameters=list_of_parameters;
if isempty(varargin(1))
    replay_bin_width = parameters.replay_bin_width;
    run_bin_width = parameters.run_bin_width;
else
    replay_bin_width = varargin{1};
    run_bin_width = varargin{1};
end 

BAYSESIAN_NORMALIZED_ACROSS_TRACKS=1;  %normalized across all tracks if set to 1, otherwise individual tracks


% load extracted_position
place_field_index = [];
place_field_index = place_fields_BAYESIAN.good_place_cells_LIBERAL;

if strcmp(bayesian_spike_count, 'replayEvents_bayesian_spike_count')
    load replayEvents_bayesian_spike_count  %if decoding replay events (which will have different time bins)
    
    if isempty(place_fields_BAYESIAN)
        load extracted_place_fields_BAYESIAN
    end
    
%     bayesian_spike_count = replayEvents_bayesian_spike_count;
%     place_field_index = place_fields_BAYESIAN.good_place_cells_LIBERAL;
elseif strcmp(bayesian_spike_count, 'replayEvents_common_good_cell_bayesian_spike_count')
    load replayEvents_common_good_cell_bayesian_spike_count
    bayesian_spike_count = replayEvents_common_good_cell_bayesian_spike_count;
    load 'rate_remapping_analysis_TRACK_PAIRS_wcorr'
    place_field_index = remapping(1).common_good_cells  
elseif strcmp(bayesian_spike_count, 'replayEvents_not_common_cell_bayesian_spike_count')
    load replayEvents_not_common_cell_bayesian_spike_count
    bayesian_spike_count = replayEvents_not_common_cell_bayesian_spike_count;
    load 'rate_remapping_analysis_TRACK_PAIRS_wcorr'
    load 'extracted_place_fields_BAYESIAN';
    place_field_index = setdiff(place_fields_BAYESIAN.good_place_cells,remapping(1).common_good_cells);% Find cells that are not good place cells

elseif isempty(bayesian_spike_count) %if decoding the whole session
    load bayesian_spike_count
    place_field_index = []
end

% If not planning to use only place fields from a specific period of time, load place fields from the whole session
if isempty(place_fields_BAYESIAN)
    load extracted_place_fields_BAYESIAN
end

% Find ratemaps
if isempty(place_field_index)
    if isfield(place_fields_BAYESIAN,'good_place_cells')
        place_field_index = place_fields_BAYESIAN.good_place_cells;  %use all place cells
    elseif ~isfield(place_fields_BAYESIAN,'good_place_cells') && isfield(place_fields_BAYESIAN.track,'good_cells')
        place_field_index = place_fields_BAYESIAN.track.good_cells;  % when using output from get_place_fields_laps
    else
        disp('ERROR- field good_place_cells missing');
        place_field_index = place_fields_BAYESIAN.pyramidal_cells;
    end
end

% Original place cell ratemap
all_place_fields = [];
all_place_fields_ratemap_shuffled = [];
for track_id = 1:length(place_fields_BAYESIAN.track)
    if isempty(modify)

        for k=1:length(place_field_index)
            single_place_field = place_fields_BAYESIAN.track(track_id).raw{place_field_index(k)}; %get raw place field
            single_place_field(find(isnan(single_place_field))) = 0; % remove NaNs in place field and replace by 0
            % single_place_field = smooth(single_place_field,parameters.smoothing_number_of_bins); %smooth place field
            if min(single_place_field)<0
                disp('error- spike rate of place field less than zero')
            end
            all_place_fields{track_id}(k,:) = single_place_field;
        end
    elseif strcmp(modify,'global_remapped')
%     else
        % for global remapping, get place fields for each event
        for event = 1:length(bayesian_spike_count.replay_events)
            for k=1:length(place_field_index)
                single_place_field = place_fields_BAYESIAN.track(track_id).(sprintf('%s',modify)){event}{place_field_index(k)}; %get raw place field
                single_place_field(find(isnan(single_place_field))) = 0; % remove NaNs in place field and replace by 0
                % single_place_field = smooth(single_place_field,parameters.smoothing_number_of_bins); %smooth place field
                if min(single_place_field)<0
                    disp('error- spike rate of place field less than zero')
                end
                all_place_fields{track_id}{event}(k,:) = single_place_field;
            end
        end
    end
end

if isempty(nseed)
    nseed = 1;
end

% place cell ratemap track id shuffled (randomly swapped)
if strcmp(modify,'global_remapped')
    for event = 1:length(all_place_fields{1})
        for k = 1:length(place_field_index) % Turning curve shuffled for each cell (across two tracks)
            s = RandStream('mrg32k3a','Seed',k+1000*nseed); % Set random seed for resampling
            track_id = randperm(s,2);
            %     track_id = randperm(2);
            all_place_fields_ratemap_shuffled{1}{event}(k,:) = all_place_fields{track_id(1)}{event}(k,:);
            all_place_fields_ratemap_shuffled{2}{event}(k,:) = all_place_fields{track_id(2)}{event}(k,:);
        end
    end
else
    for k = 1:length(place_field_index) % Turning curve shuffled for each cell (across two tracks)
        s = RandStream('mrg32k3a','Seed',k+1000*nseed); % Set random seed for resampling
        track_id = randperm(s,2);
        %     track_id = randperm(2);
        all_place_fields_ratemap_shuffled{1}(k,:) = all_place_fields{track_id(1)}(k,:);
        all_place_fields_ratemap_shuffled{2}(k,:) = all_place_fields{track_id(2)}(k,:);
    end
end


%%%%% BAYESIAN DECODING  %%%%%%
for track_id=1:length(place_fields_BAYESIAN.track)
    % Creates a vector of position bins and finds centre
    position_bin_edges = place_fields_BAYESIAN.track(track_id).x_bin_edges;
    estimated_position(track_id).position_bin_centres = place_fields_BAYESIAN.track(track_id).x_bin_centres;
    
    
    if isfield(bayesian_spike_count,'replay_events')   % When running replay events separately
        estimated_position(track_id).replay_events = bayesian_spike_count.replay_events;
        for i = 1 : length(bayesian_spike_count.replay_events)
            estimated_position(track_id).replay_events(i).replay = zeros(length(estimated_position(track_id).position_bin_centres),length(estimated_position(track_id).replay_events(i).replay_time_centered));
        end
        n.replay = bayesian_spike_count.n.replay;
        estimated_position(track_id).replay = zeros(length(estimated_position(track_id).position_bin_centres),length(n.replay));
    else    %when computing decoded position for entire experiment
        %time bins for replay
        estimated_position(track_id).replay_time_edges = bayesian_spike_count.replay_time_edges;
        estimated_position(track_id).replay_time_centered = bayesian_spike_count.replay_time_centered; %centres of bins
        estimated_position(track_id).replay = zeros(length(estimated_position(track_id).position_bin_centres),length(estimated_position(track_id).replay_time_centered));
        n.replay = bayesian_spike_count.n.replay;
    end
    
    % Apply formula of bayesian decoding
    if ~isfield(bayesian_spike_count,'run_time_edges')
        replay_id = bayesian_spike_count.replay_events_indices;
        if ~isempty(modify)
            if isempty(option)
                estimated_position(track_id).replay = modified_reconstruct(n.replay,all_place_fields{track_id},replay_id,parameters.replay_bin_width,'per event');
            elseif strcmp(option,'ratemap shuffle')
                estimated_position(track_id).replay = modified_reconstruct(n.replay,all_place_fields_ratemap_shuffled{track_id},replay_id,parameters.replay_bin_width,'per event');
            end
        else
            if isempty(option)
                estimated_position(track_id).replay = reconstruct(n.replay,all_place_fields{track_id},parameters.replay_bin_width);
            elseif strcmp(option,'ratemap shuffle') % Track 1 and Track 2 ratemap label randomly shuffled within each cell.
                estimated_position(track_id).replay = reconstruct(n.replay,all_place_fields_ratemap_shuffled{track_id},parameters.replay_bin_width);
            end
        end
    end
end

%%%%%% NORMALIZING  %%%%%%%
%       columns need to sum to 1 (total probability across positions.
%       options are normalizing across tracks or just within a track
%
summed_probability_replay = zeros(1,size(estimated_position(1).replay,2));
if isfield(bayesian_spike_count,'run_time_edges')
    summed_probability_run = zeros(1,size(estimated_position(1).run,2));
end

for track_id=1:length(place_fields_BAYESIAN.track)     %normalize probability to sum to 1
    % Sum probabilties across rows (cells)
    estimated_position(track_id).replay_OneTrack = NaN(size(estimated_position(track_id).replay));
    estimated_position(track_id).replay_Normalized = NaN(size(estimated_position(track_id).replay));
    summed_probability_replay=summed_probability_replay+sum(estimated_position(track_id).replay,1);
    if isfield(bayesian_spike_count,'run_time_edges')
        estimated_position(track_id).run_OneTrack = NaN(size(estimated_position(track_id).run));
        summed_probability_run=summed_probability_run+sum(estimated_position(track_id).run,1);
    end
end
for track_id=1:length(place_fields_BAYESIAN.track)
    summed_probability(track_id).replay = summed_probability_replay;
    if isfield(bayesian_spike_count,'run_time_edges')
        summed_probability(track_id).run = summed_probability_run;
    end
end

% Divide decoded position by summed probability  (normalize)
for track_id=1:length(place_fields_BAYESIAN.track)
    for j=1:size(estimated_position(track_id).replay,2)
        estimated_position(track_id).replay_OneTrack(:,j) = estimated_position(track_id).replay(:,j)./sum(estimated_position(track_id).replay(:,j)); %normalized by one track only
        estimated_position(track_id).replay_Normalized(:,j) = estimated_position(track_id).replay(:,j)./summed_probability(track_id).replay(j); % normalized by the sum of prob of all tracks
        if BAYSESIAN_NORMALIZED_ACROSS_TRACKS==1
       estimated_position(track_id).replay(:,j) = estimated_position(track_id).replay_Normalized(:,j);
        else
        estimated_position(track_id).replay(:,j) = estimated_position(track_id).replay_OneTrack(:,j);
        end
    end
    % Calculate replay bias  - measures which track has higher probability values for estimated positions
    estimated_position(track_id).replay_bias=sum(estimated_position(track_id).replay_Normalized,1);
   
end

% If running replay events separately, also extract individual events from the estimated position matrix
if isfield(bayesian_spike_count,'replay_events')
    for track_id = 1:length(place_fields_BAYESIAN.track)
        for event = 1 : length(bayesian_spike_count.replay_events)
            thisReplay_indxs = find(bayesian_spike_count.replay_events_indices == event);
            estimated_position(track_id).replay_events(event).replay = estimated_position(track_id).replay(:,thisReplay_indxs);
%             estimated_position(track_id).replay_events(event).replay_actual_position = interp1(position.t, estimated_position(track_id).discrete_position, estimated_position(track_id).replay_events(event).replay_time_centered, 'nearest');
            if track_id == 1
                estimated_position(track_id).replay_events(event).probability_ratio = sum(sum(estimated_position(1).replay(:,thisReplay_indxs)))/sum(sum(estimated_position(2).replay(:,thisReplay_indxs)));
            elseif track_id == 2
                estimated_position(track_id).replay_events(event).probability_ratio = sum(sum(estimated_position(2).replay(:,thisReplay_indxs)))/sum(sum(estimated_position(1).replay(:,thisReplay_indxs)));
            end
        end
    end
end

% if strcmp(save_option, 'Y')
%     save estimated_position estimated_position
% end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function estimated_position = reconstruct(n,all_place_fields,bin_width)
% Creates matrix where rows are cells and columns are position bins
bin_length = size(all_place_fields,2); %columns
number_of_cells = size(all_place_fields,1); %rows
parameters.bayesian_threshold=10.^(log2(number_of_cells)-log2(400)); % small value multiplied to all values to get rid of zeros
all_place_fields(find(all_place_fields<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
sum_of_place_fields = sum(all_place_fields,1);  % adds up spikes per bin (used later for exponential)
for j = 1: size(n,2)
    n_spikes = n(:,j)*ones(1,bin_length); %number of spikes in time bin
    pre_product = all_place_fields.^n_spikes; % pl field values raised to num of spikes
    pre_product(find(pre_product<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
    product_of_place_fields = prod(pre_product,1); %product of pl fields
    estimated_position(:,j) = product_of_place_fields.*(exp(-bin_width*sum_of_place_fields)); % bayesian formula
    %NOTE- columns do not sum to 1.  this is done at a later stage to allow normalization within a track or across tracks
end
end



function estimated_position = modified_reconstruct(n,all_place_fields,replay_id,bin_width,option)

if strcmp(option,'per event')

    for event = 1 : length(all_place_fields)
        all_place_fields_this_event = all_place_fields{event};
        thisReplay_indxs = find(replay_id == event);

        % Creates matrix where rows are cells and columns are position bins
        bin_length = size(all_place_fields_this_event,2); %columns (position bins)
        number_of_cells = size(all_place_fields_this_event,1); %rows (cells)
        parameters.bayesian_threshold=10.^(log2(number_of_cells)-log2(400)); % small value multiplied to all values to get rid of zeros
        all_place_fields_this_event(find(all_place_fields_this_event<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
        sum_of_place_fields = sum(all_place_fields_this_event,1);  % adds up spikes per posiiton bin (used later for exponential)

        for j = 1: length(thisReplay_indxs)
            n_spikes = n(:,thisReplay_indxs(j))*ones(1,bin_length); %number of spikes in time bin
            %             n_spikes = n_spikes(active_cell_index,:);
            pre_product = all_place_fields_this_event.^n_spikes; % pl field values raised to num of spikes
            pre_product(find(pre_product<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
            product_of_place_fields = prod(pre_product,1); %product of pl fields
            if length(bin_width) == 1 % if specifying only one universla time bin (i.e 20ms) for all replay events
                estimated_position(:,thisReplay_indxs(j)) = product_of_place_fields.*(exp(-bin_width*sum_of_place_fields)); % bayesian formula
                %NOTE- columns do not sum to 1.  this is done at a later stage to allow normalization within a track or across tracks
            else
                for m = 1: length(bin_width{event}) % if different time bin for each replay event
                    estimated_position(:,thisReplay_indxs(j)) = product_of_place_fields.*(exp(-bin_width{event}(m)*sum_of_place_fields)); % bayesian formula
                end
            end
        end

        thisReplay_indxs = [];
    end

else
    % Creates matrix where rows are cells and columns are position bins
    bin_length = size(all_place_fields,2); %columns (position bins)
    number_of_cells = size(all_place_fields,1); %rows (cells)
    parameters.bayesian_threshold=10.^(log2(number_of_cells)-log2(400)); % small value multiplied to all values to get rid of zeros
    all_place_fields(find(all_place_fields<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
    sum_of_place_fields = sum(all_place_fields,1);  % adds up spikes per posiiton bin (used later for exponential)

    for event = 1 : length(decoded_replay_events(1).replay_events)
        thisReplay_indxs = find(replay_id == event);


        for j = 1: length(thisReplay_indxs)
            n_spikes = n(:,thisReplay_indxs(j))*ones(1,bin_length); %number of spikes in time bin
            %             n_spikes = n_spikes(active_cell_index,:);
            pre_product = all_place_fields.^n_spikes; % pl field values raised to num of spikes
            pre_product(find(pre_product<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
            product_of_place_fields = prod(pre_product,1); %product of pl fields
            if length(bin_width) == 1 % if specifying only one universla time bin (i.e 20ms) for all replay events
                estimated_position(:,thisReplay_indxs(j)) = product_of_place_fields.*(exp(-bin_width*sum_of_place_fields)); % bayesian formula
                %NOTE- columns do not sum to 1.  this is done at a later stage to allow normalization within a track or across tracks
            else
                for m = 1: length(bin_width{event}) % if different time bin for each replay event
                    estimated_position(:,thisReplay_indxs(j)) = product_of_place_fields.*(exp(-bin_width{event}(m)*sum_of_place_fields)); % bayesian formula
                end
            end
        end

        thisReplay_indxs = [];
    end
end
end