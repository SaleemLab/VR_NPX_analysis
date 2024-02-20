function estimated_position = bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,Behaviour,option,modify,nseed,varargin)
% INPUTS:
   % place fields: matrix or []. Place fields that want to be used as a template in the decoding. If empty, it will load place fields for the whole session (extracted_place_fields_BAYESIAN.mat)
   % Save_option: enter 'Y' for saving. Else, won't save.
% OUTPUTS:
    % estimated_position structure. 
% Loads: 'extracted_position','extracted_place_fields_BAYESIAN','bayesian_spike_count'
% This code is modified for NPX 2 probe analysis. Things will change
% 26/09/2023 Masa place_field_index changed to good place cells LIBERAL

parameters = list_of_parameters;
if isempty(varargin{1})
    replay_bin_width = parameters.replay_bin_width;
    run_bin_width = parameters.run_bin_width;
else
    replay_bin_width = varargin{1};
    run_bin_width = varargin{1};
end

BAYSESIAN_NORMALIZED_ACROSS_TRACKS=1;  %normalized across all tracks if set to 1, otherwise individual tracks


% load extracted_position
place_field_index = [];
place_field_index = place_fields_BAYESIAN(1).good_cells;

% Original place cell ratemap
all_place_fields = [];
all_place_fields_ratemap_shuffled = [];
for track_id = 1:length(place_fields_BAYESIAN)
    if isempty(modify)

        all_place_fields{track_id} = place_fields_BAYESIAN(track_id).template;
        
%        [~,cell_id] = max(place_fields_BAYESIAN(track_id).template');
%        [~,cell_id] = sort(cell_id);
% imagesc(place_fields_BAYESIAN(track_id).template(cell_id,:))
    elseif strcmp(modify,'global_remapped')
%     else
        % for global remapping, get place fields for each event
        for event = 1:length(bayesian_spike_count.replay_events)
            s = RandStream('mrg32k3a','Seed',track_id+100*event); % Set random seed for resampling
            random_cell_index = randperm(s,length(place_field_index));
%             random_cell = place_fields_BAYESIAN.track(track_id).sorted_good_cells_LIBERAL(random_cell_index);

            all_place_fields{track_id}{event} = place_fields_BAYESIAN(track_id).template(random_cell_index,:);
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

position = []
%%%%% BAYESIAN DECODING  %%%%%%
for track_id=1:length(place_fields_BAYESIAN)
    % Creates a vector of position bins and finds centre
    position_bin_edges = place_fields_BAYESIAN(track_id).x_bin_edges;
    estimated_position(track_id).position_bin_centres = place_fields_BAYESIAN(track_id).x_bin_centres;
    
    % Bin position for decoding error
    estimated_position(track_id).discrete_position = NaN(size(Behaviour.position));
    discrete_position = Behaviour.position;

    estimated_position(track_id).discrete_position(Behaviour.track_ID==track_id) = Behaviour.position(Behaviour.track_ID==track_id);

    discrete_position = discretize(discrete_position,position_bin_edges); %group position points in bins delimited by edges
    estimated_position(track_id).discrete_position(Behaviour.track_ID==track_id) = estimated_position(track_id).position_bin_centres(discrete_position(Behaviour.track_ID==track_id)); %creates new positions based on centre of bins
    
    if isfield(bayesian_spike_count,'replay_events')   % When running replay events separately
        estimated_position(track_id).replay_events = bayesian_spike_count.replay_events;
        for i = 1 : length(bayesian_spike_count.replay_events)
            estimated_position(track_id).replay_events(i).replay = zeros(length(estimated_position(track_id).position_bin_centres),length(estimated_position(track_id).replay_events(i).replay_time_centered));
        end
        n.replay = bayesian_spike_count.n.replay;
        estimated_position(track_id).replay = zeros(length(estimated_position(track_id).position_bin_centres),length(n.replay));

    elseif isfield(bayesian_spike_count,'laps')   % When running replay events separately
        estimated_position(track_id).laps = bayesian_spike_count.laps;
        for i = 1 : length(bayesian_spike_count.laps)
            estimated_position(track_id).laps(i).run = zeros(length(estimated_position(track_id).position_bin_centres),length(estimated_position(track_id).laps(i).run_time_centered));
        end
        n.run = bayesian_spike_count.n.run;
        estimated_position(track_id).run = zeros(length(estimated_position(track_id).position_bin_centres),length(n.run));
    else    %when computing decoded position for entire experiment
        %time bins for replay
        estimated_position(track_id).replay_time_edges = bayesian_spike_count.replay_time_edges;
        estimated_position(track_id).replay_time_centered = bayesian_spike_count.replay_time_centered; %centres of bins
        estimated_position(track_id).replay = zeros(length(estimated_position(track_id).position_bin_centres),length(estimated_position(track_id).replay_time_centered));
        n.replay = bayesian_spike_count.n.replay;
        %time bins for run
        if isfield(bayesian_spike_count,'run_time_edges')
            estimated_position(track_id).run_time_edges = bayesian_spike_count.run_time_edges;
            estimated_position(track_id).run_time_centered =  bayesian_spike_count.run_time_centered;%centres of bins
            estimated_position(track_id).run = zeros(length(estimated_position(track_id).position_bin_centres),length(estimated_position(track_id).run_time_centered));
            n.run = bayesian_spike_count.n.run;
        end
    end
    
    % Apply formula of bayesian decoding
    if isfield(bayesian_spike_count,'replay_events') 
        replay_id = bayesian_spike_count.replay_events_indices;
        if ~isempty(modify)
            if isempty(option)
                estimated_position(track_id).replay = modified_reconstruct(n.replay,all_place_fields{track_id},replay_id,replay_bin_width,'per event');
            elseif strcmp(option,'ratemap shuffle')
                estimated_position(track_id).replay = modified_reconstruct(n.replay,all_place_fields_ratemap_shuffled{track_id},replay_id,replay_bin_width,'per event');
            end
        else
            if isempty(option)
                estimated_position(track_id).replay = reconstruct(n.replay,all_place_fields{track_id},replay_bin_width);
            elseif strcmp(option,'ratemap shuffle') % Track 1 and Track 2 ratemap label randomly shuffled within each cell.
                estimated_position(track_id).replay = reconstruct(n.replay,all_place_fields_ratemap_shuffled{track_id},replay_bin_width);
            end
        end
    end

    % RUN each lap
    if isfield(bayesian_spike_count,'laps') 
        lap_id = bayesian_spike_count.lap_indices;
        if ~isempty(modify)
            if isempty(option)
                estimated_position(track_id).run = modified_reconstruct(n.run,all_place_fields{track_id},lap_id,run_bin_width,'per event');
            elseif strcmp(option,'ratemap shuffle')
                estimated_position(track_id).run = modified_reconstruct(n.run,all_place_fields_ratemap_shuffled{track_id},lap_id,run_bin_width,'per event');
            end
        else
            if isempty(option)
                estimated_position(track_id).run = reconstruct(n.run,all_place_fields{track_id},run_bin_width);
            elseif strcmp(option,'ratemap shuffle') % Track 1 and Track 2 ratemap label randomly shuffled within each cell.
                estimated_position(track_id).run = reconstruct(n.run,all_place_fields_ratemap_shuffled{track_id},run_bin_width);
            end
        end
    end

    if isfield(bayesian_spike_count,'run_time_edges')
        if isempty(option)
           estimated_position(track_id).run = reconstruct(n.run,all_place_fields{track_id},run_bin_width);
        elseif strcmp(option,'ratemap shuffle') % Track 1 and Track 2 ratemap label randomly shuffled within each cell.
            estimated_position(track_id).run = reconstruct(n.run,all_place_fields_ratemap_shuffled{track_id},run_bin_width);
        end

    end
end

%%%%%% NORMALIZING  %%%%%%%
%       columns need to sum to 1 (total probability across positions.
%       options are normalizing across tracks or just within a track
%
if isfield(bayesian_spike_count,'replay_events')
    summed_probability_replay = zeros(1,size(estimated_position(1).replay,2));
end

if isfield(bayesian_spike_count,'run_time_edges') | isfield(bayesian_spike_count,'laps')
    summed_probability_run = zeros(1,size(estimated_position(1).run,2));
end

for track_id=1:length(place_fields_BAYESIAN)     %normalize probability to sum to 1
    if isfield(bayesian_spike_count,'replay_events')
        % Sum probabilties across rows (cells)
        estimated_position(track_id).replay_OneTrack = NaN(size(estimated_position(track_id).replay));
        estimated_position(track_id).replay_Normalized = NaN(size(estimated_position(track_id).replay));
        summed_probability_replay=summed_probability_replay+sum(estimated_position(track_id).replay,1);
    end
    if isfield(bayesian_spike_count,'run_time_edges') | isfield(bayesian_spike_count,'laps')
        estimated_position(track_id).run_OneTrack = NaN(size(estimated_position(track_id).run));
        estimated_position(track_id).run_OneTrack = NaN(size(estimated_position(track_id).run));
        summed_probability_run=summed_probability_run+sum(estimated_position(track_id).run,1);
    end
end

for track_id=1:length(place_fields_BAYESIAN)
    if isfield(bayesian_spike_count,'replay_events')
        summed_probability(track_id).replay = summed_probability_replay;
    end
    if isfield(bayesian_spike_count,'run_time_edges')| isfield(bayesian_spike_count,'laps')
        summed_probability(track_id).run = summed_probability_run;
    end
end

% Divide decoded position by summed probability  (normalize)
for track_id=1:length(place_fields_BAYESIAN)
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
    
    if isfield(bayesian_spike_count,'run_time_edges')|isfield(bayesian_spike_count,'laps')
        for j=1:size(estimated_position(track_id).run,2)
            estimated_position(track_id).run_OneTrack(:,j) = estimated_position(track_id).run(:,j)./sum(estimated_position(track_id).run(:,j));
            estimated_position(track_id).run(:,j) = estimated_position(track_id).run(:,j)./summed_probability(track_id).run(j);
        end
        [estimated_position(track_id).max_prob,index] = max(estimated_position(track_id).run,[],1);
        estimated_position(track_id).peak_position = NaN(size(index));
        valid_bins = find(~isnan(index));
        estimated_position(track_id).peak_position(valid_bins) = estimated_position(track_id).position_bin_centres(index(valid_bins));  %only compute estimated position with peak probability for valid bins (leave as NaN for other bins)
        estimated_position(track_id).run_bias = sum(estimated_position(track_id).run,1);
        estimated_position(track_id).run_error = abs(estimated_position(track_id).peak_position-interp1(position.t, estimated_position(track_id).discrete_position, estimated_position(track_id).run_time_centered, 'nearest'));
        estimated_position(track_id).run_actual_position = interp1(position.t, estimated_position(track_id).discrete_position, estimated_position(track_id).run_time_centered, 'nearest');
        estimated_position(track_id).actual_run_speed = interp1(position.t, position.v, estimated_position(track_id).run_time_centered, 'nearest');
        estimated_position(track_id).run_speed = interp1(position.t, [0 diff(position.linear(track_id).linear)./diff(position.t)], estimated_position(track_id).run_time_centered, 'nearest'); % Based on movement in VR space (when reaching end, speed becomes zero but can still be running)

    end
end

% Calculate probability ratio for RUN
if isfield(bayesian_spike_count,'run_time_edges')
    for track_id=1:length(place_fields_BAYESIAN)
        % Only calculate probability ratio when animal is moving through the track;
        speed_thresholded = estimated_position(track_id).run_speed > 5;
        if track_id == 1
            estimated_position(track_id).probability_ratio = sum(sum(estimated_position(1).run(speed_thresholded)))/sum(sum(estimated_position(2).run(speed_thresholded)));
        elseif track_id == 2
            estimated_position(track_id).probability_ratio = sum(sum(estimated_position(2).run(speed_thresholded)))/sum(sum(estimated_position(1).run(speed_thresholded)));
        end
    end
end

% If running replay events separately, also extract individual events from the estimated position matrix
if isfield(bayesian_spike_count,'replay_events')
    for track_id = 1:length(place_fields_BAYESIAN)
        for event = 1 : length(bayesian_spike_count.replay_events)
            thisReplay_indxs = find(bayesian_spike_count.replay_events_indices == event);
            estimated_position(track_id).replay_events(event).replay = estimated_position(track_id).replay(:,thisReplay_indxs);
            estimated_position(track_id).replay_events(event).replay_actual_position = interp1(position.t, estimated_position(track_id).discrete_position, estimated_position(track_id).replay_events(event).replay_time_centered, 'nearest');
            if track_id == 1
                estimated_position(track_id).replay_events(event).probability_ratio = sum(sum(estimated_position(1).replay(:,thisReplay_indxs)))/sum(sum(estimated_position(2).replay(:,thisReplay_indxs)));
            elseif track_id == 2
                estimated_position(track_id).replay_events(event).probability_ratio = sum(sum(estimated_position(2).replay(:,thisReplay_indxs)))/sum(sum(estimated_position(1).replay(:,thisReplay_indxs)));
            end
        end
    end
end


if isfield(bayesian_spike_count,'laps')
    for track_id=1:length(place_fields_BAYESIAN)
        % Only calculate probability ratio when animal is moving through the track;
        speed_thresholded = estimated_position(track_id).run_speed > 5;
        if track_id == 1
            estimated_position(track_id).probability_ratio = sum(sum(estimated_position(1).run(speed_thresholded)))/sum(sum(estimated_position(2).run(speed_thresholded)));
        elseif track_id == 2
            estimated_position(track_id).probability_ratio = sum(sum(estimated_position(2).run(speed_thresholded)))/sum(sum(estimated_position(1).run(speed_thresholded)));
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