function HMM_track_decoding(spatial_rates,spikes_running,spikes_sleep,candidate_events)

%% HMM Pipeline for Decoding Track Representations
% This script builds an HMM from positional and spiking data and optimizes the decoding process using gradient descent.

%% Input Assumptions
% - `positions`: [1 x num_time_bins] animal's position (in cm) during running.
% - `track_ids`: [1 x num_time_bins] corresponding track ID (1 or 2) for each position.
% - `spatial_rates`: [num_bins x num_cells x 2] initial spatial firing rates for each track.
% - `spikes_running`: [num_cells x num_time_bins] spike counts during running.
% - `spikes_sleep`: [num_cells x num_time_bins_sleep] spike counts during sleep.
% - `candidate_events`: [num_events x 2] start and end indices of candidate events in sleep.

%% Parameters
num_bins = size(spatial_rates, 1);
num_cells = size(spatial_rates, 2);
num_tracks = size(spatial_rates, 3);

bin_edges = linspace(0, 140, num_bins + 1);  % Bin edges for discretizing position
learning_rate = 1e-2;                       % Gradient descent learning rate
num_iterations = 50;                        % Number of iterations for optimization
num_shuffles = 1000;                        % Number of shuffles for significance testing

%% Discretize Positions into Bins
position_bins = discretize(positions, bin_edges);

%% Calculate Transition Matrix from Data
trans_mat = zeros(num_bins * num_tracks);

% Count transitions within and between tracks
for t = 1:(length(positions) - 1)
    current_track = track_ids(t);
    next_track = track_ids(t + 1);
    
    current_bin = (current_track - 1) * num_bins + position_bins(t);
    next_bin = (next_track - 1) * num_bins + position_bins(t + 1);
    
    trans_mat(current_bin, next_bin) = trans_mat(current_bin, next_bin) + 1;
end

% Normalize rows of the transition matrix
trans_mat = trans_mat ./ sum(trans_mat, 2);
trans_mat(isnan(trans_mat)) = 0;  % Handle cases where rows sum to zero

%% Initialize Weights for Spatial Firing Rates
weights = spatial_rates;  % Start with provided spatial firing rates

%% Gradient Descent Optimization of Weights
for iter = 1:num_iterations
    log_emissions = zeros(num_bins * num_tracks, size(spikes_running, 2));
    
    % Compute log-emissions for the current weights
    for t = 1:size(spikes_running, 2)
        for track = 1:num_tracks
            for bin = 1:num_bins
                state_idx = (track - 1) * num_bins + bin;
                log_emissions(state_idx, t) = sum(spikes_running(:, t) .* log(weights(bin, :, track) + 1e-12)' - weights(bin, :, track)');
            end
        end
    end
    
    % Decode states using the HMM
    [~, decoded_states] = hmmdecode(log_emissions, trans_mat);
    decoded_bins = mod(decoded_states - 1, num_bins) + 1;
    decoded_tracks = ceil(decoded_states / num_bins);
    
    % Calculate decoding accuracy
    accuracy_track = mean(decoded_tracks == track_ids);
    accuracy_bin = mean(decoded_bins == position_bins);
    fprintf('Iteration %d: Track Accuracy = %.2f%%, Bin Accuracy = %.2f%%\\n', iter, accuracy_track * 100, accuracy_bin * 100);
    
    % Compute gradients and update weights
    gradients = zeros(size(weights));
    for t = 1:size(spikes_running, 2)
        true_bin = position_bins(t);
        true_track = track_ids(t);
        state_idx = (true_track - 1) * num_bins + true_bin;
        
        % Gradient for true state
        gradients(true_bin, :, true_track) = gradients(true_bin, :, true_track) + ...
            (spikes_running(:, t)' ./ (weights(true_bin, :, true_track) + 1e-12)) - 1;
        
        % Subtract contributions from decoded state
        decoded_bin = decoded_bins(t);
        decoded_track = decoded_tracks(t);
        if true_track ~= decoded_track || true_bin ~= decoded_bin
            gradients(decoded_bin, :, decoded_track) = gradients(decoded_bin, :, decoded_track) - ...
                (spikes_running(:, t)' ./ (weights(decoded_bin, :, decoded_track) + 1e-12)) + 1;
        end
    end
    
    % Gradient descent step
    weights = weights + learning_rate * gradients;
    weights = max(weights, 1e-12);  % Prevent weights from going to zero
end

%% Decode Sleep Data
% Compute log-emissions for sleep data
log_emissions_sleep = zeros(num_bins * num_tracks, size(spikes_sleep, 2));
for t = 1:size(spikes_sleep, 2)
    for track = 1:num_tracks
        for bin = 1:num_bins
            state_idx = (track - 1) * num_bins + bin;
            log_emissions_sleep(state_idx, t) = sum(spikes_sleep(:, t) .* log(weights(bin, :, track) + 1e-12)' - weights(bin, :, track)');
        end
    end
end

% Decode candidate events
decoded_events = cell(size(candidate_events, 1), 1);
track_probabilities = zeros(size(candidate_events, 1), num_tracks);
shuffle_track_probs = zeros(size(candidate_events, 1), num_tracks, num_shuffles);

for i = 1:size(candidate_events, 1)
    event_range = candidate_events(i, 1):candidate_events(i, 2);
    [~, decoded_states_event] = hmmdecode(log_emissions_sleep(:, event_range), trans_mat);
    decoded_events{i} = decoded_states_event;
    
    % Compute track probabilities
    state_probs_event = exp(log_emissions_sleep(:, event_range));
    state_probs_event = state_probs_event ./ sum(state_probs_event, 1);
    track_probabilities(i, :) = sum(reshape(sum(state_probs_event, 2), num_bins, num_tracks), 
