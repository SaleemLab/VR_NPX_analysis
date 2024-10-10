function lagged_vector = create_lag(vector, lag)
    % Initialize the lagged vector with NaNs
    lagged_vector = nan(size(vector));

    % If the lag is positive, shift the vector to the right
    if lag > 0
        lagged_vector((lag + 1):end) = vector(1:(end - lag));
    % If the lag is negative, shift the vector to the left
    elseif lag < 0
        lagged_vector(1:(end + lag)) = vector((-lag + 1):end);
    % If the lag is zero, the lagged vector is the same as the original vector
    else
        lagged_vector = vector;
    end
end

