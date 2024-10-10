function areas = calculate_ellipse_areas(matrix)
    % Initialize the areas array
    areas = nan(size(matrix, 1), 1);
    
    % Calculate the area of each ellipse
    for i = 1:size(matrix, 1)
        % Get the major and minor diameters
        if ~any(isnan(matrix(i,:)))
            major_diam = matrix(i, 3);
            minor_diam = matrix(i, 4);
            
            % Calculate the area
            areas(i) = pi * (major_diam / 2) * (minor_diam / 2);
        end
    end
end