function polarfill(theta, upper, lower, color, alphaVal)
    % theta in radians
    theta = [theta theta(1)];
    upper = [upper upper(1)];
    lower = [lower lower(1)];
    
    % Convert polar to cartesian
    [xUpper, yUpper] = pol2cart(theta, upper);
    [xLower, yLower] = pol2cart(fliplr(theta), fliplr(lower));

    fill([xUpper xLower], [yUpper yLower], color, ...
        'FaceAlpha', alphaVal, 'EdgeColor', 'none');
end