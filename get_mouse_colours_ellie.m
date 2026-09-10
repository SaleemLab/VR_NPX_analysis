function subject_color = get_mouse_colours_ellie()

% Assigns a consistent marker colour to each mouse across scripts
% Output:
%   subject_color - containers.Map mapping mouse IDs to RGB colours


    all_mice = {'M00013', 'M00014', 'M00069', 'M00071', 'M00087', ...
                'M00088', 'M00096', 'M00098', 'M00100', 'M00113'};

    cmap = [
        0.000 0.447 0.741
        0.850 0.325 0.098
        0.929 0.694 0.125
        0.494 0.184 0.556
        0.466 0.674 0.188
        0.301 0.745 0.933
        0.635 0.078 0.184
        0.000 0.600 0.500
        0.750 0.400 0.700
        0.400 0.400 0.400
    ];

    subject_color = containers.Map(all_mice, num2cell(cmap,2));

end

