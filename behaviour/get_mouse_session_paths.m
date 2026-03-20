function sessionPaths = get_mouse_session_paths(mouseName)
%get_mouse_session_paths Finds and sorts all training session paths for a given mouse.
%
%   sessionPaths = get_mouse_session_paths(mouseName)
%
%   INPUT:
%       mouseName - A string containing the subject's ID (e.g., 'M23032').
%
%   OUTPUT:
%       sessionPaths - A cell array of full paths to each session folder,
%                      sorted chronologically from oldest to newest.
%
    % Define the base path to the subject's data repository
    baseDataPath = 'Z:\ibn-vision\DATA\SUBJECTS';
    mouseTrainingPath = fullfile(baseDataPath, mouseName, 'training');

    % Define the potential parent directories for the date folders
    searchPaths = { ...
        mouseTrainingPath, ...                 % Handles cases without a 'room' subfolder
        fullfile(mouseTrainingPath, 'matrix'), ... % Handles the 'matrix' room
        fullfile(mouseTrainingPath, 'dome') ...    % Handles the 'dome' room
    };

    % Initialize containers to store the collected paths and their corresponding dates
    allPaths = {};
    
    % --- FIX ---
    % Initialize as an empty DATETIME array instead of a generic numeric array.
    allDates = datetime.empty(0,1);

    % Loop through each possible search path
    for i = 1:length(searchPaths)
        currentSearchPath = searchPaths{i};

        % Proceed only if the directory actually exists
        if exist(currentSearchPath, 'dir')
            
            % Get all contents and filter for directories only (excluding '.' and '..')
            dirContents = dir(currentSearchPath);
            subFolders = dirContents([dirContents.isdir]);
            subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));

            % Process each found subfolder
            for j = 1:length(subFolders)
                folderName = subFolders(j).name;
                
                dateObj = []; % Variable to hold the parsed date
                
                % Check if the folder name matches an 8-digit (YYYYMMDD) or 6-digit (YYMMDD) pattern
                if ~isempty(regexp(folderName, '^\d{8}$', 'once'))
                    try
                        dateObj = datetime(folderName, 'InputFormat', 'yyyyMMdd');
                    catch
                    end
                elseif ~isempty(regexp(folderName, '^\d{6}$', 'once'))
                    try
                        dateObj = datetime(folderName, 'InputFormat', 'yyMMdd');
                    catch
                    end
                end
                
                % If a valid date was parsed, store the full path and the datetime object
                if ~isempty(dateObj)
                    allPaths{end+1, 1} = fullfile(currentSearchPath, folderName);
                    % This assignment now works correctly
                    allDates(end+1, 1) = dateObj;
                end
            end
        end
    end

    % Sort the collected paths chronologically. This works directly on datetime arrays.
    [~, sortIdx] = sort(allDates);
    sessionPaths = allPaths(sortIdx);

end