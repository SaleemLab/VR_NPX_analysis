%% MAIN SCRIPT TO ANALYZE BEHAVIORAL DATA FOR A SINGLE MOUSE
clear; clc;

% 1. DEFINE THE MOUSE TO ANALYZE
SUBJECTS = {'M24016','M23017','M24018','M24062','M24064','M24065'};
for iSub = 1:length(SUBJECTS)
    mouseName = SUBJECTS{iSub}; % <-- Set your mouse ID here

    % 2. GET ALL POTENTIAL SESSION FOLDERS
    fprintf('Searching for session folders for mouse: %s...\n', mouseName);
    all_session_paths = get_mouse_session_paths(mouseName);
    if isempty(all_session_paths)
        fprintf('No session folders found for %s.\n', mouseName);
        return;
    end
    valid_session_paths = all_session_paths;
    fprintf('Found %d sessions to process (will infer trials if Trial_info is missing).\n\n', length(valid_session_paths));


    % 3. INITIALIZE THE FINAL DATA STRUCTURE
    behaviour_all = struct();
    fprintf('Starting data extraction...\n');


    % 4. LOOP THROUGH ALL SESSIONS, EXTRACT DATA, AND DYNAMICALLY AGGREGATE
    num_sessions = length(valid_session_paths);
    for i = 1:num_sessions
        session_path = valid_session_paths{i};
        fprintf('Processing session %d/%d: %s\n', i, num_sessions, session_path);

        session_behaviour = read_behavior_session(session_path);

        if isempty(session_behaviour)
            fprintf('  -> FAILED to extract valid data set from this session. Skipping aggregation.\n');
            % If a session fails, we must still pad any existing fields to maintain alignment
            aggregate_fields = fieldnames(behaviour_all);
            for j = 1:length(aggregate_fields)
                behaviour_all.(aggregate_fields{j}){i, 1} = [];
            end
            continue; % Skip to the next session
        end

        % --- NEW DYNAMIC AGGREGATION LOGIC ---

        % Get all field names encountered so far from both the aggregate and current session
        session_fields = fieldnames(session_behaviour);
        aggregate_fields = fieldnames(behaviour_all);
        all_known_fields = union(session_fields, aggregate_fields);

        % Loop through the complete list of all unique fields
        for j = 1:length(all_known_fields)
            fn = all_known_fields{j};

            % Check if the current session has this field
            if isfield(session_behaviour, fn)
                % If yes, get its data
                data_to_add = session_behaviour.(fn);
            else
                % If no, this session is missing a field that exists in others.
                % Use an empty placeholder to maintain alignment.
                data_to_add = [];
            end

            % Add the data (or placeholder) to the correct session index 'i'.
            % This automatically creates new fields and back-fills with empty cells if needed.
            behaviour_all.(fn){i, 1} = data_to_add;
        end
        fprintf('  -> Successfully extracted and stored data.\n');
    end

    fprintf('\nData processing complete! All data is stored in the "behaviour_all" structure.\n');
    mkdir(fullfile('Z:\ibn-vision\DATA\SUBJECTS\',mouseName,'analysis'))
    save(fullfile('Z:\ibn-vision\DATA\SUBJECTS\',mouseName,'analysis','training_behaviour.mat'),'behaviour_all','-v7.3')
end