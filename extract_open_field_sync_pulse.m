% Base path to stimuli folder
base_path = 'Z:\ibn-vision\DATA\SUBJECTS\M24072\stimuli';

% Get a list of all date-based subfolders
date_folders = dir(base_path);
date_folders = date_folders([date_folders.isdir]); % Only keep directories

for i = 1:length(date_folders)
    date_name = date_folders(i).name;
    

    
    % Construct the full path to the date folder
    date_path = fullfile(base_path, date_name);
    
    % Find all OpenField videos in this folder
    video_files = dir(fullfile(date_path, '*OpenField*.avi'));
    
    for j = 1:length(video_files)
        video_name = video_files(j).name;
        video_path = fullfile(date_path, video_name);
        
        % Create a VideoReader object

            videoObj = VideoReader(video_path);

        
        fprintf('Processing video: %s\n', video_name);
        
        % Preallocate an array to store sum intensities
        sum_intensities = zeros(videoObj.NumFrames, 1);
        
        % Define crop area (adjust as needed)
        crop_x = 55; % X-coordinate of the top-left corner
        crop_y = 470; % Y-coordinate of the top-left corner
        crop_width = 130;  % Width of the cropped area
        crop_height = 40; % Height of the cropped area
        
        frame_idx = 1;
        
        % Loop through each frame
        while hasFrame(videoObj)
            frame = readFrame(videoObj);
            
            % Convert to grayscale if the frame is in color
            if size(frame, 3) == 3
                frame = rgb2gray(frame);
            end
            
            % Crop the specified area
            cropped_area = frame(crop_y:crop_y+crop_height-1, crop_x:crop_x+crop_width-1);
            
            % Calculate the sum intensity of the cropped area
            sum_intensities(frame_idx) = sum(cropped_area(:));
            frame_idx = frame_idx + 1;
        end
        
        % Save the results as a CSV file
        csv_name = sprintf('%s_sync_pulse.csv', video_name(1:end-4));
        csv_path = fullfile(date_path, csv_name);
        try
            writematrix(sum_intensities, csv_path);
            fprintf('Saved sum intensities to: %s\n', csv_path);
        catch ME
            fprintf('Error saving CSV: %s\n', csv_name);
            disp(ME.message);
        end
    end
end

disp('Processing complete.');
