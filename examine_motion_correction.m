%% This script shows you how to look the motion estimated by spikeinterface in your recordings
mouse = 'M25026';
base_folder = 'Z:\ibn-vision\DATA\SUBJECTS';
date = '20250306';


% Define recording folder
recording_folder = fullfile(base_folder, mouse, 'ephys', date);

% Find all subfolders containing 'motion'
motion_dirs = dir(fullfile(recording_folder, '*motion*'));
motion_dirs = motion_dirs([motion_dirs.isdir]);  % Only keep folders

% Loop through each motion folder
for i = 1:length(motion_dirs)
    displacement = readNPY(fullfile(recording_folder, motion_dirs(i).name,'motion/displacement_seg0.npy'));
    spatial_bins = readNPY(fullfile(recording_folder, motion_dirs(i).name,'motion/spatial_bins_um.npy'));
    temporal_bins =  readNPY(fullfile(recording_folder, motion_dirs(i).name,'motion/temporal_bins_s_seg0.npy'));
    

    

    % 3D plot
    figure('Name', motion_dirs(i).name);
    hold on;
    for j = 1:length(spatial_bins)
        plot3(temporal_bins, ...
              spatial_bins(j) * ones(size(temporal_bins)), ...
              displacement(:, j));
    end
    xlabel('Time');
    ylabel('Depth');
    zlabel('Drift');
    title(['Drift trace - ' motion_dirs(i).name]);
    view(3);
    grid on;

end
