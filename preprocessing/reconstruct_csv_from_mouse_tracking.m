% MATLAB Code: Mouse Tracking and Data Extraction in case video is saved
% but csv is lost
clear all
% Define the video file
videoFile = 'SleepChronic_2024-11-28T10_26_49.avi';
% videoFile = 'SleepChronic_PRE_2024-11-28T08_53_19.avi';
video = VideoReader(videoFile);

% Define the regions of interest
rectRegion = [50, 0, video.Width-50, video.Height]; % Region for pixel change
ledRegion = [0, 420, 30, 60]; % Region for LED tracking

% Initialize variables
prevFrame = []; % To store the previous frame for pixel difference
frameCount = 0; % Frame counter
data = []; % Array to store results

% Frame rate and timestamps
frameRate = 60; % Video object frame rate is wrong as it assumes 30Hz even it should be 60Hz
timeStep = 1 / frameRate; % Time between frames (in seconds)
background = [];
% Loop through all frames in the video
tic
while hasFrame(video)
    % Read the current frame
    %     for c = 1:20000
    frame = readFrame(video);
    frameCount = frameCount + 1;

    % Convert the frame to grayscale
    grayFrame = rgb2gray(frame);
    if frameCount == 1
        background(:,:,frameCount) = imcrop(grayFrame, rectRegion);
        meanBackground = uint8(squeeze(background));
    elseif frameCount <=30
        background(:,:,frameCount) = imcrop(grayFrame, rectRegion);
        meanBackground = uint8(round(mean(background,3)));
    end

    % Extract the regions of interest
    rectRegionFrame = imcrop(grayFrame, rectRegion)-meanBackground;
    ledRegionFrame = imcrop(grayFrame, ledRegion);

    % Calculate the sum of pixel change (relative to previous frame)
    if frameCount > 30
        foregroundMask = abs(rectRegionFrame - prevFrame)>40;
        pixelChangeSum = 255*sum(sum(foregroundMask));
    else
        pixelChangeSum = 0; % No change for the first frame
    end
    prevFrame = rectRegionFrame;

    % Calculate the sum of pixel values in the LED region
    ledSum = sum(ledRegionFrame(:));

    % Threshold the frame to create a binary image for mouse detection
    %     binaryFrame = imbinarize(rectRegionFrame, 'adaptive');
    binaryFrame = rectRegionFrame<80;
    % Find the largest connected component (assumed to be the mouse)
    stats = regionprops(binaryFrame, 'Centroid', 'Area');

    if ~isempty(stats)
        % Find the region with the largest area
        [~, idx] = max([stats.Area]);
        if stats(idx).Area > 2000
            centroid = stats(idx).Centroid;
            xCoord = centroid(1);
            yCoord = centroid(2);
        else
            xCoord = NaN;
            yCoord = NaN;
        end
    else
        % If no mouse is detected, set coordinates to NaN
        xCoord = NaN;
        yCoord = NaN;
    end

    % Calculate the timestamp
    timestamp = (frameCount - 1) * timeStep;

    % Store results in a data array
    data = [data; xCoord, yCoord, pixelChangeSum, ledSum, timestamp];
    %     end
end
% mobility = movmean(data(:,3),1);
mobility = data(:,3);
figure
plot(mobility);hold on;plot(100000*(abs([0; diff(movmean(mobility,60))])>2000))

toc
% Save the results to a file
outputFile = 'SleepChronic_2024-11-28T10_26_49_reconstructed.csv';
csvwrite(outputFile, data);