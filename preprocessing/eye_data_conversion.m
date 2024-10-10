function [pupil_ellipse,new_tracking_points_coordinates] = eye_data_conversion(eye_data)

% This function inputs the eye tracking data analysed by DeepLabCut
% The eye_data should be time x tracking coordinates (24 columns) - see
% below comments for more exlanation of the matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Outputs:
% pupil_ellipse() which contains the properties(1-6 columns:
% x coordinate of the centre of the ellipse, y coordinate, major ax,
% minor ax, angle of rotation, area of the pupil) of the fitted ellipse on 
% the DLC tracking points, movement of the pupil (7-8 columns: distance and
% angle),and the eye blink detection(9th column - logical)

% new_tracking_points_coordinates():
% are the DLC tracking points transformed into a new coordinate system on the
% basis of the eye lids rather than camera pixel coordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is an illustration of how this function works:
% 1. Tracking points with likelihood value < 0.9 are turned to NaNs
% 2. Construct the new coordinate system from 4 tracking points on the
% eye lids using the mean value - 0 degree at nasal iris, 90 degree at 
% top eye lid, 180 degree at temporal iris, 270 degree at bottom eye lid.
% 3. Calulate the intersection point between 0-180 and 90-270 lines
% 4. Transformation matrix of a 2nd order polynomial transformation 
% by assuming the original two lines are perpendicular
% 5. Reverse the coordinates of all points using transoformation matrix
% 6. Fit an ellipse of the pupil tracking points - 0 degree, 45, 90, etc
% 7. Calculate the area of the ellipse and the movement (angle and
% distance) between the centres of the fitted ellipse
% 8. Calculate the distance between top and bot eyelids and z-score it
% 9. Eye blink is thresholded at 3 std below the mean

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% notes of the table elements:
% First column is the frame number - remember it starts from 0 (not 1)
% 2nd, 3rd columns are x,y coordinates of pupil0 tracking point
% 4th column is the likelihood value of pupil0 tracking point
% 5-7 columns: pupil45, 8-10: pupil90, 11-13: pupil135, 14-16: pupil180,
% 17-19: pupil225, 20-22: pupil270, 23-25: pupil315, 26-28: eyelid0, 29-31:
% eyelid90, 32-34: eyelid180, 35-37: eyelid270


%Consider convert coordinate values to NaN when the likelihood value is low
y_coordinate_columns = 3:3:36;
%invert back the y coordinate - the DLC y pixel counts from top to bottom
eye_data_new = eye_data;
eye_data_new(:,y_coordinate_columns) = 960 - eye_data_new(:,y_coordinate_columns);
%columns in the matrix that are x,y coordinates of the tracking points
tracking_points = 1:36;
tracking_points(mod(tracking_points, 3) == 0) = [];
tracking_points_columns = tracking_points + 1;
%columns in the matrix that are likelihood values for each tracking point
tracking_points_coordinates = eye_data_new(:,tracking_points_columns);
likelihood_columns_pupil = 4:3:25;
likelihood_threshold_index_pupil =eye_data_new(:,likelihood_columns_pupil) < 0.85;
likelihood_index_XY_pupil = zeros(size(likelihood_threshold_index_pupil,1),size(likelihood_threshold_index_pupil,2)*2);
for i=1:length(likelihood_columns_pupil)
    likelihood_index_XY_pupil(:,i*2-1) = likelihood_threshold_index_pupil(:,i);
    likelihood_index_XY_pupil(:,i*2) = likelihood_threshold_index_pupil(:,i);
end
tracking_points_coordinates(logical(likelihood_index_XY_pupil)) = NaN;

likelihood_columns_eyelid = 28:3:37;
likelihood_threshold_index_eyelid =eye_data_new(:,likelihood_columns_eyelid) < 0.8;
likelihood_index_XY_eyelid = zeros(size(likelihood_threshold_index_eyelid,1),size(likelihood_threshold_index_eyelid,2)*2);
for i=1:length(likelihood_columns_eyelid)
    likelihood_index_XY_eyelid(:,i*2-1) = likelihood_threshold_index_eyelid(:,i);
    likelihood_index_XY_eyelid(:,i*2) = likelihood_threshold_index_eyelid(:,i);
end
tracking_points_coordinates(logical(likelihood_index_XY_eyelid)) = NaN;
%change basis of coordinate from eyelid labels
eyelid0 = tracking_points_coordinates(:, 17:18);
eyelid90 = tracking_points_coordinates(:, 19:20);
eyelid180 = tracking_points_coordinates(:, 21:22);
eyelid270 = tracking_points_coordinates(:, 23:24);
mean_eyelid0 = mean(eyelid0,1,'omitnan');
mean_eyelid90 = mean(eyelid90,1,'omitnan');
mean_eyelid180 = mean(eyelid180,1,'omitnan');
mean_eyelid270 = mean(eyelid270,1,'omitnan');


% Calculate the intersection
denominator = (mean_eyelid0(1)-mean_eyelid180(1))*(mean_eyelid90(2)-mean_eyelid270(2)) - (mean_eyelid0(2)-mean_eyelid180(2))*(mean_eyelid90(1)-mean_eyelid270(1));
if denominator == 0
    disp('The lines are parallel.')
else
    t = ((mean_eyelid0(1)-mean_eyelid90(1))*(mean_eyelid90(2)-mean_eyelid270(2)) - (mean_eyelid0(2)-mean_eyelid90(2))*(mean_eyelid90(1)-mean_eyelid270(1))) / denominator;
    u = -((mean_eyelid0(1)-mean_eyelid180(1))*(mean_eyelid0(2)-mean_eyelid90(2)) - (mean_eyelid0(2)-mean_eyelid180(2))*(mean_eyelid0(1)-mean_eyelid90(1))) / denominator;
    if t >= 0 && t <= 1 && u >= 0 && u <= 1 %%intersection point lies within the 4 points
        intersection = mean_eyelid0 + t * (mean_eyelid180 - mean_eyelid0);
    end
end
eyelid_intersection = intersection;
% Assume the 4 original eye lid points are perpendicular to the
% intersection point along the original axises (i.e. eyelid0 lies on (+,0),
% eyelid90 lies on (0,+), eyelid180 lies on (-,0), eyelid270 lies on (0,-).
% On the image we recorded, the original eye lid points (everything inside
% the eyes) went through polynomial geometric transform (which is the camera angle) 
moving_points= [mean_eyelid0(1,1) - eyelid_intersection(1,1), 0; %original eyelid0 point
                0, mean_eyelid90(1,2) - eyelid_intersection(1,2); % original eyelid90 point
                mean_eyelid180(1,1) - eyelid_intersection(1,1),0; % originial eyelid180 point
                0, mean_eyelid270(1,2) - eyelid_intersection(1,2); % original eyelid270 point
                0,0;  % intersection point
                (mean_eyelid0(1,1) - eyelid_intersection(1,1))/2, (mean_eyelid90(1,2) - eyelid_intersection(1,2))/2 %mid point between eyelid0 and eyelid90
];
fixed_points = [mean_eyelid0;mean_eyelid90;mean_eyelid180;mean_eyelid270; eyelid_intersection;(mean_eyelid0+mean_eyelid90)/2];
transformation_matrix = fitgeotrans(moving_points,fixed_points,'polynomial',2); %annoying thing is that polynomial only allows inverse transformation - so a bit confusing to read the code

% allocate new tracking point coordinates one column at once
new_tracking_points_coordinates = nan(size(tracking_points_coordinates));
for i = 1:12
    temp_tracking_points_coordinate_X = tracking_points_coordinates(:,i*2-1);
    temp_tracking_points_coordinate_Y = tracking_points_coordinates(:,i*2);
    NaN_index = isnan(temp_tracking_points_coordinate_X);
    temp_tracking_points_coordinate_X = temp_tracking_points_coordinate_X(~NaN_index);
    temp_tracking_points_coordinate_Y = temp_tracking_points_coordinate_Y(~NaN_index);
    [new_X, new_Y] = transformPointsInverse(transformation_matrix, temp_tracking_points_coordinate_X,temp_tracking_points_coordinate_Y); 
    new_tracking_points_coordinates(~NaN_index,i*2-1) = new_X;
    new_tracking_points_coordinates(~NaN_index,i*2) = new_Y;
end



pupil_ellipse = nan(size(new_tracking_points_coordinates,1),5);

%first column centre x coordinate, 2nd column centre y coordinate, 3rd column major diam, 4th
%column minor diam, 5th column angle, 6th column area, 7-8th column movement
for iFrame = 1:size(new_tracking_points_coordinates,1)
    tracking_points_XY_frame = nan(2,8);
    tracking_points_coordinates_frame = new_tracking_points_coordinates(iFrame,:);
    tracking_points_XY_frame(1,:) = tracking_points_coordinates_frame(1,1:2:15);
    tracking_points_XY_frame(2,:) = tracking_points_coordinates_frame(1,2:2:16);
    if ~any(isnan(tracking_points_coordinates_frame))
        fobj = ellipticalFit(tracking_points_XY_frame);
        pupil_ellipse(iFrame,1:2) = fobj.center;
        pupil_ellipse(iFrame,3) = fobj.a;
        pupil_ellipse(iFrame,4) = fobj.b;
        pupil_ellipse(iFrame,5) = fobj.angle;
    end
end
%calculate the pupil area using Area = pi * major ax * minor ax
pupil_ellipse(:,6) = calculate_ellipse_areas(pupil_ellipse(:,1:5));

% calculate the movement of the ellipse centre as the pupil movement
% Calculate the difference between each pair of points
diff_matrix = diff(pupil_ellipse(:,1:2));

% Calculate the Euclidean distance (magnitude of displacement) at each time step
distances = sqrt(sum(diff_matrix.^2, 2));

% Calculate the angle of movement at each time step
angles = atan2d(diff_matrix(:, 2), diff_matrix(:, 1));

% To handle negative angles, convert them into positive
angles = mod(angles + 360, 360);


pupil_ellipse(2:end,7) = distances;
pupil_ellipse(2:end,8) = angles;


% eye blink detection using eyelid90 and eyelid270 distance
eyelid90 = new_tracking_points_coordinates(:,19:20);
eyelid270 = new_tracking_points_coordinates(:,23:24);
eyelid_distance =  sqrt(sum((eyelid90 - eyelid270).^2, 2));
mean_distance = nanmean(eyelid_distance);
std_distance = nanstd(eyelid_distance);
% Z-score the matrix
z_eyelid_distance = (eyelid_distance - mean_distance) / std_distance;
eye_blink = z_eyelid_distance < -5; % 5 standard deviations below the mean
pupil_ellipse(:,9) = eye_blink;

% % code to check eye_blinks
% new_tracking_points_X = new_tracking_points_coordinates(:,1:2:23);
% new_tracking_points_Y = new_tracking_points_coordinates(:,2:2:24);
% figure;
% eye_blink_frames = find(eye_blink == 1);
% eye_blink_frames = [eye_blink_frames; (500:510)']; 
% % eye_blink_frames_expanded = [];
% % 
% % for i = 1:length(eye_blink_frames)
% %     % Get the current index
% %     curr_index = eye_blink_frames(i);
% %     
% %     % Calculate the previous and following 2 indices
% %     prev_indices = max(1, curr_index-2):curr_index-1;
% %     next_indices = curr_index+1:min(length(eye_blink), curr_index+2);
% %     
% %     % Append these indices to the expanded list
% %     eye_blink_frames_expanded = [eye_blink_frames_expanded, prev_indices, curr_index, next_indices];
% % end
% 
% % Remove duplicates and sort the indices
% % eye_blink_frames_expanded = unique(eye_blink_frames_expanded);
% frame_num = length(eye_blink_frames);
% % set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
% % ellipse_x = zeros(100,frame_num);
% % ellipse_y = zeros(100,frame_num);
% % ellipse_time = repmat(1:frame_num,[100,1]);
% % iFrame = 1;
% % for iTime = eye_blink_frames
% %     x_centre = pupil_ellipse(iTime,1);
% %     y_centre = pupil_ellipse(iTime,2);
% %     major_diam = pupil_ellipse(iTime,3);
% %     minor_diam = pupil_ellipse(iTime,4);
% %     angle = pupil_ellipse(iTime,5);
% %     theta=linspace(0,2*pi,100);
% %     ellipse_x(:,iFrame)=major_diam*cos(theta)*cosd(angle)-sind(angle)*minor_diam*sin(theta)+x_centre;
% %     ellipse_y(:,iFrame)=major_diam*cos(theta)*sind(angle)+cosd(angle)*minor_diam*sin(theta)+y_centre;
% %     iFrame = iFrame+1;
% % end
% % plot3(ellipse_time,ellipse_x,ellipse_y,'Color','k')
% hold on;
% for i = 1:12
%     scatter3(1:frame_num,new_tracking_points_X(eye_blink_frames,i)',new_tracking_points_Y(eye_blink_frames,i)')
% end
