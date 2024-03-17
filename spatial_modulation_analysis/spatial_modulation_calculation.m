function [place_fields] = spatial_modulation_calculation(place_fields,Task_info)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%take odd lap average response to find peak location
no_unit = length(place_fields(1).raw);
bin_size = 140/length(place_fields(1).x_bin_centres);
odd_avg_resp = nan(no_unit,2,length(place_fields(1).x_bin_centres));
even_avg_resp = nan(no_unit,2,length(place_fields(1).x_bin_centres));

for track_id = 1:2
    if mod(size(place_fields(1).raw{1},1),2) == 0
        odd_index{track_id}= 1:2:size(place_fields(track_id).raw{1},1)-1;
        even_index{track_id} = 2:2:size(place_fields(track_id).raw{1},1);
    else
        odd_index{track_id}= 1:2:size(place_fields(track_id).raw{1},1);
        even_index{track_id} = 2:2:size(place_fields(track_id).raw{1},1)-1;
    end
end


w = gausswin(21);
w = w / sum(w);

for i = 1:no_unit
    smooth_t1 = filtfilt(w,1,place_fields(1).raw{i}')';
    smooth_t2 = filtfilt(w,1,place_fields(2).raw{i}')';
    odd_avg_resp(i,1,:) = mean(smooth_t1(odd_index{track_id},:),1);
    odd_avg_resp(i,2,:) = mean(smooth_t2(odd_index{track_id},:),1);
    even_avg_resp(i,1,:) = mean(smooth_t1(even_index{track_id},:),1);
    even_avg_resp(i,2,:) = mean(smooth_t2(even_index{track_id},:),1);
end
[~,odd_peak_location] = max(odd_avg_resp(:,:,20:60),[],3);
% odd_peak_location(odd_peak_location<=10) = nan;
odd_peak_location = odd_peak_location+19;
odd_peak_location_sym = odd_peak_location - 80/bin_size;%make symmetrical location values - symmetry at 80cm
non_prefer_location = nan(size(odd_peak_location_sym));
non_prefer_location(odd_peak_location_sym<0) = odd_peak_location_sym(odd_peak_location_sym<0)+120/bin_size;
non_prefer_location(odd_peak_location_sym>=0) = odd_peak_location_sym(odd_peak_location_sym>=0)+40/bin_size;

even_R_prefer = nan([no_unit 2]);
even_R_non_prefer = nan([no_unit 2]);
for i =1:no_unit
    if ~isnan(odd_peak_location(i,1))
        even_R_prefer(i,1) = squeeze(even_avg_resp(i,1,odd_peak_location(i,1)));
        even_R_non_prefer(i,1) = max(even_avg_resp(i,1,non_prefer_location(i,1)-5:non_prefer_location(i,1)+5));

    end
    if ~isnan(odd_peak_location(i,2))
        even_R_prefer(i,2) = squeeze(even_avg_resp(i,2,odd_peak_location(i,2)));
        even_R_non_prefer(i,2) = max(even_avg_resp(i,2,non_prefer_location(i,2)-5:non_prefer_location(i,2)+5));
    end
end


SMI = (even_R_prefer - even_R_non_prefer)./(even_R_prefer + even_R_non_prefer);
place_fields(1).SMI = SMI(:,1);
place_fields(2).SMI = SMI(:,2);

% Lap


Task_info.lap_ID_all(find(abs(diff(Task_info.block_ID_all))>0))
Task_info.track_ID_all(find(abs(diff(Task_info.lap_ID_all))>1))
Task_info.track_ID_all([1; find(abs(diff(Task_info.lap_ID_all))>1)])

for i = 1:no_unit
    for track_id = 1:2
        resp = filtfilt(w,1,place_fields(track_id).raw{i}')';

        [~,peak_location] = max(resp(:,20:60),[],2);
        % odd_peak_location(odd_peak_location<=10) = nan;
        peak_location = peak_location+19;
        peak_location_sym = peak_location - 80/bin_size;%make symmetrical location values - symmetry at 80cm
        non_prefer_location = nan(size(peak_location_sym));
        non_prefer_location(peak_location_sym<0) = peak_location_sym(peak_location_sym<0)+120/bin_size;
        non_prefer_location(peak_location_sym>=0) = peak_location_sym(peak_location_sym>=0)+40/bin_size;

        for nlap = 1:size(resp,1)
            R_prefer = resp(nlap,peak_location(nlap));
            R_non_prefer = max(resp(nlap,non_prefer_location(nlap)-5:non_prefer_location(nlap)+5));
            place_fields(track_id).lap_SMI(i,nlap) = (R_prefer - R_non_prefer)./(R_prefer + R_non_prefer);
        end
        



    end
end

place_fields(1).lap_SMI = SMI(:,1);
place_fields(2).lap_SMI = SMI(:,2);




end
