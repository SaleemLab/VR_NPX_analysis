%% COMPARISON OF PLACE FIELDS BETWEEN THE TWO HALVES OF A SESSION
% 29/11/17_MH
% This functions is to compare the place fields between the first and
% second half of the session. Each variable will be named X_1 or X_2
% depending on whether it contains data from the first or second half.
% OUTCOME: 'comp' structure containing all place field information for each
% half of session


function comp = place_field_quality(clusters,place_fields,position,lap_times,tracks_compared,number_of_laps)

(t_position,x2,spikes,cm_per_pixel,velocity_limit)

%%%%%%%%%%% Variable to change

nbins=101; % change to increase resolution 
%%%%%%%%%%%

%% Divides timestamps and position data taken across the task in two (we are dividing the session in two)

x2_half= round(length(x2)/2);
x2_1= x2(1:x2_half);
x2_2= x2((x2_half+1):end);

t2_half= round(length(t_position)/2);
t2_1= t_position(1:t2_half);
t2_2= t_position((t2_half+1):end);

spikes_half= round(length(spikes)/2);
spikes_1= spikes(1:spikes_half);
spikes_2= spikes((spikes_half+1):end);

%% Filter velocity - filter depends on animal velocity

filtered_x_cm1=filtfilt(ones(1,25)/25,1,x2_1);     %filter to 1 second
filtered_x_cm2=filtfilt(ones(1,25)/25,1,x2_2);     %filter to 1 second
    
% Calculate velocity

v_1=abs(diff(filtered_x_cm1)./median(diff(t2_1)));  % absolute to get rid of negative velocity/ median(diff(t))= median time between to timestamps=0.04 (from 25Hz)  
v_1(end+1)=v_1(end);  %makes v the same size as t

v_2=abs(diff(filtered_x_cm2)./median(diff(t2_2)));  % absolute to get rid of negative velocity/ median(diff(t))= median time between to timestamps=0.04 (from 25Hz)  
v_2(end+1)=v_2(end);  %makes v the same size as t
    
% Converts to cm per s
v_cm_1= v_1/cm_per_pixel;
v_cm_2= v_2/cm_per_pixel;

index_1=find(v_cm_1>velocity_limit);
index_2=find(v_cm_2>velocity_limit);

x_cm_1=[filtered_x_cm1]*cm_per_pixel; %changes linearized position (in pixels) to cm
x_cm_2=[filtered_x_cm2]*cm_per_pixel; %changes linearized position (in pixels) to cm


%% Time spent per position (in bins) 

% Creates bins for histogram
x_bins1= linspace(min(x_cm_1), max(x_cm_1), nbins);  %where x_cm are the linearized positions
x_bins2= linspace(min(x_cm_2), max(x_cm_2), nbins);  %where x_cm are the linearized positions

%x_hist_raw= histcounts(x_cm, x_bins);         %hist from raw data
x_hist_1= histcounts(x_cm_1(index_1), x_bins1);      %hist from velocity filtered data 
x_hist_2= histcounts(x_cm_2(index_2), x_bins2);      %hist from velocity filtered data 

total_time_at_x_1=x_hist_1.*median(diff(t2_1));   %time spent per position
total_time_at_x_2=x_hist_2.*median(diff(t_position-2));   %time spent per position
 
%% Spikes per position
spikes_1(spikes_1< min(t2_1))=[];
spikes_1(spikes_1> max(t2_1))=[];

spikes_2(spikes_2< min(t2_2))=[];
spikes_2(spikes_2> max(t2_2))=[];

% Choose nearest or linear: nearest- if you want to fill the gap with the closest point/ linear- if you want to fill the gap with an estimate middle value between two points
velocity_during_spike_1= interp1(t2_1,v_cm_1,spikes_1,'nearest');  
spike_index_1= find(velocity_during_spike_1>velocity_limit); 
velocity_during_spike_2= interp1(t2_2,v_cm_2,spikes_2,'nearest');  
spike_index_2= find(velocity_during_spike_2>velocity_limit); 
 
x_position_during_spike_1= interp1(t2_1,x_cm_1,spikes_1(spike_index_1),'linear');  %spikes per position  
spike_hist_1= histcounts(x_position_during_spike_1, x_bins1);  %spikes per position in each bin
x_position_during_spike_2= interp1(t2_2,x_cm_2,spikes_2(spike_index_2),'linear');  %spikes per position  
spike_hist_2= histcounts(x_position_during_spike_2, x_bins2);  %spikes per position in each bin

%% Place fields 

place_field_1=spike_hist_1./(total_time_at_x_1+0.001);
filtered_place_field_1=filtfilt([1 1 1 1 1]/5,1,place_field_1); %better use odd number for filter. Can be changed if necessary (e.g. 3 points)

place_field_2=spike_hist_2./(total_time_at_x_2+0.001);
filtered_place_field_2=filtfilt([1 1 1 1 1]/5,1,place_field_2); %better use odd number for filter. Can be changed if necessary (e.g. 3 points)

%% FUNCTION OUTCOME

% First half of session
comp.filtered_place_field_1= filtered_place_field_1; %filtered place field
comp.place_field_1= place_field_1;   % place field
comp.spike_hist_1= spike_hist_1;  % spikes per position
comp.total_time_at_x_1= total_time_at_x_1;  %spent time per position

comp.x_bins_1= x_bins1(1:end-1)+median(diff(x_bins1))/2; %to find the middle point of the bin

% Second half of session
comp.filtered_place_field_2= filtered_place_field_2; %filtered place field
comp.place_field_2= place_field_2;   % place field
comp.spike_hist_2= spike_hist_2;  % spikes per position
comp.total_time_at_x_2= total_time_at_x_2;  %spent time per position

comp.x_bins_2= x_bins2(1:end-1)+median(diff(x_bins2))/2; %to find the middle point of the bin


end