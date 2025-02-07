function [behaviour] = import_and_align_Bonsai_OpenField(StimulusName,options)
% Program to import Open field bonsai logs and eye data and align to ePhys data and 
% align the ephys data based on the delay between quad and photodioide signal


% Find Bonsai files (For now just one bonsai file with all behavioural and photodiode signals)
bonsai_files_names = options(1).bonsai_files_names;

peripheral_path = bonsai_files_names(contains(bonsai_files_names,'.csv'));
if length(peripheral_path)>1 % possibly DLC related csv
    for nfile = 1:length(peripheral_path)
        if contains(peripheral_path{nfile},'DLC')
            file_type(nfile) = 0;
        else
            file_type(nfile) = 1;
        end
    end
    peripheral_path = peripheral_path{file_type};
end

BONSAI_DATAPATH = options(1).BONSAI_DATAPATH;
% DIR = dir(fullfile(options(1).EPHYS_DATAPATH,'*syncpulseTimes.mat'))
Nidq = [];
td = dir(fullfile(options(1).EPHYS_DATAPATH,'..','*syncpulseTimes*'));
if ~isempty(td)
    load(fullfile(td.folder,td.name));
    Nidq.on = syncTimes_ephys.on;% use upswings currently
    Nidq.off = syncTimes_ephys.off;% use upswings currently
else
    error('nidq and ephys sync pulse extraction and alignment not done!')
    return
end
frame_sampling_rate = 30; % 30Hz videos


% Import Sleap and sync pulse
% and synchronise them to async pulse
this_table = readmatrix(fullfile([BONSAI_DATAPATH,'\',char(peripheral_path)]));
behaviour.frame = this_table(:,16);
tracking = this_table(:,2:15);
Sync_pulse = this_table(:,1); %
behaviour.Sync_pulse = Sync_pulse;
% peripherals = import_bonsai_peripherals(fullfile([BONSAI_DATAPATH,'\',char(peripheral_path)]));


% photodiode_smoothed = smoothdata(table2array(this_table(:,4)),'movmedian',5); %
peripherals = [];
peripherals.Sync_pulse = Sync_pulse;
baseline_corrected_sync = Sync_pulse-movmean(Sync_pulse,120);
initial_centroids = [ prctile(baseline_corrected_sync,10);prctile(baseline_corrected_sync,90)];
peripherals.Sync = kmeans(baseline_corrected_sync,2,'Start', initial_centroids)-1;

peripherals.Time =  behaviour.frame ./ frame_sampling_rate .* 1000;
idx_trial_start = find(diff(peripherals.Sync)==1) + 1;
idx_trial_end = find(diff(peripherals.Sync)==-1) + 1;

% figure;
% plot(mean(Sync_pulse)/3+mean(Sync_pulse)/2*kmeans(Sync_pulse-movmean(Sync_pulse,120),2)-1)
% hold on;
% plot(Sync_pulse)

if idx_trial_end(1)<idx_trial_start(1)
    idx_trial_end(1)= [];
end

if length(idx_trial_start)==length(idx_trial_end)
    if sum(idx_trial_end-idx_trial_start<0)>length(idx_trial_end)
        idx_trial_start1=idx_trial_start;
        idx_trial_start=idx_trial_end;
        idx_trial_end=idx_trial_start1;
        idx_trial_start1 = [];
    end
end

if length(idx_trial_start)==length(idx_trial_end)+1
    idx_trial_start(end)=[];
end

if length(idx_trial_start)==length(idx_trial_end)
    if sum(idx_trial_end-idx_trial_start<0)>length(idx_trial_end)
        idx_trial_start1=idx_trial_start;
        idx_trial_start=idx_trial_end;
        idx_trial_end=idx_trial_start1;
        idx_trial_start1 = [];
    end

    idx_to_remove = find(idx_trial_end-idx_trial_start<10);
    
    idx_trial_end(idx_to_remove)=[];
    idx_trial_start(idx_to_remove)=[];
end

false_index = find(abs(peripherals.Time(idx_trial_end)-peripherals.Time(idx_trial_start))<0.4 * 1000); %if the pulse was less than 0.3 second
for n = 1:length(false_index)
    peripherals.Sync(idx_trial_start(false_index(n))-1:idx_trial_end(false_index(n))+1)=0; % makes it zero
end

if length(Nidq.on) == length(Nidq.off)+1
    Nidq.bonsai_sync_on(end)=[];
end


[peripherals] = alignBonsaiToEphysSyncTimes(peripherals,Nidq.on); % use upswings currently

temp = peripherals;
temp.Time= peripherals.sglxTime * 1000;
temp.Sync(peripherals.Sync==0)=1;
temp.Sync(peripherals.Sync==1)=0;
[peripherals] = alignBonsaiToEphysSyncTimes(temp,Nidq.off); % use upswings currently

[unique_time,ia,~]=unique(peripherals.sglxTime);


behaviour.sglxTime_original = unique_time;


behaviour.sglxTime = interp1(behaviour.sglxTime_original,behaviour.sglxTime_original,behaviour.sglxTime_original(1):1/60:behaviour.sglxTime_original(end),'linear');            
behaviour.tvec = behaviour.sglxTime;

behaviour.frame = interp1(behaviour.sglxTime_original,behaviour.frame(ia),behaviour.sglxTime_original(1):1/60:behaviour.sglxTime_original(end),'linear');
for i = 1:size(tracking,2)
new_tracking(:,i) = interp1(behaviour.sglxTime_original,tracking(ia,i),behaviour.sglxTime_original(1):1/60:behaviour.sglxTime_original(end),'linear');
end
behaviour.Sync_pulse = interp1(behaviour.sglxTime_original,Sync_pulse(ia)',behaviour.sglxTime_original(1):1/60:behaviour.sglxTime_original(end),'linear');
tracking_part_list = {'nose','neck','implant','earR','earL','spine','tail_start'};
no_parts = length(tracking_part_list);
for i = 1:no_parts
behaviour.(tracking_part_list{i}) = new_tracking(:,(i-1)*2+1:i*2)';
end

end