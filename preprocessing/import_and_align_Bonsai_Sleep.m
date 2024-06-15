function [behaviour] = import_and_align_Bonsai_Sleep(StimulusName,options)
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
td = dir(fullfile(options(1).EPHYS_DATAPATH,'..','*NidqTimes*'));
if ~isempty(td)
    load(fullfile(options(1).EPHYS_DATAPATH,'..',td(1).name),'Nidq');
    syncTimes_ephys = Nidq.bonsai_sync_on;% use upswings currently
else
    error('nidq and ephys sync pulse extraction and alignment not done!')
    return
end


% Import wheel data, eye data, position data and photodiodide data
% and synchronise them to async pulse
% peripherals = readmatrix(fullfile([BONSAI_DATAPATH,'\',char(peripheral_path)]));
this_table = readtable(fullfile([BONSAI_DATAPATH,'\',char(peripheral_path)]),'Delimiter',{',','=','{','}','(',')'});

computer_timevec = datevec(table2cell(this_table(:,7)),'yyyy-mm-ddTHH:MM:SS.FFF');
behaviour.computer_time_original = (computer_timevec(:,4)*60*60 + computer_timevec(:,5)*60 + computer_timevec(:,6))';
behaviour.wheel_time_original =  behaviour.computer_time_original- behaviour.computer_time_original(1);
behaviour.Sync_original = table2array(this_table(:,4))';
behaviour.mobility_original = table2array(this_table(:,5))';
behaviour.X_original = table2array(this_table(:,1))';
behaviour.Y_original = table2array(this_table(:,3))';

% peripherals = import_bonsai_peripherals(fullfile([BONSAI_DATAPATH,'\',char(peripheral_path)]));
Sync_pulse = table2array(this_table(:,4)); %

% photodiode_smoothed = smoothdata(table2array(this_table(:,4)),'movmedian',5); %
peripherals = [];
peripherals.Sync = kmeans(Sync_pulse-movmean(Sync_pulse,120),2)-1;
peripherals.Time =  behaviour.computer_time_original*1000;
idx_trial_start = find(diff(peripherals.Sync)==1) + 1;
idx_trial_end = find(diff(peripherals.Sync)==-1) + 1;

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

if idx_trial_end(1)<idx_trial_start(1)
    idx_trial_end(1)= [];
end


false_index = find(abs(peripherals.Time(idx_trial_end)-peripherals.Time(idx_trial_start))<0.4*1000); %if the pulse was less than 0.3 second
for n = 1:length(false_index)
    peripherals.Sync(idx_trial_start(false_index(n))-1:idx_trial_end(false_index(n))+1)=0; % makes it zero
end

if length(Nidq.bonsai_sync_on) == length(Nidq.bonsai_sync_off)+1
    Nidq.bonsai_sync_on(end)=[];
end

% [peripherals] = alignBonsaiToNidqSyncTimes(peripherals,sort([Nidq.bonsai_sync_on Nidq.bonsai_sync_off])); % use upswings currently
[peripherals] = alignBonsaiToEphysSyncTimes(peripherals,Nidq.bonsai_sync_on); % use upswings currently
behaviour.sglxTime_original = peripherals.sglxTime;

behaviour.sglxTime = interp1(behaviour.sglxTime_original,behaviour.sglxTime_original,behaviour.sglxTime_original(1):1/60:behaviour.sglxTime_original(end),'linear');
            
behaviour.tvec = behaviour.sglxTime;

behaviour.mobility = interp1(behaviour.sglxTime_original,behaviour.mobility_original,behaviour.sglxTime_original(1):1/60:behaviour.sglxTime_original(end),'linear');
behaviour.X = interp1(behaviour.sglxTime_original,behaviour.mobility_original,behaviour.sglxTime_original(1):1/60:behaviour.sglxTime_original(end),'linear');
behaviour.Y = interp1(behaviour.sglxTime_original,behaviour.mobility_original,behaviour.sglxTime_original(1):1/60:behaviour.sglxTime_original(end),'linear');
behaviour.Sync_pulse = interp1(behaviour.sglxTime_original,Sync_pulse',behaviour.sglxTime_original(1):1/60:behaviour.sglxTime_original(end),'linear');

% tvec = resample(Behaviour.tvec,Behaviour.tvec,60);

% plot(peripherals.sglxTime(1:300000),peripherals.Sync(1:300000)*2000+10000);hold on;plot(Nidq.sglxTime(1:3000000),Nidq.bonsai_sync(1:3000000))
% plot(behaviour.sglxTime(1:300000),behaviour.Sync_pulse(1:300000)-mean(behaviour.Sync_pulse));hold on;plot(Nidq.sglxTime(1:3000000),Nidq.bonsai_sync(1:3000000)-mean(Nidq.bonsai_sync))



end
