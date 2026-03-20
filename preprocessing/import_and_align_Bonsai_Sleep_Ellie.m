function [behaviour] = import_and_align_Bonsai_Sleep_Ellie(StimulusName,options)
% Program to import Open field bonsai logs and eye data and align to ePhys data and 
% align the ephys data based on the delay between quad and photodioide signal


% Find Bonsai files (For now just one bonsai file with all behavioural and photodiode signals)
bonsai_files_names = options(1).bonsai_files_names;

peripheral_path = bonsai_files_names(contains(bonsai_files_names,'.csv'));
if length(peripheral_path)>1 % possibly DLC related csv
    for nfile = 1:length(peripheral_path)
        if contains(peripheral_path{nfile},'DLC')
            file_type(nfile) = false;
        elseif contains(peripheral_path{nfile},'_Async') % use the Video logging csv file, not the Async pulse csv file
            file_type(nfile) = false;
        else
            file_type(nfile) = true;
        end
    end
    file_type = find(file_type==1);
    peripheral_path = peripheral_path{file_type};
end

BONSAI_DATAPATH = options(1).BONSAI_DATAPATH;
% DIR = dir(fullfile(options(1).EPHYS_DATAPATH,'*syncpulseTimes.mat'))
Nidq = [];
td = dir(fullfile(options(1).EPHYS_DATAPATH,'..','*NidqTimes*'));
if ~isempty(td)
    load(fullfile(options(1).EPHYS_DATAPATH,'..',td(1).name),'Nidq');
    syncTimes_ephys_on = Nidq.bonsai_sync_on;% use upswings currently
    syncPulse_ephys = Nidq.bonsai_sync;
    syncPulse_ephysTimes = Nidq.sglxTime;
else
    error('nidq and ephys sync pulse extraction and alignment not done!')
    return
end


% Import wheel data, eye data, position data and photodiodide data
% and synchronise them to async pulse
% peripherals = readmatrix(fullfile([BONSAI_DATAPATH,'\',char(peripheral_path)]));
if contains(peripheral_path,'reconstruct')
    this_table = readtable(fullfile([BONSAI_DATAPATH,'\',char(peripheral_path)]));

    %behaviour.computer_time_original = table2array(this_table(:,7))'; %This is actually Arduino time for Ellie's data
    behaviour.computer_time_original = table2array(this_table(:,2))'; %This is actually Bonsai time for Ellie's data
    %behaviour.wheel_time_original =  table2array(this_table(:,5))';
    behaviour.Sync_original = table2array(this_table(:,3))';
    behaviour.mobility_original = table2array(this_table(:,4))';
    behaviour.X_original = table2array(this_table(:,5))';
    behaviour.Y_original = table2array(this_table(:,6))';

    Sync_pulse = table2array(this_table(:,3)); %
else
    this_table = readtable(fullfile([BONSAI_DATAPATH,'\',char(peripheral_path)]),'Delimiter',{',','=','{','}','(',')'});

    % computer_timevec = datevec(table2cell(this_table(:,7)),'yyyy-mm-ddTHH:MM:SS.FFF');
    %behaviour.computer_time_original = table2array(this_table(:,7))';  %This is actually Arduino time for Ellie's data
    behaviour.computer_time_original = table2array(this_table(:,2))';  %This is actually Bonsai time for Ellie's data
    %behaviour.wheel_time_original =  behaviour.computer_time_original- behaviour.computer_time_original(1);
    behaviour.Sync_original = table2array(this_table(:,3))';
    behaviour.mobility_original = table2array(this_table(:,4))';
    behaviour.X_original = table2array(this_table(:,5))';
    behaviour.Y_original = table2array(this_table(:,6))';

    Sync_pulse = table2array(this_table(:,3)); %
end
% peripherals = import_bonsai_peripherals(fullfile([BONSAI_DATAPATH,'\',char(peripheral_path)]));


% photodiode_smoothed = smoothdata(table2array(this_table(:,4)),'movmedian',5); %
peripherals = [];
peripherals.Sync_pulse = Sync_pulse;
baseline_corrected_sync = Sync_pulse-movmean(Sync_pulse,120);
initial_centroids = [ prctile(baseline_corrected_sync,10);prctile(baseline_corrected_sync,90)];
peripherals.Sync = kmeans(baseline_corrected_sync,2,'Start', initial_centroids)-1;
% if mean(peripherals.Sync_pulse(peripherals.Sync==1))<mean(peripherals.Sync_pulse(peripherals.Sync==0))
%     % if 0 is when Sync goes up, swap 1 and 0
%     temp = peripherals.Sync;
%     peripherals.Sync(temp==0)=1;
%     peripherals.Sync(temp==1)=0;
% end

peripherals.Time =  behaviour.computer_time_original; %This is actually Bonsai time for Ellie's data
idx_trial_start = find(diff(peripherals.Sync)==1) + 1;
idx_trial_end = find(diff(peripherals.Sync)==-1) + 1;

figure;
plot(mean(Sync_pulse)/3+mean(Sync_pulse)/2*kmeans(Sync_pulse-movmean(Sync_pulse,120),2)-1)
hold on;
plot(Sync_pulse)

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

false_index = find(abs(peripherals.Time(idx_trial_end)-peripherals.Time(idx_trial_start))<0.4*1000); %if the pulse was less than 0.3 second
for n = 1:length(false_index)
    peripherals.Sync(idx_trial_start(false_index(n))-1:idx_trial_end(false_index(n))+1)=0; % makes it zero
end

if length(Nidq.bonsai_sync_on) == length(Nidq.bonsai_sync_off)+1
    Nidq.bonsai_sync_on(end)=[];
end


% [peripherals] = alignBonsaiToEphysSyncTimes(peripherals,Nidq.bonsai_sync_on); % use upswings currently
 [peripherals] = alignBonsaiToEphysSyncPulse(peripherals,syncPulse_ephys,syncPulse_ephysTimes);

% temp = peripherals;
% temp.Time= peripherals.sglxTime*1000;
% temp.Sync(peripherals.Sync==0)=1;
% temp.Sync(peripherals.Sync==1)=0;
% [peripherals] = alignBonsaiToEphysSyncTimes(temp,Nidq.bonsai_sync_off); % use upswings currently

disp(['minimumn value of peripherals.sglxTimes: ', mat2str(min(peripherals.sglxTime))])
disp(['first 10 peripherals.sglxTimes: ', mat2str(peripherals.sglxTime(1:10))]);
[unique_time, ia, ~] = unique(peripherals.sglxTime);
disp(['first 10 unique times: ', mat2str(unique_time(1:10))]);

behaviour.sglxTime_original = unique_time;

disp(size(behaviour.sglxTime_original))
disp(size(behaviour.mobility_original(ia)))

behaviour.sglxTime = interp1(behaviour.sglxTime_original,behaviour.sglxTime_original,behaviour.sglxTime_original(1):1/60:behaviour.sglxTime_original(end),'linear');            
behaviour.tvec = behaviour.sglxTime;

behaviour.mobility = interp1(behaviour.sglxTime_original,behaviour.mobility_original(ia),behaviour.sglxTime_original(1):1/60:behaviour.sglxTime_original(end),'linear');
behaviour.X = interp1(behaviour.sglxTime_original,behaviour.X_original(ia),behaviour.sglxTime_original(1):1/60:behaviour.sglxTime_original(end),'linear');
behaviour.Y = interp1(behaviour.sglxTime_original,behaviour.Y_original(ia),behaviour.sglxTime_original(1):1/60:behaviour.sglxTime_original(end),'linear');
behaviour.Sync_pulse = interp1(behaviour.sglxTime_original,Sync_pulse(ia)',behaviour.sglxTime_original(1):1/60:behaviour.sglxTime_original(end),'linear');

%tvec = resample(behaviour.tvec,behaviour.tvec,60);
%figure
%plot(peripherals.sglxTime(1:200000),peripherals.Sync(1:200000)*2000+10000);
%hold on;plot(Nidq.sglxTime(1:2000000),Nidq.bonsai_sync(1:2000000));...
%hold on;scatter(Nidq.bonsai_sync_on(1:50),20000*ones(1,50));
%hold on;scatter(Nidq.bonsai_sync_off(1:50),10000*ones(1,50));
% 
%figure
%plot(peripherals.sglxTime,peripherals.Sync*2000+10000);
%hold on;plot(Nidq.sglxTime,Nidq.bonsai_sync);...
%    hold on;scatter(Nidq.bonsai_sync_on(1:500),20000*ones(1,500));
%hold on;scatter(Nidq.bonsai_sync_off(1:500),10000*ones(1,500));

%figure
% plot(Nidq.sglxTime(1:3000000),Nidq.bonsai_sync(1:3000000));hold on;
% plot(peripherals.sglxTime(1:300000),Sync_pulse(1:300000)./10);
%plot(peripherals.sglxTime(1:300000),peripherals.Sync(1:300000)*2000+10000);
%hold on
%plot(behaviour.sglxTime_original(1:300000),Sync_pulse(1:300000)./10);
% 

%plot(behaviour.sglxTime(1:300000),behaviour.Sync_pulse(1:300000)-mean(behaviour.Sync_pulse));hold on;plot(Nidq.sglxTime(1:3000000),Nidq.bonsai_sync(1:3000000)-mean(Nidq.bonsai_sync))

% TEST plots for Sync_pulse timings (as a check that the clock didn't restart during the recording)
% and effect of interpolation: original timepoints and values
figure
plot(behaviour.sglxTime_original, Sync_pulse(ia)', 'o')
xlabel('Time (s)')
ylabel('Sync pulse')
hold off
figure
plot(behaviour.sglxTime_original, behaviour.X_original(ia), 'o')
%plot(behaviour.X_original, 'o')
xlabel('Time (s)')
ylabel('behaviour.X original')
hold off
% interpolated timepoints and values
figure
plot(behaviour.sglxTime, behaviour.X, 'o')
%plot(behaviour.Y_original, 'o')
xlabel('Time (s)')
ylabel('behaviour.X interpolated')
hold off

%plot(behaviour.sglxTime, behaviour.Sync_pulse, '-')

%ylabel('Sync pulse')
%legend('Original', 'Interpolated')
%title('Comparison of Original vs Interpolated Sync Pulse')

% plot(peripherals.Time(1:300000)./1000-peripherals.Time(1),peripherals.Sync(1:300000)*2000+10000);hold on;plot(Nidq.sglxTime(1:3000000),Nidq.bonsai_sync(1:3000000))
end
