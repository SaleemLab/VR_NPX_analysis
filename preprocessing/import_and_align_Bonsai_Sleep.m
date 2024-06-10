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
DIR = dir(fullfile(options(1).EPHYS_DATAPATH,'*syncpulseTimes.mat'))

if ~isempty(DIR) % everything sync to probe 1
    load(fullfile(options(1).EPHYS_DATAPATH,DIR.name))
%     syncTimes_ephys = syncTimes_ephys.on; % use upswings currently
else
    [AP_FILE,LF_FILE] = findImecBinFile(options(1).EPHYS_DATAPATH);
    %     FILE_TO_USE = AP_FILE;
    binpath = fullfile(options(1).EPHYS_DATAPATH,AP_FILE);
    %     binpath = fullfile(options(1).EPHYS_DATAPATH,LF_FILE);
    syncTimes_ephys = SGLXextractSyncPulseFromBinFile(binpath);
    parseDate = date;
    [~,fname] = fileparts(options(1).EPHYS_DATAPATH);
    save(fullfile(options(1).EPHYS_DATAPATH,[fname,'_syncpulseTimes.mat']),'syncTimes_ephys','parseDate');
%     syncTimes_ephys = syncTimes_ephys.on; % use upswings currently
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
photodiode = table2array(this_table(:,4)); %

% photodiode_smoothed = smoothdata(table2array(this_table(:,4)),'movmedian',5); %
peripherals = [];
peripherals.Sync = kmeans(photodiode-movmean(photodiode,90),2)-1;
peripherals.Time =  behaviour.computer_time_original*1000;
% idx_trial_start = peripherals.Time(find(diff(peripherals.Sync)==1) + 1);
% idx_trial_end = peripherals.Time(find(diff(peripherals.Sync)==-1) + 1);
[peripherals] = alignBonsaiToEphysSyncTimes(peripherals,syncTimes_ephys.on); % use upswings currently
behaviour.sglxTime_original = peripherals.sglxTime;

behaviour.sglxTime = resample(behaviour.sglxTime_original,behaviour.sglxTime_original,60);
behaviour.tvec = behaviour.sglxTime;

behaviour.mobility = resample(behaviour.mobility_original,behaviour.sglxTime_original,60);
behaviour.X = resample(behaviour.X_original,behaviour.sglxTime_original,60);
behaviour.Y = resample(behaviour.Y_original,behaviour.sglxTime_original,60);

% tvec = resample(Behaviour.tvec,Behaviour.tvec,60);

% plot(peripherals.sglxTime,peripherals.Sync*20000+100000);hold on;plot(0:1/30000:(length(syncTimes_ephys.Sync)-1)/30000,syncTimes_ephys.Sync*20000+110000)


end