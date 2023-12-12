% My analysis pathway
verbose = 1;
MouseName = 'M19164';
ExpDate = '191213';
%List the recording for each session excluding the SparseNoise trials
Recordings = { ...
    'M19164_2019-12-13_15-11-16_6','101_CH*.continuous','101_ADC*.continuous'; ...
    };

% First convert to kilosort compatiable dat file
ePhysRecordingChannelRoot='116_CH*.continuous'; %Specify this only for those recordings done in 2020
ePhysAnalogChannelRoot='116_ADC*.continuous';
Convert_OpenEphys_2_dat(MouseName, ExpDate, Recordings(:,1), 'SC_1', Recordings(:,2));

% Now get the other analogie channels signals
global DIRS
DIRS = SetDefaults('animal',MouseName);
DS_SamplingRate=60; %TO BE CHANGED....Set the DS freq to align the ePhys and the Bonvision signal

%% Import, Align and Save data from Bonvision with ePhys.

for thisRecording = 1:length(Recordings(:,1))
    
    % Get analogue channel stuff for each recording
    [stimON, stimOFF, ACInfo] = getPhotoDiodeSignal(MouseName,ExpDate,Recordings{thisRecording,1}, 4, verbose,Recordings{thisRecording,2}); % This will create a AC_Info.mat file as well
    ACInfo = getSyncPulseSignal(MouseName,ExpDate,Recordings{thisRecording},3,Recordings{thisRecording,3});
    % Save the syncpulse
    TimestampsDownsampled_Orig = ACInfo.TimestampsDownsampled;
    SyncPulse_Orig = ACInfo.SyncPulse;
    if iscell(SyncPulse_Orig) %some recordings might contains two consecutive recordings (mostly because of error in starting the next recording which is now merged with the previous one)
        for i=1:size(SyncPulse_Orig,2)
            reclength(i,1)=length(SyncPulse_Orig{1,i});
        end
        sessionToKeep=find(max(reclength));
        SyncPulse=SyncPulse_Orig{1,sessionToKeep};
        TimestampsDownsampled=TimestampsDownsampled_Orig{1,sessionToKeep};
        ACInfo.TimestampsDownsampled=TimestampsDownsampled;
        ACInfo.SyncPulse=SyncPulse;
    else
        SyncPulse=SyncPulse_Orig;
        TimestampsDownsampled=TimestampsDownsampled_Orig;
    end
    SyncPulse(SyncPulse<1)=0;
    SyncPulse(SyncPulse>1)=1;
    %%%%%%%%%%%%
    save(fullfile(DIRS.ePhys,ExpDate,Recordings{thisRecording},'AC_Info.mat'),'TimestampsDownsampled','SyncPulse','-append');
    %%%%%%%%%%%%
    AllTrials=loadCSVforVR(MouseName,ExpDate,Recordings{thisRecording},SparseNoiseFileCount);
    [Delay, ind_M]= alignSyncPulseBonsai(AllTrials,ACInfo,DS_SamplingRate);
    DS_AllTrials=downsamplingBonsai(AllTrials, DS_SamplingRate);
    
    %First we need to subtract the minimum timestamp value
    ACInfo.Timestamps=ACInfo.Timestamps-min(ACInfo.Timestamps);
    
    % We need to define the index interval for the ePhys data which matches with bonsai VR behavioural data
    ePhys_VR_Start=min(find(ACInfo.Timestamps>=Delay));
    ePhys.Data=ACInfo.Data(ePhys_VR_Start:end,:);
    ePhys.Timestamps=ACInfo.Timestamps(ePhys_VR_Start:end,1);
    ePhys.Timestamps=ePhys.Timestamps-min(ePhys.Timestamps);
    
    %Now we can identify for each session the start of each single trial based on stimOFF and stimON variables
    %the first element of stimOFF represent the onset of the first trial
    %The second last element of stimON represent the onset of the last trial
    %The last element of stimOFF represent the end of the VR session
    % AllOnsets=sort([stimOFF;stimON(2:end)])-Delay;
    AllOnsets=sort([stimOFF;stimON]);
    
    % AllOnsets(find(AllOnsets>DS_AllTrials.Timestamp(end)))=[]; %Remove onsets after the end of the VR
    
    save(fullfile(DIRS.ePhys,ExpDate,Recordings{thisRecording},[Recordings{thisRecording}, '_processed']),'Delay','DS_AllTrials','AllOnsets');

end

%Now load for each recording the excel file containing the trial structure and align it with the ePhys data
function [AllTrials]=loadCSVforVR(MouseName,ExpDate,Recording,SparseNoiseFileCount)
global DIRS
foldername=fullfile(DIRS.VRBehaviour,ExpDate,'Behaviour');
fstruct = dir(fullfile(foldername,'*.csv')); %List all the csv file in the folder
%we need to order the file list based on the date/time of acquisition
%first we transform the table into a structure
T=struct2table(fstruct);
%we now sort it based on date
sortedT=sortrows(T,'date');
%we go back to the structure
fstruct=table2struct(sortedT);
if str2num(Recording(end))-SparseNoiseFileCount<=0
    fileNumber=str2num(Recording(end-1:end))-SparseNoiseFileCount; %this is because we need to deal with numbers higher than 9
else
    fileNumber=str2num(Recording(end))-SparseNoiseFileCount; %Check for the number of the session to load
end
filename=fullfile(DIRS.VRBehaviour,ExpDate,'Behaviour',fstruct(fileNumber).name);
AllTrials=readtable(filename);
AllTrials.Timestamp(1,1)=AllTrials.Time(1,1);
for  ii=2:size(AllTrials.Time,1)
    AllTrials.Timestamp(ii,1)=AllTrials.Timestamp(ii-1,1)+AllTrials.Time(ii,1);
end
end

%Now align the Bonsai SyncPulse with the ePhys SyncPulse to find the delay
%between the two signals
function [Delay,ind_M]=alignSyncPulseBonsai(AllTrials, ACInfo,DS_SamplingRate)
ComparisonTime = 20; % seconds ( not the number of samples!!! )
ChunkStart = 0.5;
clear OE_chunk
clear OE_chunkTime
OE_SamplingRate=ACInfo.SamplingRateOE/ACInfo.SamplingFactor;
OE_SyncPulse=ACInfo.SyncPulse;
OE_TimeStamps=ACInfo.TimestampsDownsampled;

OE_chunk = OE_SyncPulse(round(ChunkStart*length(OE_SyncPulse)):round(ChunkStart*length(OE_SyncPulse)+OE_SamplingRate*ComparisonTime));
OE_chunkTime = OE_TimeStamps(round(ChunkStart*length(OE_TimeStamps)):round(ChunkStart*length(OE_TimeStamps)+OE_SamplingRate*ComparisonTime)) - ...
    min(OE_TimeStamps(round(ChunkStart*length(OE_TimeStamps)):round(ChunkStart*length(OE_TimeStamps)+OE_SamplingRate*ComparisonTime)));

OE_chunk_DS = interp1(OE_chunkTime,...
    OE_chunk, ...
    1/DS_SamplingRate:1/DS_SamplingRate:ComparisonTime);
OE_chunkTime_DS = interp1(OE_chunkTime,...
    OE_chunkTime, ...
    1/DS_SamplingRate:1/DS_SamplingRate:ComparisonTime);

Bonsai_SyncPulse_DS=interp1(AllTrials.Timestamp,AllTrials.SyncPulse,1/DS_SamplingRate:1/DS_SamplingRate:AllTrials.Timestamp(end));
Bonsai_DS_Timestamps=1/DS_SamplingRate:1/DS_SamplingRate:AllTrials.Timestamp(end);

for i = 1:length(Bonsai_SyncPulse_DS)-(DS_SamplingRate*ComparisonTime)-1 % add 200 samples of buffering for the end
    R = corrcoef(Bonsai_SyncPulse_DS(i+1:i+ComparisonTime*DS_SamplingRate),OE_chunk_DS); % 500 is duration of pulse
    k(i) = R(1,2);
    clear R;
end
figure; plot(k,'.')
[val_M, ind_M] = max(k);

%Plot to check that the two signals overlap so the correlation is correct
plot(OE_chunkTime_DS,OE_chunk_DS,'k')
hold on
plot(Bonsai_DS_Timestamps(ind_M:ind_M+ComparisonTime*DS_SamplingRate)-min(Bonsai_DS_Timestamps(ind_M:ind_M+ComparisonTime*DS_SamplingRate)),...
    Bonsai_SyncPulse_DS(ind_M:ind_M+ComparisonTime*DS_SamplingRate),'r.')

%To align all the signals we need to know the temporal displacement between the OE signal and the Bonsai signal
%To do so we need to find the timestamp at the beginning of the chunk OE
%data we use to find the correlation. This timestamp will be align with the
%timestamp from bonsai signal downsampled where the correlation is max
%(ind_M)
OE_start=ACInfo.TimestampsDownsampled(round(ChunkStart*length(OE_SyncPulse))); 
Bonsai_start=Bonsai_DS_Timestamps(ind_M);
Delay=OE_start-Bonsai_start; %We calculate the delay to align the two signals.
%Bonsai_DS_Timestamps_corrected=Bonsai_DS_Timestamps+Delay;
%plot(ACInfo.TimestampsDownsampled,ACInfo.SyncPulse)
% hold on
% plot(Bonsai_DS_Timestamps_corrected,Bonsai_SyncPulse_DS)
end

function DS_AllTrials=downsamplingBonsai(AllTrials, DS_SamplingRate)
%we downsample the bonsai table at 60Hz
Bonsai_DS_Timestamps=1/DS_SamplingRate:1/DS_SamplingRate:AllTrials.Timestamp(end);

var_names=AllTrials.Properties.VariableNames;
for ii=1:size(AllTrials,2)
    var_DS=interp1(AllTrials.Timestamp,AllTrials.(ii),Bonsai_DS_Timestamps);
    AllTrials_cell_DS(:,ii)=var_DS';
end
DS_AllTrials=array2table(AllTrials_cell_DS,'VariableNames',var_names);
end

