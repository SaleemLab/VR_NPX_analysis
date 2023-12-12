%% Tomaso Muzzu - UCL - 13/08/2019

%% Main script to generate a file with synchronised data of visual stimuli, behaviour and ephys
% 1) Get data from csv and bin files saved from Bonsai
% 2) Get photodiode signal
% 3) Sync it with the timestamps of Bonsai and save these info
% 4) Save spikes, single unit info and compute FR
% 5) Chunk spikes per frame --MM added 2019-09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function ES = SparseNoiseSyncR(mouseName)
% 1) Get data from csv files saved from Bonsai
% addpath(genpath('X:\CODE\DEV\general\SparseNoiseSyncR'))
% clear all
close all
CurrentFolder = cd;
%DataFolder = [CurrentFolder(1:end-(length(CurrentFolder)-3)) 'Archive - saleemlab' filesep 'Data' ];
[~,hostName] = system('hostname'); hostName = hostName(1:end-1);
if ~strcmp(hostName, 'saleem12')
    DataFolder = 'X:\DATA\SUBJECTS';
else
    if ismac
        DataFolder = '/Users/s.solomon/Desktop/TEMP';
    else
        DataFolder = 'X:\ibn-vision\DATA\SUBJECTS';
    end
end

if nargin>=1
    FileName = uigetfile_n_dir([DataFolder filesep mouseName],'Select BIN and CSV files of Sparse Noise');
else
    FileName = uigetfile_n_dir(DataFolder,'Select BIN and CSV files of Sparse Noise');
end
FileName
if size(FileName,2)==2
    % get info from CSV file
    csv_pos = find([strcmp(FileName{1,1}(end-2:end),'csv') strcmp(FileName{1,2}(end-2:end),'csv')]);
    FileName{csv_pos,1}
    fid = fopen(FileName{1,csv_pos});
    Columns = textscan(fid,'%s',1);
    fclose(fid);
    VisStimLog_Header = strsplit(Columns{1,1}{1,1},',');
    wheel_pos = find(strcmp(VisStimLog_Header,'Wheel'));
    VisStimLog = csvread(FileName{1,csv_pos},1,0);
    if sum(sign(diff(VisStimLog(:,wheel_pos))))<0 % adjustment for jan-march '19 when the rotary encoder was moved on the right.
        VisStimLog(:,wheel_pos) = -VisStimLog(:,wheel_pos)+max(VisStimLog(:,wheel_pos));
    end
    clear ES
    clear BehavData
    BehavData = array2table(VisStimLog,'VariableNames',VisStimLog_Header(1:end));
    % get info from bin file
    bin_pos = find([~strcmp(FileName{1,1}(end-2:end),'csv') ~strcmp(FileName{1,2}(end-2:end),'csv')]);
else
    bin_pos = 1;
end
fid = fopen(FileName{1,bin_pos});
SparseNoise_lin = fread(fid,inf,'uint8');
fclose(fid);
SN_dimensions = input('Please insert dimensions of Sparse Noise grid [rows columns] \n');
%SN_dimensions = [8 8];
stim_dims = input('Please insert angular span of stimulus in degrees [azimuth longitude] \n');
%stim_dims = [120 120];
screen_pos = input('Please insert centre position of display in degrees [azimuth longitude] \n');
%stim_centre_pos = [60 30];
stim_centre_pos = screen_pos;
SparseNoise = reshape(SparseNoise_lin,...
                     [SN_dimensions(1),...
                      SN_dimensions(2),...
                      length(SparseNoise_lin)/prod(SN_dimensions)]);
SparseNoise = flip(SparseNoise,1);
SparseNoise = flip(SparseNoise,2);

% figure
% colormap('gray')
% for i=1:size(SparseNoise,3)
%     figure
%     title(num2str(i))
%     colormap('gray')
%     imagesc(SparseNoise(:,:,i))
%     pause(1)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 2) Get analog signals
FileSep_i = strfind(FileName{1},filesep);
EphysDataFolder = [FileName{1}(1:FileSep_i(end-2)) 'ePhys'];
ACInfo = OE_get_AC_Signals(EphysDataFolder);
ACInfo = OE_Analyse_PD_Signal_SN(ACInfo,0.050);
Rows2Del = input('How many ACInfo.PD_changes do you want to remove from the end?\n');
ACInfo.PD_changes = ACInfo.PD_changes(1:end-Rows2Del);
%ACInfo = OE_Compute_Speed_Signal(ACInfo);
% figure
% plot(diff(ACInfo.PD_changes),'.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 3) Synchronise Bonsai data with the timing of OpenEphys
clear ES
ES.ACInfo = ACInfo;

if exist('BehavData','var')
    ES.Behav = BehavData;
else
    ES.Behav = [];
end

ES.SN_sequence = SparseNoise;
ES.SN_onsets = ACInfo.PD_changes(1:end-1);
ES.SN_singleI_dur = mean(diff(ACInfo.PD_changes));
ES.SN_gridDims = SN_dimensions;
ES.stim_dims = stim_dims;
ES.stim_centre_pos = stim_centre_pos;

if ES.SN_gridDims == 16 % just for Tomaso's recordings in July 2019
    ES.SN_square_crop=[3:10;1:8];
    ES.SN_square_dims = (ES.stim_dims*2)./ES.SN_gridDims;
    ES.ViewPortAdj = [3/16, 10/16; 1/16, 8/16] ;%[1/2, 1; 3/8, 7/8];
else
    ES.SN_square_dims = (ES.stim_dims)./ES.SN_gridDims;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 4) Select neural signal to sync (LFP or spikes)
clear SpikeData
kilosort = 1;
if kilosort == 1
    %     RecordingFolder = [EphysDataFolder filesep 'kilosort'];
    RecordingFolder = [EphysDataFolder];
    SpikeData = get_ks_Spikes(RecordingFolder);
else
    SpikeData = getSpikes(EphysDataFolder);
end
ES.MetaData = SpikeData.MetaData;

r=0; i=1;
while r==0
    if strcmp(ES.MetaData.FoldersList{i}(end-length(ES.ACInfo.ExpDate)+1:end),ES.ACInfo.ExpDate)
        r=i;
    end
    i=i+1;
end
ES.ROI = r;

ES.SpikeInfo = SpikeData.SpikeInfo;

if r==1
    stim_lims=[0 ES.MetaData.lims(r)-1]/ES.ACInfo.SamplingRateOE;
else
    stim_lims=[sum(ES.MetaData.lims(1:r-1)) sum(ES.MetaData.lims(1:r))]/ES.ACInfo.SamplingRateOE;
end

k=1;
for icell=2:size(ES.SpikeInfo,2)
    st=ES.SpikeInfo{1,icell};
    ES.StimSpiketimes{k}=st(st>stim_lims(1) & st<=stim_lims(2))-stim_lims(1);
    k=k+1;
end

%%
% 5) Count spikes per frame
numframes=length(ES.SN_onsets);
delays=-10:10:120; % -10 to 90ms delays

for icell = 1:length(ES.StimSpiketimes)
    st = ES.StimSpiketimes{icell};
    
    for delayInd=1:length(delays)
        for iframe = 1:numframes
            startTime=ES.SN_onsets(iframe);
            if iframe==numframes
                endTime=ES.SN_onsets(iframe)+ES.SN_singleI_dur;
            else
                endTime=ES.SN_onsets(iframe+1);
            end
            ES.frameSpikeCount{icell}(iframe, delayInd) = sum(st>(startTime + delays(delayInd)) & st<=(endTime + delays(delayInd))); %spike count per frame, per delay
        end
    end
end

%%
% 6) Save altogether

% Save all these in a file (this will be added with the spikes later)
% SavingDir = uigetfile_n_dir([FileName{1}(1:FileSep_i(length(FileSep_i)-1))],'Choose save directory');
if nargin>=1
    SavingDir = [DataFolder filesep mouseName filesep 'Processed'];
else
    SavingDir = uigetfile_n_dir([FileName{1}(1:FileSep_i(length(FileSep_i)-2))],'Choose save directory');
end
if kilosort == 1
    if iscell(SavingDir)
        save([SavingDir{1,1} filesep FileName{1}(FileSep_i(end-3)+1:FileSep_i(end-2)-1) '_SN_AC_Spikes_ks_' ES.ACInfo.ExpDate '.mat'], 'ES', '-v7.3');
    else
        save([SavingDir filesep FileName{1}(FileSep_i(end-3)+1:FileSep_i(end-2)-1) '_SN_AC_Spikes_ks_' ES.ACInfo.ExpDate '.mat'], 'ES', '-v7.3');
    end
else
    save([SavingDir{1,1} filesep SavingDir{1,1}(end-15:end-10) '_SN_AC_Spikes_' ES.ACInfo.ExpDate '.mat'], 'ES', '-v7.3');
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

