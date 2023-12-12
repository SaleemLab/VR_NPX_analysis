% function getRFfromLFP(thisAnimal,thisSession,thisFileNum,options)

% Function to get receptive field map from an LFP electrode

% AP 25/07/19

if ~exist('thisAnimal', 'var')
    thisAnimal = 'M19132';    
end
if ~exist('thisSession', 'var')
    thisSession = '20191130';    
end
if ~exist('thisFileNum', 'var')
    thisFileNum = 1;
end
if ~exist('options', 'var') || ~isfield(options,'serverPath')
    if ismac
        serverPath='/Users/s.solomon/Filestore/Research2/ibn-vision';
    elseif ispc
        if exist('X:\ibn-vision','dir')==7
            serverPath = 'X:\ibn-vision\';
        else
            serverPath = 'X:\';
        end
    end
end

%%%%%
% Set up paradigm
options.paradigm = 'sparsenoise';

%%%%%
% Set up paths
    options.BonsaiPath = fullfile(serverPath,'DATA','SUBJECTS',upper(thisAnimal),'Extras',upper(thisSession));
    options.conToppath = fullfile(serverPath,'DATA','SUBJECTS',upper(thisAnimal),'Ephys',upper(thisSession));
    options.matToppath = fullfile(serverPath,'DATA','SUBJECTS',upper(thisAnimal),'Processed');

if ~isfield(options,'sync_input')
    options.sync_input = 'adc';  % adc for photodiote, ttl for TTL events
    options.sync_channel = 1;
%     options.sync_input = 'ch';  % adc for photodiote, ttl for TTL events
%     options.sync_channel = 20;
    options.photo_br = 1; % Use 1 if photodiode signal is continious, otherwise 0 if it's flickering (not guaranteed to work)
    options.OE_pD_thresholds = [0.4 0.6];
end

% Include the wheel signal
options.includeWheel = 0; % 1
options.UseOneSignal = 1;
if ~isfield(options,'wheel_input') && options.includeWheel==1
    options.wheel_input = 'adc';  %
    options.wheel_channels = [3:4];
%     options.wheel_input = 'ch';  %
%     options.wheel_channels = [22:23];
end

% Properties of the sparse noise stimulus
grid_size = 10; % 10x10
stim_duration = 0.5;
options.stim_dur = stim_duration;
options.probeMap = [1:32];
%LFP channel
options.ignoreElectrodes = []; % if you want to keep some electrodes out of the average used to reference, add that here
% options.chanNumbersToInclude = [9];
options.chanNumbersToInclude = [17];
% Signal processing options (built in)
options.multiplierDat = 1;      % invert output so that negative is positive etc by making this -1
options.referenceDat = 0;       % if you want to substract the average of the other electrodes, note flag here
options.preStime = 120;

% Ephys folder matching with thisSession
thisOEfile = fullfile(options.conToppath, [ thisAnimal,'_*_', num2str(thisFileNum)]);
folder_dir = dir(thisOEfile);
folder_name = folder_dir.name;
thisOEfileName = fullfile(options.conToppath,folder_name);

%%%%%%%
% Parse Ephys file
fprintf('\nParsing Open Ephys data files...\n')
[thisOEData,ttl,wheelData, fileinfo] = readCONandPDandWheelNew(thisOEfileName,options);
ttltimestamps = union(ttl.upPhases, ttl.downPhases);

% Bonsai file
%thisBonsaiFile = [options.BonsaiPath, '\', thisAnimal,'_SparseNoise_Log_*'];
%thisBonsaiFile = [options.BonsaiPath, '\', thisAnimal,'_SparseNoise_Log_' num2str(thisFileNum), '*'];
% thisBonsaiFile = fullfile(options.BonsaiPath, [ thisAnimal,'_SparseNoise_Log_', '*']);
thisBonsaiFile = fullfile(options.BonsaiPath, [ thisAnimal,'_quaddata', '*']);
folder_dir = dir(thisBonsaiFile);
folder_name = folder_dir.name;
thisBonsaiFileName = fullfile(options.BonsaiPath,folder_name);

% Read the bin file
fileID=fopen(thisBonsaiFileName);
thisBinFile=fread(fileID);
fclose(fileID);

% Using a 10x10 grid
stim_matrix = zeros(1,length(thisBinFile));
stim_matrix(thisBinFile==0)=1;
stim_matrix(thisBinFile==255)=1;
stim_matrix(thisBinFile==128)=0;
stim_matrix = reshape(stim_matrix, [grid_size, grid_size, length(thisBinFile)/grid_size/grid_size]);
stim_matrix = stim_matrix(:,:,1:end-1);
if length(ttltimestamps)<size(stim_matrix,3)
    warning('Number of ttls (%d) LESS than number of stimulus instanses (%d).',length(ttltimestamps),length(thisBinFile)/grid_size/grid_size);
elseif length(ttltimestamps)>size(stim_matrix,3)
    warning('Number of ttls (%d) MORE than number of stimulus instanses (%d).',length(ttltimestamps),length(thisBinFile)/grid_size/grid_size);
end

% Black squares
bl_matrix = zeros(1,length(thisBinFile));
bl_matrix(thisBinFile==0)=1;
bl_matrix(thisBinFile~=0)=0;
bl_matrix = reshape(bl_matrix, [grid_size, grid_size, length(thisBinFile)/grid_size/grid_size]);
bl_matrix = bl_matrix(:,:,1:end-1);

% White squares
wh_matrix = zeros(1,length(thisBinFile));
wh_matrix(thisBinFile==255)=1;
wh_matrix(thisBinFile~=255)=0;
wh_matrix = reshape(wh_matrix, [grid_size, grid_size, length(thisBinFile)/grid_size/grid_size]);
wh_matrix = wh_matrix(:,:,1:end-1);


% Get the lfp signal for each stimulus presentation
%samples to keep after tll pulse (equals to stimulus duration)
samples_to_keep = stim_duration*fileinfo.data_info.header.sampleRate;
lfp_data = nan(length(ttltimestamps),samples_to_keep);
for i=1:length(ttltimestamps)
   lfp_data(i,:) = thisOEData(1,ttltimestamps(i):ttltimestamps(i)+samples_to_keep-1);
end

% Plot average lfp for each grid location
av_lfp = zeros(grid_size,grid_size,samples_to_keep);
%figure;
for i=1:grid_size
    for j=1:grid_size
        av_lfp(i,j,:) = nanmean(lfp_data(stim_matrix(i,j,:)==0 | stim_matrix(i,j,:)==255,:),1);
        %subplot(grid_size,grid_size,j+(i-1)*grid_size); plot(squeeze(av_lfp(i,j,:))); ylim([-300 200]); axis off; box off;
    end
end
av_lfp_reshaped = reshape(av_lfp,[grid_size*grid_size samples_to_keep]);
min_lfp = min(min(av_lfp_reshaped));
max_lfp = max(max(av_lfp_reshaped));
av_lfp_norm = (av_lfp_reshaped-min_lfp)./(max_lfp-min_lfp)+0.5;
av_lfp_norm = reshape(av_lfp_norm,[grid_size,grid_size, samples_to_keep]);
[U, S, V] = svd(av_lfp_reshaped);
[m i] = max(abs(U(:,1))); % Find the max to set max value to positive
figure; imagesc(reshape(U(:,1)*sign(U(i,1)),[grid_size,grid_size]));
step = 0.5:1/samples_to_keep:1.5;
step = step(1:end-1);
for i=1:grid_size
    for j=1:grid_size
        hold on; plot(step+j-1, squeeze(av_lfp_norm(i,j,:)+i-1),'color',[0.5 .5 .5]);
    end
end
title([thisAnimal,' ON+OFF Responses'])
set(gca,'YDir','normal')
%sdf('StandardPlots');

% ON responses
av_lfp = zeros(grid_size,grid_size,samples_to_keep);
%figure;
for i=1:grid_size
    for j=1:grid_size
        av_lfp(i,j,:) = nanmean(lfp_data(wh_matrix(i,j,:)==1,:),1);
        %subplot(grid_size,grid_size,j+(i-1)*grid_size); plot(squeeze(av_lfp(i,j,:))); ylim([-300 200]); axis off; box off;
    end
end
av_lfp_reshaped = reshape(av_lfp,[grid_size*grid_size samples_to_keep]);
av_lfp_norm = (av_lfp_reshaped-min_lfp)./(max_lfp-min_lfp)+0.5;
av_lfp_norm = reshape(av_lfp_norm,[grid_size,grid_size, samples_to_keep]);
[U, S, V] = svd(av_lfp_reshaped);
[m i] = max(abs(U(:,1))); % Find the max to set max value to positive
figure; imagesc(reshape(U(:,1)*sign(U(i,1)),[grid_size,grid_size]));
step = 0.5:1/samples_to_keep:1.5;
step = step(1:end-1);
for i=1:grid_size
    for j=1:grid_size
        hold on; plot(step+j-1, squeeze(av_lfp_norm(i,j,:)+i-1),'color',[0.5 .5 .5]);
    end
end
title([thisAnimal,' ON Responses'])
set(gca,'YDir','normal')
%sdf('StandardPlots');

% OFF responses
av_lfp = zeros(grid_size,grid_size,samples_to_keep);
%figure;
for i=1:grid_size
    for j=1:grid_size
        av_lfp(i,j,:) = nanmean(lfp_data(bl_matrix(i,j,:)==1,:),1);
        %subplot(grid_size,grid_size,j+(i-1)*grid_size); plot(squeeze(av_lfp(i,j,:))); ylim([-300 200]); axis off; box off;
    end
end
av_lfp_reshaped = reshape(av_lfp,[grid_size*grid_size samples_to_keep]);
av_lfp_norm = (av_lfp_reshaped-min_lfp)./(max_lfp-min_lfp)+0.5;
av_lfp_norm = reshape(av_lfp_norm,[grid_size,grid_size, samples_to_keep]);
[U, S, V] = svd(av_lfp_reshaped);
[m i] = max(abs(U(:,1))); % Find the max to set max value to positive
figure; imagesc(reshape(U(:,1)*sign(U(i,1)),[grid_size,grid_size]));
step = 0.5:1/samples_to_keep:1.5;
step = step(1:end-1);
for i=1:grid_size
    for j=1:grid_size
        hold on; plot(step+j-1, squeeze(av_lfp_norm(i,j,:)+i-1),'color',[0.5 .5 .5]);
    end
end
title([thisAnimal,' OFF Responses'])
set(gca,'YDir','normal')
%sdf('StandardPlots');

