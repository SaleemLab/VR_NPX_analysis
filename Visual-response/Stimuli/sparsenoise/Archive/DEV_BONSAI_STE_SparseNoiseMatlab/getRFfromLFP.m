function getRFfromLFP(thisAnimal,thisSession,thisFileNum,options)

% Function to get receptive field map from an LFP electrode

% AP 25/07/19

if ~exist('thisAnimal', 'var')
    thisAnimal = 'M19090';    
end
if ~exist('thisSession', 'var')
    thisSession = '190823';    
end
if ~exist('thisFileNum', 'var')
    thisFileNum = 6;
end

if ~exist('options', 'var')
   
%     serverPath = 'X:\ibn-vision\';
    serverPath = 'X:\';
    options.BonsaiPath = fullfile(serverPath,'DATA\SUBJECTS\',upper(thisAnimal),'Extras',upper(thisSession));
    %options.conToppath = fullfile(serverPath,'DATA\SUBJECTS\',upper(thisAnimal),'Ephys');
    options.conToppath = fullfile(serverPath,'DATA\SUBJECTS\',upper(thisAnimal),'Ephys',upper(thisSession));
    %options.matToppath = fullfile(serverPath,'DATA\PROJECTS\Tg4510\SRP\Acquisition\',upper(thisAnimal));
    options.matToppath = fullfile(serverPath,'DATA\SUBJECTS\',upper(thisAnimal),'Processed');
end


if ~isfield(options,'sync_input')
    options.sync_input = 'adc';  % adc for photodiote, ttl for TTL events
    options.sync_channel = 4;
    options.photo_br = 0; % Use 1 if photodiode signal is continious, otherwise 0 if it's flickering (not guaranteed to work)
end

% Include the wheel signal
options.includeWheel = 0; % 1
options.UseOneSignal = 1;
if ~isfield(options,'wheel_input') && options.includeWheel==1
    options.wheel_input = 'adc';  %
    options.wheel_channels = [3:4];
end

undersample = 15;

% Properties of the sparse noise stimulus
grid_size = 8; % 10x10
stim_duration = 0.1;
scal_f = 10; %Parameter for the smoothing
sigma = 5;%Parameter for the smoothing
options.probeMap = [1:32];
%LFP channel
options.ignoreElectrodes = []; % if you want to keep some electrodes out of the average used to reference, add that here
%options.chanNumbersToInclude = [9];
options.chanNumbersToInclude = [1:16];
% Signal processing options (built in)
options.multiplierDat = 1;      % invert output so that negative is positive etc by making this -1
options.referenceDat = 0;       % if you want to substract the average of the other electrodes, note flag here

% Ephys folder matching with thisSession
thisOEfile = [options.conToppath, '\', thisAnimal,'_*_', num2str(thisFileNum)];
folder_dir = dir(thisOEfile);
folder_name = folder_dir.name;
thisOEfileName = fullfile(options.conToppath,folder_name);

%%%%%%%
% Parse Ephys file
fprintf('\nParsing Open Ephys data files...\n')
tic;
[thisOEData,ttl,wheelData, fileinfo] = readCONandPDandWheel(thisOEfileName,options);
timetaken = toc;
fprintf('\n...took %3.2fs...\n',timetaken)


ttltimestamps = union(ttl.upPhases, ttl.downPhases);
ttltimestamps = ttltimestamps(1:end-1);
ttltimestamps = [997037; ttltimestamps];


% Peak detection/ get "spiketimes"
% findpeaks options (built in)
% options.findpeaksParams.MinPeakHeightSD = 3;     % At least 3SD of waveform
% options.findpeaksParams.MinPeakWidth = 10;         % At least 10 samples in peak (1/4.8 ms)
% options.findpeaksParams.MaxPeakWidth = 40;         % No more than 40 samples in peak
% options.findpeaksParams.MinPeakDistance = 10;      % Don't find peaks before 2 ms
% options.findpeaksParams.SortStr = 'descend';           % Sort the peaks
% 
% tsd = std(double(thisOEData));
% [pks,locs] = findpeaks(double(thisOEData), ...
%     'MinPeakHeight',  tsd*options.findpeaksParams.MinPeakHeightSD, ...
%     'MinPeakDistance', options.findpeaksParams.MinPeakDistance, ...
%     'SortStr', options.findpeaksParams.SortStr);
% %'MinPeakWidth', options.findpeaksParams.MinPeakWidth, ...
% %'MaxPeakWidth', options.findpeaksParams.MaxPeakWidth, ...
% 
% ts = locs/fileinfo.data_info.header.sampleRate;
                         
% Bonsai file
%thisBonsaiFile = [options.BonsaiPath, '\', thisAnimal,'_SparseNoise_Log_*'];
thisBonsaiFile = [options.BonsaiPath, '\', thisAnimal,'_SparseNoise_Log_' num2str(thisFileNum), '*'];
folder_dir = dir(thisBonsaiFile);
folder_name = folder_dir.name;
thisBonsaiFileName = fullfile(options.BonsaiPath,folder_name);

% Read the bin file
fileID=fopen(thisBonsaiFileName);
thisBinFile=fread(fileID);
fclose(fileID);

% Using a 10x10 grid
stim_matrix = reshape(thisBinFile, [grid_size, grid_size, length(thisBinFile)/grid_size/grid_size]);
stim_matrix = stim_matrix(:,:,1:end-1);
stim_matrix=flip(stim_matrix,1);
stim_matrix=flip(stim_matrix,2);
if length(ttltimestamps)<length(thisBinFile)/grid_size/grid_size
    warning('Number of ttls (%d) and stimulus instanses (%d) is not the same. Recording ended earlier.',length(ttltimestamps),length(thisBinFile)/grid_size/grid_size);
    % Assuming recording ended before the stimulus
    %stim_matrix = stim_matrix(:,:,1:length(ttltimestamps));
elseif length(ttltimestamps)<length(thisBinFile)/grid_size/grid_size
    warning('Number of ttls (%d) and stimulus instances (%d) not the same',length(ttltimestamps),length(thisBinFile)/grid_size/grid_size);
end

% Black squares
bl_matrix = zeros(1,length(thisBinFile));
bl_matrix(thisBinFile==0)=1;
bl_matrix(thisBinFile~=0)=0;
bl_matrix = reshape(bl_matrix,[grid_size, grid_size, length(thisBinFile)/grid_size/grid_size]);
bl_matrix=flip(bl_matrix,1);
bl_matrix=flip(bl_matrix,2);
if length(ttltimestamps)<length(thisBinFile)/grid_size/grid_size
    bl_matrix = bl_matrix(:,:,1:length(ttltimestamps));
end
% figure; imagesc(sum(bl_matrix,3)); axis off; axis square; colorbar
% title('black squares'); colormap hot

% White squares
wh_matrix = zeros(1,length(thisBinFile));
wh_matrix(thisBinFile==255)=1;
wh_matrix(thisBinFile~=255)=0;
wh_matrix = reshape(wh_matrix,[grid_size, grid_size, length(thisBinFile)/grid_size/grid_size]);
wh_matrix=flip(wh_matrix,1);
wh_matrix=flip(wh_matrix,2);
if length(ttltimestamps)<length(thisBinFile)/grid_size/grid_size
    wh_matrix = wh_matrix(:,:,1:length(ttltimestamps));
end
% figure; imagesc(sum(wh_matrix,3)); axis off; axis square; colorbar
% title('white squares'); colormap hot


% Get the lfp signal for each stimulus presentation
%samples to keep after tll pulse (equals to stimulus duration)
samples_to_keep = stim_duration*fileinfo.data_info.header.sampleRate;
lfp_data = nan(length(ttltimestamps),samples_to_keep);
for i=1:length(ttltimestamps)
   lfp_data(i,:) = thisOEData(1,ttltimestamps(i):ttltimestamps(i)+samples_to_keep-1);
end

%% All data: ON & OFF
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
av_lfp_norm = (av_lfp_reshaped-min(min(av_lfp_reshaped)))./(max(max(av_lfp_reshaped))-min(min(av_lfp_reshaped)));
av_lfp_norm = reshape(av_lfp_norm,[grid_size,grid_size, samples_to_keep]);
[U, S, V] = svd(av_lfp_reshaped);

figure('Name', 'All responses');
imLow = min(U(:,1));
imHigh = max(U(:,1));
[maxF,maxI] = max(abs(av_lfp_reshaped(:)));
signData = sign(av_lfp_reshaped(maxI));

[maxF,maxI] = max(abs(U(:)));
signU = sign(U(maxI));

if signU ~= signData
    U = -1*U;
end

I=reshape(U(:,1),[grid_size,grid_size]); %This is to smooth the image
J = imresize(I,scal_f);
A= imgaussfilt(J,sigma);
imagesc(A,[-0.5 0.5]); RedWhiteBlue;
axis xy

step = 0:1/samples_to_keep:1;
step = step(1:end-1);
hold on;
for i=1:grid_size
    for j=1:grid_size
        plot((step(1:undersample:end)+j-1)*10, squeeze((av_lfp_norm(i,j,1:undersample:end)+i-1)*10),'color','k');
    end
end
set(gca,'YDir','normal')
drawnow
setAxes

%% ON responses

av_lfp = zeros(grid_size,grid_size,samples_to_keep);
%figure;
for i=1:grid_size
    for j=1:grid_size
        av_lfp(i,j,:) = nanmean(lfp_data(wh_matrix(i,j,:)==1,:),1);
        %subplot(grid_size,grid_size,j+(i-1)*grid_size); plot(squeeze(av_lfp(i,j,:))); ylim([-300 200]); axis off; box off;
    end
end
av_lfp_reshaped = reshape(av_lfp,[grid_size*grid_size samples_to_keep]);
av_lfp_norm = (av_lfp_reshaped-min(min(av_lfp_reshaped)))./(max(max(av_lfp_reshaped))-min(min(av_lfp_reshaped)));
av_lfp_norm = reshape(av_lfp_norm,[grid_size,grid_size, samples_to_keep]);
[U, S, V] = svd(av_lfp_reshaped);

figure('Name', 'White responses');
[maxF,maxI] = max(abs(av_lfp_reshaped(:)));
signData = sign(av_lfp_reshaped(maxI));

[maxF,maxI] = max(abs(U(:)));
signU = sign(U(maxI));

if signU ~= signData
    U = -1*U;
end

I=reshape(U(:,1),[grid_size,grid_size]); %This is to smooth the image
J = imresize(I,scal_f);
A= imgaussfilt(J,sigma);

imagesc(A,[-0.5 0.5]); RedWhiteBlue;
axis xy

step = 0:1/samples_to_keep:1;
step = step(1:end-1);
hold on;
for i=1:grid_size
    for j=1:grid_size
        plot((step(1:undersample:end)+j-1)*10, squeeze((av_lfp_norm(i,j,1:undersample:end)+i-1)*10),'color','k');
    end
end
set(gca,'YDir','normal')
drawnow
setAxes

%% OFF responses
av_lfp = zeros(grid_size,grid_size,samples_to_keep);
%figure;
for i=1:grid_size
    for j=1:grid_size
        av_lfp(i,j,:) = nanmean(lfp_data(bl_matrix(i,j,:)==1,:),1);
        %subplot(grid_size,grid_size,j+(i-1)*grid_size); plot(squeeze(av_lfp(i,j,:))); ylim([-300 200]); axis off; box off;
    end
end
av_lfp_reshaped = reshape(av_lfp,[grid_size*grid_size samples_to_keep]);
av_lfp_norm = (av_lfp_reshaped-min(min(av_lfp_reshaped)))./(max(max(av_lfp_reshaped))-min(min(av_lfp_reshaped)));
av_lfp_norm = reshape(av_lfp_norm,[grid_size,grid_size, samples_to_keep]);
[U, S, V] = svd(av_lfp_reshaped);

%
figure('Name', 'Black responses');

[maxF,maxI] = max(abs(av_lfp_reshaped(:)));
signData = sign(av_lfp_reshaped(maxI));

[maxF,maxI] = max(abs(U(:)));
signU = sign(U(maxI));

if signU ~= signData
    U = -1*U;
end
I=reshape(U(:,1),[grid_size,grid_size]); %This is to smooth the image
J = imresize(I,scal_f);
A= imgaussfilt(J,sigma);
imagesc(-A,[-0.5 0.5]); RedWhiteBlue;
axis xy

step = 0:1/samples_to_keep:1;
step = step(1:end-1);
hold on;
for i=1:grid_size
    for j=1:grid_size
        plot((step(1:undersample:end)+j-1)*10, squeeze((av_lfp_norm(i,j,1:undersample:end)+i-1)*10),'color','k');
    end
end
set(gca,'YDir','normal')
drawnow
setAxes

% % %%%%%%%
% % % DUMMY for debugging
%Adjust for the flip due to the projector

% % % fileID=fopen('D:\Data\BV_SparseNoise\SparseNoise_Log_2019-08-16T16_41_10.1663744+01_00');
% % % thisBinFile=fread(fileID);
% % % fclose(fileID);
% % % Using a 10x10 grid
% % % stim_matrix = reshape(thisBinFile, [grid_size, grid_size, length(thisBinFile)/grid_size/grid_size]);
% imagesc(fliplr(stim_matrix(:,:,1)),[0 255]);


% set(gca,'YDir','normal')
% % 
%%
    function setAxes
        gridCentreX = 60;
        gridWidth=120;
        gridHeight= 120;
        gridCentreY = 30;
        squareSizeX = 15;
        squareSizeY = 15;
        centreLimsX = gridCentreX + (((gridWidth/2)-squareSizeX/2)*[-1 1]);
        xCentres = (centreLimsX(1):squareSizeX:centreLimsX(2));
        centreLimsY = gridCentreY + (((gridHeight/2)-squareSizeY/2)*[-1 1]);
        yCentres = (centreLimsY(1):squareSizeY:centreLimsY(2));
        
%         line(gridCentreX-[0 0], ylim);
%         line(xlim, gridCentreY-[0 0]);
        
        set(gca,'XTickLabels',xCentres)
        set(gca,'YTickLabels',yCentres)
        axis xy;
        set(gca, 'fontsize', 14, 'color','none', 'TickDir','out');
        % colormap gray
        % %
        % % % imagesc(fliplr(flipud(stim_matrix(:,:,1))),[0 255]);
        % % % stim_matrix = fliplr(flipud(stim_matrix(:,:,2)));
        % %
        % % set(gca,'YDir','normal')
        % % %sdf('StandardPlots');
        %Create
        start_i=1;
        for b=1:grid_size
            for a=1:grid_size
                matrix(a,b)=start_i;
                start_i=start_i+1;
            end
        end
        xlabel('Azimuth (^o)')
        ylabel('Elevation (^o)')
        colorbar
    end
end
