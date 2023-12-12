% Programme to estiamte receptive fields from sparse noise in NP1 / NP2
% data
% Dependencies: NPXAnalysis2022;
%               Visual-response-analysis/Stimuli/sparsenoise
% cd('/research/USERS/Masa/code')

% First load the data
if ismac
    ROOTPATH = '/Users/s.solomon/Filestore/Research2/ibn-vision';
else
    ROOTPATH = 'X:\ibn-vision';
%     ROOTPATH = '/research';
 
end

% Specify file
SUBJECT = 'M22069';%'M21200';
SESSION = '20221201';%'20211206';

% SESSION = '20221130';%'20211206';
% StimulusName = 'SparseNoise';
StimulusName = 'SparseNoise_fullscreen';

gs = [3];
nChannelsToBin = 24;    % Bin firing rates across channels for RF maps
channelRange = [234 339]; % Look at these channels 

% Some defaults
options.importMode = 'KS'; % LF or MUA or KS
% options.importMode = 'LF'; % LF or MUA or KS
options.BinWidth = 1/60; % resolution (in s) of output resps (e.g. 1/60)
options.stim_dur = 0.1;
options.AnalysisTimeWindow = [0 1/60*7];% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])
options.ks_unitType = 'good'; % 'mua', 'good' or ''

% Get correct paths
DATAPATH = fullfile(ROOTPATH,'DATA','SUBJECTS');
EPHYS_DATAPATH = fullfile(DATAPATH,SUBJECT,'ephys',SESSION);
BONSAI_DATAPATH = fullfile(DATAPATH,SUBJECT,'stimuli',SESSION);
% If Kilosort 
switch(options.importMode)
    case 'KS'
        options.KS_DATAPATH = fullfile(EPHYS_DATAPATH,'kilosort');
end

% Step 2: Find csv files associated with desired stimulus
[BehaviourDataFiles,EyeDataFiles,TrialParamsFiles,PDFiles] = getBonsaiFileNames(StimulusName,BONSAI_DATAPATH);

% Set bonsai data paths
PERIPHERALS_DATAPATH = fullfile(BONSAI_DATAPATH, BehaviourDataFiles{1});
EYEDATA_DATAPATH = fullfile(BONSAI_DATAPATH, EyeDataFiles{1});
TRIALDATA_DATAPATH = fullfile(BONSAI_DATAPATH,TrialParamsFiles{1});
PHOTODIODE_DATAPATH = fullfile(BONSAI_DATAPATH,PDFiles{1});

options.PERIPHERALS_DATAPATH = PERIPHERALS_DATAPATH; % full path to BehaviourData bonsai logfile
options.EYEDATA_DATAPATH = EYEDATA_DATAPATH; % full path to EyeTrackingData bonsai logfile
options.TRIALDATA_DATAPATH = TRIALDATA_DATAPATH; % full path to TrialData bonsai logfile
options.PHOTODIODE_DATAPATH = PHOTODIODE_DATAPATH; % full path to PDFiles bonsai logfile (Both sync pulse for spike and photodioide for quad)
options.PD_FLAG = 1;
options.paradigm = 'masa';



% Set ephys data path
options.gFileNum = gs(1);
folderName = findGFolder(EPHYS_DATAPATH,options.gFileNum);
T_EPHYS_DATAPATH = fullfile(EPHYS_DATAPATH,folderName,[folderName,'_imec0']);
options.EPHYS_DATAPATH = T_EPHYS_DATAPATH; % full path to ephys recording
options.MAP_FILE = fullfile(options.KS_DATAPATH,[SUBJECT,'_',SESSION,'_g0','_tcat.imec0.ap_kilosortChanMap.mat'])
options.KS_CATGT_FNAME = fullfile(['CatGT_',SUBJECT,'_',SESSION,'.log']);

% options.KS_CATGT_FNAME = 'CatGT_M22069_20221130.log';
% 'M22008_20220408_g0_tcat.imec0.ap_kilosortChanMap.mat'
%  
% Extract data
[resps,otherData,stimData,~,wheelData,photodiodeData,timeVector,options] = extractAndCollateNPData(options);
% [resps,~,stimData,~,~,~,timeVector,options] = extractAndCollateNPData(options);
%%
onDuration = photodiodeData.stim_off.sglxTime-photodiodeData.stim_on.sglxTime;
offDuration = photodiodeData.stim_on.sglxTime(2:end)-photodiodeData.stim_off.sglxTime(1:end-1);
% Calculate sparsenoise RFs
stim_matrix = cat(3,stimData.stim_matrix{:}); % N rows x M cols x nFrames
sn_options.grid_size = [size(stim_matrix,1) size(stim_matrix,2)];
sn_options.mapSampleRate = 60; % Hz
sn_options.mapsToShow = {'linear','black','white','contrast'};
sn_options.mapMethod = 'fitlm'; % fitlm mean
sn_options.framesToShow = [1 4 6 8 10 12];  % at 60 Hz would be 8 ms, 40, 72 etc
sn_options.plotflag = 1;

% Extract wheel and eye data
wheel_data = squeeze(nanmean(otherData(1,:,:),2))';
eye_data = squeeze(nanmean(otherData(2,:,:),2))';
% Choose a channel...or two
channelBins = channelRange(1):nChannelsToBin:channelRange(2);
channelBins = 200:30:290;%190-210 are about CA1, 100-140 are about CA3
channelBins = 260:40:300;%~260-320 are about V1

channelBins = 234:105:339;%~260-320 are about V1
SUA_data = [];
unit = 1;
On_response = [];
Off_response = [];
average_On = [];
average_Off = [];
average_contrast = [];
contrast_response = [];
On_grid_distribution = zeros(size(stim_matrix,[1 2]));
Off_grid_distribution = zeros(size(stim_matrix,[1 2]));

% x 240 degree [-120 120]
% y 120 degree [-30 90]
for x = 1:size(stim_matrix,1)
                    for y = 1:size(stim_matrix,2)
                        this_location = squeeze(stim_matrix(x,y,:));
                        On_trials = find(this_location==1); % find white quad
                        On_grid_distribution(x,y) = sum(this_location == 1);
                        Off_grid_distribution(x,y) = sum(this_location == -1);
                        
                    end
end
Overall_grid_distribution = On_grid_distribution + Off_grid_distribution;
figure;
subplot(3,1,1);imagesc(On_grid_distribution);colorbar;caxis([0 max(max(On_grid_distribution))])
subplot(3,1,2);imagesc(Off_grid_distribution);colorbar;caxis([0 max(max(Off_grid_distribution))])
subplot(3,1,3);imagesc(Overall_grid_distribution);colorbar;caxis([0 max(max(Overall_grid_distribution))])


for theseChannels = 1:length(channelBins)-1
    % Get average across channels (zscore?)
    switch(options.importMode)
        case 'KS'
            tp = options.peakChannel >= channelBins(theseChannels) & options.peakChannel < channelBins(theseChannels+1);
            unit_id = find(tp~=0);
            for n = 1:length(unit_id)
                SUA_data = squeeze(resps(unit_id(n),:,:))'; % returns 1 x N time bins x nFrames;
                sn_options.figName = sprintf('%s :: %s (file %01d) Channels %01d:%01d (unit %01d)',    SUBJECT, SESSION,  gs(1), channelBins(theseChannels),channelBins(theseChannels+1),unit_id(n));
               
                
                for x = 1:size(stim_matrix,1)
                    for y = 1:size(stim_matrix,2)
                        this_location = squeeze(stim_matrix(x,y,:));
                        On_trials = find(this_location==1); % find white quad
                        On_response(unit,x,y,:) = mean(SUA_data(On_trials,:)); % On response
                        average_On(unit,x,y,:) = mean(On_response(unit,x,y,:));
                        
                        
                        Off_trials = find(this_location==-1); % find black quad
                        Off_response(unit,x,y,:) = mean(SUA_data(Off_trials,:)); % Off response
                        average_Off(unit,x,y,:) = mean(Off_response(unit,x,y,:));
                        
                        
                        contrast_response(unit,x,y,:) = On_response(unit,x,y,:) - Off_response(unit,x,y,:);
                        average_contrast(unit,x,y,:) = mean(contrast_response(unit,x,y,:));
                    end
                end
                figure(unit)
                subplot(2,2,1)
                imagesc(squeeze(average_On(unit,:,:)))
                colorbar
                subplot(2,2,2)
                imagesc(squeeze(average_Off(unit,:,:)))
                colorbar
                subplot(2,2,3)
                imagesc(squeeze(average_contrast(unit,:,:)))
                colorbar
%                 
%                 
%                 initMap = sparseNoiseAnalysis(stim_matrix,SUA_data,wheel_data,eye_data,sn_options); % check to make sure this is the correct one - should be in Visua;-response-analysis/Stimuli/sparsenoise
%                 drawnow;
                unit = unit + 1;
            end
           
%             
%             MUA_data = mean(resps(tp,:,:),1); % returns 1 x N time bins x nFrames
%             MUA_data = squeeze(MUA_data)';   % sparsenoise needs frames x N time bins
%             
            sn_options.figName = sprintf('%s :: %s (file %01d) Channels %01d:%01d',    SUBJECT, SESSION,  gs(1), channelBins(theseChannels),channelBins(theseChannels+1));
%             initMap = sparseNoiseAnalysis(stim_matrix,SUA_data,wheel_data,eye_data,sn_options); % check to make sure this is the correct one - should be in Visua;-response-analysis/Stimuli/sparsenoise
%             drawnow;
            
        otherwise
            tp = channelBins(theseChannels):channelBins(theseChannels+1);
            lfp_data = mean(resps(tp,:,:),1); % returns 1 x N time bins x nFrames
            lfp_data = squeeze(lfp_data)';   % sparsenoise needs frames x N time bins
            sn_options.figName = sprintf('%s :: %s (file %01d) Channels %01d:%01d',    SUBJECT, SESSION,  gs(1), channelBins(theseChannels),channelBins(theseChannels+1));
%             initMap = sparseNoiseAnalysis(stim_matrix,lfp_data,wheel_data,eye_data,sn_options); % check to make sure this is the correct one - should be in Visua;-response-analysis/Stimuli/sparsenoise
            drawnow;
    end


end



%% Output sparsenoise maps for specific channels

onDuration = photodiodeData.stim_off.sglxTime-photodiodeData.stim_on.sglxTime;
offDuration = photodiodeData.stim_on.sglxTime(2:end)-photodiodeData.stim_off.sglxTime(1:end-1);
% Calculate sparsenoise RFs
stim_matrix = cat(3,stimData.stim_matrix{:}); % N rows x M cols x nFrames

sn_options.grid_size = [size(stim_matrix,1) size(stim_matrix,2)];
sn_options.mapSampleRate = 60; % Hz
sn_options.mapsToShow = {'linear','black','white','contrast'};
sn_options.mapMethod = 'fitlm'; % fitlm mean
sn_options.framesToShow = 1:7;  % at 60 Hz would be 8 ms, 40, 72 etc
sn_options.plotflag = 1; % 1 if you want to see details of the receptive fields

horizontalAngles = -30:(30+90)/8:90;
verticalAngles = -120:240/16:120;
% Extract wheel and eye data
wheel_data = squeeze(nanmean(otherData(1,:,:),2))';
eye_data = squeeze(nanmean(otherData(2,:,:),2))';
% Choose a channel...or two
channelBins = channelRange(1):nChannelsToBin:channelRange(2);
channelBins = 200:30:290;%190-210 are about CA1, 100-140 are about CA3
channelBins = 260:40:300;%~260-320 are about V1

channelBins = 260:50:310;%~260-320 are about V1
SUA_data = [];
unit = 1;
On_response = [];
Off_response = [];
average_On = [];
average_Off = [];
average_contrast = [];
contrast_response = [];
On_grid_distribution = zeros(size(stim_matrix,[1 2]));
Off_grid_distribution = zeros(size(stim_matrix,[1 2]));

% x 240 degree [-120 120]
% y 120 degree [-30 90]
for x = 1:size(stim_matrix,1)
                    for y = 1:size(stim_matrix,2)
                        this_location = squeeze(stim_matrix(x,y,:));
                        On_trials = find(this_location==1); % find white quad
                        On_grid_distribution(x,y) = sum(this_location == 1);
                        Off_grid_distribution(x,y) = sum(this_location == -1);
                        
                    end
end
Overall_grid_distribution = On_grid_distribution + Off_grid_distribution;
figure;
subplot(3,1,1);imagesc(On_grid_distribution');colorbar;caxis([0 max(max(On_grid_distribution))])
subplot(3,1,2);imagesc(Off_grid_distribution');colorbar;caxis([0 max(max(Off_grid_distribution))])
subplot(3,1,3);imagesc(Overall_grid_distribution');colorbar;caxis([0 max(max(Overall_grid_distribution))])

unit_no = 0;
for theseChannels = 1:length(channelBins)-1
    % Get average across channels (zscore?)
    switch(options.importMode)
        case 'KS'
            tp = options.peakChannel >= channelBins(theseChannels) & options.peakChannel < channelBins(theseChannels+1);
            unit_id = find(tp~=0);
            for n = 1:length(unit_id)
                unit_no = unit_no + 1;
                SUA_data = squeeze(resps(unit_id(n),:,:))'; % returns 1 x N time bins x nFrames;
                sn_options.figName = sprintf('%s :: %s (file %01d) Channels %01d:%01d (unit %01d)',    SUBJECT, SESSION,  gs(1), channelBins(theseChannels),channelBins(theseChannels+1),unit_id(n));
                initMap_temp = sparseNoiseAnalysis(stim_matrix,SUA_data,[],[],sn_options);
%                 initMap{unit_no,1} = initMap_temp;
%                 initMap_temp_black_max = max(initMap(:,:,:,2),[],3);
%                 initMap_temp_white_max = max(initMap(:,:,:,3),[],3);
%                 initMap_temp_contrast_max = max(initMap(:,:,:,4),[],3);
%                 figure;
%                 subplot(3,1,1)
%                 imagesc([-120 120],[-30 90],initMap_temp_black_max)
%                 colorbar; 
%                 subtitle('black')
%                 subplot(3,1,2)
%                 imagesc([-120 120],[-30 90],initMap_temp_white_max)
%                 colorbar;
%                 subtitle('white')               
%                 subplot(3,1,3)
%                 imagesc([-120 120],[-30 90],initMap_temp_contrast_max)
%                 colorbar;
%                 subtitle('contrast')
            end
    end
end


%% Output all kilosort units

onDuration = photodiodeData.stim_off.sglxTime-photodiodeData.stim_on.sglxTime;
offDuration = photodiodeData.stim_on.sglxTime(2:end)-photodiodeData.stim_off.sglxTime(1:end-1);
% Calculate sparsenoise RFs
stim_matrix = cat(3,stimData.stim_matrix{:}); % N rows x M cols x nFrames

sn_options.grid_size = [size(stim_matrix,1) size(stim_matrix,2)];
sn_options.mapSampleRate = 60; % Hz
sn_options.mapsToShow = {'linear','black','white','contrast'};
sn_options.mapMethod = 'fitlm'; % fitlm mean
sn_options.framesToShow = 1:7;  % at 60 Hz would be 8 ms, 40, 72 etc
sn_options.plotflag = 0; % 1 if you want to see details of the receptive fields

horizontalAngles = -30:(30+90)/8:90;
verticalAngles = -120:240/16:120;
% Extract wheel and eye data
wheel_data = squeeze(nanmean(otherData(1,:,:),2))';
eye_data = squeeze(nanmean(otherData(2,:,:),2))';


channelBins = 234:105:339;%~260-320 are about V1
SUA_data = [];
unit = 1;
On_response = [];
Off_response = [];
average_On = [];
average_Off = [];
average_contrast = [];
contrast_response = [];
On_grid_distribution = zeros(size(stim_matrix,[1 2]));
Off_grid_distribution = zeros(size(stim_matrix,[1 2]));

% x 240 degree [-120 120]
% y 120 degree [-30 90]
for x = 1:size(stim_matrix,1)
                    for y = 1:size(stim_matrix,2)
                        this_location = squeeze(stim_matrix(x,y,:));
                        On_trials = find(this_location==1); % find white quad
                        On_grid_distribution(x,y) = sum(this_location == 1);
                        Off_grid_distribution(x,y) = sum(this_location == -1);
                        
                    end
end
Overall_grid_distribution = On_grid_distribution + Off_grid_distribution;
% check on distribution of stimuli on the grid
figure;
subplot(3,1,1);imagesc(On_grid_distribution');colorbar;caxis([0 max(max(On_grid_distribution))])
subplot(3,1,2);imagesc(Off_grid_distribution');colorbar;caxis([0 max(max(Off_grid_distribution))])
subplot(3,1,3);imagesc(Overall_grid_distribution');colorbar;caxis([0 max(max(Overall_grid_distribution))])

unit_no = 0;
for unit_id = 1:size(resps,1)


                SUA_data = squeeze(resps(unit_id,:,:))'; % returns 1 x N time bins x nFrames;
                initMap_temp = sparseNoiseAnalysis(stim_matrix,SUA_data,[],[],sn_options);
                initMap{unit_id,1} = initMap_temp;


end

save(fullfile(EPHYS_DATAPATH,'analysis','receptiveFields.mat'),'initMap')