%% sparse noise analysis

% Get correct paths
ROOTPATH = 'X:\ibn-vision';
DATAPATH = fullfile(ROOTPATH,'DATA','SUBJECTS');
EPHYS_DATAPATH = fullfile(DATAPATH,SUBJECT,'ephys',SESSION);
BONSAI_DATAPATH = fullfile(DATAPATH,SUBJECT,'stimuli',SESSION);

BonsaiDir = BONSAI_DATAPATH;
sparseNoiseBinFile = dir(fullfile(BonsaiDir,['SparseNoise_fullscreen', '*' , '.bin']));

bonsai_files = dir(fullfile(BonsaiDir,['*', 'SparseNoise_fullscreen', '*']));
bonsai_files_names = {bonsai_files.name};
BehaviourDataFile = bonsai_files_names(contains(bonsai_files_names,'Wheel'));
EyeDataFile = bonsai_files_names(contains(bonsai_files_names,'Eye'));


% [BehaviourDataFiles,EyeDataFiles,TrialParamsFiles] = getBonsaiFileNames(StimulusName,BONSAI_DATAPATH);

fid = fopen([BonsaiDir, '\', sparseNoiseBinFile.name]);
SparseNoise_lin = fread(fid,inf,'uint8');

SN_dimensions = [16 8];
stim_dims = [120 120];
stim_centre_pos = [60 30];
stim_dur = 0.1;
SparseNoise = reshape(SparseNoise_lin,...
    [SN_dimensions(1),...
    SN_dimensions(2),...
    length(SparseNoise_lin)/prod(SN_dimensions)]);
SparseNoise = flip(SparseNoise,1);
SparseNoise = flip(SparseNoise,2);


% frames with a black or white square presented in this xy position
for xpos=1:SN_dimensions(1)
    for ypos=1:SN_dimensions(2)
        BStimFrames{xpos, ypos}=find(SparseNoise(xpos,ypos,:)==0);
        WStimFrames{xpos, ypos}=find(SparseNoise(xpos,ypos,:)==255);
    end
end

peripherals = import_bonsai_peripherals_SN([BonsaiDir, '\' BehaviourDataFile{:}]);
gkey_idx = find(contains(gkey.StimulusName, 'SparseNoise'));


%% use async pulse to align timestamps from sglx and bonsai
glxSyncTimes = syncTimes_sglx;

idx = glxSyncTimes>=gkey.TimeStart(gkey_idx) &...
    glxSyncTimes<=gkey.TimeEnd(gkey_idx);
relevent_g_syncTimes_sglx = glxSyncTimes(idx);

%% convert peripherals time to sglx time

idx = find(diff(peripherals.Sync)==1);
syncTimes_bonsai = peripherals.Time(idx+1)./1000; % convert to seconds

[r, lags] = xcorr(diff(syncTimes_bonsai), diff(relevent_g_syncTimes_sglx));
[val, idx] = max(r);
best_lag = lags(idx);

nSyncOffset = -best_lag+1;
t_npix = relevent_g_syncTimes_sglx(nSyncOffset:end); % sync sglx times
t_bonsai = syncTimes_bonsai(1:end); %sync bonsai times
t_npix = t_npix(1:numel(t_bonsai));
peripherals.sglxTime = interp1(t_bonsai, t_npix, peripherals.Time./1000);

%% get start times of each frame
RepCounter_temp = peripherals.RepCounter;
first_idx = find(peripherals.RepCounter==0,1,'first'); % first trial onset
RepCounter_temp(1:first_idx-1) = nan;
nTrials = max(RepCounter_temp)+1; %0 indexing;

for itrial = 1:nTrials
    start_idx(itrial) = find(RepCounter_temp==itrial-1,1,'first');
end

trial_start_times = peripherals.sglxTime(start_idx);


%% get on and off maps
tic
for iunit = 1:numel(units)
    delay = 0;
    delays=-10:10:120; % -10 to 90ms delays
    
    OFFMap = nan(SN_dimensions(1), SN_dimensions(2));
    for idelay = 1:numel(delays)
        for xpos = 1:SN_dimensions(1)
            for ypos = 1:SN_dimensions(2)
                temp_counts = [];
                for itrial = 1:numel(BStimFrames{xpos,ypos})
                    temp_counts = [temp_counts,...
                        sum(units(iunit).spiketimes > trial_start_times(BStimFrames{xpos,ypos}(itrial)) + delays(idelay) &...
                        units(iunit).spiketimes < trial_start_times(BStimFrames{xpos,ypos}(itrial)) + delays(idelay) +  0.1)];
                    
                end
                OFFMap(xpos,ypos,idelay) = mean(temp_counts)./0.1;
            end
        end
    end
    
    [maxVar, maxInd, vars] = maxVarMap(OFFMap);
    units(iunit).OFFMap = OFFMap(:,:,maxInd);
    
    ONMap = nan(SN_dimensions(1), SN_dimensions(2));
    for idelay = 1:numel(delays)
        for xpos = 1:SN_dimensions(1)
            for ypos = 1:SN_dimensions(2)
                temp_counts = [];
                for itrial = 1:numel(WStimFrames{xpos,ypos})
                    temp_counts = [temp_counts,...
                        sum(units(iunit).spiketimes > trial_start_times(WStimFrames{xpos,ypos}(itrial)) + delays(idelay) &...
                        units(iunit).spiketimes < trial_start_times(WStimFrames{xpos,ypos}(itrial)) + delays(idelay) +  0.1)];
                    
                end
                ONMap(xpos,ypos,idelay) = mean(temp_counts)./0.1;
            end
        end
    end
        [maxVar, maxInd, vars] = maxVarMap(ONMap);

        units(iunit).ONMap = ONMap(:,:,maxInd);
    
    
    
    % plot best OFF map
%     subplot(121)
%     [maxVar, maxInd, vars] = maxVarMap(OFFMap);
%     scal_f = 10; sigma = 3;
%     I = OFFMap(:,:,maxInd);
%     J = imresize(I,scal_f);
%     J = imgaussfilt(J,sigma);
%     imagesc(J)
%     colormap(gray), colorbar
%     title(maxInd)
%     % plot ON map
%     subplot(122)
%     [maxVar, maxInd, vars] = maxVarMap(ONMap);
%     scal_f = 10; sigma = 3;
%     I = ONMap(:,:,maxInd);
%     J = imresize(I,scal_f);
%     J = imgaussfilt(J,sigma);
%     imagesc(J)
%     colormap(gray), colorbar
%     title(maxInd)
%     pause
%     close
end
toc

