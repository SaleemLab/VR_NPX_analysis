% Ripple detection and check MUA activity during ripple events
% Programme to estiamte receptive fields from sparse noise in NP1 / NP2
% data
% Dependencies: NPXAnalysis2022;
%               Visual-response-analysis/Stimuli/sparsenoise
%cd('/research/USERS/Masa/code')

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
% StimulusName = 'SparseNoise_fullscreen';
StimulusName = 'replay_Masa2tracks'
% StimulusName = 'Masa2tracks'
gs = [2];

% Get correct paths
DATAPATH = fullfile(ROOTPATH,'DATA','SUBJECTS');
EPHYS_DATAPATH = fullfile(DATAPATH,SUBJECT,'ephys',SESSION);
options.ANALYSIS_DATAPATH = fullfile(EPHYS_DATAPATH,'/analysis');

BONSAI_DATAPATH = fullfile(DATAPATH,SUBJECT,'stimuli',SESSION);
options.BONSAI_DATAPATH = BONSAI_DATAPATH;

options.KS_DATAPATH = fullfile(EPHYS_DATAPATH,'kilosort');
options.gFileNum = gs(1);
folderName = findGFolder(EPHYS_DATAPATH,options.gFileNum);
T_EPHYS_DATAPATH = fullfile(EPHYS_DATAPATH,folderName,[folderName,'_imec0']);
options.EPHYS_DATAPATH = T_EPHYS_DATAPATH; % full path to ephys recording
options.MAP_FILE = fullfile(options.KS_DATAPATH,[SUBJECT,'_',SESSION,'_g0','_tcat.imec0.ap_kilosortChanMap.mat'])
options.KS_CATGT_FNAME = fullfile(['CatGT_',SUBJECT,'_',SESSION,'.log']);

% Step 2: Find csv files associated with desired stimulus
% [BehaviourDataFiles,EyeDataFiles,TrialParamsFiles,PDFiles] = getBonsaiFileNames(StimulusName,BONSAI_DATAPATH);

import_and_align_Masa_VR_Bonsai(StimulusName,BONSAI_DATAPATH,options);

% 
% % Set bonsai data paths
% PERIPHERALS_DATAPATH = fullfile(BONSAI_DATAPATH, BehaviourDataFiles{1});
% EYEDATA_DATAPATH = fullfile(BONSAI_DATAPATH, EyeDataFiles{1});
% TRIALDATA_DATAPATH = fullfile(BONSAI_DATAPATH,TrialParamsFiles{1});
% PHOTODIODE_DATAPATH = fullfile(BONSAI_DATAPATH,PDFiles{1});
% 
% options.PERIPHERALS_DATAPATH = PERIPHERALS_DATAPATH; % full path to BehaviourData bonsai logfile
% options.EYEDATA_DATAPATH = EYEDATA_DATAPATH; % full path to EyeTrackingData bonsai logfile
% options.TRIALDATA_DATAPATH = TRIALDATA_DATAPATH; % full path to TrialData bonsai logfile
% options.PHOTODIODE_DATAPATH = PHOTODIODE_DATAPATH; % full path to PDFiles bonsai logfile (Both sync pulse for spike and photodioide for quad)
% options.PD_FLAG = 1;
% options.paradigm = 'masa';

% Set ephys data path

% options.KS_CATGT_FNAME = 'CatGT_M22069_20221130.log';
% 'M22008_20220408_g0_tcat.imec0.ap_kilosortChanMap.mat'
%  
% Extract data
% [resps,otherData,stimData,~,~,photodiodeData,timeVector,options] = extractAndCollateNPData(options);
% [resps,~,stimData,~,~,~,timeVector,options] = extractAndCollateNPData(options);


%% Masa LFP

[AP_FILE,LF_FILE] = findImecBinFile(options.EPHYS_DATAPATH);
% imecMeta = ReadMeta(fullfile(options.EPHYS_DATAPATH,LF_FILE));
imecMeta = NPadmin.ReadMeta(LF_FILE,options.EPHYS_DATAPATH);
chan_config = NPadmin.getNPChannelConfig(imecMeta);
switch str2double(imecMeta.imDatPrb_type)
    case 0 
    electrode_spacing_um = 20;
    otherwise
    electrode_spacing_um = 15;
end

% Set probe columns (4 columns for NPX1)
% (x coordinate = 11, 27, 43 or 59 micron)
shanks_available = unique(chan_config.Shank);
for n=1:size(shanks_available)
    cols_available = unique(chan_config.Ks_xcoord(chan_config.Shank==shanks_available(n)));
end

% Get channels from one column (x coordinate = 11, 27, 43 or 59 micron)
col_ID = cols_available(1); % (1) is 11
sorted_config = sortrows(chan_config,'Ks_ycoord','descend');
col_idx = sorted_config.Ks_xcoord == col_ID;
sorted_config = sorted_config(col_idx,:);

% Load Data
option.BinWidth = 1/1250;
BinWidth = option.BinWidth;
start_sec = 1;
start_samp = round(start_sec*imecMeta.imSampRate);
duration_sec = imecMeta.fileTimeSecs-2;
nSamp     = round(duration_sec*imecMeta.imSampRate);
chanTypes = str2double(strsplit(imecMeta.acqApLfSy, ','));
nEPhysChan = chanTypes(1);

downSampleRate = fix(imecMeta.imSampRate*BinWidth); 
% the rate at which we downsample depends on the acquisition rate and target binwidth
% Design low pass filter (with corner frequency determined by the
% desired output binwidth)
d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',(1/BinWidth)/2,'DesignMethod','butter','SampleRate',imecMeta.imSampRate);


% Read the data
clipDur = 100; % seconds
nClipSamps = round(imecMeta.imSampRate*clipDur);
nClips = floor(nSamp/nClipSamps);
samples_to_pass = 0;
% best channel for SW, best channel for spindle, best channel for ripple,
% best channel for cortical ripple
sorted_config.Channel
selected_channels = sorted_config.Channel;

raw_LFP = [];
for clip = 1:nClips
    tic
    disp(sprintf('loading clip %i',clip))
    temp_LFP = ReadNPXBin(start_samp+samples_to_pass, nClipSamps, imecMeta, LF_FILE, options.EPHYS_DATAPATH);
    temp_LFP(nEPhysChan+1:end,:) = []; % Get rid of sync channel
    temp_LFP = GainCorrectIM(temp_LFP, 1:nEPhysChan, imecMeta);

    % % Downsample the data
    raw_LFP = [raw_LFP downsample(temp_LFP(selected_channels,:)',downSampleRate)'];
    samples_to_pass = samples_to_pass + nClipSamps + 1;
    toc
end


% raw_LFP = downsample(raw_LFP',downSampleRate)';
new_SR = imecMeta.imSampRate/downSampleRate;
% raw_LFP1 = raw_LFP;

%% MUA

% [AP_FILE,LF_FILE] = findImecBinFile('X:\ibn-vision\DATA\SUBJECTS\M22069\ephys\20221201\M22069_20221201_g0\M22069_20221201_g0_imec0');

[AP_FILE,LF_FILE] = findImecBinFile(options.EPHYS_DATAPATH);

imecMeta = NPadmin.ReadMeta(AP_FILE,options.EPHYS_DATAPATH);
chan_config = NPadmin.getNPChannelConfig(imecMeta);

KS_metrics = readtable(fullfile(options.KS_DATAPATH,'metrics.csv'));
[these_spike_times,nominal_KSLabel,cluster_id,peakChannel,maxSpkTime] = import_ks_spiketimes(options.KS_DATAPATH,options.gFileNum,options.KS_CATGT_FNAME,imecMeta.imSampRate)
mean_waveforms = readNPY([options.KS_DATAPATH,'/mean_waveforms.npy']); % Information about mean waveform

% Sort channel for spike time data
col_ID = cols_available(1); % (1) is 11
sorted_config_spikes = sortrows(chan_config,'Ks_ycoord','descend');
col_idx = sorted_config_spikes.Ks_xcoord == col_ID;
sorted_config_spikes = sorted_config_spikes(col_idx,:);

% 
% good_idx = find(nominal_KSLabel=='good');
% mua_idx  = find(nominal_KSLabel=='mua');

MUA_filter_length = 41;
MUA_filter_alpha = 4;
time_step=1/new_SR; %match timestep to LFP time
time_bins_edges= tvec(1):time_step:max(tvec);
MUA = [];
SUA = [];
tic
for nchannel = 1:size(chan_config,1)
    clusters_this_channel = find(peakChannel == nchannel)-1;
    [~,index,]= intersect(cluster_id,clusters_this_channel);
    good_units_index = index(find(nominal_KSLabel(index)=='good'));
%     good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(index)=='good'));
   
    if ~isempty(index) % If any clusters
        if ~isempty(good_units_index) % If any SUA
             SUA(nchannel).spike_times= [];
             SUA(nchannel).spike_id = [];
            for unit = 1:length(good_units_index)

                SUA(nchannel).spike_times = [SUA(nchannel).spike_times; these_spike_times{good_units_index(unit)}];
                SUA(nchannel).spike_id = [SUA(nchannel).spike_id; good_units_index(unit)*ones(length(these_spike_times{good_units_index(unit)}),1)];
            end
        end
    end
end
toc


nfft_seconds = 2;
% nfft = 2^(nextpow2(SampRate(imecMeta)*nfft_seconds));
nfft = 2^(nextpow2(new_SR*nfft_seconds));

win  = hanning(nfft);
F  = [0.5 3;4 12;9 17;30 60;60 100;125 300];

for nchannel = 1:size(sorted_config,1)

    [pxx,fxx] = pwelch(raw_LFP(nchannel,:),win,[],nfft,new_SR);
    PSD(nchannel).power = pxx;
    PSD(nchannel).powerdB = 10*log10(pxx);
    PSD(nchannel).frequency = fxx;

    P = nan(size(pxx,2),size(F,1));
    f = nan(size(F,1),2);
    for n = 1:size(F,1)
        [~,f1idx] = min(abs(fxx-F(n,1)));
        [~,f2idx] = min(abs(fxx-F(n,2)));

        P(:,n) = mean(pxx(f1idx:f2idx,:));
        f(n,:) = [fxx(f1idx),fxx(f2idx)];
    end

    PSD(nchannel).mean_power = P;
    PSD(nchannel).frequency_range = f;
end

power = [];
for nchannel = 1:size(sorted_config,1)
    power(nchannel,:) = PSD(nchannel).mean_power;
end

% Find first peak that suddenly increases the power
power_differnece = [0; diff(power(:,6)./max(power(:,6)))]; % Use ripple or theta
[~,first_in_brain_channel] = findpeaks(power_differnece,'MinPeakHeight',0.1);
first_in_brain_channel = first_in_brain_channel(1);

% Quick plotting of PSD
colour_line= {'k','r','m','b','c','g','y'};
freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','350 - 625 Hz'};
selected_powers = [1 2 3 4 5 6];

figure
subplot(1,4,2)

for n = 1:length(cluster_boundary)
    plot([0 1],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([0 1],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for n = selected_powers
        p(n) = plot(power(:,n)./max(power(:,n)),sorted_config.Ks_ycoord',colour_line{n})
    
%     p(n) = plot(power_differnece,1:96,colour_line{n})
    hold on
end
% legend('1-3','4-12','30-60','60-100','125-300')
legend([p(selected_powers)],{freq_legends{selected_powers}});

% Distribution of cell types



%% Bonsai behaviour data

% Quick and dirty loading wheel data
cd(options.EPHYS_DATAPATH)
DIR = dir('*syncpulseTimes.mat')
if ~isempty(DIR)
    if exist(DIR.name) == 2
        load(DIR.name)
        syncTimes_ephys = syncTimes_ephys.on; % use upswings currently
    end
else
    [AP_FILE,LF_FILE] = findImecBinFile(options.EPHYS_DATAPATH);
    FILE_TO_USE = AP_FILE;
    binpath = fullfile(options.EPHYS_DATAPATH,AP_FILE);
%     binpath = fullfile(options.EPHYS_DATAPATH,LF_FILE);
    syncTimes_ephys = SGLXextractSyncPulseFromBinFile(binpath);
    parseDate = date;
    [~,fname] = fileparts(EPHYS_DATAPATH);
    save(fullfile(EPHYS_DATAPATH,[fname,'_syncpulseTimes.mat']),'syncTimes_ephys','parseDate');
    syncTimes_ephys = syncTimes_ephys.on; % use upswings currently
end

cd(BONSAI_DATAPATH)
% PDAasync = readtable('Masa2tracks_saveMousePos2022-11-30T13_50_24.csv');
% temp = readtable(wheel_file);
DIR = dir('Masa2tracks*')
fName={DIR(:).name}';
searchResult=cellfun(@findstr,fName,repmat({'_replay'},length(fName),1),'UniformOutput',false);
chosen=fName(cellfun(@isempty,searchResult));

stimuli_file = 'Masa2tracks_replay_saveMousePos2022-12-01T14_58_01.csv';
% stimuli_file = 'Masa2tracks_saveMousePos2022-12-01T13_41_28.csv';
temp = readtable(stimuli_file);
% temp = import_bonsai_peripherals(fullfile([BONSAI_DATAPATH,'\',stimuli_file]));
% Load bonsai data and asynch pulse
tic
b=regexp(table2array(temp(:,4)),'\d+(\.)?(\d+)?','match');
MousePos.pos = str2double([b{:}]);
MousePos.quad=table2array(temp(:,6));
MousePos.Sync=table2array(temp(:,5));
b=regexp(table2array(temp(:,7)),'\d+(\.)?(\d+)?','match');
MousePos.Time = str2double([b{:}])';
MousePos.pos = MousePos.pos(1:length(MousePos.Time))';
toc
[MousePos] = alignBonsaiToEphysSyncTimes(MousePos,syncTimes_ephys);

% eye_data_file = 'Masa2tracks_replay_EyeLog2022-12-01T14_58_01.csv';
% eyeData = import_bonsai_eyedata(fullfile([BONSAI_DATAPATH,'\',eye_data_file]));
% [peripherals] = alignBonsaiToEphysSyncTimes(eyeData,syncTimes_ephys);

wheel_file = 'Masa2tracks_replay_WheelLog2022-12-01T14_58_01.csv';
% wheel_file = 'Masa2tracks_WheelLog2022-12-01T13_41_28.csv'
peripherals = import_bonsai_peripherals(fullfile([BONSAI_DATAPATH,'\',wheel_file]));
[peripherals] = alignBonsaiToEphysSyncTimes(peripherals,syncTimes_ephys);

photodiode_file = 'Masa2tracks_replay_PDAsync2022-12-01T14_58_01.csv';
% photodiode_file = 'Masa2tracks_PDAsync2022-12-01T13_41_28.csv';
[photodiode, photodiode_sync] = import_bonsai_photodiode(fullfile([BONSAI_DATAPATH,'\',photodiode_file]));
[photodiode] = alignBonsaiToEphysSyncTimes(photodiode,syncTimes_ephys);

tick_to_cm_conversion = 3.1415*20/1024; % radius 20 cm and one circle is 1024 ticks;
speed = [0; diff(peripherals.Wheel*tick_to_cm_conversion)];
speed(speed<-100) = 0;
speed(speed>100) = 0;

% speed_interp = interp1(peripherals.sglxTime,speed,tvec','linear');
peripherals.speed = speed;

idx_trial_start = find(diff(peripherals.QuadState)==1) + 1;
idx_trial_end = find(diff(peripherals.QuadState)==-1) + 1;
frame_trial_start = peripherals.FrameNumber(idx_trial_start);
frame_trial_end = peripherals.FrameNumber(idx_trial_end);
temp_tbl.start_time = peripherals.sglxTime(idx_trial_start);
temp_tbl.end_time = peripherals.sglxTime(idx_trial_end);
nTrials = length(temp_tbl.start_time);
% get trial start/end times from photodiode signal
pd_thresh_up = 400;
pd_thresh_down = 100;
for itrial =1:length(temp_tbl.start_time);
    % find frame where Bonsai asked to change quad
    frameidx_start = find(photodiode.FrameNumber <= frame_trial_start(itrial),1,'last');
    % find next index where the photodiode detected a quad change
    temp_idx = find(photodiode.Photodiode(frameidx_start:end)>pd_thresh_up,1,'first');
    idx_start = frameidx_start+temp_idx-1;
    % convert to spike GLX time
    pdstart(itrial) = photodiode.sglxTime(idx_start);
end

for itrial =1:length(temp_tbl.end_time)
    % to get photodiode down we need to do some median filtering to
    % remove some awkward down swings when the quad is white
    frameidx_end = find(photodiode.FrameNumber >= frame_trial_end(itrial),1,'first');
    temp_idx = find(smoothdata(photodiode.Photodiode(frameidx_end:end),'movmedian', 5)<pd_thresh_down,1,'first');
    idx_end = frameidx_end+temp_idx-1;
    pdend(itrial) = photodiode.sglxTime(idx_end);
end

photodiodeData.stim_on.sglxTime = pdstart';
photodiodeData.stim_off.sglxTime = pdend';
[MousePos] = alignBonsaiToPhotodiode(MousePos,sort([photodiodeData.stim_on.sglxTime; photodiodeData.stim_off.sglxTime]),'replay');
save MousePos MousePos
save peripherals peripherals
% [pks,locs]= findpeaks(abs(MousePos.pos-200));
% 
% % figure
% % plot(MousePos.sglxTime,MousePos.pos)
% % hold on;
% % scatter(MousePos.sglxTime(locs),400*ones(1,length(locs)))
% % 
% % MousePos.stimuli_onset = MousePos.sglxTime_corrected(locs);
% % MousePos.stimuli_id = MousePos.pos(locs);
% % MousePos.stimuli_track = MousePos.pos(locs);
% % MousePos.stimuli_track(MousePos.pos(locs)<200) = 2;
% % MousePos.stimuli_track(MousePos.pos(locs)>200) = 1;


%% Detect Sleep

[AP_FILE,LF_FILE] = findImecBinFile(options.EPHYS_DATAPATH);
% imecMeta = ReadMeta(fullfile(options.EPHYS_DATAPATH,LF_FILE));
imecMeta = NPadmin.ReadMeta(LF_FILE,options.EPHYS_DATAPATH);
chan_config = NPadmin.getNPChannelConfig(imecMeta);

% spindle
parameters = list_of_parameters;
filter_type  = 'bandpass';
filter_width = [30 100];                 % range of frequencies in Hz you want to filter between
filter_order = round(6*new_SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
norm_freq_range = filter_width/(new_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_gamma = fir1(filter_order, norm_freq_range,filter_type);

% for nchannel = 1:size(raw_LFP,1)
%     LFP.gamma(nchannel,:) = filtfilt(b_gamma,1,raw_LFP(nchannel,:));
%     LFP.gamma_zscore(nchannel,:)  = zscore(abs(hilbert(LFP.gamma(nchannel,:))));
% end

tvec = start_sec:1/new_SR:start_sec+(length(raw_LFP(1,:))-1)/new_SR;
find(peripherals.sglxTime >= start_sec);
speed = peripherals.speed;
speed_interp = interp1(peripherals.sglxTime,speed,tvec','linear');
speedTreshold = 1;

[freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
    [tvec' raw_LFP(find(sorted_config.Channel == best_channels.spindle_channel),:)'],[tvec' raw_LFP(find(sorted_config.Channel == best_channels.ripple_channel),:)'],...
    [tvec' speed_interp],speedTreshold);


behavioural_state.freezing = freezing;
behavioural_state.quietWake = quietWake;
behavioural_state.SWS = SWS;
behavioural_state.REM = REM;
behavioural_state.movement = movement;

cd(options.EPHYS_DATAPATH)
save behavioural_state behavioural_state

%% group Spikes by regioon and Detect SWR events
% Determine channel with high ripple low theta
% F  = [0.5 3;4 12;9 17;30 60;60 100;125 300;350 625];

[AP_FILE,LF_FILE] = findImecBinFile(options.EPHYS_DATAPATH);
% imecMeta = ReadMeta(fullfile(options.EPHYS_DATAPATH,LF_FILE));
imecMeta = NPadmin.ReadMeta(AP_FILE,options.EPHYS_DATAPATH);
chan_config = NPadmin.getNPChannelConfig(imecMeta);

tic
tvec = start_sec:1/new_SR:start_sec+(length(raw_LFP(1,:))-1)/new_SR;
time_bins_edges= tvec(1):time_step:max(tvec);

spike_times = [];
spike_ID = [];
CA1_spike_times = [];

channels_selected = best_channels.CA1_channel-10 : best_channels.CA1_channel+10; % Putative CA3


for nchannel = channels_selected(1):channels_selected(end)
    clusters_this_channel = find(peakChannel == nchannel)-1; % Minus one because clutser id is 0 based.
    [~,index,]= intersect(cluster_id,clusters_this_channel);
    good_units_index = index(find(nominal_KSLabel(index)=='good'));
    %     if ~isempty(clusters_this_channel) % If any clusters
    if ~isempty(good_units_index) % If any SUA

        for unit = 1:length(good_units_index)
            spike_times = [spike_times; these_spike_times{good_units_index(unit)}];
            spike_ID = [spike_ID; good_units_index(unit)*ones(length(these_spike_times{good_units_index(unit)}),1)];
        end

    end
end

%smooth SUA activity with a gaussian kernel
MUA_filter_length = 41;
MUA_filter_alpha = 4;
time_step=1/new_SR; %match timestep to LFP time

w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
CA1_MUA_zscore =zscore(filtfilt(w,1,histcounts(spike_times,time_bins_edges)));
MUA_zscore = interp1(time_bins_edges(1:end-1)+time_step/2,CA1_MUA_zscore,LFP_tvec,'linear');
[~,index] = sort(spike_times);
CA1_spike_times(:,1) = spike_ID(index);
CA1_spike_times(:,2) = spike_times(index);
MUA_time = time_bins_edges(1:end-1)+time_step/2;
toc


tic
MUA_filter_length = 41;
MUA_filter_alpha = 4;
time_step=1/new_SR; %match timestep to LFP time
time_bins_edges= tvec(1):time_step:max(tvec);
spike_times = [];
spike_ID = [];
CA3_spike_times = [];

channels_selected = find(chan_config.Ks_ycoord >1100 & chan_config.Ks_ycoord < 1450); % Putative CA3
for nchannel = channels_selected(1):channels_selected(end)
    clusters_this_channel = find(peakChannel == nchannel)-1; % Minus one because clutser id is 0 based.
    [~,index,]= intersect(cluster_id,clusters_this_channel);
    good_units_index = index(find(nominal_KSLabel(index)=='good'));
    %     if ~isempty(clusters_this_channel) % If any clusters
    if ~isempty(good_units_index) % If any SUA

        for unit = 1:length(good_units_index)
            spike_times = [spike_times; these_spike_times{good_units_index(unit)}];
            spike_ID = [spike_ID; good_units_index(unit)*ones(length(these_spike_times{good_units_index(unit)}),1)];
        end

    end
end

%smooth SUA activity with a gaussian kernel
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
CA3_MUA_zscore=zscore(filtfilt(w,1,histcounts(spike_times,time_bins_edges)));
[~,index] = sort(spike_times);
CA3_spike_times(:,1) = spike_ID(index);
CA3_spike_times(:,2) = spike_times(index);
MUA_time = time_bins_edges(1:end-1)+time_step/2;
toc


tic
MUA_filter_length = 41;
MUA_filter_alpha = 4;
time_step=1/new_SR; %match timestep to LFP time
time_bins_edges= tvec(1):time_step:max(tvec);
spike_times = [];
spike_ID = [];
DG_spike_times = [];

channels_selected = find(chan_config.Ks_ycoord >1450 & chan_config.Ks_ycoord < 1740); % Putative CA3
for nchannel = channels_selected(1):channels_selected(end)
    clusters_this_channel = find(peakChannel == nchannel)-1; % Minus one because clutser id is 0 based.
    [~,index,]= intersect(cluster_id,clusters_this_channel);
    good_units_index = index(find(nominal_KSLabel(index)=='good'));
    %     if ~isempty(clusters_this_channel) % If any clusters
    if ~isempty(good_units_index) % If any SUA

        for unit = 1:length(good_units_index)
            spike_times = [spike_times; these_spike_times{good_units_index(unit)}];
            spike_ID = [spike_ID; good_units_index(unit)*ones(length(these_spike_times{good_units_index(unit)}),1)];
        end

    end
end


%smooth SUA activity with a gaussian kernel
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
DG_MUA_zscore=zscore(filtfilt(w,1,histcounts(spike_times,time_bins_edges)));
[~,index] = sort(spike_times);
DG_spike_times(:,1) = spike_ID(index);
DG_spike_times(:,2) = spike_times(index);
MUA_time = time_bins_edges(1:end-1)+time_step/2;
toc


tic
MUA_filter_length = 41;
MUA_filter_alpha = 4;
time_step=1/new_SR; %match timestep to LFP time
time_bins_edges= tvec(1):time_step:max(tvec);
spike_times = [];
spike_ID = [];
V1_spike_times = [];
% Select channels based on depth
channels_selected = ...
    chan_config.SpikeGLXchan0(find(chan_config.Ks_ycoord < chan_config.Ks_ycoord(best_channels.SW_channel)...
    & chan_config.Ks_ycoord >chan_config.Ks_ycoord(best_channels.SW_channel)-1000));

for nchannel = channels_selected(1):channels_selected(end)
    clusters_this_channel = find(peakChannel == nchannel)-1; % Minus one because clutser id is 0 based.
    [~,index,]= intersect(cluster_id,clusters_this_channel);
    good_units_index = index(find(nominal_KSLabel(index)=='good'));
    %     if ~isempty(clusters_this_channel) % If any clusters
    if ~isempty(good_units_index) % If any SUA

        for unit = 1:length(good_units_index)
            spike_times = [spike_times; these_spike_times{good_units_index(unit)}];
            spike_ID = [spike_ID; good_units_index(unit)*ones(length(these_spike_times{good_units_index(unit)}),1)];
        end

    end
end

%smooth SUA activity with a gaussian kernel
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
V1_MUA_zscore=zscore(filtfilt(w,1,histcounts(spike_times,time_bins_edges)));
[~,index] = sort(spike_times);
V1_spike_times(:,1) = spike_ID(index);
V1_spike_times(:,2) = spike_times(index);
MUA_time = time_bins_edges(1:end-1)+time_step/2;
toc



zscore_min = 0;
zscore_max = 3;

tvec = start_sec:1/new_SR:start_sec+(length(raw_LFP(1,:))-1)/new_SR;
% LFP_tvec =  start_sec:1/new_SR:start_sec+duration_sec;
speed_interp = interp1(peripherals.sglxTime,speed,tvec','linear');

cd(options.EPHYS_DATAPATH)
channel_to_use = find(sorted_config.Channel == best_channels.ripple_channel);
[replay,reactivations] = detect_candidate_events_masa(MUA_time,tvec,raw_LFP(channel_to_use,:),...
    CA1_MUA_zscore,CA1_spike_times,zscore_min,zscore_max,options) % Need to fix SUA spike
save extracted_candidate_events replay reactivations

% bz_FindPopBursts

% durations = [30 100]; inter-ripple 30 ms and max duration 100 ms 
% 'noise', 1 here is the top noisy channel (usually in dura gel)
% [ripples] = FindRipples_masa(raw_LFP(3,:)',tvec','minDuration',20,'durations',[30 100],'frequency',new_SR,...
%     'noise',raw_LFP(5,:)','passband',[125 300],'thresholds',[2 5],'show','off')
% save extracted_ripples_events ripples

channel_to_use = find(sorted_config.Channel == best_channels.ripple_channel);
[ripples] = FindRipples_masa(raw_LFP(channel_to_use,:)',tvec','minDuration',20,'durations',[30 200],'frequency',new_SR,...
    'noise',raw_LFP(2,:)','passband',[125 300],'thresholds',[3 5])

figure
[ripples.SWS_offset,ripples.SWS_index] = RestrictInts(ripples.offset,SWS);
ripples.SWS_onset = ripples.onset(ripples.SWS_index);
ripples.SWS_peaktimes = ripples.peaktimes(ripples.SWS_index);

histogram(abs(ripples.SWS_onset-ripples.SWS_offset)',30,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r','Normalization','probability')

[ripples.awake_offset,ripples.awake_index] = RestrictInts(ripples.offset,quietWake);
ripples.awake_onset = ripples.onset(ripples.awake_index);
ripples.awake_peaktimes = ripples.peaktimes(ripples.awake_index);
hold on
histogram(abs(ripples.awake_onset-ripples.awake_offset)',30,'FaceColor','b','FaceAlpha',0.5,'EdgeColor','b','Normalization','probability')
legend('NREM Ripples','awake Ripples')
ylabel('Probability')
xlabel('Duration (sec)')
save extracted_ripples_events ripples

% need to extract best ripple channel for V1
channel_to_use = find(sorted_config.Channel == best_channels.spindle_channel);
[V1_ripples] = FindRipples_masa(raw_LFP(channel_to_use,:)',tvec','minDuration',20,'durations',[30 200],'frequency',new_SR,'noise',raw_LFP(1,:)','passband',[125 300])
[V1_ripples.SWS_offset,V1_ripples.SWS_index] = RestrictInts(V1_ripples.offset,SWS);
V1_ripples.SWS_onset = V1_ripples.onset(V1_ripples.SWS_index);
V1_ripples.SWS_peaktimes = V1_ripples.peaktimes(V1_ripples.SWS_index);

[V1_ripples.awake_offset,V1_ripples.awake_index] = RestrictInts(V1_ripples.offset,quietWake);
V1_ripples.awake_onset = V1_ripples.onset(V1_ripples.awake_index);
V1_ripples.awake_peaktimes = V1_ripples.peaktimes(V1_ripples.awake_index);
save extracted_V1_ripples_events V1_ripples

figure
histogram(abs(V1_ripples.SWS_onset-V1_ripples.SWS_offset)',20,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r','Normalization','probability')
hold on
histogram(abs(V1_ripples.awake_onset-V1_ripples.awake_offset)',20,'FaceColor','b','FaceAlpha',0.5,'EdgeColor','b','Normalization','probability')
legend('NREM cortical Ripples','awake cortical Ripples')
ylabel('Probability')
xlabel('Duration (sec)')


[spindles] = FindSpindles_masa(raw_LFP(channel_to_use,:)',tvec','durations',[400 3000],'frequency',new_SR,'noise',raw_LFP(1,:)','passband',[9 17],'thresholds',[1.5 3])
[spindles.SWS_offset,spindles.SWS_index] = RestrictInts(spindles.offset,SWS);
spindles.SWS_onset = spindles.onset(spindles.SWS_index);
spindles.SWS_peaktimes = spindles.peaktimes(spindles.SWS_index);

save extracted_spindles_events spindles
% bz_FindRipples



%% Detect Up and Down states

% figure
% hold on
% plot(normalised_slow_wave,sorted_config.Ks_ycoord)
% plot(normalised_spindle,sorted_config.Ks_ycoord)
% 
% plot(normalised_ripple,sorted_config.Ks_ycoord)
% plot(normalised_theta,sorted_config.Ks_ycoord)
% plot(normalised_gamma_power,sorted_config.Ks_ycoord)
% plot([0.5 0.5],[0 4000],'--')
% legend('Normalised slow wave','Normalised theta','Normalised ripple')
% xlabel('Normalised power')
% ylabel('Depth (um)')

% Select channels based on depth
[AP_FILE,LF_FILE] = findImecBinFile(options.EPHYS_DATAPATH);
% imecMeta = ReadMeta(fullfile(options.EPHYS_DATAPATH,LF_FILE));
imecMeta = NPadmin.ReadMeta(AP_FILE,options.EPHYS_DATAPATH);
chan_config = NPadmin.getNPChannelConfig(imecMeta);

channels_selected = ...
    chan_config.SpikeGLXchan0(find(chan_config.Ks_ycoord < chan_config.Ks_ycoord(best_channels.SW_channel)...
    & chan_config.Ks_ycoord >chan_config.Ks_ycoord(best_channels.SW_channel)-1000));

% channels_selected = ...
%     find(chan_config.Ks_ycoord < chan_config.Ks_ycoord(best_channels.SW_channel)...
%     & chan_config.Ks_ycoord >chan_config.Ks_ycoord(best_channels.SW_channel)-1000);


% Create Buz-style spike variables
spikes = [];
spikes.UID = [];
spikes.times = [];
count = 1;
for nchannel = channels_selected(1):channels_selected(end)
    clusters_this_channel = find(peakChannel == nchannel)-1; % Minus one because clutser id is 0 based.
    [~,index,]= intersect(cluster_id,clusters_this_channel);
    good_units_index = index(find(nominal_KSLabel(index)=='good'));
    good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(index)=='good'));

    if ~isempty(good_units_this_channel)
        for unit = 1:length(good_units_this_channel)
            if ~isempty(these_spike_times{good_units_index(unit)})
                spikes.times{count} = these_spike_times{good_units_index(unit)}; % Plus one because clutser id is 0 based.
                spikes.UID = [spikes.UID; good_units_this_channel(unit)];
                count = count + 1;
            end
        end
        
    end
end

channel_to_use = find(sorted_config.Channel == best_channels.spindle_channel);
slow_waves= DetectSlowWaves_masa('time',tvec,'lfp',raw_LFP(channel_to_use,:)','NREMInts',behavioural_state.SWS,'spikes',spikes);

% DetectSlowWaves
% fake_sws = [21 160; 190 477; 510 570];
% SWS;

% save best_channels best_spindle_channel best_spindle_channel best_ripple_channel
save extracted_slow_waves slow_waves


%% 

%% Identify SWR modulated cells in CA1 DG CA3 and V1

bin_size = 0.02;
ripple_events = ripples.onset;
cell_id = unique(CA1_spike_times(:,1));
SWR_modulation = [];
for cell = 1:length(unique(CA1_spike_times(:,1)))

    spikes_this_cell = CA1_spike_times(find(CA1_spike_times(:,1) == cell_id(cell)),2);
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, ripple_events, [-0.3 0.3], bin_size);
    tic

    % For each ripple event, circular shift the spike train
    for nshuffle = 1:1000
        parfor event = 1 : length(ripple_events)
            temp_shuffled(event,:) = circshift(binnedArray(event,:),ceil(rand*size(binnedArray,2)),2);
        end
        shuffled_psth(nshuffle,:) = mean(temp_shuffled./bin_size);
    end

    for nshuffle = 1:1000
        psth_difference_shuffle(nshuffle) = mean((shuffled_psth(nshuffle,:) - mean(shuffled_psth)).^2);
    end
    
    psth_difference(cell) = mean((psth - mean(shuffled_psth)).^2);

    if psth_difference(cell) > prctile(psth_difference_shuffle,97.5)
        SWR_modulation(cell) = 1;
    elseif psth_difference(cell) < prctile(psth_difference_shuffle,2.5)
        SWR_modulation(cell) = -1;
    else
        SWR_modulation(cell) = 0;
    end
    
    toc
end
CA1_SWR_modulation = SWR_modulation;

SWR_cells = find(CA1_SWR_modulation==1);
spike_times = [];
spike_ID = [];

cell_id = unique(CA1_spike_times(:,1));
CA1_SWR_spike_times = [];
for cell = 1:length(SWR_cells)
    spike_times = [spike_times; CA1_spike_times(find(CA1_spike_times(:,1) == cell_id(SWR_cells(cell))),2)];
    spike_ID = [spike_ID; cell_id(SWR_cells(cell))*ones(length(CA1_spike_times(find(CA1_spike_times(:,1) == cell_id(SWR_cells(cell))))),1)];
end

%smooth SUA activity with a gaussian kernel
% w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
% w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
% CA1_SWR_MUA_zscore=zscore(filtfilt(w,1,histcounts(spike_times,time_bins_edges)));
[~,index] = sort(spike_times);
CA1_SWR_spike_times(:,1) = spike_ID(index);
CA1_SWR_spike_times(:,2) = spike_times(index);
% MUA_time = time_bins_edges(1:end-1)+time_step/2;


bin_size = 0.02;
ripple_events = ripples.onset;
cell_id = unique(V1_spike_times(:,1));
length(unique(V1_spike_times(:,1)))
SWR_modulation = [];
for cell = 1:length(cell_id)
    spikes_this_cell = V1_spike_times(find(V1_spike_times(:,1) == cell_id(cell)),2);
    mean_FR(cell) = length(spikes_this_cell)/(spikes_this_cell(end)-spikes_this_cell(1));
end

for cell = 1:length(unique(V1_spike_times(:,1)))

    spikes_this_cell = V1_spike_times(find(V1_spike_times(:,1) == cell_id(cell)),2);
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, ripple_events, [-0.3 0.3], bin_size);
    tic

    % For each ripple event, circular shift the spike train
    for nshuffle = 1:1000
        parfor event = 1 : length(ripple_events)
            temp_shuffled(event,:) = circshift(binnedArray(event,:),ceil(rand*size(binnedArray,2)),2);
        end
        shuffled_psth(nshuffle,:) = mean(temp_shuffled./bin_size);
    end

    for nshuffle = 1:1000
        psth_difference_shuffle(nshuffle) = mean((shuffled_psth(nshuffle,:) - mean(shuffled_psth)).^2);
    end
    
    psth_difference(cell) = mean((psth - mean(shuffled_psth)).^2);

    if psth_difference(cell) > prctile(psth_difference_shuffle,97.5)
        SWR_modulation(cell) = 1;
    elseif psth_difference(cell) < prctile(psth_difference_shuffle,2.5)
        SWR_modulation(cell) = -1;
    else
        SWR_modulation(cell) = 0;
    end
    
    toc
end
V1_SWR_modulation = SWR_modulation;


tic
MUA_filter_length = 41;
MUA_filter_alpha = 4;
time_step=1/new_SR; %match timestep to LFP time
time_bins_edges= tvec(1):time_step:max(tvec);
spike_times = [];
spike_ID = [];

cell_id = unique(V1_spike_times(:,1));
SWR_cells = find(V1_SWR_modulation==1);
V1_SWR_spike_times = [];
for cell = 1:length(SWR_cells)
    spike_times = [spike_times; V1_spike_times(find(V1_spike_times(:,1) == cell_id(SWR_cells(cell))),2)];
    spike_ID = [spike_ID; cell_id(SWR_cells(cell))*ones(length(V1_spike_times(find(V1_spike_times(:,1) == cell_id(SWR_cells(cell))))),1)];
end

%smooth SUA activity with a gaussian kernel
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
V1_SWR_MUA_zscore=zscore(filtfilt(w,1,histcounts(spike_times,time_bins_edges)));
[~,index] = sort(spike_times);
V1_SWR_spike_times(:,1) = spike_ID(index);
V1_SWR_spike_times(:,2) = spike_times(index);
MUA_time = time_bins_edges(1:end-1)+time_step/2;
toc


bin_size = 0.02;
ripple_events = ripples.onset;
cell_id = unique(CA3_spike_times(:,1));
SWR_modulation = [];
for cell = 1:length(unique(CA3_spike_times(:,1)))

    spikes_this_cell = CA3_spike_times(find(CA3_spike_times(:,1) == cell_id(cell)),2);
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, ripple_events, [-0.3 0.3], bin_size);
    tic

    % For each ripple event, circular shift the spike train
    for nshuffle = 1:1000
        parfor event = 1 : length(ripple_events)
            temp_shuffled(event,:) = circshift(binnedArray(event,:),ceil(rand*size(binnedArray,2)),2);
        end
        shuffled_psth(nshuffle,:) = mean(temp_shuffled./bin_size);
    end

    for nshuffle = 1:1000
        psth_difference_shuffle(nshuffle) = mean((shuffled_psth(nshuffle,:) - mean(shuffled_psth)).^2);
    end
    
    psth_difference(cell) = mean((psth - mean(shuffled_psth)).^2);

    if psth_difference(cell) > prctile(psth_difference_shuffle,97.5)
        SWR_modulation(cell) = 1;
    elseif psth_difference(cell) < prctile(psth_difference_shuffle,2.5)
        SWR_modulation(cell) = -1;
    else
        SWR_modulation(cell) = 0;
    end
    
    toc
end
CA3_SWR_modulation = SWR_modulation;

SWR_cells = find(CA3_SWR_modulation==1);
spike_times = [];
spike_ID = [];

cell_id = unique(CA3_spike_times(:,1));
CA1_SWR_spike_times = [];
for cell = 1:length(SWR_cells)
    spike_times = [spike_times; CA3_spike_times(find(CA3_spike_times(:,1) == cell_id(SWR_cells(cell))),2)];
    spike_ID = [spike_ID; cell_id(SWR_cells(cell))*ones(length(CA3_spike_times(find(CA3_spike_times(:,1) == cell_id(SWR_cells(cell))))),1)];
end

%smooth SUA activity with a gaussian kernel
% w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
% w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
% CA1_SWR_MUA_zscore=zscore(filtfilt(w,1,histcounts(spike_times,time_bins_edges)));
[~,index] = sort(spike_times);
CA3_SWR_spike_times(:,1) = spike_ID(index);
CA3_SWR_spike_times(:,2) = spike_times(index);
% MUA_time = time_bins_edges(1:end-1)+time_step/2;


bin_size = 0.02;
ripple_events = ripples.peaktimes;
cell_id = unique(DG_spike_times(:,1));
SWR_modulation = [];
for cell = 1:length(unique(DG_spike_times(:,1)))

    spikes_this_cell = DG_spike_times(find(DG_spike_times(:,1) == cell_id(cell)),2);
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, ripple_events, [-0.3 0.3], bin_size);
    tic

    % For each ripple event, circular shift the spike train
    for nshuffle = 1:1000
        parfor event = 1 : length(ripple_events)
            temp_shuffled(event,:) = circshift(binnedArray(event,:),ceil(rand*size(binnedArray,2)),2);
        end
        shuffled_psth(nshuffle,:) = mean(temp_shuffled./bin_size);
    end

    for nshuffle = 1:1000
        psth_difference_shuffle(nshuffle) = mean((shuffled_psth(nshuffle,:) - mean(shuffled_psth)).^2);
    end
    
    psth_difference(cell) = mean((psth - mean(shuffled_psth)).^2);

    if psth_difference(cell) > prctile(psth_difference_shuffle,97.5)
        SWR_modulation(cell) = 1;
    elseif psth_difference(cell) < prctile(psth_difference_shuffle,2.5)
        SWR_modulation(cell) = -1;
    else
        SWR_modulation(cell) = 0;
    end
    
    toc
end
DG_SWR_modulation = SWR_modulation;


SWR_cells = find(DG_SWR_modulation==1);
spike_times = [];
spike_ID = [];

cell_id = unique(DG_spike_times(:,1));
DG_SWR_spike_times = [];
for cell = 1:length(SWR_cells)
    spike_times = [spike_times; DG_spike_times(find(DG_spike_times(:,1) == cell_id(SWR_cells(cell))),2)];
    spike_ID = [spike_ID; cell_id(SWR_cells(cell))*ones(length(DG_spike_times(find(DG_spike_times(:,1) == cell_id(SWR_cells(cell))))),1)];
end

%smooth SUA activity with a gaussian kernel
% w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
% w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
% CA1_SWR_MUA_zscore=zscore(filtfilt(w,1,histcounts(spike_times,time_bins_edges)));
[~,index] = sort(spike_times);
DG_SWR_spike_times(:,1) = spike_ID(index);
DG_SWR_spike_times(:,2) = spike_times(index);
% MUA_time = time_bins_edges(1:end-1)+time_step/2;


save SWR_modulation V1_SWR_modulation CA1_SWR_modulation CA3_SWR_modulation DG_SWR_modulation


%% All cells or just SWR-modulated PSTH relationship

% CCG of SLow Waves and spikes
spike_data = V1_spike_times(:,2);
spike_data = these_spike_times{360};
spike_data = these_spike_times{59};
% for nchannel = [best_spindle_channel:1:best_SW_channel]
%     spike_data = [spike_data; SUA(nchannel).spike_times];
% end


spike_data = CA1_spike_times(:,2);
num_cell = length(unique(CA1_spike_times(:,1)));

spike_data = CA1_SWR_spike_times(:,2);
num_cell = length(unique(CA1_SWR_spike_times(:,1)));


event = slow_waves.ints.UP(:,1);
% event = slow_waves.ints.DOWN(:,1);
event = ripples.onset;
% event = slow_waves.ints.UP(:,1);
% event = spindles.onset;
% event = V1_ripples.onset;
% event = reactivations.onset;
% num_cell = 1;
% for types = 1:4
MUA_filter_length = 50;
SD_alpha = 5; %2 std width
MUA_filter_alpha = (MUA_filter_length-1)/SD_alpha;
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,4
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width

% event = MousePos.stimuli_onset(MousePos.stimuli_track == 1)
event = MousePos.stimuli_onset(MousePos.stimuli_track == 2)


figure
bin_size = 0.005;
spike_data = CA1_spike_times(:,2);
num_cell = length(unique(CA1_spike_times(:,1)));
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data, event, [-1 1], bin_size);
psth = filtfilt(w,1,mean(binnedArray./num_cell./bin_size)); % normalize to Hz
psth_se = filtfilt(w,1,std(binnedArray./num_cell./bin_size)./sqrt(length(event)));

subplot(4,2,1)
plot(rasterX, rasterY,'b')
hold on
plot([0 0],get(gca,'ylim'),'r-')
plot([0.1 0.1],get(gca,'ylim'),'b-')
% plot([0.5 0.5],get(gca,'ylim'),'r-')
ylabel('events');
xlabel('Spike time relative to event onset')
title('CA1')

subplot(4,2,2)
hold on
plot(bins,psth,'b')
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');
%     ylim([0 600])
hold on
plot([0 0],get(gca,'ylim'),'r-')
plot([0.1 0.1],get(gca,'ylim'),'b-')
% plot([0.5 0.5],get(gca,'ylim'),'r-')
ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')

% CA3
spike_data = CA3_spike_times(:,2);
num_cell = length(unique(CA3_spike_times(:,1)));

% spike_data = CA3_SWR_spike_times(:,2);
% num_cell = length(unique(CA3_SWR_spike_times(:,1)));
% num_cell = 1;
bin_size = 0.005;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data, event, [-1 1], bin_size);
psth = filtfilt(w,1,mean(binnedArray./num_cell./bin_size)); % normalize to Hz
psth_se = filtfilt(w,1,std(binnedArray./num_cell./bin_size)./sqrt(length(event)));

subplot(4,2,3)
plot(rasterX, rasterY,'b')
hold on
plot([0 0],get(gca,'ylim'),'r-')
plot([0.1 0.1],get(gca,'ylim'),'b-')
% plot([0.5 0.5],get(gca,'ylim'),'r-')
ylabel('events');
xlabel('Spike time relative to event onset')
title('CA3')

subplot(4,2,4)
hold on
plot(bins,psth,'b')
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');
%     ylim([0 600])
hold on
plot([0 0],get(gca,'ylim'),'r-')
plot([0.1 0.1],get(gca,'ylim'),'b-')
% plot([0.5 0.5],get(gca,'ylim'),'r-')
ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')


% DG
spike_data = DG_spike_times(:,2);
num_cell = length(unique(DG_spike_times(:,1)));


% spike_data = DG_SWR_spike_times(:,2);
% num_cell = length(unique(DG_SWR_spike_times(:,1)));
% num_cell = 1;
bin_size = 0.005;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data, event, [-1 1], bin_size);
psth = filtfilt(w,1,mean(binnedArray./num_cell./bin_size)); % normalize to Hz
psth_se = filtfilt(w,1,std(binnedArray./num_cell./bin_size)./sqrt(length(event)));

subplot(4,2,5)
plot(rasterX, rasterY,'k')
hold on
plot([0 0],get(gca,'ylim'),'r-')
plot([0.1 0.1],get(gca,'ylim'),'b-')
% plot([0.5 0.5],get(gca,'ylim'),'r-')
ylabel('events');
xlabel('Spike time relative to event onset')
title('DG')

subplot(4,2,6)
hold on
plot(bins,psth,'k')
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'k','FaceAlpha','0.3','LineStyle','none');
%     ylim([0 600])
hold on
plot([0 0],get(gca,'ylim'),'r-')
plot([0.1 0.1],get(gca,'ylim'),'b-')
% plot([0.5 0.5],get(gca,'ylim'),'r-')
ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')

% 
% V1

spike_data = V1_spike_times(:,2);
num_cell = length(unique(V1_spike_times(:,1)));


% spike_data = V1_SWR_spike_times(:,2);
% num_cell = length(unique(V1_SWR_spike_times(:,1)));

% num_cell = 1;
bin_size = 0.005;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data, event, [-1 1], bin_size);
psth = filtfilt(w,1,mean(binnedArray./num_cell./bin_size)); % normalize to Hz
psth_se = filtfilt(w,1,std(binnedArray./num_cell./bin_size)./sqrt(length(event)));

subplot(4,2,7)
plot(rasterX, rasterY,'k')
hold on
plot([0 0],get(gca,'ylim'),'r-')
plot([0.1 0.1],get(gca,'ylim'),'b-')
% plot([0.5 0.5],get(gca,'ylim'),'r-')
ylabel('events');
xlabel('Spike time relative to event onset')
title('V1')

subplot(4,2,8)
hold on
plot(bins,psth,'k')
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'k','FaceAlpha','0.3','LineStyle','none');
%     ylim([0 600])
hold on
plot([0 0],get(gca,'ylim'),'r-')
plot([0.1 0.1],get(gca,'ylim'),'b-')
% plot([0.5 0.5],get(gca,'ylim'),'r-')
ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')


sgtitle('DOWN-UP Transition')
sgtitle('CA1 Ripple Onset')
sgtitle('V1 Spindle Onset')
sgtitle('V1 ripple Onset')
sgtitle('CA1 bursting')
sgtitle('Track 2 Stimuli (Right)')
sgtitle('Track 1 Stimuli (Left)')
% end


%% Oscillatory events cross correlation
windows = [-1.1 1.1];

figure
subplot(2,3,1)
bin_size = 0.01
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(reactivations.onset, ripples.peaktimes, windows, bin_size);
psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
psth_se = filtfilt(w,1,std(binnedArray./bin_size)./sqrt(length(reactivations.onset)));

hold on
plot(bins,psth,'b')
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');
hold on
plot([0 0],get(gca,'ylim'),'r-')
title('CA1 bursting events during CA1 ripple events')

% 
% subplot(2,3,1)
% bin_size = 0.1
% [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(ripples.peaktimes,MousePos.stimuli_onset(MousePos.stimuli_track == 1), windows, bin_size);
% % psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
% % psth_se = filtfilt(w,1,std(binnedArray./bin_size)./sqrt(length(ripples.peaktimes)));
% psth =mean(binnedArray./bin_size); % normalize to Hz
% psth_se = std(binnedArray./bin_size)./sqrt(length(ripples.onset));
% 
% hold on
% plot(bins,psth,'b')
% hold on
% % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
% patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');
% hold on
% plot([0 0],get(gca,'ylim'),'r-')
% title('CA1 ripples relative to track 1 stimuli onset')
% xlabel('time relative to stimulus onset (sec)')
% ylabel('events rate (events/sec)')
% 
% subplot(2,3,2)
% bin_size = 0.1
% [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(ripples.peaktimes,MousePos.stimuli_onset(MousePos.stimuli_track == 2), windows, bin_size);
% % psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
% % psth_se = filtfilt(w,1,std(binnedArray./bin_size)./sqrt(length(ripples.peaktimes)));
% psth =mean(binnedArray./bin_size); % normalize to Hz
% psth_se = std(binnedArray./bin_size)./sqrt(length(ripples.onset));
% 
% hold on
% plot(bins,psth,'b')
% hold on
% % patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
% patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');
% hold on
% plot([0 0],get(gca,'ylim'),'r-')
% title('CA1 ripples relative to track 2 stimuli onset')
% xlabel('time relative to stimulus onset (sec)')
% ylabel('events rate (events/sec)')

subplot(2,3,2)
bin_size = 0.01
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(ripples.peaktimes,V1_ripples.peaktimes, windows, bin_size);
psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
psth_se = filtfilt(w,1,std(binnedArray./bin_size)./sqrt(length(ripples.onset)));

plot(bins,psth,'b')
hold on
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');
hold on
plot([0 0],get(gca,'ylim'),'r-')
title('CA1 ripple events during V1 ripple events')

bin_size = 0.01
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(ripples.peaktimes, slow_waves.ints.UP(:,1), windows, bin_size);
psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
psth_se = filtfilt(w,1,std(binnedArray./bin_size)./sqrt(length(ripples.onset)));
subplot(2,3,3)
hold on
plot(bins,psth,'b')
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');
hold on
plot([0 0],get(gca,'ylim'),'r-')
title('CA1 ripples events at DOWN-UP transition')


bin_size = 0.01
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(V1_ripples.peaktimes, slow_waves.ints.UP(:,1), windows, bin_size);
    psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    psth_se = filtfilt(w,1,std(binnedArray./bin_size)./sqrt(length(V1_ripples.onset)));

subplot(2,3,4)
hold on
plot(bins,psth,'b')
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');

hold on
plot([0 0],get(gca,'ylim'),'r-')
title('V1 Ripple at DOWN-UP transition')



bin_size = 0.01
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spindles.peaktimes, slow_waves.ints.DOWN(:,2), windows, bin_size);
psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
psth_se = filtfilt(w,1,std(binnedArray./bin_size)./sqrt(length(spindles.onset)));

subplot(2,3,5)
hold on
plot(bins,psth,'b')
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');

hold on
plot([0 0],get(gca,'ylim'),'r-')
title('Cortical Spindles at DOWN-UP transition')

bin_size = 0.01
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(ripples.peaktimes, spindles.peaktimes, windows, bin_size);
psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
psth_se = filtfilt(w,1,std(binnedArray./bin_size)./sqrt(length(ripples.onset)));
subplot(2,3,6)
hold on
plot(bins,psth,'b')
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');
hold on
plot([0 0],get(gca,'ylim'),'r-')
title('CA1 ripples events during V1 spindles')


%% power and phase of slow waves peri ripple


channels_selected = chan_config.SpikeGLXchan0(best_channels.ripple_channel)-10: chan_config.SpikeGLXchan0(best_channels.ripple_channel)+10;

% channels_selected = ...
%     find(chan_config.Ks_ycoord < chan_config.Ks_ycoord(best_channels.SW_channel)...
%     & chan_config.Ks_ycoord >chan_config.Ks_ycoord(best_channels.SW_channel)-1000);


% Create Buz-style spike variables
spikes = [];
spikes.UID = [];
spikes.times = [];
count = 1;
for nchannel = channels_selected(1):channels_selected(end)
    clusters_this_channel = find(peakChannel == nchannel)-1; % Minus one because clutser id is 0 based.
    [~,index,]= intersect(cluster_id,clusters_this_channel);
    good_units_index = index(find(nominal_KSLabel(index)=='good'));
    good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(index)=='good'));

    if ~isempty(good_units_this_channel)
        for unit = 1:length(good_units_this_channel)
            if ~isempty(these_spike_times{good_units_index(unit)})
                spikes.times{count} = these_spike_times{good_units_index(unit)}; % Plus one because clutser id is 0 based.
                spikes.UID = [spikes.UID; good_units_this_channel(unit)];
                count = count + 1;
            end
        end
        
    end
end


% Slow wave oscilation band
parameters = list_of_parameters;
filter_type  = 'bandpass';
filter_width = [0.5 3];                 % range of frequencies in Hz you want to filter between

filter_order = round(6*new_SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
norm_freq_range = filter_width/(new_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_SO = fir1(filter_order, norm_freq_range,filter_type);
% filteredLFP = filtfilt(b_SO,1,raw_LFP(2,:));
filteredLFP = [];
filteredLFP.timestamps = tvec;
filteredLFP.amp = abs(hilbert(filtfilt(b_SO,1,raw_LFP(2,:))))';        
filteredLFP.phase = angle(hilbert(filtfilt(b_SO,1,raw_LFP(2,:))))'; 
filteredLFP.samplingRate = new_SR;


Ripple_event_SW_phase = [];
Ripple_event_SW_power = [];
for event = 1:length(ripples.onset)
    tic
    this_event_LFP_time = find(InIntervals(tvec,[ripples.onset(event)-1 ripples.onset(event)+1]) == 1);

    Ripple_event_SW_phase(:,event) = filteredLFP.phase(this_event_LFP_time);
    Ripple_event_SW_power(:,event) = filteredLFP.amp(this_event_LFP_time);
    toc
end

figure
subplot(2,2,1)
timebin = -1:1/new_SR:1;
event_id = 1:length(ripples.onset);
imagesc(timebin,event_id,Ripple_event_SW_phase')
colorbar
title('Cortical slow wave phase around CA1 ripple')

subplot(2,2,2)
imagesc(timebin,event_id,Ripple_event_SW_power')
colorbar
title('Cortical slow wave power around CA1 ripple')

subplot(2,2,3)
plot(timebin,circ_mean(Ripple_event_SW_phase,[],2))
title('Cortical slow wave phase around CA1 ripple')


subplot(2,2,4)
plot(timebin,mean(Ripple_event_SW_power,2))
%     psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
se = std(Ripple_event_SW_power')./sqrt(length(ripples.onset));
patch([timebin fliplr(timebin)],[mean(Ripple_event_SW_power')+se fliplr(mean(Ripple_event_SW_power')-se)],'b','FaceAlpha','0.3','LineStyle','none');
title('Cortical slow wave power around CA1 ripple')


% [PowerPhaseRatemap,spikebinIDs] = bz_PowerPhaseRatemap(spikes,filteredLFP,'int',[slow_waves.ints.UP(:,1)-1 slow_waves.ints.UP(:,1)+1])


%% All SUA PSTH 

windows = [-0.5 1];
events = slow_waves.ints.UP(:,1);
% events = ripples.onset;
events = V1_ripples.onset;
% events = MousePos.stimuli_onset(MousePos.stimuli_track == 1)
% events = MousePos.stimuli_onset(MousePos.stimuli_track == 2)

cell_id = unique(V1_spike_times(:,1));
% cell_id = unique(V1_SWR_spike_times(:,1));

% events = ripples.onset;
V1_UP_psth= [];
bin_size = 0.001;
for cell = 1:length(cell_id)
    spikes_this_cell = V1_spike_times(find(V1_spike_times(:,1) == cell_id(cell)),2);
    mean_FR(cell) = length(spikes_this_cell)/(spikes_this_cell(end)-spikes_this_cell(1));

    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, events, windows, bin_size);
     psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    V1_UP_psth(cell,:) = psth;
end
[~, index ] = sort(mean_FR,'descend');
V1_FR_sorted = cell_id(index);

figure
subplot(2,2,1)
imagesc(bins,1:1:length(index),zscore(V1_UP_psth(index,:),0,2))
% imagesc(bins,1:1:length(index),zscore(V1_UP_psth(index,:),0,2))
ylabel('Cell (sorted by FR)')
xlabel('Time relative to UP state (s)')
cbar = colorbar
cbar.Label.String = 'Firing Rate (Hz)';
hold on
plot([0 0],get(gca,'ylim'),'r','LineWidth',5)
clim([0 4])
title('V1')


% events = slow_waves.ints.UP(:,1);
cell_id = unique(CA1_spike_times(:,1));
% cell_id = unique(CA1_SWR_spike_times(:,1));
% events = slow_waves.ints.UP(:,1);
% events = ripples.onset;
CA1_UP_psth= [];
mean_FR = [];

bin_size = 0.001;
for cell = 1:length(cell_id)
    spikes_this_cell = CA1_spike_times(find(CA1_spike_times(:,1) == cell_id(cell)),2);
    mean_FR(cell) = length(spikes_this_cell)/(spikes_this_cell(end)-spikes_this_cell(1));

    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, events, windows, bin_size);
     psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    CA1_UP_psth(cell,:) = psth;
end
[~, index ] = sort(mean_FR,'descend');
% CA1_FR_sorted = cell_id(index);

subplot(2,2,2)
imagesc(bins,1:1:length(index),zscore(CA1_UP_psth(index,:),0,2))
ylabel('Cell (sorted by FR)')
xlabel('Time relative to UP state (s)')
cbar = colorbar
cbar.Label.String = 'Firing Rate (z scored)';
hold on
plot([0 0],get(gca,'ylim'),'r','LineWidth',5)
clim([0 4])
title('CA1')


cell_id = unique(CA3_spike_times(:,1));
% events = slow_waves.ints.UP(:,1);
% events = ripples.onset;
CA3_UP_psth= [];
mean_FR = [];

bin_size = 0.001;
for cell = 1:length(cell_id)
    spikes_this_cell = CA3_spike_times(find(CA3_spike_times(:,1) == cell_id(cell)),2);
    mean_FR(cell) = length(spikes_this_cell)/(spikes_this_cell(end)-spikes_this_cell(1));

    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, events, windows, bin_size);
     psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    CA3_UP_psth(cell,:) = psth;
end
[~, index ] = sort(mean_FR,'descend');
% V1_FR_sorted = cell_id(index);

subplot(2,2,3)
imagesc(bins,1:1:length(index),zscore(CA3_UP_psth(index,:),0,2))
ylabel('Cell (sorted by FR)')
xlabel('Time relative to UP state (s)')
cbar = colorbar
cbar.Label.String = 'Firing Rate (z scored)';
hold on
plot([0 0],get(gca,'ylim'),'r','LineWidth',5)
clim([0 4])
title('CA3')



cell_id = unique(DG_spike_times(:,1));
% events = slow_waves.ints.UP(:,1);
% events = ripples.onset;
DG_UP_psth= [];
mean_FR = [];

bin_size = 0.001;
for cell = 1:length(cell_id)
    spikes_this_cell = DG_spike_times(find(DG_spike_times(:,1) == cell_id(cell)),2);
    mean_FR(cell) = length(spikes_this_cell)/(spikes_this_cell(end)-spikes_this_cell(1));

    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, events, windows, bin_size);
     psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    DG_UP_psth(cell,:) = psth;
end
[~, index ] = sort(mean_FR,'descend');
% DG_FR_sorted = cell_id(index);

subplot(2,2,4)
imagesc(bins,1:1:length(index),zscore(DG_UP_psth(index,:),0,2))
ylabel('Cell (sorted by FR)')
xlabel('Time relative to UP state (s)')
cbar = colorbar
cbar.Label.String = 'Firing Rate (z scored)';
hold on
plot([0 0],get(gca,'ylim'),'r','LineWidth',5)
clim([0 4])



%% SWR-modulated PSTH

events = slow_waves.ints.UP(:,1);
events = ripples.onset;
% events = V1_ripples.onset;

cell_id = unique(V1_spike_times(:,1));
% cell_id = unique(V1_SWR_spike_times(:,1));

% events = ripples.onset;
V1_UP_psth= [];
bin_size = 0.001;
for cell = 1:length(cell_id)
    spikes_this_cell = V1_spike_times(find(V1_spike_times(:,1) == cell_id(cell)),2);
    mean_FR(cell) = length(spikes_this_cell)/(spikes_this_cell(end)-spikes_this_cell(1));

    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, events, [-1 1], bin_size);
     psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    V1_UP_psth(cell,:) = psth;
end
[~, index ] = sort(mean_FR,'descend');
V1_FR_sorted = cell_id(index);

figure
subplot(2,2,1)
imagesc(bins,1:1:length(index),zscore(V1_UP_psth(index,:),0,2))
% imagesc(bins,1:1:length(index),zscore(V1_UP_psth(index,:),0,2))
ylabel('Cell (sorted by FR)')
xlabel('Time relative to UP state (s)')
cbar = colorbar
cbar.Label.String = 'Firing Rate (Hz)';
hold on
plot([0 0],get(gca,'ylim'),'r','LineWidth',5)
clim([0 4])
title('V1')


% events = slow_waves.ints.UP(:,1);
cell_id = unique(CA1_spike_times(:,1));
% cell_id = unique(CA1_SWR_spike_times(:,1));
% events = slow_waves.ints.UP(:,1);
% events = ripples.onset;
CA1_UP_psth= [];
mean_FR = [];

bin_size = 0.001;
for cell = 1:length(cell_id)
    spikes_this_cell = CA1_spike_times(find(CA1_spike_times(:,1) == cell_id(cell)),2);
    mean_FR(cell) = length(spikes_this_cell)/(spikes_this_cell(end)-spikes_this_cell(1));

    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, events, [-1 1], bin_size);
     psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    CA1_UP_psth(cell,:) = psth;
end
[~, index ] = sort(mean_FR,'descend');
% CA1_FR_sorted = cell_id(index);

subplot(2,2,2)
imagesc(bins,1:1:length(index),zscore(CA1_UP_psth(index,:),0,2))
ylabel('Cell (sorted by FR)')
xlabel('Time relative to UP state (s)')
cbar = colorbar
cbar.Label.String = 'Firing Rate (z scored)';
hold on
plot([0 0],get(gca,'ylim'),'r','LineWidth',5)
clim([0 4])
title('CA1')


cell_id = unique(CA3_spike_times(:,1));
% cell_id = unique(CA3_SWR_spike_times(:,1));
% events = slow_waves.ints.UP(:,1);
% events = ripples.onset;
CA3_UP_psth= [];
mean_FR = [];

bin_size = 0.001;
for cell = 1:length(cell_id)
    spikes_this_cell = CA3_spike_times(find(CA3_spike_times(:,1) == cell_id(cell)),2);
    mean_FR(cell) = length(spikes_this_cell)/(spikes_this_cell(end)-spikes_this_cell(1));

    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, events, [-1 1], bin_size);
     psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    CA3_UP_psth(cell,:) = psth;
end
[~, index ] = sort(mean_FR,'descend');
% V1_FR_sorted = cell_id(index);

subplot(2,2,3)
imagesc(bins,1:1:length(index),zscore(CA3_UP_psth(index,:),0,2))
ylabel('Cell (sorted by FR)')
xlabel('Time relative to UP state (s)')
cbar = colorbar
cbar.Label.String = 'Firing Rate (z scored)';
hold on
plot([0 0],get(gca,'ylim'),'r','LineWidth',5)
clim([0 4])
title('CA3')



cell_id = unique(DG_spike_times(:,1));
% cell_id = unique(DG_SWR_spike_times(:,1));
% events = slow_waves.ints.UP(:,1);
% events = ripples.onset;
DG_UP_psth= [];
mean_FR = [];

bin_size = 0.001;
for cell = 1:length(cell_id)
    spikes_this_cell = DG_spike_times(find(DG_spike_times(:,1) == cell_id(cell)),2);
    mean_FR(cell) = length(spikes_this_cell)/(spikes_this_cell(end)-spikes_this_cell(1));

    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, events, [-1 1], bin_size);
     psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    DG_UP_psth(cell,:) = psth;
end
[~, index ] = sort(mean_FR,'descend');
% DG_FR_sorted = cell_id(index);

subplot(2,2,4)
imagesc(bins,1:1:length(index),zscore(DG_UP_psth(index,:),0,2))
ylabel('Cell (sorted by FR)')
xlabel('Time relative to UP state (s)')
cbar = colorbar
cbar.Label.String = 'Firing Rate (z scored)';
hold on
plot([0 0],get(gca,'ylim'),'r','LineWidth',5)
clim([0 4])

%% Phase-locking firing

bz_PowerPhaseRatemap


% Create Buz-style spike variables
spikes = [];
spikes.UID = [];
spikes.times = [];
count = 1;
for nchannel = channels_selected(1):channels_selected(end)
    clusters_this_channel = find(peakChannel == nchannel)-1; % Minus one because clutser id is 0 based.
    [~,index,]= intersect(cluster_id,clusters_this_channel);
    good_units_index = index(find(nominal_KSLabel(index)=='good'));
    good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(index)=='good'));

    if ~isempty(good_units_this_channel)
        for unit = 1:length(good_units_this_channel)
            if ~isempty(these_spike_times{good_units_index(unit)})
                spikes.times{count} = these_spike_times{good_units_index(unit)}; % Plus one because clutser id is 0 based.
                spikes.UID = [spikes.UID; good_units_this_channel(unit)];
                count = count + 1;
            end
        end
        
    end
end

[PhaseLockingData] = bz_PhaseModulation('lfp',[raw_LFP(2,:)' tvec'],...
    'intervals',slow_waves.ints.UP,'spikes',spikes,'passband',[0.5 4],'samplingRate',new_SR)


[PhaseLockingData] = bz_PhaseModulation('lfp',[raw_LFP(3,:)' tvec'],...
    'intervals',[ripples.onset ripples.offset],'spikes',spikes,'passband',[125 300],'samplingRate',new_SR)

%% power power correlation



lfp.data = raw_LFP(2:3,:)';

lfp.samplingRate = new_SR;


specparms.frange = [1 250];
specparms.nfreqs = 100;
specparms.ncyc = 5;
specparms.space = 'log';
specparms.samplingRate = new_SR;
specparms.showprogress = false;
specparms.saveMat = false;
specparms.fvector = [];
specparms.specnorm = 'log';
specparms.type = 'FFT'; %'FFT'
specparms.winsize = 1;
specparms.noverlap = 0.5;
% specparms.ints = behavioural_state.SWS;
feq_xcorr = [];
 previous= 1;
for n = 1:length(behavioural_state.SWS)
    tindx = FindInInterval(tvec,[behavioural_state.SWS(n,1) behavioural_state.SWS(n,2)],previous);
    lfp.timestamps = tvec((tindx(1):tindx(2)));
    lfp.data = raw_LFP(2:3,tindx(1):tindx(2))';
    specparms.winsize = 1;
    specparms.noverlap = 0.5;
    %     specparms.winsize = 1;
    [ comod ] = bz_Comodulogram(lfp,specparms,[]);
    feq_xcorr(n,:,:) = comod.corrs;

%     [coherogram,phase,t,f]  = bz_MTCoherogram(raw_LFP(2,tindx(1):tindx(2))',raw_LFP(3,tindx(1):tindx(2))',...
%         'window',1,'frequency',new_SR,'range',[0 300]);
%     power_xcorr(n,:) = mean(coherogram,2);

     previous= tindx(2);
end


% 
figure
imagesc(squeeze(mean(feq_xcorr,1)))
colorbar
clim([0 0.4])
xticks(1:10:length(feq_xcorr))
xticklabels(comod.freqs(1:10:length(feq_xcorr)))
yticks(1:10:length(feq_xcorr))
yticklabels(comod.freqs(1:10:length(feq_xcorr)))
colormap jet
xlabel('V1 Frequency')
ylabel('CA1 Frequency')


 previous= 1;
 ripple_feq_xcorr = [];
for n = 1:length(ripples.onset)
    tindx = FindInInterval(tvec,[ripples.onset(n) ripples.onset(n)+0.1],previous);
    lfp.timestamps = tvec((tindx(1):tindx(2)));
    lfp.data = raw_LFP(2:3,tindx(1):tindx(2))';
    specparms.winsize = 0.02;
    specparms.noverlap = 0.01;
    [ comod ] = bz_Comodulogram(lfp,specparms,[]);
    ripple_feq_xcorr(n,:,:) = comod.corrs;

     previous= tindx(2);
end

 previous= 1;
 ripple_feq_xcorr = [];
 comod.freqs = logspace(log10(specparms.frange(1)),...
     log10(specparms.frange(2)),specparms.nfreqs);


for n = 1:length(ripples.onset)
    tindx = FindInInterval(tvec,[ripples.onset(n) ripples.onset(n)+0.1],previous);
    tvec((tindx(1):tindx(2)));
    raw_LFP(2:3,tindx(1):tindx(2))';
   
    for nfrequency = 1:length(comod.freqs)-1
        parameters = list_of_parameters;
        filter_type  = 'bandpass';
        filter_width = [comod.freqs(nfrequency) comod.freqs(nfrequency+1)];                 % range of frequencies in Hz you want to filter between
        filter_order = round(6*new_SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
        norm_freq_range = filter_width/(new_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
        b_filter = fir1(filter_order, norm_freq_range,filter_type);

        hilbert(filtfilt(b_filter,1,raw_LFP(2,tindx(1):tindx(2))));

        filt_hilb1 = hilbert(LFP.gamma(sorted_config.Channel(n),:)); %calculates the Hilbert transform
        amp1 = abs(filt_hilb1);%calculates the instantaneous amplitude filtered signal
        filt_hilb1 = hilbert(LFP.gamma(sorted_config.Channel(n),:)); %calculates the Hilbert transform
        amp1 = abs(filt_hilb1);%calculates the instantaneous amplitude filtered signal
    end
    specparms.winsize = 0.02;
    specparms.noverlap = 0.01;
    [ comod ] = bz_Comodulogram(lfp,specparms,[]);
    ripple_feq_xcorr(n,:,:) = comod.corrs;

     previous= tindx(2);
end

figure
imagesc(squeeze(mean(ripple_feq_xcorr,1)))
colorbar
xticks(1:10:size(ripple_feq_xcorr,3))
xticklabels(comod.freqs(1:10:size(ripple_feq_xcorr,3)))
yticks(1:10:size(ripple_feq_xcorr,3))
yticklabels(comod.freqs(1:10:size(ripple_feq_xcorr,3)))
xlabel('V1 Frequency')
ylabel('CA1 Frequency')
colormap jet

clim([0 0.3])

% 
% imagesc(comod.freqs,comod.freqs,squeeze(mean(ripple_feq_xcorr,1)))
% colorbar
% % clim([-0.4 0.4])
% xlabel('CA1 Frequency')
% ylabel('V1 Frequency')
%  LogScale('xy',2)

%% CSD peri ripple

[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',ripples.peaktimes',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',[0.3 0.3],'filter',[0.5 4]);


taxis = csd.timestamps;
cmax = max(max(csd.data));
% Quick plotting of PSD
colour_line= {'k','r','m','b','c','g','y'};
freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','350 - 625 Hz'};
selected_powers = [1 2 3 4 5 6];

figure
subplot(1,9,1)

for n = 1:length(cluster_boundary)
    plot([0 1],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([0 1],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for n = selected_powers
        p(n) = plot(power(:,n)./max(power(:,n)),sorted_config.Ks_ycoord',colour_line{n})
    
%     p(n) = plot(power_differnece,1:96,colour_line{n})
    hold on
end
% legend('1-3','4-12','30-60','60-100','125-300')
legend([p(selected_powers)],{freq_legends{selected_powers}});

subplot(1,9,2);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Slow oscillation')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,3);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;


% Spindle
[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',ripples.peaktimes',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',[0.3 0.3],'filter',[9 17]);

taxis = csd.timestamps;
cmax = max(max(csd.data));


subplot(1,9,4);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Spindle')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,5);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;



% gamma
[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',ripples.peaktimes',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',[0.3 0.3],'filter',[30 100]);

taxis = csd.timestamps;
cmax = max(max(csd.data));


subplot(1,9,6);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Gamma')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,7);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;



% Ripple
[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',ripples.peaktimes',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',[0.3 0.3],'filter',[125 300]);

taxis = csd.timestamps;
cmax = max(max(csd.data));


subplot(1,9,8);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Ripple')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,9);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;

%%  UP state
[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',slow_waves.ints.UP(:,1)',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',[0.3 0.3],'filter',[0.5 4]);


taxis = csd.timestamps;
cmax = max(max(csd.data));

% Quick plotting of PSD
colour_line= {'k','r','m','b','c','g','y'};
freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','350 - 625 Hz'};
selected_powers = [1 2 3 4 5 6];

figure
subplot(1,9,1)

for n = 1:length(cluster_boundary)
    plot([0 1],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([0 1],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for n = selected_powers
        p(n) = plot(power(:,n)./max(power(:,n)),sorted_config.Ks_ycoord',colour_line{n})
    
%     p(n) = plot(power_differnece,1:96,colour_line{n})
    hold on
end
% legend('1-3','4-12','30-60','60-100','125-300')
legend([p(selected_powers)],{freq_legends{selected_powers}});

subplot(1,9,2);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Slow oscillation')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,3);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;


% Spindle
[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',slow_waves.ints.UP(:,1)',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',[0.3 0.3],'filter',[9 17]);

taxis = csd.timestamps;
cmax = max(max(csd.data));


subplot(1,9,4);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Spindle')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,5);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;



% gamma
[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',slow_waves.ints.UP(:,1)',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',[0.3 0.3],'filter',[30 100]);

taxis = csd.timestamps;
cmax = max(max(csd.data));


subplot(1,9,6);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Gamma')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,7);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;



% Ripple
[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',slow_waves.ints.UP(:,1)',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',[0.3 0.3],'filter',[125 300]);

taxis = csd.timestamps;
cmax = max(max(csd.data));


subplot(1,9,8);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Ripple')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,9);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;


%%  V1 Ripple
[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',V1_ripples.peaktimes',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',[0.3 0.3],'filter',[0.5 4]);


taxis = csd.timestamps;
cmax = max(max(csd.data));

% Quick plotting of PSD
colour_line= {'k','r','m','b','c','g','y'};
freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','350 - 625 Hz'};
selected_powers = [1 2 3 4 5 6];

figure
subplot(1,9,1)

for n = 1:length(cluster_boundary)
    plot([0 1],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([0 1],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for n = selected_powers
        p(n) = plot(power(:,n)./max(power(:,n)),sorted_config.Ks_ycoord',colour_line{n})
    
%     p(n) = plot(power_differnece,1:96,colour_line{n})
    hold on
end
% legend('1-3','4-12','30-60','60-100','125-300')
legend([p(selected_powers)],{freq_legends{selected_powers}});

subplot(1,9,2);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Slow oscillation')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,3);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;


% Spindle
[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',V1_ripples.peaktimes',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',[0.3 0.3],'filter',[9 17]);

taxis = csd.timestamps;
cmax = max(max(csd.data));


subplot(1,9,4);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Spindle')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,5);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;



% gamma
[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',V1_ripples.peaktimes',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',[0.3 0.3],'filter',[30 100]);

taxis = csd.timestamps;
cmax = max(max(csd.data));


subplot(1,9,6);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Gamma')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,7);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;



% Ripple
[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',V1_ripples.peaktimes',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',[0.3 0.3],'filter',[125 300]);

taxis = csd.timestamps;
cmax = max(max(csd.data));


subplot(1,9,8);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Ripple')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,9);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;

%% Track stimuli
event = MousePos.stimuli_onset(MousePos.stimuli_track == 1);
windows = [0.2 0.8];
% [ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',event',...
%     'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',windows);

[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',event',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',windows,'filter',[0.5 4]);
taxis = csd.timestamps;
cmax = max(max(csd.data));

% Quick plotting of PSD
colour_line= {'k','r','m','b','c','g','y'};
freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','350 - 625 Hz'};
selected_powers = [1 2 3 4 5 6];

figure
subplot(1,9,1)

for n = 1:length(cluster_boundary)
    plot([0 1],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([0 1],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for n = selected_powers
        p(n) = plot(power(:,n)./max(power(:,n)),sorted_config.Ks_ycoord',colour_line{n})
    
%     p(n) = plot(power_differnece,1:96,colour_line{n})
    hold on
end
% legend('1-3','4-12','30-60','60-100','125-300')
legend([p(selected_powers)],{freq_legends{selected_powers}});

subplot(1,9,2);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
caxis([-cmax/2 cmax/2]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Slow oscillation')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,3);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;


% Spindle
[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',event',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',windows,'filter',[9 17]);

taxis = csd.timestamps;
cmax = max(max(csd.data));


subplot(1,9,4);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Spindle')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,5);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;



% gamma
[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',event',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',windows,'filter',[30 100]);

taxis = csd.timestamps;
cmax = max(max(csd.data));


subplot(1,9,6);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Gamma')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,7);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;



% Ripple
[ csd, lfpAvg ]  = bz_eventCSD_masa(raw_LFP',tvec',event',...
    'channels',1:1:size(raw_LFP,1),'samplingRate',new_SR,'twin',windows,'filter',[125 300]);

taxis = csd.timestamps;
cmax = max(max(csd.data));


subplot(1,9,8);

contourf(taxis,sorted_config.Ks_ycoord,[zeros(1,size(csd.data,1)); csd.data'; zeros(1,size(csd.data,1))],40,'LineColor','none');hold on;
colormap jet; caxis([-cmax cmax]);
xlabel('time (s)');ylabel('channel');title('CSD');
plot([0 0],[1 size(csd.data,2)],'--k');hold on;
% set(gca,'YDir','reverse')
% yticks(1:length(sorted_config.Ks_ycoord))
% yticklabels(sorted_config.Ks_ycoord)
plot([0 0],ylim,'--r');hold on;
title('Ripple')

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


subplot(1,9,9);

for n = 1:length(cluster_boundary)
    plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])
hold on
plot([taxis(1) taxis(end)],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)


for ch=1:size(lfpAvg.data,2)

    sh_tmp = 1e0*(500000*lfpAvg.data(:,ch)) + sorted_config.Ks_ycoord(ch);
    plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
    clear sh_tmp
end
% ylim([-1000 offset+1000]);
xlim([taxis(1) taxis(end)]);
xlabel('time (ms)');ylabel('channel');title('LFP');
plot([0 0],ylim,'--r');hold on;
%%

MousePos.pos = MousePos.pos(MousePos.corrected_sglxTime<3735.50);
MousePos.corrected_sglxTime = MousePos.sglxTime(MousePos.corrected_sglxTime<3735.50);

data.linear(1).linear = nan(1,length(MousePos.pos));
data.linear(1).linear(find(MousePos.pos >= 999)) = abs(MousePos.pos(find(MousePos.pos >= 999))-1120);
data.linear(1).linear = data.linear(1).linear(MousePos.corrected_sglxTime>0);

data.linear(2).linear = nan(1,length(MousePos.pos));
data.linear(2).linear(find(MousePos.pos <= 121)) = abs(MousePos.pos(find(MousePos.pos <= 121))-120);
data.linear(2).linear = data.linear(2).linear(MousePos.corrected_sglxTime>0);

data.x = nan(1,length(MousePos.pos));
data.x(~isnan(data.linear(1).linear)) = data.linear(1).linear(~isnan(data.linear(1).linear));
data.x(~isnan(data.linear(2).linear)) = data.linear(2).linear(~isnan(data.linear(2).linear));

data.t = MousePos.corrected_sglxTime(MousePos.corrected_sglxTime>0)';

% data.t = MousePos.sglxTime;
[data.t index]= unique(data.t);

for track = 1:2
    data.linear(track).linear = data.linear(track).linear(index);
    data.linear(track).length = 1.2;
end
data.x = data.x(MousePos.corrected_sglxTime>0);
data.x = data.x(index);

data.v = nan(1,length(data.x));
data.v(2:end) = diff(data.x);
data.v_cm = data.v;

position = data;
save extracted_position position -v7.3
clusters = [];

extract_laps_masa(1)


cd(options.EPHYS_DATAPATH)
load extracted_position
x_bins_width = 12;
% x_bins_width = 24;
clusters.spike_times = V1_spike_times(:,2);
clusters.spike_id = V1_spike_times(:,1);
cell_id = unique(V1_spike_times(:,1));
for cell = 1:length(cell_id)
    clusters.spike_id(find(V1_spike_times(:,1) == cell_id(cell))) = cell;
end

place_fields = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,[])
 save extracted_place_fields_V1 place_fields

 % bayesian decoding
time_on_T2 = position.t(~isnan(position.linear(2).linear));
bayesian_spike_count = create_spike_count_masa(place_fields,clusters,time_on_T2(2),time_on_T2(end),[])
 estimated_position = bayesian_decoding(place_fields,bayesian_spike_count,[])

figure;
subplot(1,2,1)
plot(position.t,estimated_position(2).discrete_position,'k')
hold on
scatter(estimated_position(1).run_time_centered,estimated_position(1).peak_position,5,'r','filled')

subplot(1,2,2)
plot(position.t,estimated_position(2).discrete_position,'k')
hold on
scatter(estimated_position(2).run_time_centered,estimated_position(2).peak_position,'b','filled')
imagesc(estimated_position(2).run_time_centered,[1:12:120],estimated_position(2).run);




place_fields_even = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,'even laps')
place_fields_odd = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,'odd laps')


V1_track_1_cells =  cell_id(find(place_fields.track(1).mean_rate_track > 1 & place_fields.track(2).mean_rate_track < 0.5));
V1_track_2_cells =  cell_id(find(place_fields.track(1).mean_rate_track < 0.5 & place_fields.track(2).mean_rate_track > 1));
% cell_selected = find((place_fields.track(1).mean_rate_track  - place_fields.track(2).mean_rate_track)>5)

save track_cells V1_track_1_cells V1_track_2_cells CA1_track_1_cells CA1_track_2_cells

V1_track_1_cells =  find(place_fields.track(1).mean_rate_track > 1 & place_fields.track(2).mean_rate_track < 0.5);
figure
subplot(2,2,1)

bar(place_fields.track(1).mean_rate_track(V1_track_1_cells))
ylabel('Mean Firing Rate')
ylim([0 20])
title('Mean Firing Rate on Track 1 (when moving)')

subplot(2,2,2)
bar(place_fields.track(2).mean_rate_track(V1_track_1_cells))
ylabel('Mean Firing Rate')
ylim([0 20])
title('Mean Firing Rate on Track 2 (when moving)')

subplot(2,2,3)
bar(place_fields.track(1).mean_rate_track(V1_track_1_cells) - place_fields.track(2).mean_rate_track(V1_track_1_cells))
ylabel('Mean Firing Rate Difference')
ylim([0 20])
title('Firing rate difference (when moving)')
sgtitle('Track 1 ''preferred'' V1 cells')


% CA1
clusters.spike_times = CA1_spike_times(:,2);
clusters.spike_id = CA1_spike_times(:,1);
cell_id = unique(CA1_spike_times(:,1));
for cell = 1:length(cell_id)
    clusters.spike_id(find(CA1_spike_times(:,1) == cell_id(cell))) = cell;
end

 place_fields = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,[])
 save extracted_place_fields_CA1 place_fields
 place_fields_even = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,'even laps')
  place_fields_odd = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,'odd laps')

% track_1_cells =  cell_id(find(place_fields.track(1).mean_rate_track > 3 & place_fields.track(2).mean_rate_track < 2);
CA1_track_1_cells = find((place_fields.track(1).mean_rate_track  - place_fields.track(2).mean_rate_track)>0.5);
CA1_track_2_cells = find((place_fields.track(1).mean_rate_track  - place_fields.track(2).mean_rate_track)<-0.5);
pyramidal_cells = find(place_fields.mean_rate < 2);
CA1_track_1_cells = cell_id(intersect(CA1_track_1_cells,pyramidal_cells));
CA1_track_2_cells = intersect(CA1_track_2_cells,pyramidal_cells);


figure
subplot(2,2,1)
bar(place_fields.track(1).mean_rate_track(CA1_track_1_cells))
ylim([0 20])
title('Mean Firing Rate on Track 1 (when moving)')

subplot(2,2,2)
bar(place_fields.track(2).mean_rate_track(CA1_track_1_cells))
ylim([0 20])
title('Mean Firing Rate on Track 2 (when moving)')

subplot(2,2,3)
bar(place_fields.track(1).mean_rate_track(CA1_track_1_cells) - place_fields.track(2).mean_rate_track(CA1_track_1_cells))
ylim([0 20])
title('Firing rate difference (when moving)')
sgtitle('Track 1 ''preferred'' CA1 cells')




for test = 1:3
        figure;
       
        if test == 1
            place_fields = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,[]);
            sgtitle('Whole')
        elseif test == 2

            place_fields = place_fields_odd;
            sgtitle('Odd')
        elseif test == 3
            place_fields = place_fields_even;
            sgtitle('Even')
        end

        c=1;
    for kk=1:length(place_fields.track)
        for j=1:length(place_fields.track)
            y_vector=[];
            matrix=[];
            normalized_matrix=[];
            for ii=1:length(place_fields.track(j).sorted_good_cells)
                %plot sorted
                %             matrix=[];
                %             normalized_matrix=[];
                matrix(ii,:)=place_fields.track(kk).raw{place_fields.track(j).sorted_good_cells(ii)};
                normalized_matrix(ii,:)=(matrix(ii,:)-min(matrix(ii,:)))/(max(matrix(ii,:))-min(matrix(ii,:)));
                subplot(length(place_fields.track),length(place_fields.track),c)
                plfield_row= normalized_matrix(ii,:)+(1.5*ii-1);
                plot(1:length(plfield_row),plfield_row,'k'); hold on;
                xx = [1:length(plfield_row), fliplr(1:length(plfield_row))];
                inBetween = [(1.5*ii-1)*ones(size(plfield_row)), fliplr(plfield_row)];
                fill(xx, inBetween,[139,0,0]/255);
                y_vector= [y_vector, 1.5*ii-1];
            end
            xlim([0 size(normalized_matrix,2)+2]);
            ylim([0 max(y_vector)+1.2]);
            yt=place_fields.track(j).sorted_good_cells;
            set(gca,'ytick',y_vector);
            set(gca,'yticklabel',yt);
            ylabel('Unit ID');
            xlabel('sorted linearized position (bins)');
            c=c+1;
            title([{['place cells on track ' num2str(j)]} ; {['sorted by track ' num2str(kk)]}]);
        end
    end
end


%Create linearised position matrix
index_track = [];
 sorted_cells= [];

 for type = 1:3
     if type == 1
         place_fields = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,[]);
     elseif type == 2
         place_fields = place_fields_odd;
     elseif type == 3
            place_fields = place_fields_even;
     end

     for track_id = 1:2
         raw_matrix = cat(1,place_fields.track(track_id).raw{:});
         normalised_raw_matrix{type}{track_id} = bsxfun(@rdivide, raw_matrix, (place_fields.track(track_id).raw_peak)');
         [~,index_track{type}(track_id,:)] = max(normalised_raw_matrix{type}{track_id},[],2);
         %     unsorted_cells(track_id,:) =
         [~,sorted_cells{type}(track_id,:)] = sort(index_track{type}(track_id,:));
         %     ordered_matrix = normalised_raw_matrix(new_order,:);
     end
 end

 %plot heat map position
 subplot(2,2,1)
 ordered_matrix =  normalised_raw_matrix{2}{1}(sorted_cells{3}(1,:),:);
 imagesc(ordered_matrix);
 colorbar;
 xlabel('Position');
 ylabel('Cell ID');
 axis xy;
 h = colorbar;
 set(get(h,'label'),'string','Normalised firing rate');

  subplot(2,2,2)
 ordered_matrix =  normalised_raw_matrix{2}{2}(sorted_cells{3}(2,:),:);
 imagesc(ordered_matrix);
 colorbar;
 xlabel('Position');
 ylabel('Cell ID');
 axis xy;
 h = colorbar;
 set(get(h,'label'),'string','Normalised firing rate');

  subplot(2,2,3)
 ordered_matrix =  normalised_raw_matrix{2}{1}(sorted_cells{2}(1,:),:);
 imagesc(ordered_matrix);
 colorbar;
 xlabel('Position');
 ylabel('Cell ID');
 axis xy;
 h = colorbar;
 set(get(h,'label'),'string','Normalised firing rate');

  subplot(2,2,4)
 ordered_matrix =  normalised_raw_matrix{3}{1}(sorted_cells{3}(1,:),:);
 imagesc(ordered_matrix);
 colorbar;
 xlabel('Position');
 ylabel('Cell ID');
 axis xy;
 h = colorbar;
 set(get(h,'label'),'string','Normalised firing rate');

 %% Biasing track selective neruons in CA1 by V1
events = ripples;
events = reactivations;
CA1_track1_SWR_spikes = [];
V1_track1_SWR_spikes = [];
V1_track2_SWR_spikes = [];
active_cells_SWR = [];


 for event = 1:length(events.onset)
     onset = events.onset(event)-0.2;
     offset = events.offset(event);
     spike_index = find(V1_spike_times(:,2)>onset & V1_spike_times(:,2)<offset);
     for cell = 1:length(V1_track_1_cells)
         V1_track1_SWR_spikes(cell,event) = sum(find(V1_spike_times(spike_index,1) == V1_track_1_cells(cell)));
     end
     
     for cell = 1:length(V1_track_2_cells)
         V1_track2_SWR_spikes(cell,event) = sum(find(V1_spike_times(spike_index,1) == V1_track_2_cells(cell)));
     end

     onset = events.onset(event);
     offset = events.offset(event);
     spike_index = find(CA1_spike_times(:,2)>onset & CA1_spike_times(:,2)<offset);
     for cell = 1:length(CA1_track_1_cells)
         CA1_track1_SWR_spikes(cell,event) = sum(find(CA1_spike_times(spike_index,1) == CA1_track_1_cells(cell)));
     end

     for cell = 1:length(CA1_track_2_cells)
         CA1_track2_SWR_spikes(cell,event) = sum(find(CA1_spike_times(spike_index,1) == CA1_track_2_cells(cell)));
     end

%      active_cells_SWR(1,event) =  sum(V1_track1_SWR_spikes(:,event) > 0);
%      active_cells_SWR(2,event) =  sum(V1_track2_SWR_spikes(:,event) > 0) ;
%      active_cells_SWR(3,event) =  sum(CA1_track1_SWR_spikes(:,event) > 0);

          active_cells_SWR(1,event) =  sum(V1_track1_SWR_spikes(:,event) > 0) >=3;
     active_cells_SWR(2,event) =   sum(V1_track2_SWR_spikes(:,event) > 0) >=1;
     active_cells_SWR(3,event) =   sum(CA1_track1_SWR_spikes(:,event) > 0) >=1;
     active_cells_SWR(4,event) =   sum(CA1_track2_SWR_spikes(:,event) >= 0) >=0;
 end

find(reactivations.ripple_peak>3)


figure
% p1 = plot(events.onset(find(reactivations.ripple_peak>3)),cumsum(active_cells_SWR(1,find(reactivations.ripple_peak>3))...
%     /max(cumsum(active_cells_SWR(1,find(reactivations.ripple_peak>3))))),'r')
% subplot(2,2,1)
p1 = plot(events.onset,cumsum(active_cells_SWR(1,:))/max(cumsum(active_cells_SWR(1,:))),'r')
hold on
% p2 = plot(events.onset,cumsum(active_cells_SWR(4,:))/max(cumsum(active_cells_SWR(4,:))),'k')
p2 = plot(events.onset,cumsum(active_cells_SWR(3,:))/max(cumsum(active_cells_SWR(3,:))),'b')
p3 = plot(events.onset,cumsum(active_cells_SWR(4,:))/max(cumsum(active_cells_SWR(4,:))),'k')
plot(position.t,position.linear(1).linear/100,'k')
% plot(position.t,-position.linear(2).linear/100,'k')
% s1 = scatter(MousePos.stimuli_onset(MousePos.stimuli_track==1),0.5*ones(1,sum(MousePos.stimuli_track==1)),'r');
% hold on
% s2 = scatter(MousePos.stimuli_onset(MousePos.stimuli_track==2),-0.5*ones(1,sum(MousePos.stimuli_track==2)),'b');
% legend([p s1 s2],{'Cumulative CA1 population bursting events with activation of track 1 preferred V1 neurons ','Track 1 stimuli','Track 2 stimuli'})
% legend([p1 p2 s1 s2],{'Events with Track 1 V1 neurons co-activation','All ripples','Track 1 stimuli','Track 2 stimuli'})
legend([p1 p2 p3],{'Events with Track 1 V1 neurons co-activation','Events with Track 1 CA1 cells activation','All Ripples'})
xlabel('time (s)')
% ylabel('Cumulative CA1 population bursting events events with activation of track 1 preferred V1 neurons (more than 4 neurons)')
ylabel('Cumulative ripple events')
% scatter(MousePos.sglxTime(locs),400*ones(1,length(locs)))





 [rho,pval] = corr(active_cells_SWR(1,find(active_cells_SWR(1,:)>1 & active_cells_SWR(2,:)==0))'...
     ,active_cells_SWR(3,find(active_cells_SWR(1,:)>1 & active_cells_SWR(2,:)==0))','Type','Spearman')

 [rho,pval] = corr(active_cells_SWR(1,find(active_cells_SWR(2,:)>1 & active_cells_SWR(1,:)==0))'...
     ,active_cells_SWR(3,find(active_cells_SWR(1,:)>1 & active_cells_SWR(2,:)==0))','Type','Spearman')


[rho,pval] = corr(active_cells_SWR(1,:)',active_cells_SWR(3,:)','Type','Spearman')
[rho,pval] = corr(active_cells_SWR(2,:)',active_cells_SWR(3,:)','Type','Spearman')

[rho,pval] = corr(sum(V1_track1_SWR_spikes,1)',sum(CA1_track1_SWR_spikes,1)','Type','Spearman')
[rho,pval] = corr(sum(V1_track2_SWR_spikes,1)',sum(CA1_track1_SWR_spikes,1)','Type','Spearman')


figure
subplot(2,2,1)
 hold on
 arrayfun(@(x) scatter(active_cells_SWR(1,x)',active_cells_SWR(3,x)',86,'k','filled','o'),1:size(active_cells_SWR,2))
subplot(2,2,2)
 hold on
 arrayfun(@(x) scatter(active_cells_SWR(2,x)',active_cells_SWR(3,x)',86,'k','filled','o'),1:size(active_cells_SWR,2))

% figure
% subplot(2,2,1)
%  hold on
%  arrayfun(@(x) scatter(sum(V1_track1_SWR_spikes(:,x))',sum(CA1_track1_SWR_spikes(:,x))',30,'k','filled','o'),1:size(active_cells_SWR,2))
% subplot(2,2,2)
%  hold on
%  arrayfun(@(x) scatter(sum(V1_track2_SWR_spikes(:,x))',sum(CA1_track1_SWR_spikes(:,x))',30,'k','filled','o'),1:size(active_cells_SWR,2))


 mdl2 = fitlm(active_cells_SWR(2,:)',active_cells_SWR(3,:)');
 mdl1 = fitlm(active_cells_SWR(1,:)',active_cells_SWR(3,:)');
 [pval,F_stat,~] = coefTest(mdl1);
  [pval,F_stat,~] = coefTest(mdl2);
R2 = mdl1.Rsquared.Adjusted;
R2 = mdl2.Rsquared.Adjusted;

 x =[min(awake_rate) max(awake_rate)];
 b = mdl.Coefficients.Estimate';
 y_est = polyval(fliplr(b),x);
 plot(x,y_est,':','Color','k','LineWidth',3)
 xlabel('Rate of awake replay rate (log2)')
 ylabel('Number of awake replay number (log2)')

V1_track1_SWR_spikes(V1_track1_SWR_spikes > 0)

CA1_track_1_cells
V1_track_1_cells


%% NREM -> REM -> NREM
behavioural_state.SWS

for event = 1:size(behavioural_state.REM,2)
    behavioural_state.REM(event,:)

end



