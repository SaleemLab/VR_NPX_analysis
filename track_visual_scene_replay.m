% Track visual scene replay
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

% SESSION = '20221130';%'20211206';
% StimulusName = 'SparseNoise';
% StimulusName = 'SparseNoise_fullscreen';
StimulusName = 'Masa2tracks_replay'
% StimulusName = 'Masa2tracks'

% gs = [3]; % Sparse Noise full screen
gs = [0];
nChannelsToBin = 24;    % Bin firing rates across channels for RF maps
channelRange = [180 310]; % Look at these channels 

% Some defaults
% options.importMode = 'KS'; % LF or MUA or KS
options.importMode = 'KS'; % LF or MUA or KS
options.BinWidth = 1/60; % resolution (in s) of output resps (e.g. 1/60)
options.stim_dur = 0.1;
options.AnalysisTimeWindow = [0 0.1];% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])
options.ks_unitType = 'good'; % 'mua', 'good' or ''

% Get correct paths
DATAPATH = fullfile(ROOTPATH,'DATA','SUBJECTS');
EPHYS_DATAPATH = fullfile(DATAPATH,SUBJECT,'ephys',SESSION);
BONSAI_DATAPATH = fullfile(DATAPATH,SUBJECT,'stimuli',SESSION);
options.KS_DATAPATH = fullfile(EPHYS_DATAPATH,'kilosort');

% Step 2: Find csv files associated with desired stimulus
% [BehaviourDataFiles,EyeDataFiles,TrialParamsFiles,PDFiles] = getBonsaiFileNames(StimulusName,BONSAI_DATAPATH);

% import_and_align_Masa_VR_Bonsai(StimulusName,BONSAI_DATAPATH,options);


% Set ephys data path
options.gFileNum = gs(1);
folderName = findGFolder(EPHYS_DATAPATH,options.gFileNum);
T_EPHYS_DATAPATH = fullfile(EPHYS_DATAPATH,folderName,[folderName,'_imec0']);
options.EPHYS_DATAPATH = T_EPHYS_DATAPATH; % full path to ephys recording
options.MAP_FILE = fullfile(options.KS_DATAPATH,[SUBJECT,'_',SESSION,'_g0','_tcat.imec0.ap_kilosortChanMap.mat'])
options.KS_CATGT_FNAME = fullfile(['CatGT_',SUBJECT,'_',SESSION,'.log']);


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
raw_LFP = ReadNPXBin(start_samp, nSamp, imecMeta, LF_FILE, options.EPHYS_DATAPATH);
raw_LFP(nEPhysChan+1:end,:) = []; % Get rid of sync channel
raw_LFP = GainCorrectIM(raw_LFP, 1:nEPhysChan, imecMeta);

% % Downsample the data 
raw_LFP = downsample(raw_LFP',downSampleRate)';
new_SR = imecMeta.imSampRate/downSampleRate;


%% Calculate PSD down the probe (All channels)
% Including mean psd power in selected frequency ranges
% Slow wave 0.5-3 Hz, Theta 4-12 Hz, spindle 9 - 17; Slow Gamma 30 - 60 Hz, 
% High Gamma 60 - 100 Hz and Ripple 125 - 300 Hz

nfft_seconds = 2;
% nfft = 2^(nextpow2(SampRate(imecMeta)*nfft_seconds));
nfft = 2^(nextpow2(new_SR*nfft_seconds));

win  = hanning(nfft);
F  = [0.5 3;4 12;9 17;30 60;60 100;125 300;350 625];

for nchannel = 1:size(chan_config,1)
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


% Quick plotting of raw LFP traces
% tvec = start_sec:1/new_SR:start_sec+duration_sec;
tvec = start_sec:1/new_SR:start_sec+duration_sec;

figure(1)
subplot(1,4,1)
for nchannel = 1:size(sorted_config,1)
    plot(tvec,(raw_LFP(sorted_config.Channel(nchannel),:)*10000 + sorted_config.Ks_ycoord(nchannel)),'k')
    hold on
end
xlabel('Time (s)')
ylabel('Distance (um)')

subplot(1,4,2)
power = [];
for nchannel = 1:size(sorted_config,1)
    power(nchannel,:) = PSD(sorted_config.Channel(nchannel)).mean_power;
end

% Quick plotting of PSD
colour_line= {'k','r','m','b','c','g','y'};
freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','350 - 625 Hz'};
selected_powers = [1 2 3 6];
subplot(1,4,2)
for n = selected_powers
        p(n) = plot(power(:,n)./max(power(:,n)),sorted_config.Ks_ycoord',colour_line{n})
    
%     p(n) = plot(power_differnece,1:96,colour_line{n})
    hold on
end
% legend('1-3','4-12','30-60','60-100','125-300')
legend([p(selected_powers)],{freq_legends{selected_powers}});

% legend([p(1),p(2),p(5)],{freq_legends{1},freq_legends{2},freq_legends{5}});


%% Bonsai behaviour data

% Quick and dirty loading wheel data
[AP_FILE,LF_FILE] = findImecBinFile(options.EPHYS_DATAPATH);
FILE_TO_USE = AP_FILE;
binpath = fullfile(options.EPHYS_DATAPATH,AP_FILE);
syncTimes_ephys = SGLXextractSyncPulseFromBinFile(binpath);
parseDate = date;
[~,fname] = fileparts(EPHYS_DATAPATH);
save(fullfile(EPHYS_DATAPATH,[fname,'_syncpulseTimes.mat']),'syncTimes_ephys','parseDate');
syncTimes_ephys = syncTimes_ephys.on; % use upswings currently

cd(BONSAI_DATAPATH)
% PDAasync = readtable('Masa2tracks_saveMousePos2022-11-30T13_50_24.csv');
% temp = readtable(wheel_file);
stimuli_file = 'Masa2tracks_replay_saveMousePos2022-12-01T14_58_01.csv';
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
[photodiode, photodiode_sync] = import_bonsai_photodiode(fullfile([BONSAI_DATAPATH,'\',photodiode_file]));
[photodiode] = alignBonsaiToEphysSyncTimes(photodiode,syncTimes_ephys);

tick_to_cm_conversion = 3.1415*20/1024; % radius 20 cm and one circle is 1024 ticks;
speed = [0; diff(peripherals.Wheel*tick_to_cm_conversion)];
speed(speed<-100) = 0;
speed(speed>100) = 0;
speedTreshold = 1;

% speed_interp = interp1(peripherals.sglxTime,speed,tvec','linear');

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
for itrial =1:nTrials
    % find frame where Bonsai asked to change quad
    frameidx_start = find(photodiode.FrameNumber <= frame_trial_start(itrial),1,'last');
    % find next index where the photodiode detected a quad change
    temp_idx = find(photodiode.Photodiode(frameidx_start:end)>pd_thresh_up,1,'first');
    idx_start = frameidx_start+temp_idx-1;
    % convert to spike GLX time
    pdstart(itrial) = photodiode.sglxTime(idx_start);
    % to get photodiode down we need to do some median filtering to
    % remove some awkward down swings when the quad is white
    frameidx_end = find(photodiode.FrameNumber >= frame_trial_end(itrial),1,'first');
    temp_idx = find(smoothdata(photodiode.Photodiode(frameidx_end:end),'movmedian', 5)<pd_thresh_down,1,'first');
    idx_end = frameidx_end+temp_idx-1;
    pdend(itrial) = photodiode.sglxTime(idx_end);
end
photodiodeData.stim_on.sglxTime = pdstart';
photodiodeData.stim_off.sglxTime = pdend';
[MousePos] = alignBonsaiToPhotodiode(peripherals,sort([photodiodeData.stim_on.sglxTime; photodiodeData.stim_off.sglxTime]),[]);

% [pks,locs]= findpeaks(abs(MousePos.pos-200));

% figure
% plot(MousePos.sglxTime,MousePos.pos)
% hold on;
% scatter(MousePos.sglxTime(locs),400*ones(1,length(locs)))
% % 
% % MousePos.stimuli_onset = MousePos.sglxTime_corrected(locs);
% % MousePos.stimuli_id = MousePos.pos(locs);
% % MousePos.stimuli_track = MousePos.pos(locs);
% % MousePos.stimuli_track(MousePos.pos(locs)<200) = 2;
% % MousePos.stimuli_track(MousePos.pos(locs)>200) = 1;


%% Extract SUA data


[AP_FILE,LF_FILE] = findImecBinFile(options.EPHYS_DATAPATH);
imecMeta = NPadmin.ReadMeta(AP_FILE,options.EPHYS_DATAPATH);
chan_config = NPadmin.getNPChannelConfig(imecMeta);

KS_metrics = readtable(fullfile(options.KS_DATAPATH,'metrics.csv'));
[these_spike_times,nominal_KSLabel,cluster_id,peakChannel,maxSpkTime] = import_ks_spiketimes(options.KS_DATAPATH,options.gFileNum,options.KS_CATGT_FNAME,imecMeta.imSampRate)
mean_waveforms = readNPY([options.KS_DATAPATH,'/mean_waveforms.npy']); % Information about mean waveform

tic
MUA_filter_length = 41;
MUA_filter_alpha = 4;
time_step=1/new_SR; %match timestep to LFP time
time_bins_edges= tvec(1):time_step:max(tvec);
spike_times = [];
spike_ID = [];
CA1_spike_times = [];
count = 1;
% spikes.times = [];
% spikes.UID = [];

for nchannel = best_ripple_channel-15:best_ripple_channel+15
    clusters_this_channel = cluster_id(find(peakChannel == nchannel));
    good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(clusters_this_channel+1)=='good'));% Plus one because clutser id is 0 based.

    %     if ~isempty(clusters_this_channel) % If any clusters
    if ~isempty(good_units_this_channel) % If any SUA

        for unit = 1:length(good_units_this_channel)
            spike_times = [spike_times; these_spike_times{good_units_this_channel(unit)+1}]; % Plus one because clutser id is 0 based.
            spike_ID = [spike_ID; good_units_this_channel(unit)*ones(length(these_spike_times{good_units_this_channel(unit)+1}),1)];

%             if ~isempty(these_spike_times{good_units_this_channel(unit)+1})
%                 spikes.times{count} = these_spike_times{good_units_this_channel(unit)+1}; % Plus one because clutser id is 0 based.
%                 spikes.UID = [spikes.UID; good_units_this_channel(unit)];
%                 count = count + 1;
%             end
        end

%         for unit = 1:length(good_units_this_channel)
%             if ~isempty(these_spike_times{good_units_this_channel(unit)+1})
%                 spikes.times{count} = these_spike_times{good_units_this_channel(unit)+1}; % Plus one because clutser id is 0 based.
%                 spikes.UID = [spikes.UID; good_units_this_channel(unit)];
%                 count = count + 1;
%             end
%         end

    end
end


%smooth SUA activity with a gaussian kernel
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
CA1_MUA_zscore=zscore(filtfilt(w,1,histcounts(spike_times,time_bins_edges)));
[~,index] = sort(spike_times);
CA1_spike_times(:,1) = spike_ID(index);
CA1_spike_times(:,2) = spike_times(index);
MUA_time = time_bins_edges(1:end-1)+time_step/2;
% CA1_SUA_spikes = spikes;
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
    clusters_this_channel = cluster_id(find(peakChannel == nchannel));
    good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(clusters_this_channel+1)=='good'));% Plus one because clutser id is 0 based.

    %     if ~isempty(clusters_this_channel) % If any clusters
    if ~isempty(good_units_this_channel) % If any SUA

        for unit = 1:length(good_units_this_channel)
            spike_times = [spike_times; these_spike_times{good_units_this_channel(unit)+1}]; % Plus one because clutser id is 0 based.
            spike_ID = [spike_ID; good_units_this_channel(unit)*ones(length(these_spike_times{good_units_this_channel(unit)+1}),1)];
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
    clusters_this_channel = cluster_id(find(peakChannel == nchannel));
    good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(clusters_this_channel+1)=='good'));% Plus one because clutser id is 0 based.

    %     if ~isempty(clusters_this_channel) % If any clusters
    if ~isempty(good_units_this_channel) % If any SUA

        for unit = 1:length(good_units_this_channel)
            spike_times = [spike_times; these_spike_times{good_units_this_channel(unit)+1}]; % Plus one because clutser id is 0 based.
            spike_ID = [spike_ID; good_units_this_channel(unit)*ones(length(these_spike_times{good_units_this_channel(unit)+1}),1)];
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
    chan_config.SpikeGLXchan0(find(chan_config.Ks_ycoord < chan_config.Ks_ycoord(best_SW_channel)...
    & chan_config.Ks_ycoord >chan_config.Ks_ycoord(best_SW_channel)-1000));

for nchannel = channels_selected(1):channels_selected(end)
    clusters_this_channel = cluster_id(find(peakChannel == nchannel));
    good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(clusters_this_channel+1)=='good'));% Plus one because clutser id is 0 based.

    %     if ~isempty(clusters_this_channel) % If any clusters
    if ~isempty(good_units_this_channel) % If any SUA

        for unit = 1:length(good_units_this_channel)
            spike_times = [spike_times; these_spike_times{good_units_this_channel(unit)+1}]; % Plus one because clutser id is 0 based.
            spike_ID = [spike_ID; good_units_this_channel(unit)*ones(length(these_spike_times{good_units_this_channel(unit)+1}),1)];
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



[PhaseLockingData] = bz_PhaseModulation('lfp',[raw_LFP(best_spindle_channel,:)' tvec'],...
    'intervals',slow_waves.ints.UP,'spikes',spikes,'passband',[0.5 4],'samplingRate',new_SR)


[PhaseLockingData] = bz_PhaseModulation('lfp',[raw_LFP(best_spindle_channel,:)' tvec'],...
    'intervals',[ripples.onset ripples.offset],'spikes',spikes,'passband',[125 300],'samplingRate',new_SR)

%% Identify SWR modulated cells in CA1 DG CA3 and V1

bin_size = 0.02;
ripple_events = ripples.onset;
cell_id = unique(CA1_spike_times(:,1));
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


bin_size = 0.02;
ripple_events = ripples.onset;
cell_id = unique(V1_spike_times(:,1));

length(unique(V1_spike_times(:,1)))

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



bin_size = 0.02;
ripple_events = ripples.onset;
cell_id = unique(CA3_spike_times(:,1));
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



bin_size = 0.02;
ripple_events = ripples.onset;
cell_id = unique(DG_spike_times(:,1));
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

figure
subplot(2,2,1)
bin_size = 0.001
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(reactivations.onset, ripples.onset, [-0.5 0.5], bin_size);
    psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    psth_se = filtfilt(w,1,std(binnedArray./bin_size)./sqrt(length(ripples.onset)));

hold on
plot(bins,psth,'b')
hold on
patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');
    ylim([0 600])
hold on
plot([0 0],get(gca,'ylim'),'r-')
title('CA1 bursting events during CA1 ripple events')

subplot(2,2,2)
bin_size = 0.001
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(V1_ripples.onset, ripples.onset, [-0.5 0.5], bin_size);
    psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    psth_se = filtfilt(w,1,std(binnedArray./bin_size)./sqrt(length(ripples.onset)));

plot(rasterX,rasterY)
plot(bins,psth,'b')
hold on
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');
    ylim([0 600])
hold on
plot([0 0],get(gca,'ylim'),'r-')
title('V1 ripple events during ripple events')


[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(reactivations.onset, slow_waves.ints.UP(:,1), [-1 1], bin_size);
    psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    psth_se = filtfilt(w,1,std(binnedArray./bin_size)./sqrt(length(ripples.onset)));
subplot(2,2,3)
hold on
plot(bins,psth,'b')
hold on
patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');
    ylim([0 600])
hold on
plot([0 0],get(gca,'ylim'),'r-')
title('CA1 bursting during Up state')

bin_size = 0.001
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(ripples.onset, slow_waves.ints.UP(:,1), [-1 1], bin_size);
    psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    psth_se = filtfilt(w,1,std(binnedArray./bin_size)./sqrt(length(ripples.onset)));

    
subplot(2,2,4)
hold on
plot(bins,psth,'b')
hold on
patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');
    ylim([0 600])
hold on
plot([0 0],get(gca,'ylim'),'r-')
title('Ripple during Up state')
  [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(reactivations.onset, ripples.onset, [-0.5 0.5], bin_size);
subplot(2,2,4)
 plot(rasterX,rasterY)

%%

MUA_filter_length = 100;
SD_alpha = 5; %2 std width
MUA_filter_alpha = (MUA_filter_length-1)/SD_alpha;
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,4
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width

figure
bin_size = 0.001;
% spikes_this_cell = V1_spike_times(find(V1_spike_times(:,1) == cell_id(find(V1_SWR_modulation==1))),2);
SWR_V1_spikes = [];
SWR_cells = find(V1_SWR_modulation==1);
events = reactivations.onset;
% events = slow_waves.ints.UP(:,1);
cell_id = unique(V1_spike_times(:,1));

for cell = 1:length(SWR_cells)
    spikes_this_cell = V1_spike_times(find(V1_spike_times(:,1) == cell_id(SWR_cells(cell))),2);
%     SWR_V1_spikes = [SWR_V1_spikes; spikes_this_cell];
% ripples.onset
    % num_cell = length(unique(CA1_spike_times(:,1)));
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, events, [-0.5 0.5], bin_size);
    psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    psth_se = filtfilt(w,1,std(binnedArray./bin_size)./sqrt(length(ripples.onset)));
%     psth = filtfilt(w,1,mean(zscore(binnedArray./bin_size,0,2))); % normalize to Hz
%     psth_se = filtfilt(w,1,std(zscore(binnedArray./bin_size,0,2))./sqrt(length(ripples.onset)));

    figure
    subplot(4,2,1)
    plot(rasterX, rasterY,'b')
    hold on
    plot([0 0],get(gca,'ylim'),'r-')
    % plot([0.1 0.1],get(gca,'ylim'),'b-')
    % plot([0.5 0.5],get(gca,'ylim'),'r-')
    ylabel('events');
    xlabel('Spike time relative to event onset')
    title('V1')

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
    % ylabel('Spike Rate (spk/s)');
    xlabel('Spike time relative to event onset')
    sgtitle(sprintf('V1 cell %i',cell_id(SWR_cells(cell))))
end

SWR_V1_spikes = sort(SWR_V1_spikes)

% num_cell = length(unique(CA1_spike_times(:,1)));
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(SWR_V1_spikes, ripples.onset, [-0.5 0.5], bin_size);
psth = filtfilt(w,1,mean(zscore(binnedArray./bin_size,0,2))); % normalize to Hz
psth_se = filtfilt(w,1,std(zscore(binnedArray./bin_size,0,2))./sqrt(length(ripples.onset)));

subplot(4,2,1)
plot(rasterX, rasterY,'b')
hold on
plot([0 0],get(gca,'ylim'),'r-')
% plot([0.1 0.1],get(gca,'ylim'),'b-')
% plot([0.5 0.5],get(gca,'ylim'),'r-')
ylabel('events');
xlabel('Spike time relative to event onset')
title('V1')

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
% ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')


%% All SUA PSTH 

events = slow_waves.ints.UP(:,1);
events = ripples.onset;
cell_id = unique(V1_spike_times(:,1));

% events = ripples.onset;
V1_UP_psth= [];
bin_size = 0.001;
for cell = 1:length(cell_id)
    spikes_this_cell = V1_spike_times(find(V1_spike_times(:,1) == cell_id(cell)),2);
    mean_FR(cell) = length(spikes_this_cell)/(spikes_this_cell(end)-spikes_this_cell(1));

    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, events, [-0.5 0.5], bin_size);
     psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
    V1_UP_psth(cell,:) = psth;
end
[~, index ] = sort(mean_FR,'descend');
V1_FR_sorted = cell_id(index);

figure
subplot(2,2,1)
imagesc(bins,1:1:length(index),zscore(V1_UP_psth(index,:),0,2))
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
% events = slow_waves.ints.UP(:,1);
% events = ripples.onset;
CA1_UP_psth= [];
mean_FR = [];

bin_size = 0.001;
for cell = 1:length(cell_id)
    spikes_this_cell = CA1_spike_times(find(CA1_spike_times(:,1) == cell_id(cell)),2);
    mean_FR(cell) = length(spikes_this_cell)/(spikes_this_cell(end)-spikes_this_cell(1));

    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, events, [-0.5 0.5], bin_size);
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

%% PSTH
spike_data = V1_spike_times(:,2);
spike_data = these_spike_times{360};
spike_data = these_spike_times{59};
% for nchannel = [best_spindle_channel:1:best_SW_channel]
%     spike_data = [spike_data; SUA(nchannel).spike_times];
% end


figure
spike_data = CA1_spike_times(:,2);
num_cell = length(unique(CA1_spike_times(:,1)));
% event = slow_waves.ints.UP(:,1);
% event = slow_waves.ints.DOWN(:,1);
event = ripples.onset;
% event = slow_waves.ints.UP(:,1);
% event = spindles.onset;

% event = reactivations.onset(reactivations.zscore>5);
% num_cell = 1;
% for types = 1:4
MUA_filter_length = 50;
SD_alpha = 5; %2 std width
MUA_filter_alpha = (MUA_filter_length-1)/SD_alpha;
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,4
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width

event = MousePos.stimuli_onset(MousePos.stimuli_track == 1)

figure
bin_size = 0.001;
spike_data = CA1_spike_times(:,2);
num_cell = length(unique(CA1_spike_times(:,1)));
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data, event, [-0.3 1], bin_size);
psth = filtfilt(w,1,mean(binnedArray./num_cell./bin_size)); % normalize to Hz
psth_se = filtfilt(w,1,std(binnedArray./num_cell./bin_size)./sqrt(length(event)));

subplot(4,2,1)
plot(rasterX, rasterY,'b')
hold on
plot([0 0],get(gca,'ylim'),'r-')
plot([0.1 0.1],get(gca,'ylim'),'b-')
plot([0.5 0.5],get(gca,'ylim'),'r-')
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
plot([0.5 0.5],get(gca,'ylim'),'r-')
ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')


spike_data = CA3_spike_times(:,2);
num_cell = length(unique(CA3_spike_times(:,1)));
% num_cell = 1;
bin_size = 0.001;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data, event, [-0.3 1], bin_size);
psth = filtfilt(w,1,mean(binnedArray./num_cell./bin_size)); % normalize to Hz
psth_se = filtfilt(w,1,std(binnedArray./num_cell./bin_size)./sqrt(length(event)));

subplot(4,2,3)
plot(rasterX, rasterY,'b')
hold on
plot([0 0],get(gca,'ylim'),'r-')
plot([0.1 0.1],get(gca,'ylim'),'b-')
plot([0.5 0.5],get(gca,'ylim'),'r-')
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
plot([0.5 0.5],get(gca,'ylim'),'r-')
ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')


spike_data = DG_spike_times(:,2);
num_cell = length(unique(DG_spike_times(:,1)));
% num_cell = 1;
bin_size = 0.001;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data, event, [-0.3 1], bin_size);
psth = filtfilt(w,1,mean(binnedArray./num_cell./bin_size)); % normalize to Hz
psth_se = filtfilt(w,1,std(binnedArray./num_cell./bin_size)./sqrt(length(event)));

subplot(4,2,5)
plot(rasterX, rasterY,'k')
hold on
plot([0 0],get(gca,'ylim'),'r-')
plot([0.1 0.1],get(gca,'ylim'),'b-')
plot([0.5 0.5],get(gca,'ylim'),'r-')
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
ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')

sgtitle('CA1 bursting')
sgtitle('Track 2 Stimuli (Right)')


spike_data = V1_spike_times(:,2);
num_cell = length(unique(V1_spike_times(:,1)));
% num_cell = 1;
bin_size = 0.001;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data, event, [-0.3 1], bin_size);
psth = filtfilt(w,1,mean(binnedArray./num_cell./bin_size)); % normalize to Hz
psth_se = filtfilt(w,1,std(binnedArray./num_cell./bin_size)./sqrt(length(event)));

subplot(4,2,7)
plot(rasterX, rasterY,'k')
hold on
plot([0 0],get(gca,'ylim'),'r-')
plot([0.1 0.1],get(gca,'ylim'),'b-')
plot([0.5 0.5],get(gca,'ylim'),'r-')
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
plot([0.5 0.5],get(gca,'ylim'),'r-')


% xlim([-0.2 0.6])
ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')
sgtitle('Ripple Onset')

sgtitle('CA1 bursting')
sgtitle('Track 2 Stimuli (Right)')
sgtitle('Track 1 Stimuli (Left)')