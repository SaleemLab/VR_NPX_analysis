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
gs = [2];
nChannelsToBin = 24;    % Bin firing rates across channels for RF maps
channelRange = [180 310]; % Look at these channels 

% Some defaults
% options.importMode = 'KS'; % LF or MUA or KS
options.importMode = 'LF'; % LF or MUA or KS
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

%% MUA

[AP_FILE,LF_FILE] = findImecBinFile(options.EPHYS_DATAPATH);
imecMeta = NPadmin.ReadMeta(AP_FILE,options.EPHYS_DATAPATH);
chan_config = NPadmin.getNPChannelConfig(imecMeta);

KS_metrics = readtable(fullfile(options.KS_DATAPATH,'metrics.csv'));
[these_spike_times,nominal_KSLabel,cluster_id,peakChannel,maxSpkTime] = import_ks_spiketimes(options.KS_DATAPATH,options.gFileNum,options.KS_CATGT_FNAME,imecMeta.imSampRate)
mean_waveforms = readNPY([options.KS_DATAPATH,'/mean_waveforms.npy']); % Information about mean waveform
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
    clusters_this_channel = cluster_id(find(peakChannel == nchannel));
    good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(clusters_this_channel+1)=='good'));% Plus one because clutser id is 0 based.
    
    MUA(nchannel).zscore = zeros(1,length(time_bins_edges));
    MUA(nchannel).time_bins_edges= time_bins_edges;
    MUA(nchannel).time_bins =(MUA(nchannel).time_bins_edges(1:end-1))+time_step/2;


    SUA(nchannel).zscore = zeros(1,length(time_bins_edges));
    SUA(nchannel).time_bins_edges= time_bins_edges;
    SUA(nchannel).time_bins =(SUA(nchannel).time_bins_edges(1:end-1))+time_step/2;
    
    MUA(nchannel).cluster_ID = clusters_this_channel;
    if ~isempty(clusters_this_channel) % If any clusters
        MUA(nchannel).spike_times = [];
        for unit = 1:length(clusters_this_channel)
            MUA(nchannel).spike_times = [MUA(nchannel).spike_times; these_spike_times{clusters_this_channel(unit)+1}]; % Plus one because clutser id is 0 based.
        end
        %smooth MUA activity with a gaussian kernel

        %         time_step=0.001; %1 ms timestep

        w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
        w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
        %         MUA(nchannel).zscored=zscore(histcounts(MUA(nchannel).spike_times,MUA(nchannel).time_bins_edges));
        MUA(nchannel).zscore=zscore(filtfilt(w,1,histcounts(MUA(nchannel).spike_times,MUA(nchannel).time_bins_edges)));


        if ~isempty(good_units_this_channel) % If any SUA
             SUA(nchannel).spike_times= [];
            for unit = 1:length(good_units_this_channel)

                SUA(nchannel).spike_times = [SUA(nchannel).spike_times; these_spike_times{good_units_this_channel(unit)+1}]; % Plus one because clutser id is 0 based.
            end
            SUA(nchannel).unit_ID = good_units_this_channel;
            %smooth SUA activity with a gaussian kernel
            w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
            w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
            %         MUA(nchannel).zscored=zscore(histcounts(MUA(nchannel).spike_times,MUA(nchannel).time_bins_edges));
            SUA(nchannel).zscore=zscore(filtfilt(w,1,histcounts(SUA(nchannel).spike_times,SUA(nchannel).time_bins_edges)));
        end
    end
end
toc

% Sort channel for spike time data
col_ID = cols_available(1); % (1) is 11
sorted_config_spikes = sortrows(chan_config,'Ks_ycoord','descend');
col_idx = sorted_config_spikes.Ks_xcoord == col_ID;
sorted_config_spikes = sorted_config_spikes(col_idx,:);

%% LFP gamma coherence -> Gradient descent

% Gamma band
parameters = list_of_parameters;
filter_type  = 'bandpass';
filter_width = [30 100];                 % range of frequencies in Hz you want to filter between
filter_order = round(6*new_SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
norm_freq_range = filter_width/(new_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_gamma = fir1(filter_order, norm_freq_range,filter_type);

for nchannel = 1:size(chan_config,1)
    LFP.gamma(nchannel,:) = filtfilt(b_gamma,1,raw_LFP(nchannel,:));
%     LFP(nchannel).gamma_zscore = zscore(abs(hilbert(LFP(nchannel).gamma)));
end


for n = 1:size(sorted_config,1)
    tic
    for m = 1:size(sorted_config,1)
%         [coherence,f] = mscohere(raw_LFP(sorted_config.Channel(n),:),raw_LFP(sorted_config.Channel(m),:),[],[],[],new_SR);

         filt_hilb1 = hilbert(LFP.gamma(sorted_config.Channel(n),:)); %calculates the Hilbert transform of eeg1
         amp1 = abs(filt_hilb1);%calculates the instantaneous amplitude of eeg1 filtered between low_freq and high_freq
         amp1=amp1-mean(amp1); %removes mean of the signal because the DC component of a signal does not change the correlation
         filt_hilb2 = hilbert(LFP.gamma(sorted_config.Channel(m),:));%calculates the Hilbert transform of eeg2
         amp2 = abs(filt_hilb2);%calculates the instantaneous amplitude of eeg2 filtered between low_freq and high_freq
         amp2=amp2-mean(amp2);
         [crosscorr,lags]=xcorr(amp1, amp2,round(new_SR/10),'coeff'); %calculates crosscorrelations between amplitude vectors
         lags=(lags./new_SR)*1000; %converts lags to miliseconds
         g=find(crosscorr==max(crosscorr));%identifies index where the crosscorrelation peaks
         max_crosscorr_lag=lags(g);%identifies the lag at which the crosscorrelation peaks
         gamma_phase_coherence(n,m) = unwrap(circ_mean(circ_dist(angle(filt_hilb1),angle(filt_hilb2)),[],2));
         gamma_coherence(n,m) = max(crosscorr);
%          gamma_coherence_ms(n,m) = mean(coherence(find(f <= 100 & f>= 30)));
%          [coherogram,phase,t,f]  = bz_MTCoherogram(raw_LFP(sorted_config.Channel(n),:)',raw_LFP(sorted_config.Channel(m),:)','window',5,'frequency',new_SR,'range',[30 100]);
%          gamma_coherence_bz(n,m) = mean(mean(coherogram));
%          gamma_phase_coherence_bz(n,m) = circ_mean(reshape(phase,size(phase,1)*size(phase,2),1));
% %          gamma_coherence2(n,m) = mean(coherence( find(f>=30 & f<=100)));
%          figure('color',[1 1 1])
%          plot(lags, crosscorr,'color',[0 0 1],'linewidth',2),hold on %plots crosscorrelations
%          plot(lags(g),crosscorr(g),'rp','markerfacecolor',[1 0 0],'markersize',10)%plots marker at the peak of the cross correlation
%          plot([0 0],[1.05*max(crosscorr) 0.95*min(crosscorr)],'color',[0 0 0],'linestyle',':', 'linewidth',2) %plots dashed line at zero lag
% %          set(gca,'xtick',[-100 -50 0 50 100])
% %          axis tight, box off, xlim([-101 100])
%          xlabel('Lag (ms)','fontsize',14)
%          ylabel('Crosscorrelation','fontsize',14)

    end
    toc
end


save gamma_coherence gamma_coherence gamma_phase_coherence

save gamma_coherence_bz gamma_coherence_bz gamma_phase_coherence_bz
% [coherogram,phase,~,~,~,t,f] = cohgramc(LFP(n).gamma_zscore',LFP(m).gamma_zscore',[window window-window/2],parameters);
% 
%           [coherogram,phase,t,f]  = bz_MTCoherogram(LFP(n).gamma_zscore',LFP(m).gamma_zscore','frequency',new_SR,'range',[30 100]);
% [coherogram,phase,t,f]  = bz_MTCoherogram(raw_LFP(sorted_config.Channel(n),:)',raw_LFP(sorted_config.Channel(m),:)','frequency',new_SR,'range',[30 100]);
hold on
cutoffs = [0 1];
figure;hold on;
subplot(2,1,1);
PlotColorMap(coherogram,'x',t,'y',f,'cutoffs',cutoffs,'newfig','off');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Coherogram Amplitude');
colorbar
% 
subplot(2,1,2);
PlotColorMap(phase,'x',t,'y',f,'cutoffs',[-pi pi],'newfig','off');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Coherogram Phase');
colorbar
% 
% parameters.Fs = new_SR;
% parameters.tapers = [3 5];
% parameters.pad = 0;

% -11 here because it is in dura gel
% brainmap{1} = 1:1:96;
% brainmap{2} = [96 96];
% [ cluass,cluNRG ] = bz_GradDescCluster(gamma_coherence,'brainmap',brainmap)

% Find where does the probe go into the brain (only brain signal for gamma coherence clustering analysis)
% Find first peak that suddenly increases the power
power_differnece = [0; diff(power(:,6)./max(power(:,6)))]; % Use ripple or theta
[~,first_in_brain_channel] = findpeaks(power_differnece,'MinPeakHeight',0.1);
first_in_brain_channel = first_in_brain_channel(1);
% power_differnece = [0; diff(power(:,2)./max(power(:,2)))]; % Use ripple and theta
% [~,theta_loc] = findpeaks(power_differnece,'MinPeakHeight',0.1);

probe_length_in_brain = sorted_config.Ks_ycoord(first_in_brain_channel(1));

brainmap{1} = 1:1:96-first_in_brain_channel(1)+1;
brainmap{2} = [96-first_in_brain_channel(1)+1 96-first_in_brain_channel(1)+1];
[ cluass,cluNRG ] = bz_GradDescCluster(gamma_coherence(first_in_brain_channel(1):end,first_in_brain_channel(1):end),'brainmap',brainmap)

% Plot cross-coherence matrix
figure
subplot(1,2,1)
imagesc(gamma_coherence)
colorbar
xticks(1:2:96)
xticklabels(sorted_config.Ks_ycoord(1:2:96))
yticks(1:2:96)
yticklabels(sorted_config.Ks_ycoord(1:2:96))
hold on
cluster_boundary = find(diff(cluass)~=0);
plot([0 96],[first_in_brain_channel-1 first_in_brain_channel-1],'--k','LineWidth',2)
plot([first_in_brain_channel-1 first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
hold on

for n = 1:length(cluster_boundary)
    plot([0 96],[cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],'--k','LineWidth',2)
    plot([cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
%         plot([0 96],[cluster_boundary(n)+11 cluster_boundary(n)+11],'--k','LineWidth',2)
%     plot([cluster_boundary(n)+11 cluster_boundary(n)+11],[0 96],'--k','LineWidth',2)
    hold on
end
title('LFP Gamma Coherence (30 - 100Hz)')

subplot(1,2,2)
imagesc(gamma_phase_coherence)
colorbar

xticks(1:2:96)
xticklabels(sorted_config.Ks_ycoord(1:2:96))
yticks(1:2:96)
yticklabels(sorted_config.Ks_ycoord(1:2:96))
hold on
cluster_boundary = find(diff(cluass)~=0);
plot([0 96],[11 11],'--k','LineWidth',2)
plot([11 11],[0 96],'--k','LineWidth',2)
hold on

for n = 1:length(cluster_boundary)
    plot([0 96],[cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],'--k','LineWidth',2)
    plot([cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
    hold on
end
colorbar
title('LFP Gamma Phase Coherence (30 - 100Hz)')


subplot(1,4,3)
for nchannel = 1:size(sorted_config,1)
    plot(tvec(sample_to_view),(MUA(sorted_config_spikes.SpikeGLXchan0(nchannel)).zscore(sample_to_view)*5 + sorted_config.Ks_ycoord(nchannel)),'k')
    hold on
end

ylim([0 4000])
xlabel('Time (s)')
ylabel('Distance (um)')
title('Smoothed MUA zsocre activity')

subplot(1,4,4)
cluster_boundary = find(diff(cluass)~=0);
plot([0 1],[sorted_config.Ks_ycoord(11) sorted_config.Ks_ycoord(11)],'--k','LineWidth',2)
hold on

for n = 1:length(cluster_boundary)
    plot([0 1],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])

% Quick plotting of PSD
colour_line= {'k','r','m','b','c','g','y'};
freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','350 - 625 Hz'};
selected_powers = [1 2 3 6];

for n = selected_powers
        p(n) = plot(power(:,n)./max(power(:,n)),sorted_config.Ks_ycoord',colour_line{n})
    
%     p(n) = plot(power_differnece,1:96,colour_line{n})
    hold on
end
% legend('1-3','4-12','30-60','60-100','125-300')
legend([p(selected_powers)],{freq_legends{selected_powers}});
ylim([0 4000])
title('Gradient descent clusters (layer boundary)')

% chan_config.SpikeGLXchan0(find(chan_config.Ks_ycoord < 3420 & chan_config.Ks_ycoord >2340))
%% CSD gamma band
 lfp = [];
for nchannel = 1:size(sorted_config,1)
    lfp.data(:,nchannel) = LFP(sorted_config.Channel(nchannel)).gamma;% switch to data X nchannel
end
lfp.timestamps = tvec;
lfp.samplingRate = new_SR;
[ csd ] = bz_CSD (lfp);

% [ csd ] = bz_eventCSD (lfp);
gamma_coherence_CSD = zeros(size(sorted_config,1),size(sorted_config,1));
gamma_phase_coherence_CSD = zeros(size(sorted_config,1),size(sorted_config,1));
for n = 1:size(sorted_config,1)-2
    tic
    for m = 1:size(sorted_config,1)-2
%          filt_hilb1 = hilbert(csd.data(:,n))'; %calculates the Hilbert transform of eeg1
%          filt_hilb2 = hilbert(csd.data(:,m))';%calculates the Hilbert transform of eeg2

         [coherogram,phase,t,f]  = bz_MTCoherogram(csd.data(:,n),csd.data(:,m),'window',5,'frequency',new_SR,'range',[30 100]);
         gamma_coherence_CSD(n+1,m+1) = mean(mean(coherogram));
         gamma_phase_coherence_CSD(n+1,m+1) = circ_mean(reshape(phase,size(phase,1)*size(phase,2),1));
%          gamma_phase_coherence_CSD(n+1,m+1) = unwrap(circ_mean(circ_dist(angle(filt_hilb1),angle(filt_hilb2)),[],2));
        
    end
    toc
end

save gamma_coherence_CSD gamma_coherence_CSD gamma_phase_coherence_CSD

figure
subplot(1,2,1)
imagesc(gamma_coherence_CSD)
colorbar
xticks(1:2:96)
xticklabels(sorted_config.Ks_ycoord(1:2:96))
yticks(1:2:96)
yticklabels(sorted_config.Ks_ycoord(1:2:96))
hold on
cluster_boundary = find(diff(cluass)~=0);
plot([0 96],[first_in_brain_channel-1 first_in_brain_channel-1],'--k','LineWidth',2)
plot([first_in_brain_channel-1 first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
hold on

for n = 1:length(cluster_boundary)
    plot([0 96],[cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],'--k','LineWidth',2)
    plot([cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
%         plot([0 96],[cluster_boundary(n)+11 cluster_boundary(n)+11],'--k','LineWidth',2)
%     plot([cluster_boundary(n)+11 cluster_boundary(n)+11],[0 96],'--k','LineWidth',2)
    hold on
end
title('CSD Gamma Coherence (30 - 100Hz)')
subplot(1,2,2)
imagesc(gamma_phase_coherence_CSD)
colorbar
xticks(1:2:96)
xticklabels(sorted_config.Ks_ycoord(1:2:96))
yticks(1:2:96)
yticklabels(sorted_config.Ks_ycoord(1:2:96))
hold on
cluster_boundary = find(diff(cluass)~=0);
plot([0 96],[first_in_brain_channel-1 first_in_brain_channel-1],'--k','LineWidth',2)
plot([first_in_brain_channel-1 first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
hold on

for n = 1:length(cluster_boundary)
    plot([0 96],[cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],'--k','LineWidth',2)
    plot([cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
%         plot([0 96],[cluster_boundary(n)+11 cluster_boundary(n)+11],'--k','LineWidth',2)
%     plot([cluster_boundary(n)+11 cluster_boundary(n)+11],[0 96],'--k','LineWidth',2)
    hold on
end
title('CSD Phase Gamma Coherence (30 - 100Hz)')



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

speed_interp = interp1(peripherals.sglxTime,speed,tvec','linear');

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
[MousePos] = alignBonsaiToPhotodiode(MousePos,photodiodeData.stim_on.sglxTime);

[pks,locs]= findpeaks(abs(MousePos.pos-200));

figure
plot(MousePos.sglxTime,MousePos.pos)
hold on;
scatter(MousePos.sglxTime(locs),400*ones(1,length(locs)))

MousePos.stimuli_onset = MousePos.sglxTime_corrected(locs);
MousePos.stimuli_id = MousePos.pos(locs);
MousePos.stimuli_track = MousePos.pos(locs);
MousePos.stimuli_track(MousePos.pos(locs)<200) = 2;
MousePos.stimuli_track(MousePos.pos(locs)>200) = 1;


%% Sleep detection

normalised_ripple = power(:,6)/max(power(:,6)); %power normalized by max power across channel
normalised_theta = power(:,2)/max(power(:,2)); %power normalized by max power across channel
normalised_slow_wave = power(:,1)./max(power(:,1));
normalised_spindle = power(:,3)./max(power(:,3));

[value, candidate_index] = findpeaks(normalised_slow_wave,'MinPeakHeight',0.3); % Find channels with normalised Slow wave power bigger than 0.3
[~,index]  = max(normalised_slow_wave(candidate_index)-normalised_ripple(candidate_index));  %find channel with maximum power difference between gamma and slow wave
best_SW_channel_this_column = candidate_index(1); % First peak (superficial layer of slow wave)
best_SW_channel = sorted_config.Channel(best_SW_channel_this_column+1);

[value, candidate_index] = findpeaks(normalised_spindle,'MinPeakHeight',0.3); % Find channels with normalised Slow wave power bigger than 0.3
[~,index]  = max(normalised_spindle(candidate_index)-normalised_ripple(candidate_index));  %find channel with maximum power difference between theta and slow wave
best_spindle_channel_this_column = candidate_index(1); % First peak 
best_spindle_channel = sorted_config.Channel(best_spindle_channel_this_column+1);


% spindle
parameters = list_of_parameters;
filter_type  = 'bandpass';
filter_width = [30 100];                 % range of frequencies in Hz you want to filter between
filter_order = round(6*new_SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
norm_freq_range = filter_width/(new_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_gamma = fir1(filter_order, norm_freq_range,filter_type);

for nchannel = 1:size(chan_config,1)
    LFP.gamma(nchannel,:) = filtfilt(b_gamma,1,raw_LFP(nchannel,:));
    LFP.gamma_zscore(nchannel,:)  = zscore(abs(hilbert(LFP.gamma(nchannel,:))));
end


[freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
    [tvec' raw_LFP(best_SW_channel,:)'],[tvec' raw_LFP(best_ripple_channel,:)'],...
    [tvec' speed_interp],speedTreshold);


behavioural_state.freezing = freezing;
behavioural_state.quietWake = quietWake;
behavioural_state.SWS = SWS;
behavioural_state.REM = REM;
behavioural_state.movement = movement;

save behavioural_state behavioural_state

best_SW_channel;
best_spindle_channel;
% High spindle power as quietwake
% Low spindle High theta as freezing

figure
subplot(1,2,1)
imagesc(gamma_phase_coherence_CSD)
colorbar
xticks(1:2:96)
xticklabels(sorted_config.Ks_ycoord(1:2:96))
yticks(1:2:96)
yticklabels(sorted_config.Ks_ycoord(1:2:96))
hold on
cluster_boundary = find(diff(cluass)~=0);
plot([0 96],[first_in_brain_channel-1 first_in_brain_channel-1],'--k','LineWidth',2)
plot([first_in_brain_channel-1 first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
hold on
for n = 1:length(cluster_boundary)
    plot([0 96],[cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],'--k','LineWidth',2)
    plot([cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
%         plot([0 96],[cluster_boundary(n)+11 cluster_boundary(n)+11],'--k','LineWidth',2)
%     plot([cluster_boundary(n)+11 cluster_boundary(n)+11],[0 96],'--k','LineWidth',2)
    hold on
end
title('CSD Phase Gamma Coherence (30 - 100Hz)')
% xlim([find(sorted_config.Channel == best_ripple_channel) 96])
% ylim([find(sorted_config.Channel == best_ripple_channel) 96])

subplot(1,2,2)
imagesc(gamma_coherence)
colorbar
xticks(1:2:96)
xticklabels(sorted_config.Ks_ycoord(1:2:96))
yticks(1:2:96)
yticklabels(sorted_config.Ks_ycoord(1:2:96))
hold on
cluster_boundary = find(diff(cluass)~=0);
plot([0 96],[first_in_brain_channel-1 first_in_brain_channel-1],'--k','LineWidth',2)
plot([first_in_brain_channel-1 first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
hold on
for n = 1:length(cluster_boundary)
    plot([0 96],[cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],'--k','LineWidth',2)
    plot([cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
%         plot([0 96],[cluster_boundary(n)+11 cluster_boundary(n)+11],'--k','LineWidth',2)
%     plot([cluster_boundary(n)+11 cluster_boundary(n)+11],[0 96],'--k','LineWidth',2)
    hold on
end
title('LFP Gamma Coherence (30 - 100Hz)')
% xlim([find(sorted_config.Channel == best_ripple_channel) 96])
% ylim([find(sorted_config.Channel == best_ripple_channel) 96])


%% Detect SWR events
% Determine channel with high ripple low theta
% F  = [0.5 3;4 12;9 17;30 60;60 100;125 300;350 625];

for n = 1:length(sorted_config.Ks_ycoord)
    spike_counts(n) = length(MUA(n).spike_times);
end
figure
plot(spike_counts,sorted_config.Ks_ycoord)

%calculate channel with biggest difference between theta and ripple normalized power.
normalised_ripple = power(:,6)/max(power(:,6)); %power normalized by max power across channel
normalised_theta = power(:,2)/max(power(:,2)); %power normalized by max power across channel
% normalised_low_gamma = power(:,3)/max(power(:,3)); %power normalized by max power across channel
% normalised_high_gamma = power(:,4)/max(power(:,4)); %power normalized by max power across channel
% normalised_ripple(find(normalised_ripple<0))=0; %if normalized ripple power is less than zero, treat it as 0
% normalised_gamma_power = gamma_power'/max(gamma_power);

[value, candidate_index] = findpeaks(normalised_ripple,'MinPeakHeight',0.5); % Find channels with normalised ripple power bigger than 0.5
[~,index]  = max(normalised_ripple(candidate_index)-normalised_theta(candidate_index));  %find channel with maximum power difference between theta and ripple
best_ripple_channel_this_column = candidate_index(index);
best_ripple_channel = sorted_config.Channel(best_ripple_channel_this_column);

figure
plot(normalised_ripple,sorted_config.Ks_ycoord)
hold on
plot(normalised_slow_wave,sorted_config.Ks_ycoord)

plot(normalised_theta,sorted_config.Ks_ycoord)
plot(normalised_gamma_power,sorted_config.Ks_ycoord)
plot([0.5 0.5],[0 4000],'--')
legend('Normalised ripple','Normalised theta','Normalised gamma')
xlabel('Normalised power')
ylabel('Depth (um)')


tic
MUA_filter_length = 41;
MUA_filter_alpha = 4;
time_step=1/new_SR; %match timestep to LFP time
time_bins_edges= tvec(1):time_step:max(tvec);
spike_times = [];
spike_ID = [];
CA1_spike_times = [];

for nchannel = best_ripple_channel-10:best_ripple_channel+10
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
(SUA(nchannel).time_bins_edges(1:end-1))+time_step/2;

%smooth SUA activity with a gaussian kernel
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
CA1_MUA_zscore=zscore(filtfilt(w,1,histcounts(spike_times,time_bins_edges)));
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


zscore_min = 0;
zscore_max = 3;

% tvec = start_sec:1/new_SR:start_sec+duration_sec;
LFP_tvec =  start_sec:1/new_SR:start_sec+duration_sec;
[replay,reactivations] = detect_candidate_events_masa(MUA_time,LFP_tvec,raw_LFP(best_ripple_channel,:),...
    CA1_MUA_zscore,CA1_spike_times,zscore_min,zscore_max,options) % Need to fix SUA spike
save extracted_candidate_events replay reactivations
% bz_FindPopBursts

% durations = [30 100]; inter-ripple 30 ms and max duration 100 ms 
% 'noise', 1 here is the top noisy channel (usually in dura gel)
[ripples] = FindRipples_masa(raw_LFP(best_ripple_channel,:)',tvec','minDuration',20,'durations',[30 100],'frequency',new_SR,'noise',raw_LFP(1,:)','passband',[125 300],'show','on')
save extracted_ripples_events ripples

[spindles] = FindSpindles_masa(raw_LFP(best_spindle_channel,:)',tvec','durations',[400 3000],'frequency',new_SR,'noise',raw_LFP(1,:)','passband',[9 17],'thresholds',1.25)
save extracted_spindles_events spindles
% bz_FindRipples

%% CSD peri-event
csd = bz_eventCSD(raw_LFP',ripples.onset',...
    'channels',sorted_config.Channel,'samplingRate',new_SR,'twin',[2 2])

csd = bz_eventCSD(raw_LFP',slow_waves.ints.UP(:,1),...
    'channels',sorted_config.Channel,'samplingRate',new_SR,'twin',[2 2])

csd = bz_eventCSD(raw_LFP',replay.onset,...
    'channels',sorted_config.Channel,'samplingRate',new_SR,'twin',[2 2])

csd = bz_eventCSD(raw_LFP',spindles.onset,...
    'channels',sorted_config.Channel,'samplingRate',new_SR,'twin',[3 3])

csd = bz_eventCSD(raw_LFP',MousePos.stimuli_onset(MousePos.stimuli_track == 1),...
    'channels',sorted_config.Channel,'samplingRate',new_SR,'twin',[1 1])

csd = bz_eventCSD(raw_LFP',MousePos.stimuli_onset(MousePos.stimuli_track == 2),...
    'channels',sorted_config.Channel,'samplingRate',new_SR,'twin',[1 1])

figure
plot(MUA_time(80000:100000),CA1_MUA_zscore(80000:100000))
hold on
xline(reactivations.onset,'r')
xlim([MUA_time(80000) MUA_time(100000)])
%% Detect Up and Down states

figure
hold on
plot(normalised_slow_wave,sorted_config.Ks_ycoord)
plot(normalised_spindle,sorted_config.Ks_ycoord)

plot(normalised_ripple,sorted_config.Ks_ycoord)
plot(normalised_theta,sorted_config.Ks_ycoord)
plot(normalised_gamma_power,sorted_config.Ks_ycoord)
plot([0.5 0.5],[0 4000],'--')
legend('Normalised slow wave','Normalised theta','Normalised ripple')
xlabel('Normalised power')
ylabel('Depth (um)')

% Select channels based on depth
channels_selected = ...
    chan_config.SpikeGLXchan0(find(chan_config.Ks_ycoord < chan_config.Ks_ycoord(best_SW_channel)...
    & chan_config.Ks_ycoord >chan_config.Ks_ycoord(best_SW_channel)-1000));

% Create Buz-style spike variables
spikes = [];
spikes.UID = [];
count = 1;
for nchannel = channels_selected(1):channels_selected(end)
    clusters_this_channel = cluster_id(find(peakChannel == nchannel))';
    good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(clusters_this_channel+1)=='good'))';% Plus one because clutser id is 0 based.

    if ~isempty(good_units_this_channel)
        for unit = 1:length(good_units_this_channel)
            if ~isempty(these_spike_times{good_units_this_channel(unit)+1})
                spikes.times{count} = these_spike_times{good_units_this_channel(unit)+1}; % Plus one because clutser id is 0 based.
                spikes.UID = [spikes.UID; good_units_this_channel(unit)];
                count = count + 1;
            end
        end
        
    end
end

slow_waves= DetectSlowWaves_masa('time',tvec,'lfp',raw_LFP(best_spindle_channel,:)','NREMInts',behavioural_state.SWS,'spikes',spikes);
% DetectSlowWaves
% fake_sws = [21 160; 190 477; 510 570];
% SWS;

save best_channels best_spindle_channel best_spindle_channel best_ripple_channel
save slow_waves slow_waves

%% Spike LFP relationship

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

% CCG of SLow Waves and spikes
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
event = MousePos.stimuli_onset(MousePos.stimuli_track == 1)
% event = reactivations.onset(reactivations.zscore>5);
% num_cell = 1;
% for types = 1:4
MUA_filter_length = 10;
SD_alpha = 2; %2 std width
MUA_filter_alpha = (10-1)/SD_alpha
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,4
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width


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
ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')


spike_data = CA3_spike_times(:,2);
num_cell = length(unique(CA3_spike_times(:,1)));
% num_cell = 1;
bin_size = 0.001;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data, event, [-0.3 1], bin_size);
psth = mean(binnedArray./num_cell./bin_size); % normalize to Hz
psth_se = std(binnedArray./num_cell./bin_size)./sqrt(length(event));

subplot(4,2,3)
plot(rasterX, rasterY,'b')
hold on
plot([0 0],get(gca,'ylim'),'r-')
ylabel('events');
xlabel('Spike time relative to event onset')
title('CA3')

subplot(4,2,4)
hold on
plot(bins,psth,'b')
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'g','FaceAlpha','0.3','LineStyle','none');
%     ylim([0 600])
hold on
plot([0 0],get(gca,'ylim'),'r-')
ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')


spike_data = DG_spike_times(:,2);
num_cell = length(unique(DG_spike_times(:,1)));
% num_cell = 1;
bin_size = 0.001;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data, event, [-0.3 1], bin_size);
psth = mean(binnedArray./num_cell./bin_size); % normalize to Hz
psth_se = std(binnedArray./num_cell./bin_size)./sqrt(length(event));

subplot(4,2,5)
plot(rasterX, rasterY,'k')
hold on
plot([0 0],get(gca,'ylim'),'r-')
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
psth = mean(binnedArray./num_cell./bin_size); % normalize to Hz
psth_se = std(binnedArray./num_cell./bin_size)./sqrt(length(event));
CA3_spike_times

subplot(4,2,7)
plot(rasterX, rasterY,'k')
hold on
plot([0 0],get(gca,'ylim'),'r-')
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
ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')
sgtitle('Ripple Onset')

sgtitle('CA1 bursting')
sgtitle('Track 2 Stimuli (Right)')
sgtitle('Track 1 Stimuli (Left)')
% end

%% CCG method (do't use)
CCGvec = [slow_waves.ints.UP(:,1),ones(size(slow_waves.ints.UP(:,1)));...
    spike_data,2.*ones(size(spike_data))];

[SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02,'norm','rate');
figure
subplot(1,2,1)
bar(t_CCG,SWspike_CCG(:,1,2)./length(slow_waves.ints.UP(:,1)),'facecolor',[0.5 0.5 0.5])
xlim([-1 1])
set(gca,'xticklabel',[]);
ylabel('Spike Rate');
xlabel('Relative to Up State')
title('CA1')
hold on
plot([0 0],get(gca,'ylim'),'r-')

spike_data = V1_spike_times(:,2);
CCGvec = [slow_waves.ints.UP(:,1),ones(size(slow_waves.ints.UP(:,1)));...
    spike_data,2.*ones(size(spike_data))];

[SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02,'norm','rate');
subplot(1,2,2)
bar(t_CCG,SWspike_CCG(:,1,2)./length(slow_waves.ints.UP(:,1)),'facecolor',[0.5 0.5 0.5])
xlim([-1 1])
set(gca,'xticklabel',[]);
ylabel('Spike Rate');
xlabel('Relative to Up State')
sgtitle('CA1 and V1 Populational spiking relative to Up State')
title('V1')
hold on
plot([0 0],get(gca,'ylim'),'r-')



spike_data = CA1_spike_times(:,2);
CCGvec = [ripples.onset,ones(size(ripples.peaktimes));...
    spike_data,2.*ones(size(spike_data))];

[SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02,'norm','rate');
figure
subplot(1,2,1)
bar(t_CCG,SWspike_CCG(:,1,2)./length(ripples.peaktimes),'facecolor',[0.5 0.5 0.5])
xlim([-0.3 0.3])
set(gca,'xticklabel',[]);
ylabel('Spike Rate');
xlabel('Relative to ripple onset')
title('CA1')
hold on
plot([0 0],get(gca,'ylim'),'r-')

spike_data = V1_spike_times(:,2);
CCGvec = [ripples.onset,ones(size(ripples.peaktimes));...
    spike_data,2.*ones(size(spike_data))];

[SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02,'norm','rate');
subplot(1,2,2)
bar(t_CCG,SWspike_CCG(:,1,2)./length(ripples.peaktimes),'facecolor',[0.5 0.5 0.5])
xlim([-0.3 0.3])
set(gca,'xticklabel',[]);
ylabel('Spike Rate');
xlabel('Relative to ripple onset')
sgtitle('CA1 and V1 Populational spiking relative to ripple onset')
title('V1')
hold on
plot([0 0],get(gca,'ylim'),'r-')


spike_data = CA1_spike_times(:,2);
CCGvec = [reactivations.onset',ones(size(reactivations.onset'));...
    spike_data,2.*ones(size(spike_data))];

[SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02,'norm','rate');
figure
subplot(1,2,1)
bar(t_CCG,SWspike_CCG(:,1,2)./length(reactivations.onset),'facecolor',[0.5 0.5 0.5])
xlim([-1 1])
set(gca,'xticklabel',[]);
ylabel('Spike Rate');
xlabel('Relative to populational bursting')
title('CA1')
hold on
plot([0 0],get(gca,'ylim'),'r-')

spike_data = V1_spike_times(:,2);
CCGvec = [reactivations.onset',ones(size(reactivations.onset'));...
    spike_data,2.*ones(size(spike_data))];

[SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02,'norm','rate');
subplot(1,2,2)
bar(t_CCG,SWspike_CCG(:,1,2)./length(reactivations.onset),'facecolor',[0.5 0.5 0.5])
xlim([-1 1])
set(gca,'xticklabel',[]);
ylabel('Spike Rate');
xlabel('Relative to populational bursting')
sgtitle('CA1 and V1 Populational spiking relative to CA1 Populational bursting')
title('V1')
hold on
plot([0 0],get(gca,'ylim'),'r-')





figure
spike_data = CA1_spike_times(:,2);

for n = 1:2
    CCGvec = [MousePos.stimuli_onset(MousePos.stimuli_track == n),ones(size(MousePos.stimuli_onset(MousePos.stimuli_track == n)));...
        spike_data,2.*ones(size(spike_data))];
    [SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02,'norm','rate');

    subplot(1,2,n)
    bar(t_CCG,SWspike_CCG(:,1,2)./length(MousePos.stimuli_onset(MousePos.stimuli_track == n)),'facecolor',[0.5 0.5 0.5])
    xlim([-1 1])
    hold on

    set(gca,'xticklabel',[]);
    ylabel('Spike Rate');
    xlabel('Relative to track stimulus')
    title(sprintf('Track %i stimulus',n))
    ylim([0 3])
    sgtitle('CA1')

    plot([0 0],get(gca,'ylim'),'r-')
    plot([0.04 0.04],get(gca,'ylim'),'k-')
    plot([0.1 0.1],get(gca,'ylim'),'b-')

end

figure
spike_data = V1_spike_times(:,2);

for n = 1:2
    CCGvec = [MousePos.stimuli_onset(MousePos.stimuli_track == n),ones(size(MousePos.stimuli_onset(MousePos.stimuli_track == n)));...
        spike_data,2.*ones(size(spike_data))];

    [t_spkmat,inNREMidx] = RestrictInts(spike_data,[MousePos.stimuli_onset(MousePos.stimuli_track == n)-1 MousePos.stimuli_onset(MousePos.stimuli_track == n)+1]); %Replace with InInterval
    
    
    [SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02,'norm','rate');

    subplot(1,2,n)
    bar(t_CCG,SWspike_CCG(:,1,2)./length(MousePos.stimuli_onset(MousePos.stimuli_track == n)),'facecolor',[0.5 0.5 0.5])
    xlim([-1 1])
    hold on

    set(gca,'xticklabel',[]);
    ylabel('Spike Rate');
    xlabel('Relative to track stimulus')
    title(sprintf('Track %i stimulus',n))
    ylim([0 40])
    sgtitle('V1')
    plot([0 0],get(gca,'ylim'),'r-')
    plot([0.04 0.04],get(gca,'ylim'),'k-')
    plot([0.1 0.1],get(gca,'ylim'),'b-')
end


MousePos.stimuli_onset(MousePos.stimuli_track == 1)
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
gs = [2];
nChannelsToBin = 24;    % Bin firing rates across channels for RF maps
channelRange = [180 310]; % Look at these channels 

% Some defaults
% options.importMode = 'KS'; % LF or MUA or KS
options.importMode = 'LF'; % LF or MUA or KS
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
% [resps,otherData,stimData,~,~,photodiodeData,timeVector,options] = extractAndCollateNPData(options);
% [resps,~,stimData,~,~,~,timeVector,options] = extractAndCollateNPData(options);

%% LFP putative theta and ripple visualisation
% https://github.com/cortex-lab/neuropixels/issues/40
% How to read specific files
%

% 300 clips of non-overlapping 1 seconds ()
lfpFilename = 'D:\Neuropixel_recording\M22069\20221130\M22069_20221130_g0\M22069_20221130_g0_imec0\M22069_20221130_g0_t0.imec0.lf.bin';
lfpFilename = '/research/DATA/SUBJECTS/M22069/ephys/20221201/M22069_20221201_g0/M22069_20221201_g0_imec0/M22069_20221201_g0_t0.imec0.lf.bin';

lfpFilename = 'X:\ibn-vision\DATA\SUBJECTS\M22069\ephys\20221201\M22069_20221201_g0\M22069_20221201_g0_imec0/M22069_20221201_g0_t0.imec0.lf.bin'
lfpFs = 2500;
freqBand = [125 300];
nChansInFile = 385;
% [lfpByChannel, allPowerEst, allPowerEstByBand,F, allPowerVar] = lfpBandPower(lfpFilename, lfpFs, nChansInFile, freqBand);
% plot(lfpByChannel)

[lfpByChannel_ripple, allPowerEst_ripple, allPowerEstByBand_ripple,F, allPowerVar_ripple] = lfpBandPower(lfpFilename, lfpFs, nChansInFile, freqBand);

% caxis([0 200])

freqBand = [4 12];
[lfpByChannel_theta, allPowerEst_theta, allPowerEstByBand_theta,F, allPowerVar_theta] = lfpBandPower(lfpFilename, lfpFs, nChansInFile, freqBand);
% plot(lfpByChannel_theta)

subplot(2,1,1)
imagesc(flip(allPowerEstByBand_theta,2)')
colorbar
title('Theta (4-12Hz) Power')

subplot(2,1,2)
imagesc(flip(allPowerEstByBand_ripple,2)')
colorbar
title('Ripple (125-300Hz) Power')
caxis([0 0.5])
% caxis([0 200])

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

%% Theta phase shift
parameters = list_of_parameters;
filter_type  = 'bandpass';
filter_width = [4 12];                 % range of frequencies in Hz you want to filter between
filter_order = round(6*new_SR/(max(filter_width)-min(filter_width)));  % creates filter with length 0.417 s for theta
norm_freq_range = filter_width/(new_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_theta = fir1(filter_order, norm_freq_range,filter_type);

for nchannel = 1:size(chan_config,1)
    LFP.theta(nchannel,:) = filtfilt(b_theta,1,raw_LFP(nchannel,:));
    LFP.theta_zscore(nchannel,:)  = zscore(abs(hilbert(LFP.theta(nchannel,:))));
    LFP.theta_phase(nchannel,:)  = angle(hilbert(LFP.theta(nchannel,:)));
end

col_ID = cols_available(1); % (1) is 11
sorted_config = sortrows(chan_config,'Ks_ycoord','descend');
col_idx = sorted_config.Ks_xcoord == col_ID;
sorted_config = sorted_config(col_idx,:);

columns_available = unique(chan_config.Ks_xcoord);
phase_diff_mean  = NaN(size(sorted_config,1),numel(columns_available));
phase_diff_std   = NaN(size(phase_diff_mean));
s = RandStream('mrg32k3a','Seed',1);
idx = datasample(s,1:length(LFP.theta(1,:)),2000,'Replace',false);% Distribution of theta phase from each channel
sample_phases = [];

% numel(columns_available)
for n = 1:1
    col_ID = cols_available(n); % (1) is 11
    sorted_config = sortrows(chan_config,'Ks_ycoord','descend');
    col_idx = sorted_config.Ks_xcoord == col_ID;
    sorted_config = sorted_config(col_idx,:);

    ref_chan = interp1(sorted_config.Electrode,1:size(sorted_config,1),200,'nearest'); % edit field gives electrode, find nearest recorded electrode, this gives index for mapping applied below
    
    for nchannel = 1:size(sorted_config,1)
        sample_phases(:,nchannel) = LFP.theta_phase(sorted_config.Channel(nchannel),idx);
    end

    % circ_ functions below from CircStats toolbox
    phase_diff      = circ_dist(sample_phases,repmat(sample_phases(:,ref_chan),1,size(sample_phases,2)));
    phase_diff_mean(1:size(phase_diff,2),n) = unwrap(circ_mean(phase_diff));
    phase_diff_std(1:size(phase_diff,2),n)  = circ_std(phase_diff);
end
sorted_config.Channel
subplot(1,4,3)
p(n) = plot(phase_diff_mean(:,1),sorted_config.Ks_ycoord','k')
hold on
plot(phase_diff_mean(:,1)-phase_diff_std(:,1),sorted_config.Ks_ycoord','color',[0 0 0 0.5],'LineWidth',1,'LineStyle',":")
plot(phase_diff_mean(:,1)+phase_diff_std(:,1),sorted_config.Ks_ycoord','color',[0 0 0 0.5],'LineWidth',1,'LineStyle',":")
xlabel('Radian')
ylabel('Distance (um)')


%% Filter LFP

% Ripple band
parameters = list_of_parameters;
filter_type  = 'bandpass';
filter_width = [125 300];                 % range of frequencies in Hz you want to filter between
filter_order = round(6*new_SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
norm_freq_range = filter_width/(new_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_ripple = fir1(filter_order, norm_freq_range,filter_type);

for nchannel = 1:size(chan_config,1)
    LFP.ripple(nchannel,:) = filtfilt(b_ripple,1,raw_LFP(nchannel,:));
%     LFP.ripple_zscore(nchannel,:)  = zscore(abs(hilbert(LFP.ripple(nchannel,:))));
end

% spindle
parameters = list_of_parameters;
filter_type  = 'bandpass';
filter_width = [9 17];                 % range of frequencies in Hz you want to filter between
filter_order = round(6*new_SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
norm_freq_range = filter_width/(new_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_spindle = fir1(filter_order, norm_freq_range,filter_type);

for nchannel = 1:size(chan_config,1)
    LFP.spindle(nchannel,:) = filtfilt(b_spindle,1,raw_LFP(nchannel,:));
%     LFP.spindle_zscore(nchannel,:)  = zscore(abs(hilbert(LFP.spindle(nchannel,:))));
end


% Slow wave oscilation band
parameters = list_of_parameters;
filter_type  = 'bandpass';
filter_width = [0.5 3];                 % range of frequencies in Hz you want to filter between
filter_order = round(6*new_SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
norm_freq_range = filter_width/(new_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_SO = fir1(filter_order, norm_freq_range,filter_type);

for nchannel = 1:size(chan_config,1)
    LFP.slow_oscillation(nchannel,:) = filtfilt(b_SO,1,raw_LFP(nchannel,:));
%     LFP.slow_oscillation_zscore(nchannel,:)  = zscore(abs(hilbert(LFP.slow_oscillation(nchannel,:))));
end


% Gamma band
parameters = list_of_parameters;
filter_type  = 'bandpass';
filter_width = [30 100];                 % range of frequencies in Hz you want to filter between
filter_order = round(6*new_SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
norm_freq_range = filter_width/(new_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_gamma = fir1(filter_order, norm_freq_range,filter_type);

for nchannel = 1:size(chan_config,1)
    LFP.gamma(nchannel,:) = filtfilt(b_gamma,1,raw_LFP(nchannel,:));
%     LFP.gamma_zscore(nchannel,:)  = zscore(abs(hilbert(LFP.gamma(nchannel,:))));
end

LFP.timevec = tvec(1:end-1);
save LFP LFP -v7.3

%% MUA

[AP_FILE,LF_FILE] = findImecBinFile(options.EPHYS_DATAPATH);
imecMeta = NPadmin.ReadMeta(AP_FILE,options.EPHYS_DATAPATH);
chan_config = NPadmin.getNPChannelConfig(imecMeta);

KS_metrics = readtable(fullfile(options.KS_DATAPATH,'metrics.csv'));
[these_spike_times,nominal_KSLabel,cluster_id,peakChannel,maxSpkTime] = import_ks_spiketimes(options.KS_DATAPATH,options.gFileNum,options.KS_CATGT_FNAME,imecMeta.imSampRate)
mean_waveforms = readNPY([options.KS_DATAPATH,'/mean_waveforms.npy']); % Information about mean waveform
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
    clusters_this_channel = cluster_id(find(peakChannel == nchannel));
    good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(clusters_this_channel+1)=='good'));% Plus one because clutser id is 0 based.
    
    MUA(nchannel).zscore = zeros(1,length(time_bins_edges));
    MUA(nchannel).time_bins_edges= time_bins_edges;
    MUA(nchannel).time_bins =(MUA(nchannel).time_bins_edges(1:end-1))+time_step/2;


    SUA(nchannel).zscore = zeros(1,length(time_bins_edges));
    SUA(nchannel).time_bins_edges= time_bins_edges;
    SUA(nchannel).time_bins =(SUA(nchannel).time_bins_edges(1:end-1))+time_step/2;
    
    MUA(nchannel).cluster_ID = clusters_this_channel;
    if ~isempty(clusters_this_channel) % If any clusters
        MUA(nchannel).spike_times = [];
        for unit = 1:length(clusters_this_channel)
            MUA(nchannel).spike_times = [MUA(nchannel).spike_times; these_spike_times{clusters_this_channel(unit)+1}]; % Plus one because clutser id is 0 based.
        end
        %smooth MUA activity with a gaussian kernel

        %         time_step=0.001; %1 ms timestep

        w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
        w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
        %         MUA(nchannel).zscored=zscore(histcounts(MUA(nchannel).spike_times,MUA(nchannel).time_bins_edges));
        MUA(nchannel).zscore=zscore(filtfilt(w,1,histcounts(MUA(nchannel).spike_times,MUA(nchannel).time_bins_edges)));


        if ~isempty(good_units_this_channel) % If any SUA
             SUA(nchannel).spike_times= [];
            for unit = 1:length(good_units_this_channel)

                SUA(nchannel).spike_times = [SUA(nchannel).spike_times; these_spike_times{good_units_this_channel(unit)+1}]; % Plus one because clutser id is 0 based.
            end
            SUA(nchannel).unit_ID = good_units_this_channel;
            %smooth SUA activity with a gaussian kernel
            w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
            w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
            %         MUA(nchannel).zscored=zscore(histcounts(MUA(nchannel).spike_times,MUA(nchannel).time_bins_edges));
            SUA(nchannel).zscore=zscore(filtfilt(w,1,histcounts(SUA(nchannel).spike_times,SUA(nchannel).time_bins_edges)));
        end
    end
end
toc

% Sort channel for spike time data
col_ID = cols_available(1); % (1) is 11
sorted_config_spikes = sortrows(chan_config,'Ks_ycoord','descend');
col_idx = sorted_config_spikes.Ks_xcoord == col_ID;
sorted_config_spikes = sorted_config_spikes(col_idx,:);

sample_to_view = 100100:167100;
figure(1)
subplot(1,4,1)
for nchannel = 1:size(sorted_config,1)
    plot(tvec(sample_to_view),(LFP.slow_oscillation(sorted_config.Channel(nchannel),sample_to_view)*30000 + sorted_config.Ks_ycoord(nchannel)),'k')
    hold on
end
xlabel('Time (s)')
ylabel('Distance (um)')
title('Slow wave oscillation (0.5 - 3 Hz)')
ylim([0 4000])

for nchannel = 1:size(sorted_config,1)
    plot(tvec(sample_to_view),(LFP.ripple(sorted_config.Channel(nchannel),sample_to_view)*100000 + sorted_config.Ks_ycoord(nchannel)),'r')
    hold on
end
xlabel('Time (s)')
ylabel('Distance (um)')
% title('Ripple (125 - 300 Hz)')
ylim([0 4000])

subplot(1,4,2)
for nchannel = 1:size(sorted_config,1)
    plot(tvec(sample_to_view),(LFP.theta(sorted_config.Channel(nchannel),sample_to_view)*30000 + sorted_config.Ks_ycoord(nchannel)),'k')
    hold on
end
xlabel('Time (s)')
ylabel('Distance (um)')
title('Theta (4 - 12 Hz)')
ylim([0 4000])

subplot(1,4,3)
for nchannel = 1:size(sorted_config,1)
    plot(tvec(sample_to_view),(LFP.ripple(sorted_config.Channel(nchannel),sample_to_view)*100000 + sorted_config.Ks_ycoord(nchannel)),'k')
    hold on
end
xlabel('Time (s)')
ylabel('Distance (um)')
title('Ripple (125 - 300 Hz)')
ylim([0 4000])


subplot(1,4,4)
for nchannel = 1:size(sorted_config,1)
    plot(tvec(sample_to_view),(MUA(sorted_config_spikes.SpikeGLXchan0(nchannel)).zscore(sample_to_view)*5 + sorted_config.Ks_ycoord(nchannel)),'k')
    hold on
end

xlabel('Time (s)')
ylabel('Distance (um)')
title('Smoothed MUA zsocre activity')
ylim([0 4000])

for nchannel = 1:size(sorted_config,1)
    plot(tvec(sample_to_view),(LFP.ripple(sorted_config.Channel(nchannel),sample_to_view)*100000 + sorted_config.Ks_ycoord(nchannel)),'k')
    hold on
end
xlabel('Time (s)')
ylabel('Distance (um)')
% title('Ripple (125 - 300 Hz)')
ylim([0 4000])

% 
% subplot(1,4,4)
% for nchannel = 1:size(sorted_config,1)
%     plot(tvec(sample_to_view),(MUA(sorted_config_spikes.Channel(nchannel)).zscored(sample_to_view)*5 + sorted_config.Ks_ycoord(nchannel)),'k')
%     hold on
% end
% 
% xlabel('Time (s)')
% ylabel('Distance (um)')
% title('Smoothed MUA zsocre activity')
% ylim([0 4000])
% 
% for nchannel = 1:size(sorted_config,1)
%     spikes(nchannel) = length(MUA(sorted_config_spikes.Channel(nchannel)).spike_times);
%     spikes_GLX_channel(nchannel) = length(MUA(sorted_config_spikes.SpikeGLXchan0(nchannel)).spike_times);
% end
% subplot(1,4,1)
% plot(spikes_GLX_channel,sorted_config.Ks_ycoord,'k')
% hold on
% plot(spikes,sorted_config.Ks_ycoord,'r')
% legend('SpikeGLXchan0','Channel')


%% LFP gamma coherence -> Gradient descent

% Gamma band
parameters = list_of_parameters;
filter_type  = 'bandpass';
filter_width = [30 100];                 % range of frequencies in Hz you want to filter between
filter_order = round(6*new_SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
norm_freq_range = filter_width/(new_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_gamma = fir1(filter_order, norm_freq_range,filter_type);

for nchannel = 1:size(chan_config,1)
    LFP.gamma(nchannel,:) = filtfilt(b_gamma,1,raw_LFP(nchannel,:));
%     LFP(nchannel).gamma_zscore = zscore(abs(hilbert(LFP(nchannel).gamma)));
end


for n = 1:size(sorted_config,1)
    tic
    for m = 1:size(sorted_config,1)
%         [coherence,f] = mscohere(raw_LFP(sorted_config.Channel(n),:),raw_LFP(sorted_config.Channel(m),:),[],[],[],new_SR);

         filt_hilb1 = hilbert(LFP.gamma(sorted_config.Channel(n),:)); %calculates the Hilbert transform of eeg1
         amp1 = abs(filt_hilb1);%calculates the instantaneous amplitude of eeg1 filtered between low_freq and high_freq
         amp1=amp1-mean(amp1); %removes mean of the signal because the DC component of a signal does not change the correlation
         filt_hilb2 = hilbert(LFP.gamma(sorted_config.Channel(m),:));%calculates the Hilbert transform of eeg2
         amp2 = abs(filt_hilb2);%calculates the instantaneous amplitude of eeg2 filtered between low_freq and high_freq
         amp2=amp2-mean(amp2);
         [crosscorr,lags]=xcorr(amp1, amp2,round(new_SR/10),'coeff'); %calculates crosscorrelations between amplitude vectors
         lags=(lags./new_SR)*1000; %converts lags to miliseconds
         g=find(crosscorr==max(crosscorr));%identifies index where the crosscorrelation peaks
         max_crosscorr_lag=lags(g);%identifies the lag at which the crosscorrelation peaks
         gamma_phase_coherence(n,m) = unwrap(circ_mean(circ_dist(angle(filt_hilb1),angle(filt_hilb2)),[],2));
         gamma_coherence(n,m) = max(crosscorr);
%          gamma_coherence_ms(n,m) = mean(coherence(find(f <= 100 & f>= 30)));
%          [coherogram,phase,t,f]  = bz_MTCoherogram(raw_LFP(sorted_config.Channel(n),:)',raw_LFP(sorted_config.Channel(m),:)','window',5,'frequency',new_SR,'range',[30 100]);
%          gamma_coherence_bz(n,m) = mean(mean(coherogram));
%          gamma_phase_coherence_bz(n,m) = circ_mean(reshape(phase,size(phase,1)*size(phase,2),1));
% %          gamma_coherence2(n,m) = mean(coherence( find(f>=30 & f<=100)));
%          figure('color',[1 1 1])
%          plot(lags, crosscorr,'color',[0 0 1],'linewidth',2),hold on %plots crosscorrelations
%          plot(lags(g),crosscorr(g),'rp','markerfacecolor',[1 0 0],'markersize',10)%plots marker at the peak of the cross correlation
%          plot([0 0],[1.05*max(crosscorr) 0.95*min(crosscorr)],'color',[0 0 0],'linestyle',':', 'linewidth',2) %plots dashed line at zero lag
% %          set(gca,'xtick',[-100 -50 0 50 100])
% %          axis tight, box off, xlim([-101 100])
%          xlabel('Lag (ms)','fontsize',14)
%          ylabel('Crosscorrelation','fontsize',14)

    end
    toc
end


save gamma_coherence gamma_coherence gamma_phase_coherence

save gamma_coherence_bz gamma_coherence_bz gamma_phase_coherence_bz
% [coherogram,phase,~,~,~,t,f] = cohgramc(LFP(n).gamma_zscore',LFP(m).gamma_zscore',[window window-window/2],parameters);
% 
%           [coherogram,phase,t,f]  = bz_MTCoherogram(LFP(n).gamma_zscore',LFP(m).gamma_zscore','frequency',new_SR,'range',[30 100]);
% [coherogram,phase,t,f]  = bz_MTCoherogram(raw_LFP(sorted_config.Channel(n),:)',raw_LFP(sorted_config.Channel(m),:)','frequency',new_SR,'range',[30 100]);
hold on
cutoffs = [0 1];
figure;hold on;
subplot(2,1,1);
PlotColorMap(coherogram,'x',t,'y',f,'cutoffs',cutoffs,'newfig','off');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Coherogram Amplitude');
colorbar
% 
subplot(2,1,2);
PlotColorMap(phase,'x',t,'y',f,'cutoffs',[-pi pi],'newfig','off');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Coherogram Phase');
colorbar
% 
% parameters.Fs = new_SR;
% parameters.tapers = [3 5];
% parameters.pad = 0;

% -11 here because it is in dura gel
% brainmap{1} = 1:1:96;
% brainmap{2} = [96 96];
% [ cluass,cluNRG ] = bz_GradDescCluster(gamma_coherence,'brainmap',brainmap)

% Find where does the probe go into the brain (only brain signal for gamma coherence clustering analysis)
% Find first peak that suddenly increases the power
power_differnece = [0; diff(power(:,6)./max(power(:,6)))]; % Use ripple or theta
[~,first_in_brain_channel] = findpeaks(power_differnece,'MinPeakHeight',0.1);
first_in_brain_channel = first_in_brain_channel(1);
% power_differnece = [0; diff(power(:,2)./max(power(:,2)))]; % Use ripple and theta
% [~,theta_loc] = findpeaks(power_differnece,'MinPeakHeight',0.1);

probe_length_in_brain = sorted_config.Ks_ycoord(first_in_brain_channel(1));

brainmap{1} = 1:1:96-first_in_brain_channel(1)+1;
brainmap{2} = [96-first_in_brain_channel(1)+1 96-first_in_brain_channel(1)+1];
[ cluass,cluNRG ] = bz_GradDescCluster(gamma_coherence(first_in_brain_channel(1):end,first_in_brain_channel(1):end),'brainmap',brainmap)

% Plot cross-coherence matrix
figure
subplot(1,2,1)
imagesc(gamma_coherence)
colorbar
xticks(1:2:96)
xticklabels(sorted_config.Ks_ycoord(1:2:96))
yticks(1:2:96)
yticklabels(sorted_config.Ks_ycoord(1:2:96))
hold on
cluster_boundary = find(diff(cluass)~=0);
plot([0 96],[first_in_brain_channel-1 first_in_brain_channel-1],'--k','LineWidth',2)
plot([first_in_brain_channel-1 first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
hold on

for n = 1:length(cluster_boundary)
    plot([0 96],[cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],'--k','LineWidth',2)
    plot([cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
%         plot([0 96],[cluster_boundary(n)+11 cluster_boundary(n)+11],'--k','LineWidth',2)
%     plot([cluster_boundary(n)+11 cluster_boundary(n)+11],[0 96],'--k','LineWidth',2)
    hold on
end
title('LFP Gamma Coherence (30 - 100Hz)')

subplot(1,2,2)
imagesc(gamma_phase_coherence)
colorbar

xticks(1:2:96)
xticklabels(sorted_config.Ks_ycoord(1:2:96))
yticks(1:2:96)
yticklabels(sorted_config.Ks_ycoord(1:2:96))
hold on
cluster_boundary = find(diff(cluass)~=0);
plot([0 96],[11 11],'--k','LineWidth',2)
plot([11 11],[0 96],'--k','LineWidth',2)
hold on

for n = 1:length(cluster_boundary)
    plot([0 96],[cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],'--k','LineWidth',2)
    plot([cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
    hold on
end
colorbar
title('LFP Gamma Phase Coherence (30 - 100Hz)')


sample_to_view = 100100:180000;
figure(1)
subplot(1,4,1)
for nchannel = 1:size(sorted_config,1)
    plot(tvec(sample_to_view),(LFP.gamma(sorted_config.Channel(nchannel),sample_to_view)*100000 + sorted_config.Ks_ycoord(nchannel)),'k')
    hold on
end
xlabel('Time (s)')
ylabel('Distance (um)')
title('Gamma (30 - 100 Hz)')
ylim([0 4000])

subplot(1,4,2)
for nchannel = 1:size(sorted_config,1)
    plot(tvec(sample_to_view),(LFP.theta(sorted_config.Channel(nchannel),sample_to_view)*100000 + sorted_config.Ks_ycoord(nchannel)),'k')
    hold on
end
xlabel('Time (s)')
ylabel('Distance (um)')
title('Theta (4 - 12 Hz)')
ylim([0 4000])

subplot(1,4,3)
for nchannel = 1:size(sorted_config,1)
    plot(tvec(sample_to_view),(MUA(sorted_config_spikes.SpikeGLXchan0(nchannel)).zscore(sample_to_view)*5 + sorted_config.Ks_ycoord(nchannel)),'k')
    hold on
end

ylim([0 4000])
xlabel('Time (s)')
ylabel('Distance (um)')
title('Smoothed MUA zsocre activity')

subplot(1,4,4)
cluster_boundary = find(diff(cluass)~=0);
plot([0 1],[sorted_config.Ks_ycoord(11) sorted_config.Ks_ycoord(11)],'--k','LineWidth',2)
hold on

for n = 1:length(cluster_boundary)
    plot([0 1],[sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1) sorted_config.Ks_ycoord(cluster_boundary(n)+first_in_brain_channel-1)],'--k','LineWidth',2)
    hold on
end
ylim([0 4000])

% Quick plotting of PSD
colour_line= {'k','r','m','b','c','g','y'};
freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','350 - 625 Hz'};
selected_powers = [1 2 3 6];

for n = selected_powers
        p(n) = plot(power(:,n)./max(power(:,n)),sorted_config.Ks_ycoord',colour_line{n})
    
%     p(n) = plot(power_differnece,1:96,colour_line{n})
    hold on
end
% legend('1-3','4-12','30-60','60-100','125-300')
legend([p(selected_powers)],{freq_legends{selected_powers}});
ylim([0 4000])
title('Gradient descent clusters (layer boundary)')

% chan_config.SpikeGLXchan0(find(chan_config.Ks_ycoord < 3420 & chan_config.Ks_ycoord >2340))
%% CSD gamma band
 lfp = [];
for nchannel = 1:size(sorted_config,1)
    lfp.data(:,nchannel) = LFP(sorted_config.Channel(nchannel)).gamma;% switch to data X nchannel
end
lfp.timestamps = tvec;
lfp.samplingRate = new_SR;
[ csd ] = bz_CSD (lfp);

% [ csd ] = bz_eventCSD (lfp);
gamma_coherence_CSD = zeros(size(sorted_config,1),size(sorted_config,1));
gamma_phase_coherence_CSD = zeros(size(sorted_config,1),size(sorted_config,1));
for n = 1:size(sorted_config,1)-2
    tic
    for m = 1:size(sorted_config,1)-2
%          filt_hilb1 = hilbert(csd.data(:,n))'; %calculates the Hilbert transform of eeg1
%          filt_hilb2 = hilbert(csd.data(:,m))';%calculates the Hilbert transform of eeg2

         [coherogram,phase,t,f]  = bz_MTCoherogram(csd.data(:,n),csd.data(:,m),'window',5,'frequency',new_SR,'range',[30 100]);
         gamma_coherence_CSD(n+1,m+1) = mean(mean(coherogram));
         gamma_phase_coherence_CSD(n+1,m+1) = circ_mean(reshape(phase,size(phase,1)*size(phase,2),1));
%          gamma_phase_coherence_CSD(n+1,m+1) = unwrap(circ_mean(circ_dist(angle(filt_hilb1),angle(filt_hilb2)),[],2));
        
    end
    toc
end

save gamma_coherence_CSD gamma_coherence_CSD gamma_phase_coherence_CSD

figure
subplot(1,2,1)
imagesc(gamma_coherence_CSD)
colorbar
xticks(1:2:96)
xticklabels(sorted_config.Ks_ycoord(1:2:96))
yticks(1:2:96)
yticklabels(sorted_config.Ks_ycoord(1:2:96))
hold on
cluster_boundary = find(diff(cluass)~=0);
plot([0 96],[first_in_brain_channel-1 first_in_brain_channel-1],'--k','LineWidth',2)
plot([first_in_brain_channel-1 first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
hold on

for n = 1:length(cluster_boundary)
    plot([0 96],[cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],'--k','LineWidth',2)
    plot([cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
%         plot([0 96],[cluster_boundary(n)+11 cluster_boundary(n)+11],'--k','LineWidth',2)
%     plot([cluster_boundary(n)+11 cluster_boundary(n)+11],[0 96],'--k','LineWidth',2)
    hold on
end
title('CSD Gamma Coherence (30 - 100Hz)')
subplot(1,2,2)
imagesc(gamma_phase_coherence_CSD)
colorbar
xticks(1:2:96)
xticklabels(sorted_config.Ks_ycoord(1:2:96))
yticks(1:2:96)
yticklabels(sorted_config.Ks_ycoord(1:2:96))
hold on
cluster_boundary = find(diff(cluass)~=0);
plot([0 96],[first_in_brain_channel-1 first_in_brain_channel-1],'--k','LineWidth',2)
plot([first_in_brain_channel-1 first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
hold on

for n = 1:length(cluster_boundary)
    plot([0 96],[cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],'--k','LineWidth',2)
    plot([cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
%         plot([0 96],[cluster_boundary(n)+11 cluster_boundary(n)+11],'--k','LineWidth',2)
%     plot([cluster_boundary(n)+11 cluster_boundary(n)+11],[0 96],'--k','LineWidth',2)
    hold on
end
title('CSD Phase Gamma Coherence (30 - 100Hz)')



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

speed_interp = interp1(peripherals.sglxTime,speed,tvec','linear');

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
[MousePos] = alignBonsaiToPhotodiode(MousePos,photodiodeData.stim_on.sglxTime);

[pks,locs]= findpeaks(abs(MousePos.pos-200));

figure
plot(MousePos.sglxTime,MousePos.pos)
hold on;
scatter(MousePos.sglxTime(locs),400*ones(1,length(locs)))

MousePos.stimuli_onset = MousePos.sglxTime_corrected(locs);
MousePos.stimuli_id = MousePos.pos(locs);
MousePos.stimuli_track = MousePos.pos(locs);
MousePos.stimuli_track(MousePos.pos(locs)<200) = 2;
MousePos.stimuli_track(MousePos.pos(locs)>200) = 1;


%% Sleep detection

normalised_ripple = power(:,6)/max(power(:,6)); %power normalized by max power across channel
normalised_theta = power(:,2)/max(power(:,2)); %power normalized by max power across channel
normalised_slow_wave = power(:,1)./max(power(:,1));
normalised_spindle = power(:,3)./max(power(:,3));

[value, candidate_index] = findpeaks(normalised_slow_wave,'MinPeakHeight',0.3); % Find channels with normalised Slow wave power bigger than 0.3
[~,index]  = max(normalised_slow_wave(candidate_index)-normalised_ripple(candidate_index));  %find channel with maximum power difference between gamma and slow wave
best_SW_channel_this_column = candidate_index(1); % First peak (superficial layer of slow wave)
best_SW_channel = sorted_config.Channel(best_SW_channel_this_column+1);

[value, candidate_index] = findpeaks(normalised_spindle,'MinPeakHeight',0.3); % Find channels with normalised Slow wave power bigger than 0.3
[~,index]  = max(normalised_spindle(candidate_index)-normalised_ripple(candidate_index));  %find channel with maximum power difference between theta and slow wave
best_spindle_channel_this_column = candidate_index(1); % First peak 
best_spindle_channel = sorted_config.Channel(best_spindle_channel_this_column+1);


% spindle
parameters = list_of_parameters;
filter_type  = 'bandpass';
filter_width = [30 100];                 % range of frequencies in Hz you want to filter between
filter_order = round(6*new_SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
norm_freq_range = filter_width/(new_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_gamma = fir1(filter_order, norm_freq_range,filter_type);

for nchannel = 1:size(chan_config,1)
    LFP.gamma(nchannel,:) = filtfilt(b_gamma,1,raw_LFP(nchannel,:));
    LFP.gamma_zscore(nchannel,:)  = zscore(abs(hilbert(LFP.gamma(nchannel,:))));
end


[freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
    [tvec' raw_LFP(best_SW_channel,:)'],[tvec' raw_LFP(best_ripple_channel,:)'],...
    [tvec' speed_interp],speedTreshold);


behavioural_state.freezing = freezing;
behavioural_state.quietWake = quietWake;
behavioural_state.SWS = SWS;
behavioural_state.REM = REM;
behavioural_state.movement = movement;

save behavioural_state behavioural_state

best_SW_channel;
best_spindle_channel;
% High spindle power as quietwake
% Low spindle High theta as freezing

figure
subplot(1,2,1)
imagesc(gamma_phase_coherence_CSD)
colorbar
xticks(1:2:96)
xticklabels(sorted_config.Ks_ycoord(1:2:96))
yticks(1:2:96)
yticklabels(sorted_config.Ks_ycoord(1:2:96))
hold on
cluster_boundary = find(diff(cluass)~=0);
plot([0 96],[first_in_brain_channel-1 first_in_brain_channel-1],'--k','LineWidth',2)
plot([first_in_brain_channel-1 first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
hold on
for n = 1:length(cluster_boundary)
    plot([0 96],[cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],'--k','LineWidth',2)
    plot([cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
%         plot([0 96],[cluster_boundary(n)+11 cluster_boundary(n)+11],'--k','LineWidth',2)
%     plot([cluster_boundary(n)+11 cluster_boundary(n)+11],[0 96],'--k','LineWidth',2)
    hold on
end
title('CSD Phase Gamma Coherence (30 - 100Hz)')
% xlim([find(sorted_config.Channel == best_ripple_channel) 96])
% ylim([find(sorted_config.Channel == best_ripple_channel) 96])

subplot(1,2,2)
imagesc(gamma_coherence)
colorbar
xticks(1:2:96)
xticklabels(sorted_config.Ks_ycoord(1:2:96))
yticks(1:2:96)
yticklabels(sorted_config.Ks_ycoord(1:2:96))
hold on
cluster_boundary = find(diff(cluass)~=0);
plot([0 96],[first_in_brain_channel-1 first_in_brain_channel-1],'--k','LineWidth',2)
plot([first_in_brain_channel-1 first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
hold on
for n = 1:length(cluster_boundary)
    plot([0 96],[cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],'--k','LineWidth',2)
    plot([cluster_boundary(n)+first_in_brain_channel-1 cluster_boundary(n)+first_in_brain_channel-1],[0 96],'--k','LineWidth',2)
%         plot([0 96],[cluster_boundary(n)+11 cluster_boundary(n)+11],'--k','LineWidth',2)
%     plot([cluster_boundary(n)+11 cluster_boundary(n)+11],[0 96],'--k','LineWidth',2)
    hold on
end
title('LFP Gamma Coherence (30 - 100Hz)')
% xlim([find(sorted_config.Channel == best_ripple_channel) 96])
% ylim([find(sorted_config.Channel == best_ripple_channel) 96])


%% Detect SWR events
% Determine channel with high ripple low theta
% F  = [0.5 3;4 12;9 17;30 60;60 100;125 300;350 625];

for n = 1:length(sorted_config.Ks_ycoord)
    spike_counts(n) = length(MUA(n).spike_times);
end
figure
plot(spike_counts,sorted_config.Ks_ycoord)

%calculate channel with biggest difference between theta and ripple normalized power.
normalised_ripple = power(:,6)/max(power(:,6)); %power normalized by max power across channel
normalised_theta = power(:,2)/max(power(:,2)); %power normalized by max power across channel
% normalised_low_gamma = power(:,3)/max(power(:,3)); %power normalized by max power across channel
% normalised_high_gamma = power(:,4)/max(power(:,4)); %power normalized by max power across channel
% normalised_ripple(find(normalised_ripple<0))=0; %if normalized ripple power is less than zero, treat it as 0
% normalised_gamma_power = gamma_power'/max(gamma_power);

[value, candidate_index] = findpeaks(normalised_ripple,'MinPeakHeight',0.5); % Find channels with normalised ripple power bigger than 0.5
[~,index]  = max(normalised_ripple(candidate_index)-normalised_theta(candidate_index));  %find channel with maximum power difference between theta and ripple
best_ripple_channel_this_column = candidate_index(index);
best_ripple_channel = sorted_config.Channel(best_ripple_channel_this_column);

figure
plot(normalised_ripple,sorted_config.Ks_ycoord)
hold on
plot(normalised_slow_wave,sorted_config.Ks_ycoord)

plot(normalised_theta,sorted_config.Ks_ycoord)
plot(normalised_gamma_power,sorted_config.Ks_ycoord)
plot([0.5 0.5],[0 4000],'--')
legend('Normalised ripple','Normalised theta','Normalised gamma')
xlabel('Normalised power')
ylabel('Depth (um)')


tic
MUA_filter_length = 41;
MUA_filter_alpha = 4;
time_step=1/new_SR; %match timestep to LFP time
time_bins_edges= tvec(1):time_step:max(tvec);
spike_times = [];
spike_ID = [];
CA1_spike_times = [];

for nchannel = best_ripple_channel-10:best_ripple_channel+10
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
(SUA(nchannel).time_bins_edges(1:end-1))+time_step/2;

%smooth SUA activity with a gaussian kernel
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
CA1_MUA_zscore=zscore(filtfilt(w,1,histcounts(spike_times,time_bins_edges)));
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


zscore_min = 0;
zscore_max = 3;

% tvec = start_sec:1/new_SR:start_sec+duration_sec;
LFP_tvec =  start_sec:1/new_SR:start_sec+duration_sec;
[replay,reactivations] = detect_candidate_events_masa(MUA_time,LFP_tvec,raw_LFP(best_ripple_channel,:),...
    CA1_MUA_zscore,CA1_spike_times,zscore_min,zscore_max,options) % Need to fix SUA spike
save extracted_candidate_events replay reactivations
% bz_FindPopBursts

% durations = [30 100]; inter-ripple 30 ms and max duration 100 ms 
% 'noise', 1 here is the top noisy channel (usually in dura gel)
[ripples] = FindRipples_masa(raw_LFP(best_ripple_channel,:)',tvec','minDuration',20,'durations',[30 100],'frequency',new_SR,'noise',raw_LFP(1,:)','passband',[125 300],'show','on')
save extracted_ripples_events ripples

[spindles] = FindSpindles_masa(raw_LFP(best_spindle_channel,:)',tvec','durations',[400 3000],'frequency',new_SR,'noise',raw_LFP(1,:)','passband',[9 17],'thresholds',1.25)
save extracted_spindles_events spindles
% bz_FindRipples

%% CSD peri-event
csd = bz_eventCSD(raw_LFP',ripples.onset',...
    'channels',sorted_config.Channel,'samplingRate',new_SR,'twin',[2 2])

csd = bz_eventCSD(raw_LFP',slow_waves.ints.UP(:,1),...
    'channels',sorted_config.Channel,'samplingRate',new_SR,'twin',[2 2])

csd = bz_eventCSD(raw_LFP',replay.onset,...
    'channels',sorted_config.Channel,'samplingRate',new_SR,'twin',[2 2])

csd = bz_eventCSD(raw_LFP',spindles.onset,...
    'channels',sorted_config.Channel,'samplingRate',new_SR,'twin',[3 3])

csd = bz_eventCSD(raw_LFP',MousePos.stimuli_onset(MousePos.stimuli_track == 1),...
    'channels',sorted_config.Channel,'samplingRate',new_SR,'twin',[1 1])

csd = bz_eventCSD(raw_LFP',MousePos.stimuli_onset(MousePos.stimuli_track == 2),...
    'channels',sorted_config.Channel,'samplingRate',new_SR,'twin',[1 1])

figure
plot(MUA_time(80000:100000),CA1_MUA_zscore(80000:100000))
hold on
xline(reactivations.onset,'r')
xlim([MUA_time(80000) MUA_time(100000)])
%% Detect Up and Down states

figure
hold on
plot(normalised_slow_wave,sorted_config.Ks_ycoord)
plot(normalised_spindle,sorted_config.Ks_ycoord)

plot(normalised_ripple,sorted_config.Ks_ycoord)
plot(normalised_theta,sorted_config.Ks_ycoord)
plot(normalised_gamma_power,sorted_config.Ks_ycoord)
plot([0.5 0.5],[0 4000],'--')
legend('Normalised slow wave','Normalised theta','Normalised ripple')
xlabel('Normalised power')
ylabel('Depth (um)')

% Select channels based on depth
channels_selected = ...
    chan_config.SpikeGLXchan0(find(chan_config.Ks_ycoord < chan_config.Ks_ycoord(best_SW_channel)...
    & chan_config.Ks_ycoord >chan_config.Ks_ycoord(best_SW_channel)-1000));

% Create Buz-style spike variables
spikes = [];
spikes.UID = [];
count = 1;
for nchannel = channels_selected(1):channels_selected(end)
    clusters_this_channel = cluster_id(find(peakChannel == nchannel))';
    good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(clusters_this_channel+1)=='good'))';% Plus one because clutser id is 0 based.

    if ~isempty(good_units_this_channel)
        for unit = 1:length(good_units_this_channel)
            if ~isempty(these_spike_times{good_units_this_channel(unit)+1})
                spikes.times{count} = these_spike_times{good_units_this_channel(unit)+1}; % Plus one because clutser id is 0 based.
                spikes.UID = [spikes.UID; good_units_this_channel(unit)];
                count = count + 1;
            end
        end
        
    end
end

slow_waves= DetectSlowWaves_masa('time',tvec,'lfp',raw_LFP(best_spindle_channel,:)','NREMInts',behavioural_state.SWS,'spikes',spikes);
% DetectSlowWaves
% fake_sws = [21 160; 190 477; 510 570];
% SWS;

save best_channels best_spindle_channel best_spindle_channel best_ripple_channel
save slow_waves slow_waves

%% Spike LFP relationship

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

% CCG of SLow Waves and spikes
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
event = MousePos.stimuli_onset(MousePos.stimuli_track == 1)
% event = reactivations.onset(reactivations.zscore>5);
% num_cell = 1;
% for types = 1:4
MUA_filter_length = 10;
SD_alpha = 2; %2 std width
MUA_filter_alpha = (10-1)/SD_alpha
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,4
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width


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
ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')


spike_data = CA3_spike_times(:,2);
num_cell = length(unique(CA3_spike_times(:,1)));
% num_cell = 1;
bin_size = 0.001;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data, event, [-0.3 1], bin_size);
psth = mean(binnedArray./num_cell./bin_size); % normalize to Hz
psth_se = std(binnedArray./num_cell./bin_size)./sqrt(length(event));

subplot(4,2,3)
plot(rasterX, rasterY,'b')
hold on
plot([0 0],get(gca,'ylim'),'r-')
ylabel('events');
xlabel('Spike time relative to event onset')
title('CA3')

subplot(4,2,4)
hold on
plot(bins,psth,'b')
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'g','FaceAlpha','0.3','LineStyle','none');
%     ylim([0 600])
hold on
plot([0 0],get(gca,'ylim'),'r-')
ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')


spike_data = DG_spike_times(:,2);
num_cell = length(unique(DG_spike_times(:,1)));
% num_cell = 1;
bin_size = 0.001;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data, event, [-0.3 1], bin_size);
psth = mean(binnedArray./num_cell./bin_size); % normalize to Hz
psth_se = std(binnedArray./num_cell./bin_size)./sqrt(length(event));

subplot(4,2,5)
plot(rasterX, rasterY,'k')
hold on
plot([0 0],get(gca,'ylim'),'r-')
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
psth = mean(binnedArray./num_cell./bin_size); % normalize to Hz
psth_se = std(binnedArray./num_cell./bin_size)./sqrt(length(event));
CA3_spike_times

subplot(4,2,7)
plot(rasterX, rasterY,'k')
hold on
plot([0 0],get(gca,'ylim'),'r-')
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
ylabel('Spike Rate (spk/s)');
xlabel('Spike time relative to event onset')
sgtitle('Ripple Onset')

sgtitle('CA1 bursting')
sgtitle('Track 2 Stimuli (Right)')
sgtitle('Track 1 Stimuli (Left)')
% end

%% CCG method (do't use)
CCGvec = [slow_waves.ints.UP(:,1),ones(size(slow_waves.ints.UP(:,1)));...
    spike_data,2.*ones(size(spike_data))];

[SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02,'norm','rate');
figure
subplot(1,2,1)
bar(t_CCG,SWspike_CCG(:,1,2)./length(slow_waves.ints.UP(:,1)),'facecolor',[0.5 0.5 0.5])
xlim([-1 1])
set(gca,'xticklabel',[]);
ylabel('Spike Rate');
xlabel('Relative to Up State')
title('CA1')
hold on
plot([0 0],get(gca,'ylim'),'r-')

spike_data = V1_spike_times(:,2);
CCGvec = [slow_waves.ints.UP(:,1),ones(size(slow_waves.ints.UP(:,1)));...
    spike_data,2.*ones(size(spike_data))];

[SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02,'norm','rate');
subplot(1,2,2)
bar(t_CCG,SWspike_CCG(:,1,2)./length(slow_waves.ints.UP(:,1)),'facecolor',[0.5 0.5 0.5])
xlim([-1 1])
set(gca,'xticklabel',[]);
ylabel('Spike Rate');
xlabel('Relative to Up State')
sgtitle('CA1 and V1 Populational spiking relative to Up State')
title('V1')
hold on
plot([0 0],get(gca,'ylim'),'r-')

