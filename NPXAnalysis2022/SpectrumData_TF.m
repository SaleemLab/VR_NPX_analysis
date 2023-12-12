% function SpectrumData_TF(tstimData,ttimeVector,resps)
TrialInfo = tstimData.StimIndex;
Trial_1 = TrialInfo==1;
Trial_2 = TrialInfo==2;
Trial_4 = TrialInfo==4;
Trial_8 = TrialInfo==8;

% TempSpectralAnalysis_NPX %
params_spectrogram_full.tapers = [4, 10];
params_spectrogram_full.Fs = 1000;
params_spectrogram_full.fpass = [0, 150];
params_spectrogram_full.err = 0;
params_spectrogram_full.trialave = 1;

movingwin_full = [0.5, 0.01];

postStimTime = ttimeVector>=3;
PostStimResps=resps(:,postStimTime,:);

StimRespsFilt_1 = bandpassSig(:,~postStimTime,Trial_1);
StimRespsFilt_2 = bandpassSig(:,~postStimTime,Trial_2);
StimRespsFilt_4 = bandpassSig(:,~postStimTime,Trial_4);
StimRespsFilt_8 = bandpassSig(:,~postStimTime,Trial_8);

PostStimRespsFilt_1 = bandpassSig(:,postStimTime,Trial_1);
PostStimRespsFilt_2 = bandpassSig(:,postStimTime,Trial_2);
PostStimRespsFilt_4 = bandpassSig(:,postStimTime,Trial_4);
PostStimRespsFilt_8 = bandpassSig(:,postStimTime,Trial_8);

[ChannelMapData] = SGLXMetaToCoords_ChannelMap_FR('X:\ibn-vision\DATA\SUBJECTS\M22008\ephys\20220408\M22008_20220408_g0\M22008_20220408_g0_imec0');
Thalamus = ChannelMapData.ycoords<1000;
Cortex = ChannelMapData.ycoords>1000;

[~,ShankOrderVT] = sort(ChannelMapData.elecInd(Thalamus),'ascend');
[~,ShankOrderV1] = sort(ChannelMapData.elecInd(Cortex),'ascend');

% for i = 1:size(PostStimResps,1)
%     LFPTemp = squeeze(PostStimResps(i,:,:));
%     [S{i,1},t{i,1}, f{i,1}]=mtspecgramc(LFPTemp,movingwin_full,params_spectrogram_full);
%     clear LFPTemp
% end
%alternatively

for i = 1:size(StimResps_8,1)
    LFPTemp = squeeze(StimResps_2(i,:,:));
    FFT = fft(mean(LFPTemp,2))/length(mean(LFPTemp,2));
    fs=1000;
    amp1 = abs(FFT2);
    hz1 = linspace( 0, fs/2, (length(mean(LFPTemp,2))/2)+1);
    posfreqs = 2:floor(length(mean(LFPTemp,2))/2);
    amp1(posfreqs) = amp1(posfreqs)*2;    
    power1 = amp1.^2;
    figure;     
    subplot(311)
    t=ttimeVector(~postStimTime);
    plot(t,mean(LFPTemp,2))
    
    subplot(312,'next','add')
    plot(hz1, amp1(1:length(hz1)))
    
    subplot(313,'next','add')
    plot(hz1, power1(1:length(hz1)))


    clear LFPTemp
end
% end