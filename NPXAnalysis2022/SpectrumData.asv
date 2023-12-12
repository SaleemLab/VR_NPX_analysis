% TempSpectralAnalysis_NPX %
params_spectrogram_full.tapers = [4, 10];
params_spectrogram_full.Fs = 1000;
params_spectrogram_full.fpass = [0, 150];
params_spectrogram_full.err = 0;
params_spectrogram_full.trialave = 1;

movingwin_full = [0.5, 0.01];

postStimTime = ttimeVector>=1;
PostStimResps=resps(:,postStimTime,:);

for i = 1:size(PostStimResps,1)
    LFPTemp = squeeze(PostStimResps(i,:,:));
    [S{i,1},t{i,1}, f{i,1}]=mtspecgramc(LFPTemp,movingwin_full,params_spectrogram_full);
    clear LFPTemp
end
 %alternatively
 
 for i = 1:size(PostStimResps,1)
     LFPTemp = squeeze(PostStimResps(1,:,:));
     FFT = fft(LFPTemp');
     fs=1000;
     f = 0:(fs/750):750-1*fs/750;     
     clear LFPTemp
 end