function [gamma_coherence gamma_phase_coherence] = gamma_coherence_analysis(raw_LFP,tvec,SR,best_channels,sorted_config)

% Gamma band
parameters = list_of_parameters;
filter_type  = 'bandpass';
filter_width = [30 100];                 % range of frequencies in Hz you want to filter between
filter_order = round(6*SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
norm_freq_range = filter_width/(SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_gamma = fir1(filter_order, norm_freq_range,filter_type);

lfp = [];
for nchannel = 1:size(sorted_config,1)
    tic
    lfp.data(:,nchannel) = filtfilt(b_gamma,1,raw_LFP(nchannel,:));
    toc
end
lfp.timestamps = tvec;
lfp.samplingRate = SR;
disp('Gamma band filtering finished!')
% [ csd ] = bz_CSD (lfp);


for n = 1:size(raw_LFP,1)
    %     s = RandStream('mrg32k3a','Seed',n);
    %     idx = datasample(s,1:length(lfp.data(n,:)),100000,'Replace',false);% Distribution of theta phase from each channel
    idx = length(lfp.data(:,n))-100000:length(lfp.data(:,n)); % Just to speed things up, use only part of the data

    tic
    for m = 1:size(raw_LFP,1)
        %         [coherence,f] = mscohere(raw_LFP(sorted_config.Channel(n),:),raw_LFP(sorted_config.Channel(m),:),[],[],[],new_SR);
        
        filt_hilb1 = hilbert(lfp.data(idx,n)'); %calculates the Hilbert transform of eeg1
        amp1 = abs(filt_hilb1);%calculates the instantaneous amplitude of eeg1 filtered between low_freq and high_freq
        amp1=amp1-mean(amp1); %removes mean of the signal because the DC component of a signal does not change the correlation
        filt_hilb2 = hilbert(lfp.data(idx,m)');%calculates the Hilbert transform of eeg2
        amp2 = abs(filt_hilb2);%calculates the instantaneous amplitude of eeg2 filtered between low_freq and high_freq
        amp2=amp2-mean(amp2);
        [crosscorr,lags]=xcorr(amp1, amp2,round(SR/10),'coeff'); %calculates crosscorrelations between amplitude vectors
        lags=(lags./SR)*1000; %converts lags to miliseconds
        g=find(crosscorr==max(crosscorr));%identifies index where the crosscorrelation peaks
        max_crosscorr_lag=lags(g);%identifies the lag at which the crosscorrelation peaks
        gamma_phase_coherence(n,m) = circ_mean(circ_dist(angle(filt_hilb1),angle(filt_hilb2)),[],2);
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
% gamma_coherence_CSD = zeros(size(sorted_config,1),size(sorted_config,1));
% gamma_phase_coherence_CSD = zeros(size(sorted_config,1),size(sorted_config,1));
% 
% for n = 1:size(sorted_config,1)-2
%     tic
%     for m = 1:size(sorted_config,1)-2
%         %          filt_hilb1 = hilbert(csd.data(:,n))'; %calculates the Hilbert transform of eeg1
%         %          filt_hilb2 = hilbert(csd.data(:,m))';%calculates the Hilbert transform of eeg2
% 
%         [coherogram,phase,t,f]  = bz_MTCoherogram(csd.data(:,n),csd.data(:,m),'window',5,'frequency',new_SR,'range',[30 100]);
%         gamma_coherence_CSD(n+1,m+1) = mean(mean(coherogram));
%         gamma_phase_coherence_CSD(n+1,m+1) = circ_mean(reshape(phase,size(phase,1)*size(phase,2),1));
%         %          gamma_phase_coherence_CSD(n+1,m+1) = unwrap(circ_mean(circ_dist(angle(filt_hilb1),angle(filt_hilb2)),[],2));
% 
%     end
%     toc
% end


%     Plot cross-coherence matrix
    figure
    subplot(1,2,1)
    imagesc(gamma_coherence)
    colorbar
    xticks(1:2:96)
    xticklabels(sorted_config.Ks_ycoord(1:2:96))
    yticks(1:2:96)
    yticklabels(sorted_config.Ks_ycoord(1:2:96))

    hold on
    plot([0 96],[find(sorted_config.Channel == best_channels.first_in_brain_channel) find(sorted_config.Channel == best_channels.first_in_brain_channel)],'--k','LineWidth',2)
    plot([0 96],[find(sorted_config.Channel == best_channels.L4_channel) find(sorted_config.Channel == best_channels.L4_channel)],'--b','LineWidth',2)
    plot([0 96],[find(sorted_config.Channel == best_channels.L5_channel) find(sorted_config.Channel == best_channels.L5_channel)],'--c','LineWidth',2)
    plot([0 96],[find(sorted_config.Channel == best_channels.CA1_channel) find(sorted_config.Channel == best_channels.CA1_channel)],'--r','LineWidth',2)

    plot([find(sorted_config.Channel == best_channels.first_in_brain_channel) find(sorted_config.Channel == best_channels.first_in_brain_channel)],[0 96],'--k','LineWidth',2)
    plot([find(sorted_config.Channel == best_channels.L4_channel) find(sorted_config.Channel == best_channels.L4_channel)],[0 96],'--b','LineWidth',2)
    plot([find(sorted_config.Channel == best_channels.L5_channel) find(sorted_config.Channel == best_channels.L5_channel)],[0 96],'--c','LineWidth',2)
    plot([find(sorted_config.Channel == best_channels.CA1_channel) find(sorted_config.Channel == best_channels.CA1_channel)],[0 96],'--r','LineWidth',2)

    title('LFP Gamma Coherence (30 - 100Hz)')

    subplot(1,2,2)
    imagesc(gamma_phase_coherence)
    colorbar

    xticks(1:2:96)
    xticklabels(sorted_config.Ks_ycoord(1:2:96))
    yticks(1:2:96)
    yticklabels(sorted_config.Ks_ycoord(1:2:96))
    hold on
   plot([0 96],[find(sorted_config.Channel == best_channels.first_in_brain_channel) find(sorted_config.Channel == best_channels.first_in_brain_channel)],'--k','LineWidth',2)
    plot([0 96],[find(sorted_config.Channel == best_channels.L4_channel) find(sorted_config.Channel == best_channels.L4_channel)],'--b','LineWidth',2)
    plot([0 96],[find(sorted_config.Channel == best_channels.L5_channel) find(sorted_config.Channel == best_channels.L5_channel)],'--c','LineWidth',2)
    plot([0 96],[find(sorted_config.Channel == best_channels.CA1_channel) find(sorted_config.Channel == best_channels.CA1_channel)],'--r','LineWidth',2)

    plot([find(sorted_config.Channel == best_channels.first_in_brain_channel) find(sorted_config.Channel == best_channels.first_in_brain_channel)],[0 96],'--k','LineWidth',2)
    plot([find(sorted_config.Channel == best_channels.L4_channel) find(sorted_config.Channel == best_channels.L4_channel)],[0 96],'--b','LineWidth',2)
    plot([find(sorted_config.Channel == best_channels.L5_channel) find(sorted_config.Channel == best_channels.L5_channel)],[0 96],'--c','LineWidth',2)
    plot([find(sorted_config.Channel == best_channels.CA1_channel) find(sorted_config.Channel == best_channels.CA1_channel)],[0 96],'--r','LineWidth',2)

    colorbar
    title('LFP Gamma Phase Coherence (30 - 100Hz)')
end

