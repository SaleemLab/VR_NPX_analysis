function [PSD power best_channels] = calculate_channel_PSD(raw_LFP,SR,sorted_config,options,varargin)

% Default values
p = inputParser;
addParameter(p,'selected_channels',[1:size(raw_LFP,1)],@isnumeric) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn))
addParameter(p,'nfft_seconds',[2],@isnumeric) % Hanning window for pwelch analysis in seconds
addParameter(p,'frequency_range',[0.5 3;4 12;9 17;30 60;60 100;125 300; 300 600],@isnumeric) % Frequency range: [freq1(1) freq1;freq2 freq2;......]
addParameter(p,'selected_frequency',[1 2 3 4 5 6 7],@isnumeric) % Powers Selected frequency for plotting
addParameter(p,'plot_option',1,@isnumeric) % Powers Selected frequency for plotting


% assign parameters (either defaults or given)
parse(p,varargin{:});
selected_channels = p.Results.selected_channels;
nfft_seconds = p.Results.nfft_seconds;
F = p.Results.frequency_range;
plot_option =  p.Results.plot_option;
selected_frequency = p.Results.selected_frequency;

nfft = 2^(nextpow2(SR*nfft_seconds));
win  = hanning(nfft);
% F  = [0.5 3;4 12;9 17;30 60;60 100;125 300];

for nchannel = 1:size(sorted_config,1)

    [pxx,fxx] = pwelch(raw_LFP(nchannel,:),win,[],nfft,SR);
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


% Find peak channel for each oscillation band (currently hard-coded)
normalised_high_frequency = power(:,7)/max(power(:,7));  %power normalized by max power across channel
normalised_ripple = power(:,6)/max(power(:,6)); 
normalised_theta = power(:,2)/max(power(:,2)); 
normalised_slow_wave = power(:,1)./max(power(:,1));
normalised_spindle = power(:,3)./max(power(:,3));

% Find first peak that suddenly increases in power as the probe enters the
% brain
high_frequency_power_differnece = [0; diff(normalised_high_frequency)];
theta_power_differnece = [0; diff(normalised_theta)]; 
ripple_power_differnece = [0; diff(normalised_ripple)]; 
candidate_index = [];
[~,candidate_index{1}] = findpeaks(high_frequency_power_differnece,'MinPeakHeight',0.05);
[~,candidate_index{2}] = findpeaks(theta_power_differnece,'MinPeakHeight',0.05);
[~,candidate_index{3}] = findpeaks(ripple_power_differnece,'MinPeakHeight',0.05);

first_in_brain_channel = min([candidate_index{1}; candidate_index{2}; candidate_index{3}]);
best_channels.first_in_brain_channel = sorted_config.Channel(first_in_brain_channel);

% Find channels with normalised high-frequency power bigger than 0.1 (above 1500 micron below the brain)
candidate_index = [];
[value, candidate_index] = findpeaks(normalised_high_frequency(1:first_in_brain_channel+40),'MinPeakHeight',0.1); 
[~,index]  = max(normalised_high_frequency(candidate_index));  %find channel with maximum high-frequency-power
best_high_frequency_channel_this_column = candidate_index(index);
best_L5_channel = sorted_config.Channel(best_high_frequency_channel_this_column);

% Find CA1 channels below layer 5 and above theta-peaked channel
% Find channel with max theta power (putative dentate region)
[value, peak_theta_channel_this_column] = max(normalised_theta); 
channel_ranges = best_high_frequency_channel_this_column:peak_theta_channel_this_column;
[value, candidate_index] = findpeaks(normalised_ripple(best_high_frequency_channel_this_column:peak_theta_channel_this_column),'MinPeakHeight',0.1); 
[~,index]  = max(normalised_ripple(channel_ranges(candidate_index))-normalised_theta(channel_ranges(candidate_index)));  %find channel with maximum difference between theta and ripple power
best_ripple_channel_this_column = channel_ranges(candidate_index(index)); % First peak to pass that threshold is usually CA1
best_CA1_ripple_channel = sorted_config.Channel(best_ripple_channel_this_column);


 % Find channels with normalised high-frequency power bigger than 0.4 ( more than 120 micron above layer 5)
channel_ranges = first_in_brain_channel:best_high_frequency_channel_this_column-3;
[value, candidate_index] = findpeaks(normalised_high_frequency(channel_ranges),'MinPeakHeight',0.1); % should be at least 120 micron (Distance between each channel in one column is 40 micron)
[~,index]  = max(normalised_high_frequency(channel_ranges(candidate_index))-normalised_high_frequency(channel_ranges(candidate_index)-1));  %find channel with biggest high-frequency-power gain (compared to the channel above)
best_high_frequency_channel_this_column = channel_ranges(candidate_index(index));
best_L4_channel = sorted_config.Channel(best_high_frequency_channel_this_column);


best_channels.L4_channel = best_L4_channel;
best_channels.L5_channel = best_L5_channel;
best_channels.CA1_channel = best_CA1_ripple_channel;


if plot_option == 1
    % Quick plotting of PSD
%     colour_line= {'k','r','m','b','c','g','y'};
    colour_line= {[177,0,38]/256,[213,62,79]/256,[252,141,89]/256,...
        [153,213,148]/256,[26,152,80]/256,[66,146,198]/256,[8,69,148]/256};
    freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','300-600 Hz'};
    selected_frequency = [2 6 7];
    figure
%     subplot(1,4,2)

    for n = selected_frequency
        pl(n) = plot(power(:,n)./max(power(:,n)),sorted_config.Ks_ycoord','Color',colour_line{n})

        %     p(n) = plot(power_differnece,1:96,colour_line{n})
        hold on
    end
    hold on
    plot([0 1],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)

    if ~isempty(best_L4_channel)
        plot([0 1],[sorted_config.Ks_ycoord(find(sorted_config.Channel == best_L4_channel)) sorted_config.Ks_ycoord(find(sorted_config.Channel == best_L4_channel))],'--b','LineWidth',2)
    end

    if ~isempty(best_L5_channel)
        plot([0 1],[sorted_config.Ks_ycoord(find(sorted_config.Channel == best_L5_channel)) sorted_config.Ks_ycoord(find(sorted_config.Channel == best_L5_channel))],'--c','LineWidth',2)
    end

    if ~isempty(best_CA1_ripple_channel)
        plot([0 1],[sorted_config.Ks_ycoord(find(sorted_config.Channel == best_CA1_ripple_channel)) sorted_config.Ks_ycoord(find(sorted_config.Channel == best_CA1_ripple_channel))],'--r','LineWidth',2)
    end

    % legend('1-3','4-12','30-60','60-100','125-300')
    legend([pl(selected_frequency)],{freq_legends{selected_frequency}});
    ylim([0 4000])
end


end
