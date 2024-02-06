function [PSD power best_channels] = calculate_channel_PSD(raw_LFP,SR,chan_config,options,varargin)

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

for nchannel = 1:size(chan_config,1)

    [pxx,fxx] = pwelch(raw_LFP(nchannel,:),win,[],nfft,SR);
    PSD(nchannel).channel = chan_config.Channel(nchannel);
    PSD(nchannel).shank = chan_config.Shank(nchannel);
    PSD(nchannel).xcoord = chan_config.Ks_xcoord(nchannel);
    PSD(nchannel).ycoord = chan_config.Ks_ycoord(nchannel);

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
for nchannel = 1:size(chan_config,1)
    power(nchannel,:) = PSD(nchannel).mean_power;
    xcoord(nchannel) = PSD(nchannel).xcoord;
    ycoord(nchannel) = PSD(nchannel).ycoord;
end
xcoord_avaliable = unique(xcoord);
% sort channel according to y coordinate
[ycoord idx] = sort(ycoord,'ascend');
xcoord = xcoord(idx);
power = power(idx,:);

if plot_option == 1


    stimulus_name = extractAfter(options.ANALYSIS_DATAPATH,[fullfile('analysis',options.SESSION),'\']);

    fig = figure;
    fig.Position = [260,50,1140,900]
    fig.Name = sprintf('%s %s %s PSD profile heatmap probe %i',options.SUBJECT,options.SESSION,stimulus_name,options.probe_id+1);

    freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','300-600 Hz'};
    % Loop over the frequency ranges
    for nfreq = 1:size(power, 2)
        % Create a subplot for this frequency range
        subplot(ceil(size(power, 2) / 2), 2, nfreq);

        % Create a scatter plot for this frequency range
        for col = 1:length(xcoord_avaliable)
            scatter(xcoord(xcoord == xcoord_avaliable(col)), ycoord(xcoord == xcoord_avaliable(col)), 24, power(xcoord == xcoord_avaliable(col),nfreq)/max(power(xcoord == xcoord_avaliable(col),nfreq)), 'filled'); hold on
        end
        xlim([0 1.25*max(xcoord)])
        ylim([0 1.25*max(ycoord)])
        % Add a colorbar
        colorbar;
        colormap(flip(gray))

        % Set the title
        title(['Frequency range ',freq_legends{nfreq}]);

        % Set the axis labels
        xlabel('X coordinate (micron)');
        ylabel('Y coordinate (micron)');
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    end
    sgtitle(sprintf('%s %s %s PSD profile heatmap probe %i',options.SUBJECT,options.SESSION,stimulus_name,options.probe_id+1),'Interpreter', 'none')

    fig = figure;
    fig.Position = [260,50,1140,900]
    fig.Name = sprintf('%s %s %s PSD profile probe %i',options.SUBJECT,options.SESSION,stimulus_name,options.probe_id+1);
    freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','300-600 Hz'};
    

    scaling_factor = ceil(max(xcoord_avaliable)/length(xcoord_avaliable)/10)*10;

    % Loop over the frequency ranges
    for nfreq = 1:size(power, 2)
        % Create a subplot for this frequency range
        subplot(ceil(size(power, 2) / 2), 2, nfreq);

        hold on
        Xticks = [];
        for col = 1:length(xcoord_avaliable)
            if col < length(xcoord_avaliable)
                plot(xcoord_avaliable(col)+(xcoord_avaliable(col+1)-xcoord_avaliable(col))/2+scaling_factor*power(xcoord == xcoord_avaliable(col),nfreq)/max(power(xcoord == xcoord_avaliable(col),nfreq)),ycoord(xcoord == xcoord_avaliable(col)))

                Xticks = [Xticks xcoord_avaliable(col)+(xcoord_avaliable(col+1)-xcoord_avaliable(col))/2+scaling_factor/2];
            else
                plot(1.1*xcoord_avaliable(col)+scaling_factor*power(xcoord == xcoord_avaliable(col),nfreq)/max(power(xcoord == xcoord_avaliable(col),nfreq)),ycoord(xcoord == xcoord_avaliable(col)))
                Xticks = [Xticks 1.1*xcoord_avaliable(col)+scaling_factor/2];
            end
        end

        xticks(Xticks)
        xticklabels(xcoord_avaliable)
        xlim([0 1.25*(scaling_factor+max(xcoord_avaliable(end)))])
        ylim([0 1.25*max(ycoord)])
        % Set the title
        title(['Frequency range ',freq_legends{nfreq}]);

        % Set the axis labels
        xlabel('X coordinate (micron)');
        ylabel('Y coordinate (micron)');
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    end
    sgtitle(sprintf('%s %s %s PSD profile probe %i',options.SUBJECT,options.SESSION,stimulus_name,options.probe_id+1),'Interpreter', 'none')
end

best_channels = [];


end
