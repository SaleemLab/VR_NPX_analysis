function PSTH = plot_perievent_spiketime_histogram(all_spike_data,events,varargin)

%INPUTS
%   all_spike_data      spike
%
%   (options)
%   'lfp'               -A buzcode-style lfp structure... if you would
%                        rather just input the lfp instead of loading from
%                        basepath
%                           Default: load from basePath with bz_GetLFP
%   'spikes'            -A buzcode-style spike structure
%                           Default: load from basePath with bz_GetSpikes
%   'NREMInts'          -Interval of times for NREM (seconds)
%                        (Default: loaded from SleepState.states.mat,
%                                   run SleepScoreMaster if not exist)
%                        use [0 Inf] to detect over all time points
%   'DetectionChannel'  -Channel with the most robust Slow Waves. (0-Indexing a la neuroscope).
%                        (Default: 'autoselect')
%                        'useold' to use channel from existing SlowWaves.events.mat
%                        If providing lfp via the 'lfp' input, make sure to
%                        give a detection channel here.
%   'noSpikes'          -true/false - set to true to not use spike information
%                        (default: false)
%   'MUAspikes'       -true/false - use MUA peaks (500-5000Hz) extracted
%                        from the .dat file instead of spikes
%   'CTXChans'          -LFP channels that are in the cortex...
%                        default: region 'CTX' from baseName.sessionInfo.mat or xml
%   'sensitivity'       -sensititivity (0-1) for determining LFP thresholds
%                        sensitivity for setting gamma/delta thresholds.
%                        lower sensitivity will result in fewer False Positives,
%                        but more Missed slow waves. (default 0.6)
%   'filterparms'       -filtering parameters structure with fields:
%           .deltafilter    [low high] bounds (default: [0.5 8]Hz)
%           .gammafilter    [low high] bounds (default: [100 400]Hz)
%           .gammasmoothwin  window for smoothing gamma power (default: 0.08 s)
%           .gammanormwin    window for normalizing gamma power (default: 20s)
%   'showFig'           -true/false show a quality control figure (default: true)
%   'saveMat'           -logical (default=true) to save in buzcode format
%   'forceReload'       -logical (default: false) to re-detect
%   'noPrompts'         -true/false disable any user prompts (default: false)
%
%
%OUTPUTS
%   SlowWaves    a buzcode structure
%   VerboseOut   extra output stuff for detection quality checks/figures
%
%
%
%DLevenstein 2016/2017
%If used, please cite: Levenstein et al 2018, currently on bioRxiv
%TO DO
%-incorporate multiple channels for detection of slow wave, which is robust
%on all (deep) lfp channels in the local cortical population
%-update input parameter list
%NOTE: requires 2017a or higher. sad. (functions: movingmad and movingmedian)
% Modified by Masahiro Takigawa 2023

% Default values
p = inputParser;
addParameter(p,'group','by region',@isstr) %
addParameter(p,'mode',1,@isnumeric) % Hanning window for pwelch analysis in seconds
addParameter(p,'group_name',{},@iscell) % Frequency range: [freq1(1) freq1;freq2 freq2;......]
addParameter(p,'event_name','event',@isstr) % Frequency range: [freq1(1) freq1;freq2 freq2;......]
addParameter(p,'twin',[-1 1],@isnumeric) % time window around the event for analysis
addParameter(p,'plot_option',1,@isnumeric) % Powers Selected frequency for plotting
addParameter(p,'sort_option','mean FR',@isstr) % Option for cell sorting (FR or latency or sorted or [])
addParameter(p,'sorted_cell_id',[],@iscell) % for sorted case, specify the id
addParameter(p,'bin_size',0.005,@isnumeric) % Powers Selected frequency for plotting
addParameter(p,'channel_map',[],@isnumeric) % Powers Selected frequency for plotting
addParameter(p,'smooth_option',1,@isnumeric) % Powers Selected frequency for plotting



% assign parameters (either defaults or given)
parse(p,varargin{:});
group_name = p.Results.group_name;
event_name = p.Results.event_name;
twin = p.Results.twin;
group = p.Results.group;
sort_option = p.Results.sort_option;
channel_map = p.Results.channel_map;
sorted_cell_id = p.Results.sorted_cell_id;
bin_size = p.Results.bin_size;
plot_option = p.Results.plot_option;
smooth_option = p.Results.smooth_option;

% MUA_filter_length = 30;
% SD_alpha = 3; %2 std width
% MUA_filter_alpha = (MUA_filter_length-1)/SD_alpha;
% w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,4
% w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
% Define Gaussian window for smoothing
w = gausswin(0.1*1/bin_size);

% Normalize to have an area of 1 (i.e., to be a probability distribution)
w = w / sum(w);

PSTH = [];
if plot_option == 1
    fig = figure()
    fig.Position = [661 431 1000 650];
end

if iscell(all_spike_data)
    switch group
        case 'by region'
            no_of_groups = length(all_spike_data);

            count = 1;
            for n = 1:size(all_spike_data,2)
                spike_data = all_spike_data{n};
                num_cell = length(unique(spike_data(:,1)));
                [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data(:,2), events, twin, bin_size);

                if smooth_option == 1
                    psth = filtfilt(w,1,mean(binnedArray./num_cell./bin_size)); % normalize to Hz
                    psth_se = filtfilt(w,1,std(binnedArray./num_cell./bin_size)./sqrt(length(events)));
                else
                    psth = mean(binnedArray./num_cell./bin_size); % normalize to Hz
                    psth_se = std(binnedArray./num_cell./bin_size)./sqrt(length(events));
                end

                if plot_option == 1
                    subplot(no_of_groups,2,count)
                    plot(rasterX, rasterY,'b')
                    hold on
                    plot([0 0],get(gca,'ylim'),'r-')
                    plot([0.1 0.1],get(gca,'ylim'),'b-')
                    % plot([0.5 0.5],get(gca,'ylim'),'r-')
                    ylabel('events');
                    xlabel('Spike time relative to event onset')
                    title(group_name{n})
                    count = count + 1;
                    subplot(no_of_groups,2,count)
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
                    count = count + 1;
                    sgtitle(event_name)
                    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
                end
            end
            PSTH = [];

        case 'by cell zscore'
            no_of_groups = length(all_spike_data);
            count = 1;
            for n = 1:size(all_spike_data,2)
                spike_data = all_spike_data{n};
                if isempty(spike_data)
                    PSTH(n).MUA_psth = [];
                    PSTH(n).cell_psth = [];
                    PSTH(n).spike_count = [];

                    PSTH(n).group = [];
                    PSTH(n).group_name = [];
                    PSTH(n).event_name = [];
                    PSTH(n).peak_latency = [];
                    PSTH(n).mean_FR = [];

                    continue
                else
                    cell_id = unique(spike_data(:,1));
                end

                switch sort_option
                    case 'sorted'
                        cell_id = sorted_cell_id{n};
                    otherwise
                        cell_id = unique(spike_data(:,1));
                end

                cell_psth= [];
                mean_FR = [];
                peak_latency = [];
                peak_FR = [];
                for cell = 1:length(cell_id)
                    spikes_this_cell = spike_data(find(spike_data(:,1) == cell_id(cell)),2);
                    if length(spikes_this_cell) < 1
                        mean_FR(cell) = 0;
                        peak_latency(cell) = 0;
                    else
                        mean_FR(cell) = length(spikes_this_cell)/(spikes_this_cell(end)-spikes_this_cell(1));
                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, events, twin, bin_size);

                        cell_spike_counts(cell,:,:) = binnedArray; % cell x trial x timebin



                        if size(binnedArray,1) ~= 1 % if multiple perievent spike trains
                            if smooth_option == 1
                                for nevent = 1:size(binnedArray,1)
                                    binnedArray(nevent,:) = conv(binnedArray(nevent,:)/bin_size,w,'same');
                                end

                                psth = mean(binnedArray./bin_size); % normalize to Hz
                            else
                                psth = mean(binnedArray./bin_size); % normalize to Hz
                            end
%                         else
%                             if smooth_option == 1
%                                 for nevent = 1:size(binnedArray,1)
%                                     binnedArray(nevent,:) = conv(binnedArray(nevent,:)/bin_size,w,'same');
%                                 end
% 
%                                 psth = mean(binnedArray./bin_size); % normalize to Hz
%                             else
%                                 psth = binnedArray./bin_size; % normalize to Hz
%                             end
                        end

                        [val,index]=min(abs(cumsum(psth)/sum(psth)-0.5));
                        [peakFR,~] = max(psth);
                        peak_latency(cell) = index;
                        peak_FR(cell) = peakFR;
                        cell_psth(cell,:) = psth;
                    end
                end

                if ~strcmp(sort_option,'sorted') & bin_size <= 0.01
                    mean_FR(peak_FR<0.5) = [];
                    peak_latency(peak_FR<0.5) = [];
                    cell_id(peak_FR<0.5) = [];
                    cell_psth(peak_FR<0.5,:) = [];
                    cell_spike_counts(peak_FR<0.5,:,:) = [];
                end

                switch sort_option
                    case 'mean FR'

                        [~, index ] = sort(mean_FR,'descend');
                        cell_sorted = cell_id(index);
                    case 'latency'

                        [~, index ] = sort(peak_latency,'ascend');
                        cell_sorted = cell_id(index);
                    otherwise
                        index = 1:length(cell_id);
                        cell_sorted = cell_id;
                end

                if smooth_option == 1

                    if size(cell_psth,1) == 1
                        psth = zscore(cell_psth,0,2); % normalize to Hz
                        psth_se = zeros(1,size(cell_psth,2));
                    else
                        psth = mean(zscore(cell_psth,0,2)); % normalize to Hz
                        psth_se = std(zscore(cell_psth,0,2))./sqrt(size(cell_psth,1));
                    end
                else
                    if size(cell_psth,1) == 1
                        psth = zscore(cell_psth,0,2); % normalize to Hz
                        psth_se = zeros(1,size(cell_psth,2));
                    else
                        psth = mean(zscore(cell_psth,0,2)); % normalize to Hz
                        psth_se = std(zscore(cell_psth,0,2))./sqrt(size(cell_psth,1));
                    end
                end

                % Save data
                PSTH(n).cell_id = cell_sorted;
                PSTH(n).sorted_option = sort_option;
                PSTH(n).twin = twin;
                PSTH(n).bin_size = bin_size;

                switch group
                    case 'by region'
                        PSTH(n).MUA_psth = psth;
                    otherwise
                        PSTH(n).MUA_psth = psth;
                        PSTH(n).cell_psth = cell_psth(index,:);
                        PSTH(n).spike_count = cell_spike_counts(index,:,:);
                end
                PSTH(n).group = group;
                PSTH(n).group_name = group_name{n};
                PSTH(n).event_name = event_name;
                PSTH(n).peak_latency = peak_latency(index);
                PSTH(n).mean_FR = mean_FR(index);

                if plot_option == 1
                    subplot(no_of_groups,2,count)
                    imagesc(bins,1:1:length(index),zscore(cell_psth(index,:),0,2))
                    %                 colormap jet
                    yticks(1:length(index))
                    yticklabels(cell_sorted)
                    % imagesc(bins,1:1:length(index),zscore(V1_UP_psth(index,:),0,2))
                    ylabel(sprintf('Cell (sorted by %s)',sort_option))
                    xlabel('Time relative to event onset (s)')
                    cbar = colorbar
                    cbar.Label.String = 'Firing Rate (zscore)';
                    hold on
                    plot([0 0],get(gca,'ylim'),'r','LineWidth',2)
                    clim([-1 1])
                    title(group_name{n})
                    count = count + 1;
%                     fontsize(fig, 14, "points")
                    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

                    subplot(no_of_groups,2,count)
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
                    ylabel('z-scored Spike Rate');
                    xlabel('Spike time relative to event onset')
                    count = count + 1;
                    sgtitle(event_name)
%                     fontsize(fig, 14, "points")
                    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
                end
            end
    end
elseif isstruct(all_spike_data) %
    switch group
        case 'by channel'
            %     no_of_channels = length(all_spike_data);
            count = 1;
            for n = 1:size(all_spike_data,2)
                if isempty(all_spike_data(n).spike_times)
                    channel_latency(n) = nan;
                elseif ~isempty(all_spike_data(n).spike_times)
                    spike_data = [all_spike_data(n).spike_id all_spike_data(n).spike_times];

                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data(:,2), events, twin, bin_size);
                    psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
                    channel_psth(n,:) = psth;

                    temp = find(binnedArray' ~= 0)/size(binnedArray,2);
                    [~,ia,~] = unique(floor(temp));
                    temp = temp(ia);
                    temp = (temp - floor(temp))*size(binnedArray,2);
                    channel_latency(n) = mean(temp); % Average first spike latency

                end
            end

            subplot(1,2,1)
            imagesc(bins,1:1:size(channel_psth,1),zscore(flip(channel_psth),0,2))
            yticks(1:10:size(channel_psth,1))
            yticklabels(flip(channel_map(1:10:end)))
            xlim([0 0.1])
            % imagesc(bins,1:1:length(index),zscore(V1_UP_psth(index,:),0,2))
            ylabel('Channel Depth')
            xlabel('Time relative to event onset (s)')
            cbar = colorbar
            cbar.Label.String = 'Firing Rate (zscore)';
            hold on
            %             plot([0 0],get(gca,'ylim'),'r','LineWidth',2)
            clim([0 2])
            title('zscored MUA PSTH by channel')

            subplot(1,2,2)
            hold on
            scatter(channel_latency,channel_map,'filled')
            ylim([min(channel_map) max(channel_map)])
            xlim([0 50])
            title('First spike latency by channel')
            ylabel('Channel Depth')
            xlabel('First spike latency by channel (ms)')

            sgtitle(event_name)
            %             fontsize(fig, 14, "points")
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
            
    end
else
    spike_data = all_spike_data;

end


end
