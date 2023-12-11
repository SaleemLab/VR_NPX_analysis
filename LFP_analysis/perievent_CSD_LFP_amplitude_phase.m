
function [ csd, lfpAvg ] = perievent_CSD_LFP_amplitude_phase (lfp, timestamps,all_events, varargin)

% Based on [ CSD ] = bz_eventCSD (lfp, events, varargin)
% Calculates event-triggered (i.e. SWRs) CSD map from a linear array of LFPs
% lfp usually = timevec X nchannel 
% If using event-evoked lfp, use tvec x nchannel x nevent
% twin [0.1 0.1] meaning timewindow that includes  0.1s before and 0.1s
% after the onset of the event.

%% Parse inputs

p = inputParser;
addParameter(p,'channels',[1 384],@isvector);
addParameter(p,'samplingRate',1250,@isnumeric);
addParameter(p,'twin',[0.1 0.1],@isnumeric);
addParameter(p,'spat_sm',11,@isnumeric);
addParameter(p,'temp_sm',11,@isnumeric);
addParameter(p,'doDetrend',false,@islogical);
addParameter(p,'plotCSD',false,@islogical);
addParameter(p,'plotLFP',false,@islogical);
addParameter(p,'filter',[],@isnumeric);
addParameter(p,'channel_boundary',[],@isnumeric);


parse(p,varargin{:});
channels = p.Results.channels;
samplingRate = p.Results.samplingRate;
spat_sm = p.Results.spat_sm;
temp_sm = p.Results.temp_sm;
doDetrend = p.Results.doDetrend;
plotCSD = p.Results.plotCSD;
plotLFP = p.Results.plotLFP;
filter = p.Results.filter;
channel_boundary = p.Results.channel_boundary;

%lfp input
if isstruct(lfp)
    data = lfp.data;
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
elseif iscell(lfp) %for multiple trials
    celllengths = cellfun(@length,lfp);
    data = vertcat(lfp{:});
elseif isnumeric(lfp)
    tic
    data = lfp;

    if ~isempty(filter)
        filter_type  = 'bandpass';
        filter_width = filter;                 % range of frequencies in Hz you want to filter between
        filter_order = round(6*samplingRate/(max(filter_width)-min(filter_width)));  % creates filter for ripple
        norm_freq_range = filter_width/(samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
        b_filter = fir1(filter_order, norm_freq_range,filter_type);

        for nchannel = 1:size(channels,2)
            data(:,nchannel) = filtfilt(b_filter,1,data(:,nchannel));
        end
    end
    %     timestamps = [1:length(lfp)]'./samplingRate;
    toc
end

twin = p.Results.twin*samplingRate;

%% Conpute event-triggered LFP average

if ~iscell(all_events) % just one event type

    if ~isempty(all_events) % If all_events is just timestamp of one event type
        events = all_events;
        events = round((events-timestamps(1))*samplingRate);
        events = events((events + twin(2) <= size(data,1)) & (events - twin(1) > 0));
        lfp_temp = nan(twin(1)+twin(2)+1,length(channels),length(events));
        lfp_raw = nan(twin(1)+twin(2)+1,length(channels),length(events));
        lfp_power = nan(twin(1)+twin(2)+1,length(channels),length(events));
        lfp_phase = nan(twin(1)+twin(2)+1,length(channels),length(events));

        for e = 1:length(events)
            lfp_temp(:,:,e) = data(events(e)-twin(1):events(e)+twin(2),channels);
            lfp_raw(:,:,e) = lfp(events(e)-twin(1):events(e)+twin(2),channels);
            lfp_power(:,:,e) = zscore(abs(hilbert(lfp_temp(:,:,e))),0,1);
            lfp_phase(:,:,e) = angle(hilbert(lfp_temp(:,:,e)));
        end
    else % If already in trial format
        lfp_temp = nan(size(lfp,1),size(lfp,2),size(lfp,3));
        lfp_raw = nan(size(lfp,1),size(lfp,2),size(lfp,3));
        lfp_power = nan(size(lfp,1),size(lfp,2),size(lfp,3));
        lfp_phase = nan(size(lfp,1),size(lfp,2),size(lfp,3));

        for e = 1:size(lfp,3)
            lfp_temp(:,:,e) = data(:,:,e);
            lfp_raw(:,:,e) = lfp(:,:,e);
            lfp_power(:,:,e) = zscore(abs(hilbert(lfp_temp(:,:,e))),0,1);
            lfp_phase(:,:,e) = angle(hilbert(lfp_temp(:,:,e)));
        end
    end

    % lfp_avg = nanmean(lfp_temp,3);
    lfp_avg = nanmean(lfp_temp,3)*-1;
    lfp_raw_avg = nanmean(lfp_raw,3);
    lfp_power_avg = nanmean(lfp_power,3);
    lfp_phase_avg = circ_mean(lfp_phase,[],3);
    %
    % figure
    % for i = 1:length(events)
    %     for n = 1:30
    %         if n <12
    %             plot(lfp_temp(:,n,i)*10000 + 40*30 - n*40,'r')
    %         else
    %             plot(lfp_temp(:,n,i)*10000 + 40*30 - n*40,'k')
    %         end
    %         hold on
    %     end
    % end
    %% Conpute CSD

    % detrend
    if doDetrend
        lfp_avg = detrend(lfp_avg')';
    end

    % temporal smoothing
    if temp_sm > 0
        for ch = 1:size(lfp_avg,2)
            lfp_avg(:,ch) = smooth(lfp_avg(:,ch),temp_sm,'sgolay');
        end
    end

    % spatial smoothing
    if spat_sm > 0
        for t = 1:size(lfp_avg,1)
            lfp_avg(t,:) = smooth(lfp_avg(t,:),spat_sm,'lowess');
        end
    end

    % calculate CSD
    CSD = diff(lfp_avg,2,2);

    % generate output structure
    csd.data = CSD;
    if ~isempty(all_events)
        csd.timestamps = (-(twin(1)/samplingRate):(1/samplingRate):(twin(2)/samplingRate))*1e3;
    else
        csd.timestamps = timestamps;
    end
    csd.samplingRate = samplingRate;
    csd.channels = channels;
    csd.params.spat_sm = spat_sm;
    csd.params.temp_sm = temp_sm;
    csd.params.detrend = doDetrend;


    lfpAvg.raw = lfp_raw_avg;
    lfpAvg.filtered = lfp_avg;
    lfpAvg.power = lfp_power_avg;
    lfpAvg.phase = lfp_phase_avg;
    lfpAvg.filter_range = filter;
    if ~isempty(all_events)
        lfpAvg.timestamps = (-(twin(1)/samplingRate):(1/samplingRate):(twin(2)/samplingRate))*1e3;
    else
        lfpAvg.timestamps = timestamps;
    end
    lfpAvg.samplingRate = samplingRate;
    lfpAvg.channels = channels;
    lfpAvg.params.spat_sm = spat_sm;
    lfpAvg.params.temp_sm = temp_sm;
    lfpAvg.params.detrend = doDetrend;
else
    for n = 1:length(all_events)
        events = all_events{n};
        events = round((events-timestamps(1))*samplingRate);
        events = events((events + twin(2) <= size(data,1)) & (events - twin(1) > 0));
        lfp_temp = nan(twin(1)+twin(2)+1,length(channels),length(events));
        lfp_raw = nan(twin(1)+twin(2)+1,length(channels),length(events));
        lfp_power = nan(twin(1)+twin(2)+1,length(channels),length(events));
        lfp_phase = nan(twin(1)+twin(2)+1,length(channels),length(events));

        for e = 1:length(events)
            lfp_temp(:,:,e) = data(events(e)-twin(1):events(e)+twin(2),channels);
            lfp_raw(:,:,e) = lfp(events(e)-twin(1):events(e)+twin(2),channels);
            lfp_power(:,:,e) = zscore(abs(hilbert(lfp_temp(:,:,e))),0,1);
            lfp_phase(:,:,e) = angle(hilbert(lfp_temp(:,:,e)));
        end

        % lfp_avg = nanmean(lfp_temp,3);
        lfp_avg = nanmean(lfp_temp,3)*-1;
        lfp_raw_avg = nanmean(lfp_raw,3);
        lfp_power_avg = nanmean(lfp_power,3);
        lfp_phase_avg = circ_mean(lfp_phase,[],3);

        %
        % figure
        % for i = 1:length(events)
        %     for n = 1:30
        %         if n <12
        %             plot(lfp_temp(:,n,i)*10000 + 40*30 - n*40,'r')
        %         else
        %             plot(lfp_temp(:,n,i)*10000 + 40*30 - n*40,'k')
        %         end
        %         hold on
        %     end
        % end
        %% Conpute CSD

        % detrend
        if doDetrend
            lfp_avg = detrend(lfp_avg')';
        end

        % temporal smoothing
        if temp_sm > 0
            for ch = 1:size(lfp_avg,2)
                lfp_avg(:,ch) = smooth(lfp_avg(:,ch),temp_sm,'sgolay');
            end
        end

        % spatial smoothing
        if spat_sm > 0
            for t = 1:size(lfp_avg,1)
                lfp_avg(t,:) = smooth(lfp_avg(t,:),spat_sm,'lowess');
            end
        end

        % calculate CSD
        CSD = diff(lfp_avg,2,2);
        
        % generate output structure
        csd(n).data = CSD;
        csd(n).timestamps = (-(twin(1)/samplingRate):(1/samplingRate):(twin(2)/samplingRate))*1e3;
        csd(n).samplingRate = samplingRate;
        csd(n).channels = channels;
        csd(n).params.spat_sm = spat_sm;
        csd(n).params.temp_sm = temp_sm;
        csd(n).params.detrend = doDetrend;


        lfpAvg(n).raw = lfp_raw_avg;
        lfpAvg(n).filtered = lfp_avg;
        lfpAvg(n).power = lfp_power_avg;
        lfpAvg(n).phase = lfp_phase_avg;
        lfpAvg(n).filter_range = filter;
        lfpAvg(n).timestamps = (-(twin(1)/samplingRate):(1/samplingRate):(twin(2)/samplingRate))*1e3;
        lfpAvg(n).samplingRate = samplingRate;
        lfpAvg(n).channels = channels;
        lfpAvg(n).params.spat_sm = spat_sm;
        lfpAvg(n).params.temp_sm = temp_sm;
        lfpAvg(n).params.detrend = doDetrend;
    end
end


%% Plot

if plotLFP

    taxis = lfpAvg.timestamps;
    cmax = max(max(CSD)); 
    figure;
    subplot(1,2,1);
    contourf(taxis,1:size(CSD,2),CSD',40,'LineColor','none');hold on;
    colormap jet; caxis([-cmax cmax]);
    set(gca,'YDir','reverse');xlabel('time (s)');ylabel('channel');title('CSD'); 
    plot([0 0],[1 size(CSD,2)],'--k');hold on;
    
    subplot(1,2,2);
    for ch=1:size(lfp_avg,2)
        offset = 400*(ch-1);
        sh_tmp = 1e0*(10000000*lfp_avg(:,ch)) + offset;
        plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
        clear sh_tmp
    end
    set(gca,'YDir','reverse','YTickLabel',[]);ylim([-1000 offset+1000]);xlim([taxis(1) taxis(end)]);
    xlabel('time (ms)');ylabel('channel');title('LFP');   
    plot([0 0],ylim,'--r');hold on;

       
elseif plotCSD  
    
     cmax = max(max(CSD)); 
   
     figure;
     contourf(taxis,1:size(CSD,2),CSD',40,'LineColor','none');hold on;
     colormap jet; caxis([-cmax cmax]);
     set(gca,'YDir','reverse','YTickLabel',[]);ylim([-1000 offset+1000]);xlim([taxis(1) taxis(end)]);
     plot([0 0],[1 size(CSD,2)],'--k');hold on;
   
end

end
