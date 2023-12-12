function [clusters chan_config sorted_config] = load_KS_NPX1(options,column,varargin)

[file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column);

% Default values
p = inputParser;
addParameter(p,'selected_channels',[1:size(chan_config,1)],@isnumeric) % Select channels for analysis (default is all the channles)
addParameter(p,'group','by channel',@isstr) % spike data grouped 'by channel' or 'by region'
addParameter(p,'plot_option',1,@isnumeric) % Powers Selected frequency for plotting
addParameter(p,'SR',1250,@isnumeric) % Frequency for MUA time course
addParameter(p,'LFP_tvec',[],@isnumeric) % LFP time vector
addParameter(p,'cell_exporer','off',@isstr) % LFP time vector

% addParameter(p,'stdev',[],@isnumeric)
% addParameter(p,'show','on',@isstr)
% addParameter(p,'noise',[],@ismatrix)
% addParameter(p,'passband',[125 300],@isnumeric)
% addParameter(p,'EMGThresh',.9,@isnumeric);
% addParameter(p,'saveMat',false,@islogical);
% addParameter(p,'minDuration',20,@isnumeric)
% addParameter(p,'plotType',1,@isnumeric)
% addParameter(p,'savepath',[],@isstr)

% assign parameters (either defaults or given)
parse(p,varargin{:});
selected_channels = p.Results.selected_channels;
group = p.Results.group;
LFP_tvec = p.Results.LFP_tvec;
SR = p.Results.SR;
cell_exporer = p.Results.cell_exporer;

% Load cluster table
load([options.KS_DATAPATH,'\cluster_table.mat'])

% Load cell explorer generated cell metrics
switch cell_exporer
    case 'on'

        if exist([options.KS_DATAPATH,'\',options.SUBJECT,'_',options.SESSION,'.cell_metrics.cellinfo.mat']) == 2
            load([options.KS_DATAPATH,'\',options.SUBJECT,'_',options.SESSION,'.cell_metrics.cellinfo.mat'])
            cell_metrics.KS_cluster_id = double(cluster_table.cluster_id)';
            cell_metrics.good_unit = double(cluster_table.cluster_id)';
        else
            cell_metrics = [];
            disp('cell metrics not found!')
        end
        
    case 'off'
        cell_metrics = [];
        disp('cell metrics not used!')
end

KS_metrics = readtable(fullfile(options.KS_DATAPATH,'metrics.csv'));

% Load spike time data
[these_spike_times,nominal_KSLabel,cluster_id,peakChannel,maxSpkTime] = import_ks_spiketimes(options,options.gFileNum,options.KS_CATGT_FNAME,imecMeta.imSampRate);
mean_waveforms = readNPY([options.KS_DATAPATH,'/mean_waveforms.npy']); % Information about mean waveform
cluster_id = cluster_id + 1; % Cluster ID transformed from 0-based to 1-based
for id = 1:length(cluster_id)
    peak_channel_waveforms(id,:) = squeeze(mean_waveforms(cluster_id(id),peakChannel(id),:)); % Here ID is post-kilosort id
end

clusters = [];

if ~isempty(LFP_tvec)
    %smooth SUA activity with a gaussian kernel
    MUA_filter_length = 41;
    MUA_filter_alpha = 4;
    time_step=1/SR; %match timestep to LFP time
    w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,2
    w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
    time_bins_edges= LFP_tvec(1):time_step:max(LFP_tvec);
end

tic
switch group
    case 'by channel'
        for nchannel = 1:size(chan_config,1)
            clusters_this_channel = find(peakChannel == nchannel); % ID for all clusters
            [~,index,]= intersect(cluster_id,clusters_this_channel);
            good_units_index = index(find(nominal_KSLabel(index)=='good'));
%             good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(index)=='good')); % ID for all clusters
            clusters(nchannel).spike_times= [];
            clusters(nchannel).spike_id = [];
            clusters(nchannel).cell_type  = [];
            clusters(nchannel).waveform  = [];

            if ~isempty(index) % If any clusters
                if ~isempty(good_units_index) % If any SUA
                    clusters(nchannel).spike_times= [];
                    clusters(nchannel).spike_id = [];
                    clusters(nchannel).cell_type  = [];
                    for unit = 1:length(good_units_index)

                        clusters(nchannel).spike_times = [clusters(nchannel).spike_times; these_spike_times{good_units_index(unit)}];
                        clusters(nchannel).spike_id = [clusters(nchannel).spike_id; good_units_index(unit)*ones(length(these_spike_times{good_units_index(unit)}),1)];
                        if ~isempty(cell_metrics)
                            switch char(cell_metrics.putativeCellType(good_units_index(unit) == cell_metrics.KS_cluster_id+1)) % +1 beause KS_cluster_id is based on 0
                                case 'Narrow Interneuron'
                                    clusters(nchannel).cell_type(unit) = 1; % Putative cell types
                                case 'Wide Interneuron'
                                    clusters(nchannel).cell_type(unit) = 2; % Putative cell types
                                case 'Pyramidal Cell'
                                    clusters(nchannel).cell_type(unit) = 3; % Putative cell types
                            end
                        end
                        clusters(nchannel).waveform(unit,:) = peak_channel_waveforms(good_units_index(unit),:);
                    end

                    if ~isempty(LFP_tvec)
                        MUA_zscore = [];
                        MUA_zscore = zscore(filtfilt(w,1,histcounts(clusters(nchannel).spike_times,time_bins_edges)));
                        clusters(nchannel).MUA_zscore= interp1(time_bins_edges(1:end-1)+time_step/2,MUA_zscore,LFP_tvec,'linear'); %MUA to match LFP time
                        clusters(nchannel).MUA_tvec = LFP_tvec;
                    end
                end
            end
        end
        toc
    case 'by region'

        count = 1;
        unit_id = [];
        unit_type = [];
        SUA_spike_time = [];
        SUA_spike_id = [];
        SUA_peak_channel_waveforms = [];
        SUA_peak_channel = [];
        SUA_cell_type = [];
        id_conversion = [];
        
        for nchannel = selected_channels(1):selected_channels(2)
            clusters_this_channel = find(peakChannel == nchannel); % ID for all clusters
            [~,index,]= intersect(cluster_id,clusters_this_channel);
            good_units_index = index(find(nominal_KSLabel(index)=='good'));
            %             good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(index)=='good')); % ID for all clusters

            if ~isempty(index) % If any clusters
                if ~isempty(good_units_index) % If any SUA
                    for unit = 1:length(good_units_index)

                        SUA_spike_time = [SUA_spike_time; these_spike_times{good_units_index(unit)}];
                        SUA_spike_id = [SUA_spike_id; cluster_id(good_units_index(unit))*ones(length(these_spike_times{good_units_index(unit)}),1)];
                        SUA_peak_channel(count) = nchannel;
                        id_conversion(count,:) = [count good_units_index(unit)];
                        SUA_peak_channel_waveforms(count,:) = peak_channel_waveforms(good_units_index(unit),:);
                        if ~isempty(cell_metrics)
                            
                            switch char(cell_metrics.putativeCellType(good_units_index(unit) == cell_metrics.KS_cluster_id+1)); % +1 beause KS_cluster_id is based on 0
                                case 'Narrow Interneuron'
                                    SUA_cell_type = [SUA_cell_type 1]; % Putative cell types
                                case 'Wide Interneuron'
                                    SUA_cell_type = [SUA_cell_type 2]; % Putative cell types
                                case 'Pyramidal Cell'
                                    SUA_cell_type = [SUA_cell_type 3]; % Putative cell types
                                otherwise
                                     SUA_cell_type = [SUA_cell_type 0]; % Putative cell types
                            end
                        end
                        count = count + 1;
                    end
                end
            end
        end

        %smooth SUA activity with a gaussian kernel
        if ~isempty(LFP_tvec)
            MUA_zscore=zscore(filtfilt(w,1,histcounts(SUA_spike_time,time_bins_edges)));
            MUA_zscore= interp1(time_bins_edges(1:end-1)+time_step/2,MUA_zscore,LFP_tvec,'linear'); %MUA to match LFP time
            MUA_tvec = LFP_tvec;
        else
            MUA_tvec = [];
            MUA_zscore = [];
        end
        toc
% 
%         for unit = 1:length(nominal_KSLabel)
%             if peakChannel(unit) > selected_channels(1) & peakChannel(unit) < selected_channels(2) % If units from certain channel range
% 
%                 if nominal_KSLabel(unit) == 'good'
%                     SUA_peak_channel(count) = peakChannel(unit);
%                     SUA_peak_channel_waveforms(count,:) = peak_channel_waveforms(unit,:);
%                     SUA_spike_time = [SUA_spike_time; these_spike_times{unit}];
%                     SUA_spike_id = [SUA_spike_id; count*ones(length(these_spike_times{unit}),1)];
%                     id_conversion(count,:) = [count cluster_id(unit)];
% 
%                     if ~isempty(cell_metrics)
%                         switch char(cell_metrics.putativeCellType(unit))
%                             case 'Narrow Interneuron'
%                                 cell_type = [cell_type 1]; % Putative cell types
%                             case 'Wide Interneuron'
%                                 cell_type = [cell_type 2]; % Putative cell types
%                             case 'Pyramidal Cell'
%                                 cell_type = [cell_type 3]; % Putative cell types
%                         end
%                     end
% 
%                     count = count + 1;
%                 end
%             end
%         end

        % Sort spike time from start to end
        [~,index] = sort(SUA_spike_time);
        clusters.spike_id = SUA_spike_id(index);
        clusters.spike_times = SUA_spike_time(index);
        clusters.peak_channel = SUA_peak_channel;
        clusters.peak_channel_waveforms = SUA_peak_channel_waveforms;
        clusters.cell_type = SUA_cell_type;
        clusters.id_conversion = id_conversion;
        clusters.MUA_zscore = MUA_zscore;
        clusters.MUA_tvec = MUA_tvec;
    case 'Buz style'
        % Create Buz-style spike variables
        spikes = [];
        spikes.UID = [];
        spikes.times = [];
        count = 1;

        % channels_selected =...
        % chan_config.SpikeGLXchan0(find(chan_config.Ks_ycoord < chan_config.Ks_ycoord(best_channels.SW_channel)...
        %     & chan_config.Ks_ycoord >chan_config.Ks_ycoord(best_channels.SW_channel)-1000));

        for nchannel = selected_channels(1):selected_channels(2)
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

        clusters = spikes;

end

end