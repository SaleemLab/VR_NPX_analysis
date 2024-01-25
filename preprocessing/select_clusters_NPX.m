function [clusters chan_config sorted_config] = select_clusters_NPX(options,varargin)

%%% Function to extract clusters based on different sorters
% Was modified/simplified from old function: 
% [clusters chan_config sorted_config] = load_KS_NPX1(options,column,varargin)
% Modified from original version of load_KS_NPX1 function


% Extract NPX channel configuration
options.importMode = 'KS';
[file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,1);

% Default values
p = inputParser;
addParameter(p,'selected_channels',[1 size(chan_config,1)],@isnumeric) % Select channels for analysis (default is all the channles)
addParameter(p,'group','good clusters',@isstr) % spike data grouped 'good clusters'
addParameter(p,'plot_option',0,@isnumeric) % option to plot
addParameter(p,'SR',60,@isnumeric) % SR for cluster spike count time course (normally 60Hz)
addParameter(p,'tvec',[],@isnumeric) % Time vector

% addParameter(p,'cell_exporer','off',@isstr) % cell explorer option
% removed (replaced by few lines of code that can classify the cell based on )
addParameter(p,'waveform','off',@isstr) % extract individual spike waveform for each cluster (not mean waveform)
addParameter(p,'sorter','off',@isstr) % What sorters to use. If off, it is the default kilosort 3.
% If 'KS3' or 'KS2', it is spike interface KS3 and KS2.


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
tvec = p.Results.tvec;
SR = p.Results.SR;
sorter = p.Results.sorter;

% Load cluster table and find good units based on quality metrics
if contains(sorter,'off')
    load([options.KS_DATAPATH,'\cluster_table.mat'])
    cluster_metrics = readtable(fullfile(options.KS_DATAPATH,'metrics.csv')); % usually not all clusters (some that does not pass Kilosort threshold would not be included here)

    % Currently unit selection criteria based allen sdk
    % https://allensdk.readthedocs.io/en/latest/_static/examples/nb/ecephys_quality_metrics.html
    amplitude_cutoff = cluster_metrics.amplitude_cutoff; % <0.1 (just in case there are cells that are selectively active)
    isi_violation = cluster_metrics.isi_viol; % < 0.5
    presence_ratio = cluster_metrics.presence_ratio; % > 0.5 %Should be at least active more than 50% of the time
    troughToPeak = cluster_metrics.duration; % Duration (in seconds) between peak and trough

    good_unit = cluster_metrics.amplitude_cutoff <= 0.1...
        &cluster_metrics.isi_viol <= 0.1...
        &cluster_metrics.presence_ratio >=0.5;

    mean_waveforms = readNPY([options.KS_DATAPATH,'/mean_waveforms.npy']); % Information about mean waveform (for all clusters)
    %     cluster_ContamPct = tdfread(fullfile(options.KS_DATAPATH,'cluster_ContamPct.tsv'));% contamination percentage

    % Load spike time data and quality metrics
    [these_spike_times,nominal_KSLabel,cluster_id,peakChannel,maxSpkTime] = import_ks_spiketimes(options,options.gFileNum,options.KS_CATGT_FNAME,imecMeta.imSampRate);
    cluster_id = cluster_id + 1; % Cluster ID transformed from 0-based to 1-based

    % scatter3(cluster_metrics.amplitude_cutoff,cluster_metrics.isi_viol,cluster_metrics.presence_ratio)
    % hold on
    % scatter3(cluster_metrics.amplitude_cutoff(cluster_metrics.amplitude_cutoff<0.1),cluster_metrics.isi_viol(cluster_metrics.amplitude_cutoff<0.1),cluster_metrics.presence_ratio(cluster_metrics.amplitude_cutoff<0.1),'filled','r','MarkerFaceAlpha','0.4')
    % xlabel('amplitude cutoff')
    % ylabel('isi ratio')
    % zlabel('presence ratio')
    tic
    parfor id = 1:length(cluster_id)
        peak_channel_waveforms(id,:) = squeeze(mean_waveforms(cluster_id(id),peakChannel(id),:)); % Here ID is pre-KS id
        [ccg, t] = CCG(    these_spike_times{id}, ones(size(these_spike_times{id})),...
            'binSize', 0.0005, 'duration', 0.1,'norm', 'rate');
        if isempty(ccg)
            cluster_acg(id,:) = zeros(size(t));
            tau_rise(id) = nan;
            troughToPeak_all(id) = nan;
            good_unit_all(id) = 0;
        else
            cluster_acg(id,:) = ccg;
            fit_params_out = [];
            fit_params_out = fit_ACG(ccg,false);% fit_ACG and isToolboxInstalled function from cell explorer repo
            tau_rise(id) = fit_params_out.acg_tau_rise;
            troughToPeak_all(id) = troughToPeak(cluster_metrics.cluster_id == cluster_id(id)-1) % -1 as cluster_id is 1 based (after adding 1) but cluster_metrics.cluster_id is 0 based (post-KS)
            good_unit_all(id) = good_unit(cluster_metrics.cluster_id == cluster_id(id)-1); % good_unit_all for pre-ks id.
        end
    end
    toc

    % 1 is pyramidal, 2 is narrow and 3 is wide
    narrow_idx = find([troughToPeak_all]<=0.425);
    wide_idx = find([troughToPeak_all]>0.425 & [tau_rise]>6);
    pyr_idx = find(~ismember(1:length(cluster_id), [narrow_idx,wide_idx]));

    [cell_type_all] = deal(nan);
    [cell_type_all(pyr_idx)] = deal(1);
    [cell_type_all(narrow_idx)] = deal(2);
    [cell_type_all(wide_idx)] = deal(3);
    % cell_type_all(good_unit_all==0) = nan;

else % place holder for other sorters


end

clusters = [];

if ~isempty(tvec)
    %smooth SUA activity with a gaussian kernel
    time_step = 1/SR;
    filter_length = round(0.2*SR);
    filter_alpha = 2;
    %     stdev = 3;
    %     filter_alpha = (filter_length - 1)/stdev/2;
    stdev = (filter_length-1)/(2*filter_alpha);
    gk=gausswin(filter_length);
    gk = gk./sum(gk);

%         clusters.spike_count_smoothed = conv2(SUA_spike_count_raw,gk,'same');

    time_bins_edges= tvec(1)-(tvec(2) - tvec(1))/2:time_step:tvec(end)+(tvec(end) - tvec(end-1))/2;
end


tic
switch group
    case 'good clusters'

        count = 1;
        unit_id = [];
        unit_type = [];
        SUA_spike_time = [];
        SUA_spike_id = [];
        SUA_peak_channel_waveforms = [];
        SUA_peak_channel = [];
        SUA_peak_depth = [];
        SUA_cell_type = [];
        id_conversion = [];
        SUA_spike_count_raw = [];
        SUA_spike_count_smoothed = [];
        SUA_zscore = [];
        good_clusters_all= [];
        SUA_good_unit = [];
        
        for nchannel = selected_channels(1):selected_channels(2)
            clusters_this_channel = find(peakChannel == nchannel); % ID for all clusters
            [~,index,]= intersect(cluster_id,clusters_this_channel);

            good_units_index = index(find(good_unit_all(index)==1));
            %             good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(index)=='good')); % ID for all clusters

            if ~isempty(index) % If any clusters
                if ~isempty(good_units_index) % If any SUA
                    for unit = 1:length(good_units_index)

                        SUA_spike_time = [SUA_spike_time; these_spike_times{good_units_index(unit)}];
                        SUA_spike_id = [SUA_spike_id; cluster_id(good_units_index(unit))*ones(length(these_spike_times{good_units_index(unit)}),1)];

                        SUA_spike_count_raw(count,:) = histcounts(these_spike_times{good_units_index(unit)},time_bins_edges);
%                         SUA_spike_count_smoothed(count,:) = filtfilt(w,1,histcounts(these_spike_times{good_units_index(unit)},time_bins_edges));
%                         SUA_zscore(count,:) = zscore(SUA_spike_count_smoothed(count,:));

                        SUA_peak_channel(count) = nchannel;
                        SUA_peak_depth(count) = chan_config.Ks_ycoord(chan_config.Channel == nchannel);

                        SUA_good_unit(count) = good_units_index(unit);
                        good_clusters_all(count) = good_unit_all(good_units_index(unit));
                        SUA_peak_channel_waveforms(count,:) = peak_channel_waveforms(good_units_index(unit),:);
                        SUA_cell_type(count) = cell_type_all(good_units_index(unit));
                        
                        count = count + 1;
                    end
                end
            end
        end

        % Continuous spike data
        clusters.timevec = tvec;
        clusters.spike_count_raw = SUA_spike_count_raw;
%         clusters.spike_count_smoothed = SUA_spike_count_smoothed; 
         clusters.spike_count_smoothed = filtfilt(gk,1,SUA_spike_count_raw')';
        clusters.zscore_smoothed = zscore(clusters.spike_count_smoothed,0,2);
        clusters.unit_id = SUA_good_unit;
        clusters.good_clusters = good_clusters_all;

        % Spike time and spike id
        [~,index] = sort(SUA_spike_time); % Sort spike time from start to end
        clusters.spike_id = SUA_spike_id(index);
        clusters.spike_times = SUA_spike_time(index);
        clusters.peak_channel = SUA_peak_channel;
        clusters.peak_depth = SUA_peak_depth;
        clusters.peak_channel_waveforms = SUA_peak_channel_waveforms;
        clusters.cell_type = SUA_cell_type;
        clusters.probe_id = options.probe_id; % 0 is probe 1 and 1 is probe 2.

        if isfield(options,'probe_hemisphere')
            clusters.probe_hemisphere = options.probe_hemisphere; % 1 is left and 2 is right;
        end

    case 'all clusters'

        count = 1;
        unit_id = [];
        unit_type = [];
        SUA_spike_time = [];
        SUA_spike_id = [];
        SUA_peak_channel_waveforms = [];
        SUA_peak_channel = [];
        SUA_peak_depth = [];
        SUA_cell_type = [];
        id_conversion = [];
        SUA_spike_count_raw = [];
        SUA_spike_count_smoothed = [];
        SUA_zscore = [];
        good_clusters_all= [];
        SUA_good_unit = [];

        for nchannel = selected_channels(1):selected_channels(2)
            clusters_this_channel = find(peakChannel == nchannel); % ID for all clusters
            [~,index,]= intersect(cluster_id,clusters_this_channel);
            %             good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(index)=='good')); % ID for all clusters

            if ~isempty(index) % If any clusters
                for unit = 1:length(index)

                    SUA_spike_time = [SUA_spike_time; these_spike_times{index(unit)}];
                    SUA_spike_id = [SUA_spike_id; cluster_id(index(unit))*ones(length(these_spike_times{index(unit)}),1)];

                    SUA_spike_count_raw(count,:) = histcounts(these_spike_times{index(unit)},time_bins_edges);
%                     SUA_spike_count_smoothed(count,:) = filtfilt(w,1,histcounts(these_spike_times{index(unit)},time_bins_edges));
%                     SUA_zscore(count,:) = zscore(SUA_spike_count_smoothed(count,:));

                    SUA_peak_channel(count) = nchannel;
                    SUA_peak_depth(count) = chan_config.Ks_ycoord(chan_config.Channel == nchannel);
                    good_clusters_all(count) = good_unit_all(index(unit));

                    SUA_good_unit(count) = index(unit);
                    SUA_peak_channel_waveforms(count,:) = peak_channel_waveforms(index(unit),:);
                    SUA_cell_type(count) = cell_type_all(index(unit));

                    count = count + 1;
                end
            end
        end

        % Continuous spike data
        clusters.timevec = tvec;
        clusters.spike_count_raw = SUA_spike_count_raw;
%         clusters.spike_count_smoothed = SUA_spike_count_smoothed; 
         clusters.spike_count_smoothed = filtfilt(gk,1,SUA_spike_count_raw')';
        clusters.zscore_smoothed = zscore(clusters.spike_count_smoothed,0,2);
        clusters.zscore_smoothed = SUA_zscore;
        clusters.unit_id = SUA_good_unit;
        clusters.good_clusters = good_clusters_all;
        

        % Spike time and spike id
        [~,index] = sort(SUA_spike_time); % Sort spike time from start to end
        clusters.spike_id = SUA_spike_id(index);
        clusters.spike_times = SUA_spike_time(index);
        clusters.peak_channel = SUA_peak_channel;
        clusters.peak_depth = SUA_peak_depth;
        clusters.peak_channel_waveforms = SUA_peak_channel_waveforms;
        clusters.cell_type = SUA_cell_type;
        clusters.probe_id = options.probe_id; % 0 is probe 1 and 1 is probe 2.
        clusters.group = group;% good clusters or all clusters or etc

        if isfield(options,'probe_hemisphere')
            clusters.probe_hemisphere = options.probe_hemisphere; % 1 is left and 2 is right;
        end


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
toc

end