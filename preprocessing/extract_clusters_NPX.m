function [clusters chan_config sorted_config] = extract_clusters_NPX(options,varargin)

%%% Function to extract clusters based on different sorters
% modified/simplified from old function:
% [clusters chan_config sorted_config] = load_KS_NPX1(options,column,varargin)
% Modified from original version of load_KS_NPX1 function


% Extract NPX channel configuration
options.importMode = 'KS';
[file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,1);

% Default values
p = inputParser;
addParameter(p,'selected_channels',[1 size(chan_config,1)],@isnumeric) % Select channels for analysis (default is all the channles)
addParameter(p,'group','all',@isstr) % spike data grouped 'all clusters'
addParameter(p,'plot_option',0,@isnumeric) % option to plot
addParameter(p,'SR',60,@isnumeric) % SR for cluster spike count time course (normally 60Hz)
addParameter(p,'tvec',[],@isnumeric) % Time vector

% addParameter(p,'cell_exporer','off',@isstr) % cell explorer option
% removed (replaced by few lines of code that can classify the cell based on CCG)

% addParameter(p,'waveform','off',@isstr) % extract individual spike waveform for each cluster (not mean waveform)
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

if ~contains(sorter,'off') % if not default ks
    sorter_folder = [];
    if contains(sorter,'KS2')
        sorter_folder = 'kilosort2_5';
%     elseif contains(sorter,'KS3_merged')
%         sorter_folder = 'kilosort3_merged';
%     elseif contains(sorter,'KS4_merged')
%         sorter_folder = 'kilosort4_merged';

    elseif contains(sorter,'KS3')
        sorter_folder = 'kilosort3';
    elseif contains(sorter,'KS4')
        sorter_folder = 'kilosort4';
    end
end


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


if contains(sorter,'off') % Load cluster table and find good units based on quality metrics
    load([options.KS_DATAPATH,'\cluster_table.mat'])
    cluster_metrics = readtable(fullfile(options.KS_DATAPATH,'metrics.csv')); % usually not all clusters (some that does not pass Kilosort threshold would not be included here)

    if sum(contains(cluster_metrics.Properties.VariableNames,'epoch_name_quality_metrics')) == 1
        cluster_metrics.epoch_name_quality_metrics = [];
    end

    if sum(contains(cluster_metrics.Properties.VariableNames,'epoch_name_waveform_metrics')) == 1
        cluster_metrics.epoch_name_waveform_metrics = [];
    end

    % Currently unit selection criteria based allen sdk
    % https://allensdk.readthedocs.io/en/latest/_static/examples/nb/ecephys_quality_metrics.html
    amplitude_cutoff = cluster_metrics.amplitude_cutoff; % <0.1 (just in case there are cells that are selectively active)
    isi_violation = cluster_metrics.isi_viol; % < 0.5
    presence_ratio = cluster_metrics.presence_ratio; % > 0.5 %Should be at least active more than 50% of the time
    troughToPeak = cluster_metrics.duration; % Duration (in seconds) between peak and trough
    amplitude = cluster_metrics.amplitude;


    if contains(group,'good') % good clusters modified from Allen metrics
        good_unit = cluster_metrics.amplitude_cutoff <= 0.1...
            &cluster_metrics.isi_viol <= 0.1...
            &cluster_metrics.amplitude >=50;
    elseif contains(group,'all') % all clusters
        good_unit = ones(size(cluster_metrics,1),1);
    end

    mean_waveforms = readNPY([options.KS_DATAPATH,'/mean_waveforms.npy']); % Information about mean waveform (for all clusters)
    %     cluster_ContamPct = tdfread(fullfile(options.KS_DATAPATH,'cluster_ContamPct.tsv'));% contamination percentage

    % Load spike time data and quality metrics
    [these_spike_times,nominal_KSLabel,cluster_id,peak_channel,maxSpkTime] = import_ks_spiketimes(options,options.gFileNum,options.KS_CATGT_FNAME,imecMeta.imSampRate);
    cluster_id = cluster_id + 1; % Cluster ID transformed from 0-based to 1-based

    % scatter3(cluster_metrics.amplitude_cutoff,cluster_metrics.isi_viol,cluster_metrics.presence_ratio)
    % hold on
    % scatter3(cluster_metrics.amplitude_cutoff(cluster_metrics.amplitude_cutoff<0.1),cluster_metrics.isi_viol(cluster_metrics.amplitude_cutoff<0.1),cluster_metrics.presence_ratio(cluster_metrics.amplitude_cutoff<0.1),'filled','r','MarkerFaceAlpha','0.4')
    % xlabel('amplitude cutoff')
    % ylabel('isi ratio')
    % zlabel('presence ratio')

    tic
    parfor id = 1:length(cluster_id)
        peak_channel_waveforms(id,:) = squeeze(mean_waveforms(cluster_id(id),peak_channel(id),:)); % Here ID is pre-KS id
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
            troughToPeak_all(id) = troughToPeak(cluster_metrics.cluster_id == cluster_id(id)-1); % -1 as cluster_id is 1 based (after adding 1) but cluster_metrics.cluster_id is 0 based (post-KS)
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


    tic

    count = 1;
    clusters = [];

    SUA_spike_time = [];
    SUA_spike_id = [];
    SUA_peak_channel_waveforms = [];
    SUA_peak_channel = [];
    SUA_peak_depth = [];
    SUA_cell_type = [];
    SUA_spike_count_raw = [];
    % SUA_spike_count_smoothed = [];
    % SUA_zscore = [];
    % good_clusters_all= [];
    SUA_good_unit = [];

    for nchannel = selected_channels(1):selected_channels(2)
        clusters_this_channel = find(peak_channel == nchannel); % ID for all clusters
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
                    %                         good_clusters_all(count) = good_unit_all(good_units_index(unit));
                    SUA_peak_channel_waveforms(count,:) = peak_channel_waveforms(good_units_index(unit),:);
                    SUA_cell_type(count) = cell_type_all(good_units_index(unit));

                    count = count + 1;
                end
            end
        end
    end

    [~,~,index] = intersect(SUA_good_unit-1,cluster_metrics.cluster_id);
    clusters = table2struct(cluster_metrics(index,:),"ToScalar",true);
    clusters.cluster_id = clusters.cluster_id+1; %1 based

    % Continuous spike data (retired)
    clusters.timevec = tvec; % for future continuous spike data at 60Hz
    %         clusters.spike_count_raw = SUA_spike_count_raw;
    %         clusters.spike_count_smoothed = SUA_spike_count_smoothed;
    %          clusters.spike_count_smoothed = filtfilt(gk,1,SUA_spike_count_raw')';
    %         clusters.zscore_smoothed = zscore(clusters.spike_count_smoothed,0,2);
    % clusters.unit_id = SUA_good_unit; % 1 based cluster id.
    %         clusters.good_clusters = good_clusters_all;
    % Spike time and spike id
    [~,index] = sort(SUA_spike_time); % Sort spike time from start to end
    clusters.spike_id = SUA_spike_id(index);
    clusters.spike_times = SUA_spike_time(index);
    clusters.peak_channel = SUA_peak_channel';
    clusters.peak_depth = SUA_peak_depth';
    clusters.peak_channel_waveforms = SUA_peak_channel_waveforms;
    clusters.cell_type = SUA_cell_type';
    clusters.probe_id = options.probe_id'; % 0 is probe 1 and 1 is probe 2.
    clusters.sorter = sorter; % what sorter is used for spike sorting and clustering

    if isfield(options,'probe_hemisphere')
        clusters.probe_hemisphere = options.probe_hemisphere; % 1 is left and 2 is right;
    end

    toc


elseif contains(sorter,'other sorter')  % place holder for other sorters



else % if KS2 and KS3
    quality_metrics = readtable(fullfile(options.SORTER_DATAPATH,'waveform',[sorter_folder,'_merged'],'extensions','quality_metrics','metrics.csv'));
    template_metrics = readtable(fullfile(options.SORTER_DATAPATH,'waveform',[sorter_folder,'_merged'],'extensions','template_metrics','metrics.csv'));


    CCG_bin = readNPY(fullfile(options.SORTER_DATAPATH,'waveform',[sorter_folder,'_merged'],'extensions','correlograms','bins.npy'));
    CCGs= readNPY(fullfile(options.SORTER_DATAPATH,'waveform',[sorter_folder,'_merged'],'extensions','correlograms','CCGs.npy'));

    cluster_metrics = join(quality_metrics,template_metrics);

    if contains(group,'good') % good clusters modified from Allen metrics
        good_unit = cluster_metrics.amplitude_cutoff <= 0.1...
            &cluster_metrics.isi_violations_ratio <= 0.1...
            &cluster_metrics.amplitude_median >=50;
    elseif contains(group,'all') % all clusters
        good_unit = ones(size(cluster_metrics,1),1);
    end

    %     cluster_ContamPct = tdfread(fullfile(options.KS_DATAPATH,'cluster_ContamPct.tsv'));% contamination percentage

    % Load spike time data and quality metrics
    options.sorter_folder = sorter_folder;
    [these_spike_times,cluster_id,peak_channel,peak_depth,peak_channel_waveforms] = import_spikeinterface_spiketimes(options,imecMeta.imSampRate);% Only clusters after spike interface post processing
    cluster_id = cluster_id + 1; % Cluster ID transformed from 0-based to 1-based

    % scatter3(cluster_metrics.amplitude_cutoff,cluster_metrics.isi_viol,cluster_metrics.presence_ratio)
    % hold on
    % scatter3(cluster_metrics.amplitude_cutoff(cluster_metrics.amplitude_cutoff<0.1),cluster_metrics.isi_viol(cluster_metrics.amplitude_cutoff<0.1),cluster_metrics.presence_ratio(cluster_metrics.amplitude_cutoff<0.1),'filled','r','MarkerFaceAlpha','0.4')
    % xlabel('amplitude cutoff')
    % ylabel('isi ratio')
    % zlabel('presence ratio')

    tic
    parfor id = 1:length(cluster_id)
        %         peak_channel_waveforms(id,:) = squeeze(mean_waveforms(cluster_id(id),peakChannel(id),:)); % Here ID is pre-KS id
        [ccg, t] = CCG(    these_spike_times{id}, ones(size(these_spike_times{id})),...
            'binSize', 0.0005, 'duration', 0.1,'norm', 'rate');
        if isempty(ccg)
            cluster_acg(id,:) = zeros(size(t));
            tau_rise(id) = nan;
            %good_unit_all(id) = 0;
        else
            cluster_acg(id,:) = ccg;
            fit_params_out = [];
            fit_params_out = fit_ACG(ccg,false);% fit_ACG and isToolboxInstalled function from cell explorer repo
            tau_rise(id) = fit_params_out.acg_tau_rise;
            %good_unit_all(id) = good_unit(id);
        end
    end
    good_unit_all = good_unit;
    toc

    % 1 is pyramidal, 2 is narrow and 3 is wide
    narrow_idx = find([1000*cluster_metrics.peak_to_valley]<=0.425);
    wide_idx = find([1000*cluster_metrics.peak_to_valley]>0.425 & [tau_rise']>6);
    pyr_idx = find(~ismember(1:length(cluster_id), [narrow_idx;wide_idx]));

    [cell_type_all] = deal(nan);
    [cell_type_all(pyr_idx)] = deal(1);
    [cell_type_all(narrow_idx)] = deal(2);
    [cell_type_all(wide_idx)] = deal(3);
    % cell_type_all(good_unit_all==0) = nan;
    % cell_type_all(good_unit_all==0) = nan;

    tic

    count = 1;
    clusters = [];
    SUA_spike_time = [];
    SUA_spike_id = [];
    SUA_peak_channel_waveforms = [];
    SUA_peak_channel = [];
    SUA_peak_depth = [];
    SUA_cell_type = [];
    SUA_spike_count_raw = [];
    % SUA_spike_count_smoothed = [];
    % SUA_zscore = [];
    % good_clusters_all= [];
    % SUA_good_unit = [];

    for nchannel = selected_channels(1):selected_channels(2)
        index = find(peak_channel == nchannel); % ID for post-processed clusters
        %         [~,index,]= intersect(cluster_id,clusters_this_channel);

        good_units_index = index(find(good_unit_all(index)==1)); % sometimes the unit may not fire during a short recording
        %             good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(index)=='good')); % ID for all clusters

        if ~isempty(index) % If any clusters
            if ~isempty(good_units_index) % If any SUA
                for unit = 1:length(good_units_index)

                    SUA_spike_time = [SUA_spike_time; these_spike_times{good_units_index(unit)}];
                    SUA_spike_id = [SUA_spike_id; cluster_id(good_units_index(unit))*ones(length(these_spike_times{good_units_index(unit)}),1)];

                    SUA_spike_count_raw(count,:) = histcounts(these_spike_times{good_units_index(unit)},time_bins_edges);
                    %                         SUA_spike_count_smoothed(count,:) = filtfilt(w,1,histcounts(these_spike_times{good_units_index(unit)},time_bins_edges));
                    %                         SUA_zscore(count,:) = zscore(SUA_spike_count_smoothed(count,:));

                    %this part below has a different way of indexing than cluster_id, hence wrong 
%                     SUA_peak_channel(count) = nchannel;
%                     SUA_peak_depth(count) = chan_config.Ks_ycoord(chan_config.Channel == nchannel);
% 
%                     SUA_good_unit(count) = good_units_index(unit);
%                     %                         good_clusters_all(count) = good_unit_all(good_units_index(unit));
%                     SUA_peak_channel_waveforms(count,:) = peak_channel_waveforms(good_units_index(unit),:);
%                     SUA_cell_type(count) = cell_type_all(good_units_index(unit));

                    count = count + 1;
                end
            end
        end
    end

    %     [~,~,index] = intersect(SUA_good_unit-1,table2array(cluster_metrics(:,1)));
    clusters = table2struct(cluster_metrics,"ToScalar",true);
    clusters.cluster_id = clusters.Var1+1; %1 based

    % Continuous spike data (retired)
    clusters.timevec = tvec; % for future continuous spike data at 60Hz
    %         clusters.spike_count_raw = SUA_spike_count_raw;
    %         clusters.spike_count_smoothed = SUA_spike_count_smoothed;
    %          clusters.spike_count_smoothed = filtfilt(gk,1,SUA_spike_count_raw')';
    %         clusters.zscore_smoothed = zscore(clusters.spike_count_smoothed,0,2);
    % clusters.unit_id = SUA_good_unit; % 1 based cluster id.
    %         clusters.good_clusters = good_clusters_all;
    % Spike time and spike id
    [~,index] = sort(SUA_spike_time); % Sort spike time from start to end
    clusters.spike_id = SUA_spike_id;
    clusters.spike_times = SUA_spike_time;
    clusters.peak_channel = peak_channel';
    clusters.peak_depth = peak_depth;
    clusters.peak_channel_waveforms = peak_channel_waveforms;
    clusters.cell_type = cell_type_all';
    clusters.probe_id = options.probe_id'; % 0 is probe 1 and 1 is probe 2.
    clusters.sorter = sorter; % what sorter is used for spike sorting and clustering
    %     clusters.CCGs = CCGs;

    if isfield(options,'probe_hemisphere')
        clusters.probe_hemisphere = options.probe_hemisphere; % 1 is left and 2 is right;
    end

    toc


end

end

