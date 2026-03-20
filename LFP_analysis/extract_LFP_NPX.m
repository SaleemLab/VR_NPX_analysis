function extract_LFP_NPX(session_info,stimulus_name,best_channels)

options = session_info.probe(1);
% load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));
if isempty(DIR)
    disp('No extracted clusters')
    return
end
DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
DIR1 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));

if ~isempty(DIR)
    load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
    session_clusters_RUN=session_clusters;
    clear session_clusters
end

if ~isempty(DIR1)
    load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
    session_clusters_RUN=session_clusters;
    clear session_clusters
end

if contains(stimulus_name,'Masa2tracks')
    load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters%s.mat',erase(stimulus_name,'Masa2tracks'))),'session_clusters');
    load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name,'Masa2tracks'))));
    load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks4%s.mat',erase(stimulus_name,'Masa2tracks'))));
    clusters = clusters_ks4;
else
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
    load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name,'Chronic'))),'session_clusters'); % Session clusters for SUA
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));% clusters for MUA
    clusters = clusters_ks4;
end

% From cell structure to all cluster spike time and spike id
session_clusters.spike_times=vertcat(session_clusters.spike_times{:});
session_clusters.spike_id=vertcat(session_clusters.spike_id{:});
[session_clusters.spike_times,index] =sort(session_clusters.spike_times);
session_clusters.spike_id=session_clusters.spike_id(index);

% check cluster region label
if ~isfield(clusters(1),'region')
    for nprobe = 1:length(session_info.probe)
        clusters(nprobe).region = strings(length(clusters(nprobe).cluster_id),1);
        if clusters(nprobe).probe_hemisphere == 1
            clusters(nprobe).region(:) = 'n.a_L';
            kClusters=kmeans(clusters(nprobe).peak_depth,2);
            if mean(clusters(nprobe).peak_depth(kClusters==1))>mean(clusters(nprobe).peak_depth(kClusters==2))
                % if mean ocation of cluster one is above cluster two, it is
                % Cortex.
                % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==1);
                V1_cell_id = find(kClusters==1);
                HPC_cell_id = find(kClusters==2);
            else
                % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==2);
                V1_cell_id = find(kClusters==2);
                HPC_cell_id = find(kClusters==1);
            end
            clusters(nprobe).region(V1_cell_id) = 'V1_L';
            clusters(nprobe).region(HPC_cell_id) = 'HPC_L';

        elseif clusters(nprobe).probe_hemisphere == 2
            clusters(nprobe).region(:) = 'n.a_R';

            kClusters=kmeans(clusters(nprobe).peak_depth,2);
            if mean(clusters(nprobe).peak_depth(kClusters==1))>mean(clusters(nprobe).peak_depth(kClusters==2))
                % if mean ocation of cluster one is above cluster two, it is
                % Cortex.
                % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==1);
                V1_cell_id = find(kClusters==1);
                HPC_cell_id = find(kClusters==2);
            else
                % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==2);
                V1_cell_id = find(kClusters==2);
                HPC_cell_id = find(kClusters==1);
            end
            clusters(nprobe).region(V1_cell_id) = 'V1_R';
            clusters(nprobe).region(HPC_cell_id) = 'HPC_R';
        end
    end
end


for nprobe = 1:length(session_info.probe)
    options = session_info.probe(nprobe);
    selected_channels = [];
    channel_regions = [];
    shank_id = [];
    options.importMode = 'KS';
    [file_to_use imecMeta chan_config ~] = extract_NPX_channel_config(options,[]);% Since it is LF

    all_fields = fieldnames(best_channels{nprobe});
    all_fields = {all_fields{contains(all_fields,'depth')}};
    all_fields = erase(all_fields,'_depth');
    all_fields{length(all_fields)+1} = 'V1_sparse'; % Add cortex sparse
    all_fields{length(all_fields)+1} = 'HPC_sparse'; % Add HPC sparse

    for nregion = 1:length(all_fields)
        region_name = all_fields{nregion};
        [region_channels,channels_temp] = determine_region_channels(best_channels{nprobe},options,'region',region_name,'group','by shank','clusters',clusters(nprobe));
        %                 noise_channel{options.probe_hemisphere} = channels_temp;
        if contains(all_fields{nregion},'HPC')
            region_channels(chan_config.Ks_ycoord(region_channels) > mean(chan_config.Ks_ycoord))=[];
        elseif contains(all_fields{nregion},'V1')
            region_channels(chan_config.Ks_ycoord(region_channels) < mean(chan_config.Ks_ycoord))=[];
        end
        
        if contains(all_fields{nregion},'sparse') % want sparsely sampled channels
            selected_channels = [selected_channels region_channels'];
            channel_regions = [channel_regions nregion*ones(1,length(region_channels))];
%             shank_id = [shank_id chan_config.Shank(region_channels)'];
        else % Just want the best channel
            selected_channels = [selected_channels channels_temp];
            channel_regions = [channel_regions nregion*ones(1,length(channels_temp))];
%             shank_id = [shank_id chan_config.Shank(region_channels)'];
        end
        %                     shank_id = [shank_id chan_config.Shank(channels_temp)'];
    end

    channel_regions(isnan(selected_channels)) = []; % remove nan channel (Missing best channels for some shanks e.g. only 3 shanks with CA1)
%     shank_id(isnan(selected_channels)) = [];
    selected_channels(isnan(selected_channels)) = [];
    [unique_selected_channels,ia,ic] = unique(selected_channels);
    shank_id = chan_config.Shank(selected_channels);
    unique_shank_id = chan_config.Shank(unique_selected_channels);
%     unique_shank_id = shank_id(ia);

    % Extract LFP
    [raw_LFP,tvec,SR,chan_config,~] = load_LFP_NPX(options,[],'selected_channels',unique_selected_channels);

    % Calculate PSD
    selected_chan_config = chan_config(unique_selected_channels,:);
    [PSD{nprobe},power{nprobe}] = calculate_channel_PSD(raw_LFP,SR,selected_chan_config,options,'nfft_seconds',2,'plot_option',0,'tvec',tvec);

    rms_values = rms(raw_LFP, 2);
    channel_std = std(raw_LFP, 0, 2);

    %%%%% Find CA1 channels with best ripple power
    HPC_channel_id = find(ismember(unique_selected_channels,selected_channels(channel_regions==find(contains(all_fields,'HPC_sparse')))));
    HPC_channels = unique_selected_channels(HPC_channel_id);
    % Putatively remove bad channels based on PSD being too
    % high for wide frequency range
    bad_channels =[];
    for nChannel = 1:length(HPC_channel_id)
        not_this_shank = unique_shank_id(HPC_channel_id(nChannel))~=unique_shank_id(HPC_channel_id);% find channels not this shanks
        this_shank = unique_shank_id(HPC_channel_id(nChannel))==unique_shank_id(HPC_channel_id);% find channels this shanks
        bad_channels(nChannel,:) = rms_values(nChannel) > mean(rms_values) + 3 * std(rms_values) | channel_std(nChannel) > 5 * median(channel_std); % remove noisy channels based RMS and STD
%         power{nprobe}(HPC_channel_id(nChannel),:)>2*mean(power{nprobe}(not_this_shank,:))
    end
%     bad_channels = sum(bad_channels,2)>3;

    good_HPC_channels=[];
    good_HPC_channels = HPC_channel_id(~bad_channels);
    [~,best_channel]=max(power{nprobe}(good_HPC_channels,6)); % Best CA1 channel with highest ripple power
    %                 good_channels = find(~bad_channels);
    best_HPC_channel = good_HPC_channels(best_channel);
    %                 CA1_LFP = raw_LFP(best_CA1_channel,:);

    best_HPC_shank_channel_id=[];
    % Best HPC channels for each shank
    unique_shanks = unique(unique_shank_id(good_HPC_channels));
    for nShank = 1:length(unique_shanks)
        this_shank = find(unique_shank_id(good_HPC_channels)==unique_shanks(nShank));% find channels not this shanks
        [~,channel_id]=max(power{nprobe}(good_HPC_channels(this_shank),6)); % Best CA1 channel with highest ripple power
        best_HPC_shank_channel_id(nShank) = good_HPC_channels(this_shank(channel_id));
    end
    
    %%%%% Find best V1 channels in terms of high freq power
    V1_channel_id = find(ismember(unique_selected_channels,selected_channels(channel_regions==find(contains(all_fields,'V1_sparse')))));
    V1_channels = unique_selected_channels(V1_channel_id);
    bad_channels =[];
    for nChannel = 1:length(V1_channel_id)
        not_this_shank = unique_shank_id(V1_channel_id(nChannel))~=unique_shank_id(V1_channel_id);% find channels not this shanks
        bad_channels(nChannel,:) = rms_values(nChannel) > mean(rms_values) + 3 * std(rms_values) | channel_std(nChannel) > 5 * median(channel_std); % remove noisy channels based RMS and STD
    end
%     bad_channels = sum(bad_channels,2)>3;

    good_V1_channels=[];
    good_V1_channels = V1_channel_id(~bad_channels);
    [~,best_channel]=max(power{nprobe}(good_V1_channels,7)); % best cortex channel based on high freq power (best putative L5)
    best_V1_channel = good_V1_channels(best_channel);
    % %                  plot(5000*power{nprobe}(good_channels,7));hold on; plot(power{nprobe}(good_channels,1))
    %                 cortex_LFP=raw_LFP(good_channels(best_V1_channel),:); % get mean V1 cortex LFP
    % Best V1 channels high freq power for each shank
    best_V1_high_power_shank_channel_id=[];
    unique_shanks = unique(unique_shank_id(good_V1_channels));
    for nShank = 1:length(unique_shanks)
        this_shank = find(unique_shank_id(good_V1_channels)==unique_shanks(nShank));% find channels not this shanks
        [~,channel_id]=max(power{nprobe}(good_V1_channels(this_shank),7)); % Best V1 channel with highest high freq power
        best_V1_high_power_shank_channel_id(nShank) = good_V1_channels(this_shank(channel_id));
    end

    %%%%% V1 MUA spikes

    metric_param =[]; % get all V1 spikes from left or right hemisphere
    metric_param = create_cluster_selection_params('sorting_option','masa');
    metric_param.firing_rate = @(x)x>=0.01&x<50;
    if options.probe_hemisphere==1
        metric_param.region = @(x) contains(x,'V1_L');
        V1_clusters = select_clusters(clusters(nprobe),metric_param);
        metric_param.region = @(x) contains(x,'HPC_L');
        HPC_clusters = select_clusters(clusters(nprobe),metric_param);
    elseif options.probe_hemisphere==2
        metric_param.region = @(x) contains(x,'V1_R');
        V1_clusters = select_clusters(clusters(nprobe),metric_param);
        metric_param.region = @(x) contains(x,'HPC_R');
        HPC_clusters = select_clusters(clusters(nprobe),metric_param);
    end

    LFP_SR = 1/mean(diff(tvec));
    tvec_edges = [tvec(1)-1/(1/mean(diff(tvec))*2) tvec+1/(1/mean(diff(tvec))*2)];


    %%%%%% Identify best LFP with SO phase 

    deltaspikecorr=[];
    gammaspikecorr=[];
    deltagammacorr=[];
    filterparms.deltafilter = [0.5 4];
    lfp.timestamps = tvec;

%     chan_config.Ks_ycoord(ismember(chan_config.Channel,V1_channels))
    n_bins = 20; % Number of phase bins
    edges = linspace(-pi, pi, n_bins+1);
    centers = (edges(1:end-1) + edges(2:end)) / 2;


    for nChannel = 1:length(good_V1_channels)
        tic
        deltaLFP=[];
        %             gammaLFP=[];
        
%         lfp.data= raw_LFP(good_V1_channels(nChannel),:)';
%         deltaLFP1 = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);

        filt_order = round(1*LFP_SR./filterparms.deltafilter(1));
        [b a] = fir1(filt_order,filterparms.deltafilter/LFP_SR*2);
        deltaLFP = filtfilt(b,a,double(raw_LFP(good_V1_channels(nChannel),:)));
        deltaLFP= angle(hilbert(deltaLFP));

%         deltaLFP.phase = interp1(deltaLFP.timestamps,deltaLFP.phase',tvec_NREM,'nearest');
%         deltaLFP.phase = interp1(deltaLFP.timestamps,deltaLFP.phase',tvec_NREM,'nearest');
        %% 1. Settings

        %% 2. Pre-processing: Filter Spikes & Unwrap Phase
        % Filter spikes to ensure they fall within the LFP recording time
        valid_spike_indices = V1_clusters.spike_times >= min(tvec) & ...
            V1_clusters.spike_times <= max(tvec);
        valid_spikes = V1_clusters.spike_times(valid_spike_indices);

        % Unwrap phase to allow for linear interpolation across the -pi/pi boundary
        phi_unwrapped = unwrap(deltaLFP);

        % Interpolate LFP phase to the exact spike times
        spike_phases_unwrapped = interp1(tvec, phi_unwrapped, valid_spikes, 'linear');

        % Re-wrap spike phases back to [-pi, pi]
        spike_phases = mod(spike_phases_unwrapped + pi, 2*pi) - pi;

        %% 3. Calculate Time Occupancy (Denominator)
        % Calculate the duration of a single LFP sample
        dt = mean(diff(tvec));

        % Count how many LFP samples fall in each phase bin
        [lfp_counts, ~] = histcounts(deltaLFP, edges);

        % Convert sample counts to seconds (Time Spent in each phase)
        time_occupancy = lfp_counts * dt;

        %% 4. Calculate Spike Distribution (Numerator)
        [spike_counts, ~] = histcounts(spike_phases, edges);

        %% 5. Compute Firing Rate (Hz)
        % Rate = Spikes / Time
        % We use max(..., eps) to prevent division by zero in empty bins
        phase_FR(nChannel,:) = spike_counts ./ max(time_occupancy, eps) ./length(V1_clusters.cluster_id);
        toc
    end

    peak_bins = mean(phase_FR(:,centers>-pi/2 & centers<pi/2),2,'omitnan');
    trough_bins =  mean(phase_FR(:,centers>pi/2 | centers<-pi/2),2,'omitnan');


    %%%% Save downsampled LFP from key channels

    options.importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP (in case of NP1)
    [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(options,[]);
    if isfield(options,'probe_hemisphere')
        LFP(nprobe).probe_hemisphere = options.probe_hemisphere;
        LFP(nprobe).probe_id = options.probe_id;
    else
        LFP(nprobe).probe_id = options.probe_id;
    end

    LFP(nprobe).tvec = tvec;


    % BEST V1 LFP channels for SO phase
    V1_channels=unique_selected_channels(good_V1_channels);
    unique_shanks = unique(chan_config.Shank(ismember(chan_config.Channel,V1_channels)));
    mean_shank_LFP=[];
    Region = 'best_SO_V1';

    best_V1_channel_index = [];
%     trough_peak_ratio=[];
    for nShank = 1:length(unique_shanks)
        channel_id = find(ismember(good_V1_channels,good_V1_channels(chan_config.Shank(ismember(chan_config.Channel,unique_selected_channels(good_V1_channels)))==unique_shanks(nShank))));
        [temp,index] = max(((trough_bins(channel_id)-peak_bins(channel_id))./(trough_bins(channel_id)+peak_bins(channel_id))));
        best_V1_channel_index(nShank) = good_V1_channels(channel_id(index));
%         trough_peak_ratio(nShank) = temp;

        LFP(nprobe).(sprintf('%s_trough_peak_ratio',Region))(nShank)=temp;
        LFP(nprobe).(sprintf('%s_phase_FR',Region))(nShank,:) =phase_FR(channel_id(index),:);
    end

    LFP(nprobe).(Region) = raw_LFP(best_V1_channel_index,:);
    LFP(nprobe).(sprintf('%s_shank_id',Region)) = unique_shank_id(best_V1_channel_index); % only avaliable shanks
    LFP(nprobe).(sprintf('%s_channel',Region))= unique_selected_channels(best_V1_channel_index);

    [xcoord,sorted_order] = sort(chan_config.Ks_xcoord(ismember(chan_config.Channel,unique_selected_channels(best_V1_channel_index))));
    ycoord = chan_config.Ks_ycoord(ismember(chan_config.Channel,unique_selected_channels(best_V1_channel_index)))';
    ycoord=ycoord(sorted_order);

    LFP(nprobe).(sprintf('%s_depth',Region))=ycoord';
    LFP(nprobe).(sprintf('%s_xcoord',Region))=xcoord';
    LFP(nprobe).(sprintf('%s_power',Region))=power{nprobe}(best_V1_channel_index,:);

    %         LFP(nprobe).(sprintf('%s_power',Region))=power{nprobe}(best_V1_shank_channel_id,:);


    % Save average V1 LFP per shank
    V1_channels=unique_selected_channels(good_V1_channels);
    unique_shanks = unique(chan_config.Shank(ismember(chan_config.Channel,V1_channels)));
    mean_shank_LFP=[];
    Region = 'average_V1';
    for nShank = 1:length(unique_shanks)
        channel_id = ismember(unique_selected_channels,V1_channels(chan_config.Shank(ismember(chan_config.Channel,V1_channels))==unique_shanks(nShank)));

        LFP(nprobe).(Region)(nShank,:) = mean(raw_LFP(channel_id,:));
        LFP(nprobe).(sprintf('%s_shank_id',Region))(nShank) = unique_shanks(nShank); % only avaliable shanks
    end
    LFP(nprobe).(sprintf('%s_channel',Region)) = V1_channels; % Record all good V1 channels;
    LFP(nprobe).(sprintf('%s_depth',Region)) = chan_config.Ks_ycoord(ismember(chan_config.Channel,V1_channels)); % Record all good V1 channels;
    LFP(nprobe).(sprintf('%s_xcoord',Region)) = chan_config.Ks_xcoord(ismember(chan_config.Channel,V1_channels)); % Record all good V1 channels;


%     if ~isempty(SWS)
%         % Save best SW V1 channels (here is equivalent of L5 based on high freq power)
%         Region = 'best_V1';
%         LFP(nprobe).(Region) = raw_LFP(best_V1_high_power_shank_channel_id,:);
%         LFP(nprobe).(sprintf('%s_shank_id',Region)) = unique_shank_id(best_V1_high_power_shank_channel_id); % only avaliable shanks
%         LFP(nprobe).(sprintf('%s_channel',Region))= unique_selected_channels(best_V1_high_power_shank_channel_id);
% 
%         [xcoord,sorted_order] = sort(chan_config.Ks_xcoord(ismember(chan_config.Channel,unique_selected_channels(best_V1_high_power_shank_channel_id))));
%         ycoord = chan_config.Ks_ycoord(ismember(chan_config.Channel,unique_selected_channels(best_V1_high_power_shank_channel_id)))';
%         ycoord=ycoord(sorted_order);
% 
%         LFP(nprobe).(sprintf('%s_depth',Region))=ycoord';
%         LFP(nprobe).(sprintf('%s_xcoord',Region))=xcoord';
%         LFP(nprobe).(sprintf('%s_power',Region))=power{nprobe}(best_V1_high_power_shank_channel_id,:);
% %         LFP(nprobe).(sprintf('%s_power',Region))=power{nprobe}(best_V1_shank_channel_id,:);
%     end

    % Save best V1 channels based on high freq power (eqivalent of L5 if channel didn't move much)
    Region = 'best_V1_high_freq';
    LFP(nprobe).(Region) = raw_LFP(best_V1_high_power_shank_channel_id,:);
    LFP(nprobe).(sprintf('%s_shank_id',Region)) = unique_shank_id(best_V1_high_power_shank_channel_id); % only avaliable shanks
    LFP(nprobe).(sprintf('%s_channel',Region))= unique_selected_channels(best_V1_high_power_shank_channel_id);

    [xcoord,sorted_order] = sort(chan_config.Ks_xcoord(ismember(chan_config.Channel,unique_selected_channels(best_V1_high_power_shank_channel_id))));
    ycoord = chan_config.Ks_ycoord(ismember(chan_config.Channel,unique_selected_channels(best_V1_high_power_shank_channel_id)))';
    ycoord=ycoord(sorted_order);

    LFP(nprobe).(sprintf('%s_depth',Region))=ycoord';
    LFP(nprobe).(sprintf('%s_xcoord',Region))=xcoord';
    LFP(nprobe).(sprintf('%s_power',Region))=power{nprobe}(best_V1_high_power_shank_channel_id,:);



    % Save best HPC channels based on ripple power
    Region = 'best_HPC';
    LFP(nprobe).(Region) = raw_LFP(best_HPC_shank_channel_id,:);
    LFP(nprobe).(sprintf('%s_shank_id',Region)) = unique_shank_id(best_HPC_shank_channel_id); % only avaliable shanks
    LFP(nprobe).(sprintf('%s_channel',Region))= unique_selected_channels(best_HPC_shank_channel_id);

    [xcoord,sorted_order] = sort(chan_config.Ks_xcoord(ismember(chan_config.Channel,unique_selected_channels(best_HPC_shank_channel_id))));
    ycoord = chan_config.Ks_ycoord(ismember(chan_config.Channel,unique_selected_channels(best_HPC_shank_channel_id)))';
    ycoord=ycoord(sorted_order);

    LFP(nprobe).(sprintf('%s_depth',Region))=ycoord';
    LFP(nprobe).(sprintf('%s_xcoord',Region))=xcoord';

    LFP(nprobe).(sprintf('%s_depth',Region))=chan_config.Ks_ycoord(ismember(chan_config.Channel,unique_selected_channels(best_HPC_shank_channel_id)))';
    LFP(nprobe).(sprintf('%s_xcoord',Region))=chan_config.Ks_xcoord(ismember(chan_config.Channel,unique_selected_channels(best_HPC_shank_channel_id)))';
    LFP(nprobe).(sprintf('%s_power',Region))=power{nprobe}(best_HPC_shank_channel_id,:);


    all_fields(contains(all_fields,'V1_sparse'))=[];%remove
    all_fields(contains(all_fields,'HPC_sparse'))=[];%remove
    % Save other region-specific channels (based on checkerboard channel characterisation)
    for nregion = 1:length(all_fields)
        if sum(channel_regions == nregion)>0
            LFP(nprobe).(all_fields{nregion}) = raw_LFP(ismember(unique_selected_channels,selected_channels(channel_regions == nregion)),:);
            LFP(nprobe).(sprintf('%s_shank_id',all_fields{nregion})) = shank_id(channel_regions == nregion); % only avaliable shanks
            LFP(nprobe).(sprintf('%s_channel',all_fields{nregion})) = selected_channels(channel_regions == nregion);
            LFP(nprobe).(sprintf('%s_depth',all_fields{nregion})) = chan_config.Ks_ycoord(selected_channels(channel_regions == nregion));
            LFP(nprobe).(sprintf('%s_xcoord',all_fields{nregion})) = chan_config.Ks_xcoord(selected_channels(channel_regions == nregion));
            LFP(nprobe).(sprintf('%s_power',all_fields{nregion})) = power{nprobe}(ismember(unique_selected_channels,selected_channels(channel_regions == nregion)),:);
        else
            LFP(nprobe).(all_fields{nregion}) = [];
            LFP(nprobe).(sprintf('%s_shank_id',all_fields{nregion})) = []; % only avaliable shanks
            LFP(nprobe).(sprintf('%s_channel',all_fields{nregion})) = [];
            LFP(nprobe).(sprintf('%s_depth',all_fields{nregion})) = [];
            LFP(nprobe).(sprintf('%s_xcoord',all_fields{nregion})) = [];
            LFP(nprobe).(sprintf('%s_power',all_fields{nregion})) = [];
        end
        %                 LFP(nprobe).(sprintf('%s_depth',all_fields{nregion})) = chan_config.Ks_ycoord(selected_channels(channel_regions == nregion));
    end

end
if contains(stimulus_name,'Masa2tracks')
    % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_slow_wave_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'slow_waves')
    % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('behavioural_state%s.mat',erase(stimulus_name,'Masa2tracks'))),'behavioural_state')
    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name,'Masa2tracks'))),'PSD')
    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name,'Masa2tracks'))),'LFP','-v7.3')
else

    % save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'),'slow_waves')
    % save(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state.mat'),'behavioural_state')
    save(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'),'PSD','power')% save PSD for the sleep session
    % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP','-v7.3')
    save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP','-v7.3')
end


close all