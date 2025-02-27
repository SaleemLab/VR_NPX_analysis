function  [region_channels,best_channels_all_shanks_avaliable] = determine_region_channels(best_channels,options,varargin)
% This is a function to quickly determine what is best channel for each
% shank for each regions

% Default values
p = inputParser;
addParameter(p,'region',[],@isstr) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn))
addParameter(p,'groups','by shank',@isstr) % selects channels based on shanks (rather than all columns)
addParameter(p,'clusters',[],@isstruct) % selects channels based on shanks (rather than all columns)
% addParameter(p,'groups','by shank',@isstr) % selects channels based on shanks (rather than all columns)

% assign parameters (either defaults or given)
parse(p,varargin{:});
region = [p.Results.region,'_depth'];
groups = p.Results.groups;
clusters = p.Results.clusters;

% grab channel info
options.importMode = 'KS';
[file_to_use imecMeta chan_config ~] = extract_NPX_channel_config(options,[]);% Since it is LF

all_fields = fieldnames(best_channels);

% shank_id_avaliable = ceil(best_channels.xcoord./250); % based on xcoord in best_channels
shank_id_avaliable = ceil(unique(chan_config.Ks_xcoord)./250)'; % based on xcoord in chan_config


region_channels = [];

if isempty(region)

elseif sum(contains(all_fields,region)) == 1

    best_depths_this_region = best_channels.(all_fields{contains(all_fields,region)});
    best_depths_this_region(best_depths_this_region==0) = nan;
    best_depths_all_shanks_avaliable = [];

    if contains(groups,'by shank') % Find depth for each shank
        for nShank = 1:length(unique(shank_id_avaliable))
            best_depths_all_shanks_avaliable = [best_depths_all_shanks_avaliable mean(best_depths_this_region(shank_id_avaliable == nShank),'omitnan')];
        end

        best_depths_all_shanks_avaliable = round(best_depths_all_shanks_avaliable);

        shanks_avaliable = unique(shank_id_avaliable);
        %     elseif contains(groups,'all')

    elseif contains(groups,'by probe') % Find median depth across all shanks
        for nShank = 1:length(unique(shank_id_avaliable))
            best_depths_all_shanks_avaliable = [best_depths_all_shanks_avaliable mean(best_depths_this_region(shank_id_avaliable == nShank),'omitnan')];
        end

        best_depths_all_shanks_avaliable = round(median(best_depths_this_region,'omitnan'));

        shanks_avaliable = 1;% just one
    end



    for nShank = shanks_avaliable
        depth_this_shank = best_depths_all_shanks_avaliable(nShank);

        if contains(groups,'by shank')
            channels_this_shank = (chan_config.Shank == nShank);
        elseif contains(groups,'by probe')
            channels_this_shank = ones(length(chan_config.Shank),1);
        end

        if contains(region,'L5')
            all_channels_this_shank = chan_config.Channel((chan_config.Ks_ycoord>depth_this_shank-100 ...
                &chan_config.Ks_ycoord<depth_this_shank+100) & channels_this_shank) ;

        elseif contains(region,'L4')
            all_channels_this_shank = chan_config.Channel((chan_config.Ks_ycoord>depth_this_shank-50 ...
                &chan_config.Ks_ycoord<depth_this_shank+50) & channels_this_shank) ;

        elseif contains(region,'CA1')
            all_channels_this_shank = chan_config.Channel((chan_config.Ks_ycoord>depth_this_shank-100 ...
                &chan_config.Ks_ycoord<depth_this_shank+100) & channels_this_shank) ;

        elseif contains(region,'surface')
            all_channels_this_shank = chan_config.Channel((chan_config.Ks_ycoord>depth_this_shank-1 ...
                &chan_config.Ks_ycoord<depth_this_shank+1) & channels_this_shank) ;
            
        elseif contains(region,'MEC_entry') % putatively just taking 2000 mcrion from entry 
            all_channels_this_shank = chan_config.Channel((chan_config.Ks_ycoord>=0 ...
                &chan_config.Ks_ycoord<depth_this_shank) & channels_this_shank) ;

        elseif contains(region,'MEC_ripple') % for now just taking few channel 
            all_channels_this_shank = chan_config.Channel((chan_config.Ks_ycoord>depth_this_shank-500 ...
                &chan_config.Ks_ycoord<depth_this_shank+1) & channels_this_shank) ;

        elseif contains(region,'MEC_theta') % for now just taking few channel
            all_channels_this_shank = chan_config.Channel((chan_config.Ks_ycoord>depth_this_shank-1 ...
                &chan_config.Ks_ycoord<depth_this_shank+1) & channels_this_shank) ;

        elseif contains(region,'HVA') % for now just taking few channel for HVA
            all_channels_this_shank = chan_config.Channel((chan_config.Ks_ycoord>depth_this_shank-1 ...
                &chan_config.Ks_ycoord<depth_this_shank+1) & channels_this_shank) ;


        end

        region_channels = [region_channels; all_channels_this_shank];
    end


elseif contains(region,'V1')


    % if NPX2 with four shanks (with V1 and HPC recording), best_depths is Middle of top channel blocks
    if length(unique(chan_config.Shank))>1 && sum(contains(all_fields,'CA1'))>0 

%         best_depths_this_region = best_channels.surface_depth;
%         best_depths_this_region(best_depths_this_region==0) = nan;
        best_depths_all_shanks_avaliable = [];
        region_channels=[];

        for nShank = 1:length(unique(shank_id_avaliable))
            channels_this_shank = find(chan_config.Shank == nShank);
            channels_depth_this_shank = chan_config.Ks_ycoord(channels_this_shank);

            kClusters=kmeans(channels_depth_this_shank,2);
            if mean(channels_depth_this_shank(kClusters==1))>mean(channels_depth_this_shank(kClusters==2))
                % if mean ocation of cluster one is above cluster two, it is
                % Cortex.
                % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==1);
                V1_channel_id = find(kClusters==1);
            else
                % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==2);
                V1_channel_id = find(kClusters==2);
            end
            if contains(region,'sparse') % Sample sparse channels from all shanks of V1 (every 45 micron)
                selected_channels_ycoord = chan_config.Ks_ycoord(channels_this_shank(V1_channel_id));
                unique_ycoords = unique(selected_channels_ycoord);
                unique_ycoords=unique_ycoords(1:3:end);
                V1_channel_id = V1_channel_id(ismember(selected_channels_ycoord,unique_ycoords));
            end

            V1_channel_depths = channels_depth_this_shank(V1_channel_id);

            best_depths_all_shanks_avaliable = [best_depths_all_shanks_avaliable round(median(V1_channel_depths,'omitnan'))]; % for now, take the median channels of the shank as 'best channels'
            region_channels = [region_channels; channels_this_shank(V1_channel_id)];
        end

        if contains(groups,'by shank') % Find depth for each shank
            shanks_avaliable = unique(shank_id_avaliable);
            %     elseif contains(groups,'all')

        elseif contains(groups,'by probe') % Find median depth across all shanks
            best_depths_all_shanks_avaliable = round(median(best_depths_all_shanks_avaliable));

            shanks_avaliable = 1;% just one
        end

    else % if NPX1, best_depths is surface channel and then get 1mm below surface as putative V1 channels

        best_depths_this_region = best_channels.surface_depth;
        best_depths_this_region(best_depths_this_region==0) = nan;
        best_depths_all_shanks_avaliable = [];

        if contains(groups,'by shank') % Find depth for each shank
            for nShank = 1:length(unique(shank_id_avaliable))
                best_depths_all_shanks_avaliable = [best_depths_all_shanks_avaliable mean(best_depths_this_region(shank_id_avaliable == nShank),'omitnan')];
            end

            best_depths_all_shanks_avaliable = round(best_depths_all_shanks_avaliable);

            shanks_avaliable = unique(shank_id_avaliable);
            %     elseif contains(groups,'all')

        elseif contains(groups,'by probe') % Find median depth across all shanks
            for nShank = 1:length(unique(shank_id_avaliable))
                best_depths_all_shanks_avaliable = [best_depths_all_shanks_avaliable mean(best_depths_this_region(shank_id_avaliable == nShank),'omitnan')];
            end

            best_depths_all_shanks_avaliable = round(median(best_depths_this_region,'omitnan'));

            shanks_avaliable = 1;% just one
        end


        for nShank = shanks_avaliable
            depth_this_shank = best_depths_all_shanks_avaliable(nShank);

            if contains(groups,'by shank')
                channels_this_shank = (chan_config.Shank == nShank);
            elseif contains(groups,'by probe')
                channels_this_shank = ones(length(chan_config.Shank),1);
            end

            all_channels_this_shank = chan_config.Channel(chan_config.Ks_ycoord>depth_this_shank-1000 ...
                & channels_this_shank);

            region_channels = [region_channels; all_channels_this_shank];
        end
    end


elseif contains(region,'HPC')
    best_depths_this_region = best_channels.CA1_depth;
    best_depths_this_region(best_depths_this_region==0) = nan;
    best_depths_all_shanks_avaliable = [];


    if length(unique(chan_config.Shank))>1 && sum(contains(all_fields,'CA1'))>0

        %         best_depths_this_region = best_channels.surface_depth;
        %         best_depths_this_region(best_depths_this_region==0) = nan;
        best_depths_all_shanks_avaliable = [];
        region_channels=[];

        for nShank = 1:length(unique(shank_id_avaliable))
            channels_this_shank = find(chan_config.Shank == nShank);
            channels_depth_this_shank = chan_config.Ks_ycoord(channels_this_shank);

            kClusters=kmeans(channels_depth_this_shank,2);
            if mean(channels_depth_this_shank(kClusters==1))>mean(channels_depth_this_shank(kClusters==2))
                % if mean ocation of cluster one is above cluster two, it is
                % Cortex.
                % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==1);
                HPC_channel_id = find(kClusters==2);
            else
                % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==2);
                HPC_channel_id = find(kClusters==1);
            end
            if contains(region,'sparse') % Sample sparse channels from all shanks of V1 (every 45 micron)
                selected_channels_ycoord = chan_config.Ks_ycoord(channels_this_shank(HPC_channel_id));
                unique_ycoords = unique(selected_channels_ycoord);
                unique_ycoords=unique_ycoords(1:3:end);
                HPC_channel_id = HPC_channel_id(ismember(selected_channels_ycoord,unique_ycoords));
            end

            HPC_channel_depths = channels_depth_this_shank(HPC_channel_id);

            best_depths_all_shanks_avaliable = [best_depths_all_shanks_avaliable round(median(HPC_channel_depths,'omitnan'))]; % for now, take the median channels of the shank as 'best channels'
            region_channels = [region_channels; channels_this_shank(HPC_channel_id)];
        end

        if contains(groups,'by shank') % Find depth for each shank
            shanks_avaliable = unique(shank_id_avaliable);
            %     elseif contains(groups,'all')

        elseif contains(groups,'by probe') % Find median depth across all shanks
            best_depths_all_shanks_avaliable = round(median(best_depths_all_shanks_avaliable));

            shanks_avaliable = 1;% just one
        end

    else

        if contains(groups,'by shank') % Find depth for each shank
            for nShank = 1:length(unique(shank_id_avaliable))
                best_depths_all_shanks_avaliable = [best_depths_all_shanks_avaliable mean(best_depths_this_region(shank_id_avaliable == nShank),'omitnan')];
            end

            best_depths_all_shanks_avaliable = round(best_depths_all_shanks_avaliable);

            shanks_avaliable = unique(shank_id_avaliable);
            %     elseif contains(groups,'all')

        elseif contains(groups,'by probe') % Find median depth across all shanks

            best_depths_all_shanks_avaliable = round(median(best_depths_this_region,'omitnan'));

            shanks_avaliable = 1;% just one
        end


        for nShank = shanks_avaliable
            depth_this_shank = best_depths_all_shanks_avaliable(nShank);

            if contains(groups,'by shank')
                channels_this_shank = (chan_config.Shank == nShank);
            elseif contains(groups,'by probe')
                channels_this_shank = ones(length(chan_config.Shank),1);
            end

            all_channels_this_shank = chan_config.Channel(chan_config.Ks_ycoord>depth_this_shank-700 ...
                & chan_config.Ks_ycoord<depth_this_shank+100 ...
                & channels_this_shank);

            region_channels = [region_channels; all_channels_this_shank];
        end
    end

    if contains(groups,'by shank') % Find depth for each shank
        shanks_avaliable = unique(shank_id_avaliable);
        %     elseif contains(groups,'all')

    elseif contains(groups,'by probe') % Find median depth across all shanks
        best_depths_all_shanks_avaliable = round(median(best_depths_all_shanks_avaliable));

        shanks_avaliable = 1;% just one
    end
end

for nShank = 1:length(shanks_avaliable)
    channels_this_shank = find(chan_config.Shank == shanks_avaliable(nShank));
    [minValue,idx] = min(abs(best_depths_all_shanks_avaliable(nShank)-  chan_config.Ks_ycoord(channels_this_shank)));
    if minValue<200 % if more than 200 micron apart from the best channel (CA1,L4 or L5) for this shank (ignore that shank)
        best_channels_all_shanks_avaliable(nShank) = channels_this_shank(idx);
    else
        best_channels_all_shanks_avaliable(nShank) = nan;
    end
end

end