function  [region_channels,median_channels] = determine_region_channels(best_channels,options,varargin)
% This is a function to quickly determine what is best channel for each
% shank for each regions

% Default values
p = inputParser;
addParameter(p,'region',[],@isstr) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn))
addParameter(p,'groups','by shank',@isstr) % selects channels based on shanks (rather than all columns)

% assign parameters (either defaults or given)
parse(p,varargin{:});
region = [p.Results.region,'_depth'];
groups = p.Results.groups;


% grab channel info
options.importMode = 'KS';
[file_to_use imecMeta chan_config ~] = extract_NPX_channel_config(options,[]);% Since it is LF

all_fields = fieldnames(best_channels);

shank_id_avaliable = ceil(best_channels.xcoord./250); % based on xcoord in best_channels

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

        elseif contains(region,'MEC')
            all_channels_this_shank = [];

        end

        region_channels = [region_channels; all_channels_this_shank];
        median_channels(nShank) = find(chan_config.Channel,median(chan_config.Ks_ycoord(all_channels_this_shank)),'first');
    end


elseif contains(region,'V1')

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
        median_channels(nShank) = find(chan_config.Channel,median(chan_config.Ks_ycoord(all_channels_this_shank)),'first');
    end

elseif contains(region,'HPC')
    best_depths_this_region = best_channels.CA1_depth;
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
        median_channels(nShank) = find(chan_config.Channel,median(chan_config.Ks_ycoord(all_channels_this_shank)),'first');
    end

end

end