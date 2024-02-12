function [updated_channels] = update_best_channels(options,chan_config)

options.importMode = 'KS';
[file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,[]);% Since it is LF

if ~contains(imecMeta.acqApLfSy,'384,0') % NPX2 only has AP but NPX1 has AP and LF
    options.importMode = 'LF';
    [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,[]);% Since it is LF
    probe_type = 1;
else
    probe_type = 2;
end

shanks_avaliable = unique(chan_config.Shank);
columns_avaliable = unique(chan_config.Ks_xcoord);
T.shanks_avaliable = shanks_avaliable;
T.columns_avaliable = columns_avaliable;

promp = sprintf('This is a NPX%i\n with shank(s):[%s]\n and unique column(s):[%s]\n',probe_type,join(string(T.shanks_avaliable),' '),join(string(T.columns_avaliable),' '))
updated_channels = [];
% shanks_to_process = input('What shank(s) do you want to analyse? e.g. If you want to analyse shank 1 and 3, input [1 3]:  ');
columns_to_process = input('What columns(s) do you want to analyse? e.g. If you want to analyse 1st and 3rd column, input [1 3]:  ');
regions_to_process = input('What regions do you want to process? (First one should be surface) e.g.{surface,L4,L5,CA1} (string inside each cell) ');


for nregion = 1:length(regions_to_process)
    updated_channels.([regions_to_process{nregion},'_channel'])(1:length(columns_avaliable)) =nan;
    updated_channels.([regions_to_process{nregion},'_depth'])(1:length(columns_avaliable)) =nan;
end
updated_channels.xcoord =columns_avaliable;

for col = columns_to_process
    nshank = unique(chan_config.Shank( chan_config.Ks_xcoord  ==  columns_avaliable(col)));
    sprintf('This is Shank %i column %i (X coord %i micron)',nshank,col,columns_avaliable(col))
    openfig(fullfile(options.ANALYSIS_DATAPATH,sprintf('%s %s Checkerboard event (all filtered) probe %i X coord %i.fig',options.SUBJECT,options.SESSION,options.probe_no,columns_avaliable(col))))
     datacursormode('on');
     openfig(fullfile(options.ANALYSIS_DATAPATH,sprintf('%s %s cluster density probe %i X coord %i.fig',options.SUBJECT,options.SESSION,options.probe_no,columns_avaliable(col))))
     datacursormode('on');

    channel_this_column = find(chan_config.Ks_xcoord == columns_avaliable(col));
    updated_channels.xcoord(col) =columns_avaliable(col);
    for nregion = 1:length(regions_to_process)
        confirm = 'N';
        while sum(strncmp(confirm,'y',1))==0
            best_depths(nregion) = input(sprintf('What is the best depth for %s (nan if not applicable):  ',regions_to_process{nregion}));
            confirm = input('Do you confirm your entry? (y/n)  ','s');
        end

        [~,idx] = min(abs(best_depths(nregion)-chan_config.Ks_ycoord(channel_this_column)));
        best_depths(nregion) = chan_config.Ks_ycoord(channel_this_column(idx));
        best_channels(nregion) =  chan_config.Channel(channel_this_column(idx));
        updated_channels.([regions_to_process{nregion},'_channel'])(col) =best_channels(nregion);
        updated_channels.([regions_to_process{nregion},'_depth'])(col) =best_depths(nregion);
    end
    close all
end


% first_in_brain_channel = find(sorted_config.Channel == best_channels.first_in_brain_channel);
% best_L4_channel =  find(sorted_config.Channel == best_channels.L4_channel);
% best_L5_channel =  find(sorted_config.Channel == best_channels.L5_channel);
% 
% best_CA1_ripple_channel= [];
% if ~isempty(best_channels.CA1_channel)
%     best_CA1_ripple_channel =  find(sorted_config.Channel == best_channels.CA1_channel);
% end
% 
% % Quick plotting of PSD
% %     colour_line= {'k','r','m','b','c','g','y'};
% colour_line= {[177,0,38]/256,[213,62,79]/256,[252,141,89]/256,...
%     [153,213,148]/256,[26,152,80]/256,[66,146,198]/256,[8,69,148]/256};
% freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','300-600 Hz'};
% selected_frequency = [2 6 7];
% figure
% %     subplot(1,4,2)
% 
% for n = selected_frequency
%     pl(n) = plot(power(:,n)./max(power(:,n)),sorted_config.Ks_ycoord','Color',colour_line{n})
% 
%     %     p(n) = plot(power_differnece,1:96,colour_line{n})
%     hold on
% end
% hold on
% plot([0 1],[sorted_config.Ks_ycoord(first_in_brain_channel) sorted_config.Ks_ycoord(first_in_brain_channel)],'--k','LineWidth',2)
% plot([0 1],[sorted_config.Ks_ycoord(best_L4_channel) sorted_config.Ks_ycoord(best_L4_channel)],'--b','LineWidth',2)
% plot([0 1],[sorted_config.Ks_ycoord(best_L5_channel) sorted_config.Ks_ycoord(best_L5_channel)],'--c','LineWidth',2)
% if ~isempty(best_CA1_ripple_channel)
%     plot([0 1],[sorted_config.Ks_ycoord(best_CA1_ripple_channel) sorted_config.Ks_ycoord(best_CA1_ripple_channel)],'--r','LineWidth',2)
% end
% 
% % legend('1-3','4-12','30-60','60-100','125-300')
% legend([pl(selected_frequency)],{freq_legends{selected_frequency}});
% ylim([0 4000])
% openfig(sprintf('Checkerboard event (all filtered) probe %i.fig',options.probe_no))
% openfig(sprintf('%s %s Checkerboard event (all filtered) probe %i.fig',options.SUBJECT,options.SESSION,options.probe_no))
% yyaxis left
% 
% disp(['Select four channels manually: 1. first channel that enters the brain 2. L4 channel (putative)...' ...
%     ' 3. L5 channel (based on high frequency power) 4. CA1 channel (based on low theta and high ripple power). ...' ...
%     'Please press ENTRE key after selecting all four points'])
% [ X , updated_channels ] = getpts;
% 
% % find the closest channels
% for nchannel = 1:length(updated_channels)
%     [~,index]= min(abs(sorted_config.Ks_ycoord- updated_channels(nchannel)));
%     updated_channels(nchannel) = sorted_config.Channel(index);
% end

% best_channels.first_in_brain_channel = updated_channels(1);
% best_channels.L4_channel = updated_channels(2);
% best_channels.L5_channel = updated_channels(3);
% best_channels.CA1_channel = updated_channels(4);

end
