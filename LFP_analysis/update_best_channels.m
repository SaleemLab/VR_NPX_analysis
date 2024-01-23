function [updated_channels] = update_best_channels(options,sorted_config)

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
openfig(sprintf('%s %s Checkerboard event (all filtered) probe %i.fig',options.SUBJECT,options.SESSION,options.probe_no))
yyaxis left

disp(['Select four channels manually: 1. first channel that enters the brain 2. L4 channel (putative)...' ...
    ' 3. L5 channel (based on high frequency power) 4. CA1 channel (based on low theta and high ripple power). ...' ...
    'Please press ENTRE key after selecting all four points'])
[ X , updated_channels ] = getpts;

% find the closest channels
for nchannel = 1:length(updated_channels)
    [~,index]= min(abs(sorted_config.Ks_ycoord- updated_channels(nchannel)));
    updated_channels(nchannel) = sorted_config.Channel(index);
end

% best_channels.first_in_brain_channel = updated_channels(1);
% best_channels.L4_channel = updated_channels(2);
% best_channels.L5_channel = updated_channels(3);
% best_channels.CA1_channel = updated_channels(4);

end
