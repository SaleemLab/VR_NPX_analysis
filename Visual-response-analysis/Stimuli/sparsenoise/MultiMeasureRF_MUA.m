%% 
function m = MultiMeasureRF_MUA(ES,m)
k=1;
%selected=[1 7 15 16 17 22 27 28];%[5,15,32,34,35,45,53];%[2,4,5,36,43,47,55];%[1,2,6,7,10,24,26];%[5,15,32,34,35,45,53];%[1,3,12,15,18,23];%[6,14, 19, 23, 37, 40];%[6,15,20,22,23,24]
% 
% if ismac
%     ROOTPATH = '/Users/s.solomon/Filestore/Research2/ibn-vision/';
% else
%     if filesep == '\'
%         [~,hostName] = system('hostname'); hostName = hostName(1:end-1);
%         if ~strcmp(hostName, 'saleem12')
%             ROOTPATH = 'X:';
%         else
%             ROOTPATH = 'X:\ibn-vision';
%         end
%     else
%         ROOTPATH = '/mnt/pfs09/';
%     end
% end
% addpath(fullfile(ROOTPATH,'Archive - saleemlab','Code','BonVision','Bonsai_analysis'))

% initialize
numframes=length(ES.SN_onsets);
delays=-10:10:120; % -10 to 90ms delays

stim_frames=ES.SN_sequence/max(ES.SN_sequence(:));
numframes=length(ES.SN_onsets);

% define whether to crop image
if isfield(ES,'SN_square_crop')
    cropX=ES.SN_square_crop(1,:);
    cropY=ES.SN_square_crop(2,:);
else
    cropX=1:ES.SN_gridDims(1);
    cropY=1:ES.SN_gridDims(2);
end

% valid units
clear ne_i th_i
for i = 1:length(ES.StimSpiketimes)
    ne_i(i) = ~isempty(ES.StimSpiketimes{1,i}); % not empty array
    th_i(i) = (sum(ES.StimSpiketimes{1,i})/(max(ES.SN_onsets)-min(ES.SN_onsets)))>0.1; % mean FR>0.1 Hz;
end

h=figure
set(gcf, 'Position',  [100, 100, 1000, 900])

subplot(2,2,1)
imagesc(sum(m.STA(cropX,cropY,ne_i & th_i),3)/sum(ne_i & th_i));
title('STA for MUA')% ' ' num2str(m.ONdelays.max(icell))
colormap(gca,'gray')
hcb1=colorbar;
box on
plotspec(ES,cropX,cropY);

subplot(2,2,2)
imagesc(sum(m.STC(cropX,cropY,ne_i & th_i),3)/sum(ne_i & th_i))
title('STC for MUA')% ' ' num2str(m.ONdelays.max(icell))
colormap(gca,'gray')
hcb3=colorbar;
box on
plotspec(ES,cropX,cropY);

subplot(2,2,3)
imagesc(sum(m.OFFMap(cropX,cropY,ne_i & th_i),3)/sum(ne_i & th_i))
title('OFF map for MUA')% ' ' num2str(m.ONdelays.max(icell))
colormap(gca,flipud(gray))
hcb3=colorbar;
box on
title(hcb3,'Hz')
plotspec(ES,cropX,cropY);

subplot(2,2,4)
imagesc(sum(m.ONMap(cropX,cropY,ne_i & th_i),3)/sum(ne_i & th_i))
title('ON map for MUA')% ' ' num2str(m.ONdelays.max(icell))
colormap(gca,'gray')
hcb4=colorbar;
box on
title(hcb4,'Hz')
plotspec(ES,cropX,cropY);

if isfield(ES.MetaData,'FileName')
    [~,name,~] = fileparts(ES.MetaData.FileName);
else
    [~,name,~] = fileparts(ES.MetaData.OutFileName);
end
annotation(h,'textbox',[0.45 0.95 0.04 0.04],'String',{name},'FitBoxToText','on', 'FontSize',14,'LineStyle','none', 'Interpreter','none');

end

function plotspec(ES,cropX,cropY)
line([0+0.5 0+0.5],[0  cropY(end)+0.5],'Color','blue','LineStyle','-','LineWidth',2)
line([0+0.5 length(cropX)+0.5],[(((ES.stim_centre_pos(2)+ES.stim_dims(2)/2))/ES.stim_dims(2))*length(cropY) (((ES.stim_centre_pos(2)+ES.stim_dims(2)/2))/ES.stim_dims(2))*length(cropY)]+0.5,'Color','blue','LineStyle','-','LineWidth',1)
set(gca,...
    'ylim',[0+0.5 length(cropX)+0.5],...
    'YTick',[0.5 6.5 length(cropX)+0.5],...
    'YTickLabel',[(ES.stim_centre_pos(2)+ES.stim_dims(2)/2)  0  (ES.stim_centre_pos(2)-ES.stim_dims(2)/2)],...
    'xlim',[0+0.5 length(cropY)+0.5],...
    'XTick',[0.5 length(cropY)+0.5],...
    'XTickLabel',[(ES.stim_centre_pos(1)-ES.stim_dims(1)/2)    (ES.stim_centre_pos(1)+ES.stim_dims(1)/2)])
xlabel('Azimuth (^o)');
ylabel('Elevation (^o)');
set(gca,'FontSize',13,'TickDir','out','FontName','Arial');

end