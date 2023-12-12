
%% 
function m = MultiMeasureRF_MUA_snooth(ES,m)
k=1;
%selected=[1 7 15 16 17 22 27 28];%[5,15,32,34,35,45,53];%[2,4,5,36,43,47,55];%[1,2,6,7,10,24,26];%[5,15,32,34,35,45,53];%[1,3,12,15,18,23];%[6,14, 19, 23, 37, 40];%[6,15,20,22,23,24]

addpath( 'X:\Archive - saleemlab\Code\BonVision\Bonsai_analysis');
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

n = 256; gamma = 1;
cm = RedWhiteBlue;
% colormap from white to red
Nmap = 100;
cMin = [1 1 1];
cMax = [1 0 0];
cMapRed = zeros(Nmap,3);
for i = 1:Nmap;
    cMapRed(i,:) = cMin*(Nmap - i)/(Nmap - 1) + cMax*(i - 1)/(Nmap - 1);
end
cMax = [0 0 1];
cMapblu = zeros(Nmap,3);
for i = 1:Nmap;
    cMapblu(i,:) = cMin*(Nmap - i)/(Nmap - 1) + cMax*(i - 1)/(Nmap - 1);
end

subplot(2,2,1)
imagesc(sum(m.OFFMap(cropX,cropY,ne_i & th_i),3)/sum(ne_i & th_i))
title('OFF map for MUA')% ' ' num2str(m.ONdelays.max(icell))
colormap(gca,cm);
hcb3=colorbar;
box on
title(hcb3,'Hz')
plotspec(ES,cropX,cropY);

subplot(2,2,2)
imagesc(sum(m.ONMap(cropX,cropY,ne_i & th_i),3)/sum(ne_i & th_i))
title('ON map for MUA')% ' ' num2str(m.ONdelays.max(icell))
colormap(gca,cm);
hcb4=colorbar;
box on
title(hcb4,'Hz')
plotspec(ES,cropX,cropY);

scal_f = 10; sigma = 5;
I = sum(m.OFFMap(cropX,cropY,ne_i & th_i),3)/sum(ne_i & th_i);
J = imresize(I,scal_f);
J = imgaussfilt(J,sigma);
J = (J/max(J(:)));
subplot(2,2,3)
imagesc(J)
title('OFF map for MUA smoothed')% ' ' num2str(m.ONdelays.max(icell))
colormap(gca,cMapblu)
%colormap(gca,jet)
hcb3=colorbar;
box on
title(hcb3,'Hz')
plotspecSmoothed(ES,cropX,cropY,scal_f);

I = sum(m.ONMap(cropX,cropY,ne_i & th_i),3)/sum(ne_i & th_i);
K = imresize(I,scal_f);
K = imgaussfilt(K,sigma);
K = (K/max(K(:)));
subplot(2,2,4)
imagesc(K)
title('ON map for MUA smoothed')% ' ' num2str(m.ONdelays.max(icell))
colormap(gca,cMapRed)
%colormap(gca,jet)
hcb4=colorbar;
box on
title(hcb4,'Hz')
plotspecSmoothed(ES,cropX,cropY,scal_f);
    
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

function plotspecSmoothed(ES,cropX,cropY,scal_f)
% line([0+0.5 0+0.5],[0  cropY(end)*scal_f+0.5],'Color','blue','LineStyle','-','LineWidth',2)
% line([0+0.5 length(cropX)*scal_f+0.5],[(((ES.stim_centre_pos(2)+ES.stim_dims(2)/2))/ES.stim_dims(2))*length(cropY)*scal_f (((ES.stim_centre_pos(2)+ES.stim_dims(2)/2))/ES.stim_dims(2))*length(cropY)*scal_f]+0.5,'Color','blue','LineStyle','-','LineWidth',1)
StepsX = ((abs(ES.stim_centre_pos(1)-ES.stim_dims(1)/2)+abs(ES.stim_centre_pos(1)+ES.stim_dims(1)/2))/4);
XTickLabels = min(ES.stim_centre_pos(1)-ES.stim_dims(1)/2,ES.stim_centre_pos(1)+ES.stim_dims(1)/2):StepsX: ...
                max(ES.stim_centre_pos(1)-ES.stim_dims(1)/2,ES.stim_centre_pos(1)+ES.stim_dims(1)/2);
StepsY = ((abs(ES.stim_centre_pos(2)-ES.stim_dims(2)/2)+abs(ES.stim_centre_pos(2)+ES.stim_dims(2)/2))/4);
YTickLabels = min(ES.stim_centre_pos(2)-ES.stim_dims(2)/2,ES.stim_centre_pos(2)+ES.stim_dims(2)/2):StepsY: ...
                max(ES.stim_centre_pos(2)-ES.stim_dims(2)/2,ES.stim_centre_pos(2)+ES.stim_dims(2)/2);
set(gca,...
    'ylim',[0 length(cropX)*scal_f+0.5],...
    'YTick',[0:length(cropX)*scal_f*0.25:length(cropX)*scal_f],...
    'YTickLabel',sort(YTickLabels,'descend'),...
    'xlim',[0 length(cropX)*scal_f],...
    'XTick',[0:length(cropX)*scal_f*0.25:length(cropX)*scal_f],...
    'XTickLabel',XTickLabels)
xlabel('Azimuth (^o)');
ylabel('Elevation (^o)');
set(gca,'FontSize',13,'TickDir','out','FontName','Arial');
grid on
end

%% plot some example sparse noise signals
% Frames = round(rand(3,1)*3000);
% 
% figure
% for i = 1:length(Frames)
%     subplot(1,3,i)
%     imagesc(squeeze(ES.SN_sequence(:,:,Frames(i))))
%     colormap(gca,'gray')
%     box off       
%     set(gca,...
%         'YTick',[0.5:(size(ES.SN_sequence,1)*0.25):size(ES.SN_sequence,2)+0.5],...
%         'YTickLabel',[90 60 30 0 -30],...
%         'XTick',[0.5:(size(ES.SN_sequence,1)*0.25):size(ES.SN_sequence,1)+0.5],...
%         'XTickLabel',[0 30 60 90 120]);
%     xlabel('Azimuth (^o)');
%     ylabel('Elevation (^o)');
%     set(gca,'FontSize',13,'TickDir','out','FontName','Arial');    
%     %grid on
% end   
%     


