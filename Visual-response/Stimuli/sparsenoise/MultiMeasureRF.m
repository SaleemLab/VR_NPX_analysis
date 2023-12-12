%% 
function m = MultiMeasureRF(ES,m,selected)
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

h=figure;
for icell=selected
    figure(icell+100)
    subplot(4,1,1)
    imagesc(m.STA(cropX,cropY,icell))
    %         imagesc(m.ONMap(:,:,icell))
    title({['clust' num2str(icell) ' ID' num2str(ES.SpikeInfo{5,icell+1})],'STA'})% ' ' num2str(m.ONdelays.max(icell))
    colormap gray
%     axis square off
    box on
    hcb=colorbar;
    %         title(hcb,'hz')
    plotspec(ES,cropX,cropY);
    k=k+1;
end

for icell=selected
    figure(icell+100)
    subplot(4,1,2)
    imagesc(m.STC(cropX,cropY,icell))
    %         imagesc(m.ONMap(:,:,icell))
    title('STC')% ' ' num2str(m.ONdelays.max(icell))
    colormap gray
%     axis square off
    box on
    hcb=colorbar;
    plotspec(ES,cropX,cropY);
    k=k+1;
end

for icell=selected
    figure(icell+100)
    subplot(4,1,3)
    imagesc(m.OFFMap(cropX,cropY,icell))
    %         imagesc(m.ONMap(:,:,icell))
    title('OFF map')% ' ' num2str(m.ONdelays.max(icell))
    colormap gray
%     axis square off
    box on
    hcb=colorbar;
    title(hcb,'hz')
    plotspec(ES,cropX,cropY);
    k=k+1;
end

for icell=selected
    figure(icell+100)
    set(gcf, 'position', [360   201   178   497]);
    subplot(4,1,4)
    imagesc(m.ONMap(cropX,cropY,icell))
    %         imagesc(m.ONMap(:,:,icell))
    title('ON map')% ' ' num2str(m.ONdelays.max(icell))
    colormap gray
%     axis square off
    box on
    hcb=colorbar;
    title(hcb,'hz')
    plotspec(ES,cropX,cropY);
    k=k+1;
end

if isfield(ES.MetaData,'FileName')
    [~,name,~] = fileparts(ES.MetaData.FileName);
else
    [~,name,~] = fileparts(ES.MetaData.OutFileName);
end
annotation(h,'textbox',[0.45 0.95 0.04 0.04],'String',{name},'FitBoxToText','on', 'FontSize',14,'LineStyle','none', 'Interpreter','none');

end

% function plotspec(ES,cropX,cropY)
% line([cropX(1)+0.5 cropX(1)+0.5],[0  cropY(end)+0.5],'Color','blue','LineStyle','-','LineWidth',1)
% line([cropX(1)+0.5 cropX(end)+0.5],[(((ES.stim_centre_pos(2)+ES.stim_dims(2)/2))/ES.stim_dims(2))*ES.SN_gridDims(2) (((ES.stim_centre_pos(2)+ES.stim_dims(2)/2))/ES.stim_dims(2))*ES.SN_gridDims(2)]+0.5,'Color','blue','LineStyle','-','LineWidth',1)
% set(gca,...
%     'ylim',[0+0.5 length(cropX)+0.5],...
%     'YTick',[1 cropX(end)],...
%     'YTickLabel',[ES.stim_centre_pos(2)+ES.stim_dims(2)/2    ES.stim_centre_pos(2)-ES.stim_dims(2)/2],...
%     'xlim',[0+0.5 length(cropY)+0.5],...
%     'XTick',[1 cropY(end)],...
%     'XTickLabel',[ES.stim_centre_pos(1)-ES.stim_dims(1)/2    ES.stim_centre_pos(1)+ES.stim_dims(1)/2])
% end