%% Reverse correlation, STC (no ON OFF stim discrimination, easy spike triggered covariance)
function m = ReverseCorrelationSTC(ES,m,plotSingleUnits)
addpath( 'X:\Archive - saleemlab\Code\BonVision\Bonsai_analysis');
% initialize
numframes=length(ES.SN_onsets);
delays=-10:10:120; % -10 to 90ms delays

STC=zeros(ES.SN_gridDims(1),ES.SN_gridDims(2),length(delays),length(ES.StimSpiketimes));
stim_frames=2*abs(ES.SN_sequence/max(ES.SN_sequence(:))-0.5);

% define whether to crop image
if isfield(ES,'SN_square_crop')
    cropX=ES.SN_square_crop(1,:);
    cropY=ES.SN_square_crop(2,:);
else
    cropX=1:ES.SN_gridDims(1);
    cropY=1:ES.SN_gridDims(2);
end

for icell=1:length(ES.StimSpiketimes)
    
    totalSpikeCount = 0;
    for delayInd=1:length(delays)
        for frame=1:numframes
            STC(:,:,delayInd,icell)= STC(:,:,delayInd,icell) + stim_frames(:,:,frame)*ES.frameSpikeCount{icell}(frame,delayInd);
            totalSpikeCount = totalSpikeCount +  ES.frameSpikeCount{icell}(frame,delayInd);
        end
    end
    
    % pick the best delay for this cluster (delay that maximizes spatial variance)
    [maxVar, maxInd, vars] = maxVarMap(STC(:,:,:,icell));
    
    m.STCdelays.delays= delays;
    m.STCdelays.max(icell)=delays(maxInd);
    m.STCvars.vars(icell,:)= vars;
    m.STCvars.max(icell)= maxVar;
    m.STC(:,:,icell) = STC(:,:,maxInd,icell)./totalSpikeCount;
    
end

if plotSingleUnits
    figure
    % define nr of subplots
    SProws = floor(length(ES.frameSpikeCount)/8);
    if mod(length(ES.frameSpikeCount),8)
        SProws = SProws+1;
    end
    k=1;
    for icell=1:length(ES.frameSpikeCount)%[6,15,20,22,23,24]
        %         subplot(1,6,k)
        subplot(SProws,8,icell)
        imagesc(m.STC(cropX,cropY,icell))
        title(['STC clust ' num2str(icell)])
        colormap gray
        %axis off
        axis square
        box on
        hcb=colorbar;
        title(hcb,'gray')
        plotspec(ES,cropX,cropY);
        k=k+1;
    end
end

end

% function plotspec(ES,cropX,cropY)
% line([0+0.5 0+0.5],[0  ES.SN_gridDims(2)+0.5],'Color','blue','LineStyle','-','LineWidth',1)
% line([0+0.5 ES.SN_gridDims(1)+0.5],[(((ES.stim_centre_pos(2)+ES.stim_dims(2)/2))/ES.stim_dims(2))*ES.SN_gridDims(2) (((ES.stim_centre_pos(2)+ES.stim_dims(2)/2))/ES.stim_dims(2))*ES.SN_gridDims(2)]+0.5,'Color','blue','LineStyle','-','LineWidth',1)
% set(gca,...
%     'ylim',[cropY(1)-0.5 cropY(end)+0.5],...
%     'YTick',[1 cropY(end)],...
%     'YTickLabel',[ES.stim_centre_pos(2)+ES.stim_dims(2)/2    ES.stim_centre_pos(2)-ES.stim_dims(2)/2],...
%     'xlim',[cropX(1)-0.5 cropX(end)+0.5],...
%     'XTick',[1 cropX(end)],...
%     'XTickLabel',[ES.stim_centre_pos(1)-ES.stim_dims(1)/2    ES.stim_centre_pos(1)+ES.stim_dims(1)/2])
% end
