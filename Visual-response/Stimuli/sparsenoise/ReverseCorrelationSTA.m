%% Reverse correlation, STA, spike triggered average
function m = ReverseCorrelationSTA(ES,plotSingleUnits)

% weighted sum STA
% initialize
%     STA=0.5*ones(p.GridNum,p.GridNum,length(chans));
%     STA=zeros(ES.SN_gridDims(1),ES.SN_gridDims(2),1,1);

delays=-10:10:120; % -10 to 90ms delays

STA=zeros(ES.SN_gridDims(1),ES.SN_gridDims(2),length(delays),length(ES.StimSpiketimes));
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

for icell=1:length(ES.StimSpiketimes)
    
    totalSpikeCount = 0;
    for delayInd=1:length(delays)
        for frame=1:numframes
                STA(:,:,delayInd,icell)= STA(:,:,delayInd,icell) + stim_frames(:,:,frame)*ES.frameSpikeCount{icell}(frame,delayInd);
            totalSpikeCount = totalSpikeCount +  ES.frameSpikeCount{icell}(frame,delayInd);
        end
    end
    
    % pick the best delay for this cluster (delay that maximizes spatial variance)
    [maxVar, maxInd, vars] = maxVarMap(STA(:,:,:,icell));
    
    m.STAdelays.delays= delays;
    m.STAdelays.max(icell)=delays(maxInd);
    m.STAvars.vars(icell,:)= vars;
    m.STAvars.max(icell)= maxVar;
    m.STA(:,:,icell) = STA(:,:,maxInd,icell)./totalSpikeCount;
    
end

if plotSingleUnits
    h=figure;
    % define nr of subplots
    SProws = floor(length(ES.frameSpikeCount)/8);
    if mod(length(ES.frameSpikeCount),8)
        SProws = SProws+1;
    end
    k=1;
    for icell=1:length(ES.frameSpikeCount)%[6,15,20,22,23,24]
        % subplot(1,6,k)
        subplot(SProws,8,icell)
        
        scal_f = 10; sigma = 3;
        I = m.STA(cropX,cropY,icell);
        J = imresize(I,scal_f);
        J = imgaussfilt(J,sigma);
        
%         imagesc(m.STA(cropX,cropY,icell))
        imagesc(J)
        % imagesc(m.STA(:,:,icell))
%         title(['STA clust ' num2str(ES.SpikeInfo{5,icell+1})])
        try
            RedWhiteBlue;
        catch
            colormap gray
        end
%         axis square off
%         box on
%         hcb=colorbar;
%         title(hcb,'gray')
%         plotspecSmoothed(ES,cropX,cropY,scal_f);
        k=k+1;
    end
    
    if isfield(ES.MetaData,'FileName')
        [~,name,~] = fileparts(ES.MetaData.FileName);
    elseif isfield(ES.MetaData,'Filename')
        [~,name,~] = fileparts(ES.MetaData.Filename);
    else
        [~,name,~] = fileparts(ES.MetaData.OutFileName);
    end
    annotation(h,'textbox',[0.45 0.95 0.04 0.04],'String',{[name ' - STA maps']},'FitBoxToText','on', 'FontSize',14,'LineStyle','none', 'Interpreter','none');
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


