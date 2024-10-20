%% Forward correlation (ON map)
function m = ForwardCorrelationONMap(ES,m,plotSingleUnits)
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


% frames with a white square presented in this xy position
for xpos=1:ES.SN_gridDims(1)
    for ypos=1:ES.SN_gridDims(2)
        WStimFrames{xpos, ypos}=find(ES.SN_sequence(xpos,ypos,1:numframes)==255);
    end
end

% Count spikes per pixel for each cluster
for icell=1:length(ES.frameSpikeCount)
    for delayInd=1:length(delays)
        for xpos=1:ES.SN_gridDims(1)
            for ypos=1:ES.SN_gridDims(2)
                SpikeCount=sum(ES.frameSpikeCount{icell}(WStimFrames{xpos,ypos},delayInd));
                ONMap(xpos,ypos,delayInd,icell)=SpikeCount/(length(WStimFrames{xpos,ypos})*ES.SN_singleI_dur); % normalize with repeats per pixel
            end
        end
    end
    % pick the best delay for this cluster (delay that maximizes spatial variance)
    [maxVar, maxInd, vars] = maxVarMap(ONMap(:,:,:,icell));
    
    m.ONdelays.delays= delays;
    m.ONdelays.max(icell)=delays(maxInd);
    m.ONvars.vars(icell,:)= vars;
    m.ONvars.max(icell)= maxVar;
    m.ONMap(:,:,icell) = ONMap(:,:,maxInd,icell);
end

if plotSingleUnits
    k=1;
    figure;
    % define nr of subplots
    SProws = floor(length(ES.frameSpikeCount)/8);
    if mod(length(ES.frameSpikeCount),8)
        SProws = SProws+1;
    end
    for icell=1:length(ES.frameSpikeCount)%[6,15,20,22,23,24]
        subplot(SProws,8,icell)
        %         subplot(1,6,k)
        scal_f = 10; sigma = 3;
        I = m.ONMap(cropX,cropY,icell);
        J = imresize(I,scal_f);
        J = imgaussfilt(J,sigma);
%         J = (J/max(J(:)));
        imagesc(J)
%         imagesc(m.ONMap(cropX,cropY,icell))
        %         imagesc(m.ONMap(:,:,icell))
        title(['ON clust ' num2str(ES.SpikeInfo{5,icell+1})])% ' ' num2str(m.ONdelays.max(icell))
        colormap ((gray))
        axis square
        box on
        hcb=colorbar;
        title(hcb,'hz')
        plotspecSmoothed(ES,cropX,cropY,scal_f);
        k=k+1;
    end  
end

end
% 
% function plotspec(ES,cropX,cropY)
% line([0+0.5 0+0.5],[0  ES.SN_gridDims(2)+0.5],'Color','blue','LineStyle','-','LineWidth',1)
% line([0+0.5 ES.SN_gridDims(1)+0.5],[(((ES.stim_centre_pos(2)+ES.stim_dims(2)/2))/ES.stim_dims(2))*ES.SN_gridDims(2) (((ES.stim_centre_pos(2)+ES.stim_dims(2)/2))/ES.stim_dims(2))*ES.SN_gridDims(2)]+0.5,'Color','blue','LineStyle','-','LineWidth',1)
% set(gca,...
%     'ylim',[cropY(1)+0.5 cropY(end)+0.5],...
%     'YTick',[1 cropY(end)],...
%     'YTickLabel',[ES.stim_centre_pos(2)+ES.stim_dims(2)/2    ES.stim_centre_pos(2)-ES.stim_dims(2)/2],...
%     'xlim',[cropX(1)+0.5 cropX(end)+0.5],...
%     'XTick',[1 cropX(1)],...
%     'XTickLabel',[ES.stim_centre_pos(1)-ES.stim_dims(1)/2    ES.stim_centre_pos(1)+ES.stim_dims(1)/2])
% end
