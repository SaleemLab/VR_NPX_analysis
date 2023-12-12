%% Reverse correlation, STC (no ON OFF stim discrimination, easy spike triggered covariance)
function m = ReverseCorrelationSTC(ES,m,plotSingleUnits)


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
    
%     totalSpikeCount = NaN(size(ES.frameSpikeCount{icell}));
    for delayInd=1:length(delays)
        
        for frame=1:numframes
            
            STC(:,:,delayInd,icell)= STC(:,:,delayInd,icell) + stim_frames(:,:,frame)*ES.frameSpikeCount{icell}(frame,delayInd);
%             totalSpikeCount(frame,delayInd) = totalSpikeCount(frame,delayInd) +  ES.frameSpikeCount{icell}(frame,delayInd);
            
        end
    end
    
    
    % pick the best delay for this cluster (delay that maximizes spatial variance)
    [maxVar, maxInd, vars] = maxVarMap(STC(:,:,:,icell));
    
    m.spikecount=ES.frameSpikeCount{icell}(1:numframes,maxInd);
    m.STCdelays.delays= delays;
    m.STCdelays.max(icell)=delays(maxInd);
    m.STCvars.vars(icell,:)= vars;
    m.STCvars.max(icell)= maxVar;
    m.STC(:,:,icell) = STC(:,:,maxInd,icell)./sum(m.spikecount);
    
%     [forfit_xval,forfit_yval]=meshgrid(1:8,1:8);
    [forfit_xval,forfit_yval]=meshgrid(7.5:15:112.5,-22.5:15:82.5);
    options.rangeBg=2;
    options.sdInit=3;
    options.interpModelOutMultiplier=10;
    if all(isnan(m.STC(:,:,icell)))
        fitparams(icell,:)=NaN(1,5);
        modelout(:,:,icell)=NaN(8,8);
        interpmodelout(:,:,icell)=NaN(80,80);
        allcorr(icell)=NaN;
        allpValCorr(icell)=NaN;
        R2(icell)=NaN;
    else
        [fitparams(icell,:), modelout(:,:,icell), options,R2(icell),initparamsBl] = ibn_fit2DGaussLsq(m.STC(:,:,icell),forfit_xval,forfit_yval,options);
        interpmodelout(:,:,icell)=options.interpmodelout;
        %%Prediction of the response
        thismodel=squeeze(modelout(:,:,icell))-fitparams(icell,end);
        thispred=NaN(1,numframes);
        for thisframe=1:numframes
            tframe=squeeze(stim_frames(:,:,thisframe));
            thispred(thisframe)=sum(thismodel(:).*tframe(:));
        end
        thispred=thispred+initparamsBl(5);
        [r,p]=corrcoef(thispred(:),m.spikecount(:));
        allcorr(icell)=r(1,2);
        allpValCorr(icell)=p(1,2);
    end
    
end

m.fitparams=fitparams;
m.modelout=modelout;
m.options=options;
m.interpmodelout=interpmodelout;
m.R2=R2;
m.allcorr=allcorr;
m.allpValCorr=allpValCorr;

if plotSingleUnits
    h = figure
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
%         title(['STC clust ' num2str(icell)])
%         colormap gray
        %axis off
%         axis square
%         box on
%         hcb=colorbar;
%         title(hcb,'gray')
%         plotspec(ES,cropX,cropY);
        k=k+1;
    end
     if isfield(ES.MetaData,'FileName')
        [~,name,~] = fileparts(ES.MetaData.FileName);
     elseif isfield(ES.MetaData,'Filename')
        [~,name,~] = fileparts(ES.MetaData.Filename);
     else
        [~,name,~] = fileparts(ES.MetaData.OutFileName);
     end
    annotation(h,'textbox',[0.45 0.95 0.04 0.04],'String',{[name ' - STC maps']},'FitBoxToText','on', 'FontSize',14,'LineStyle','none', 'Interpreter','none');
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
