% 

folderRecs =  'X:\CODE\DEV\general\SparseNoise\Tests'

ProcessedFiles = dir(folderRecs);

f_n_i=length(ProcessedFiles);


for rec  = 3:length(ProcessedFiles)
    
    load([ProcessedFiles(rec).folder filesep ProcessedFiles(rec).name]);
    
    ES_multi{rec-2} = ES;
    
    TotalDuration(rec-2) = ES.SN_onsets(end)-ES.SN_onsets(1);
    
    clear PD_changes m_d
    PD_changes = ES.SN_onsets;
    % most common values of SN frame duration
    frameDur = diff(PD_changes)*1000;
    frameDur_t = frameDur;
    m_d(1) = mode(frameDur_t);
    frameDur_t(frameDur_t>100) = NaN;
    m_d(2) = mode(frameDur_t);
    % change manually the value of each onset time
    for i = 1:length(frameDur)
        if frameDur(i)>90 && frameDur(i)<100
            PD_changes(i+1) = PD_changes(i+1)+abs(diff(m_d)/2)/1000;
        end
    end   
    
    edges = (0:0.2:200);
    
    Distr2plot(rec-2,:) = histcounts(diff(PD_changes)*1000,edges);
    SN_dur(rec-2,1:2998) = diff(PD_changes(1:2999))*1000;
 
end

% find
frameDur = 16.667;
edges<(100-16.667)
% shorter
sum(sum(Distr2plot(:,edges<(100-16.667)))) / size(Distr2plot,1)
sum(SN_dur'<(100-16.667))
% longer
sum(sum(Distr2plot(:,edges(1:end-1)>(100+16.667)))) / size(Distr2plot,1)
sum(SN_dur'>(100+16.667))
% right
sum(sum(Distr2plot(:,(edges(1:end-1)>(100-16.667) & edges(1:end-1)<(100+16.667)) ))) / size(Distr2plot,1)
std(sum(SN_dur'>(100-16.667) & SN_dur'<(100+16.667)))


0.1003*3000


1/(100.3/6)*1000

median(diff(BonVisionData.Time))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




figure
RedWhiteBlue
for icell = 1:33
subplot(1,3,1)
imagesc(ES.m.ONMap(:,:,icell))
colorbar
title([num2str(icell) '  ' num2str(ES.m.ONdelays.max(icell))])
caxis([min(min(ES.m.ONMap(:,:,icell))) max(max(ES.m.ONMap(:,:,icell)))])
subplot(1,3,2)
imagesc(ES_1.m_1.ONMap(:,:,icell))
colorbar
caxis([min(min(ES.m.ONMap(:,:,icell))) max(max(ES.m.ONMap(:,:,icell)))])
title(num2str(ES_1.m_1.ONdelays.max(icell)));
subplot(1,3,3)
imagesc(ES_1.m_1.ONMap(:,:,icell)-ES.m.ONMap(:,:,icell))
colorbar
caxis([min(min(ES.m.ONMap(:,:,icell))) max(max(ES.m.ONMap(:,:,icell)))])
% title(num2str(ES_1.m_1.ONdelays.max(icell)));

pause
end

ES_1.m_1.ONdelays.max - ES.m.ONdelays.max    


figure
imagesc(sum(m.OFFMap(cropX,cropY,ne_i & th_i),3) - ...
    sum(m_1.OFFMap(cropX,cropY,ne_i & th_i),3))




folderRecs = 'X:\CODE\DEV\general\SparseNoise\Tests';

