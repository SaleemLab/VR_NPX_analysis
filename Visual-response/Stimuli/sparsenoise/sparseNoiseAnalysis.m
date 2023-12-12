% Program to take in Stimulus, Physiological variable and & Wheel / Eye data and organise it appropriately to make RF maps
% Inputs -
%            - options, including
%           - responses (time x frames psth; ie each row a different time
%           point in the perievent histogram, each column a different
%           sample of noise)
%           - wheel_data (optional; 1 x frames)
%           - eye_data (optional; 1 x frames)
% Outputs
function [initMap,options] = sparseNoiseAnalysis(stim_matrix,responses,wheel_data,eye_data,options)
% Some defaults
if ~exist('options','var')
    % Following defaults really should be defined elsewhere
    options.grid_size = [10 10]; % 10x10
end
%%%
% Other defulats that its OK to leave defined here if you like
if ~isfield(options,'mapsToShow') || isempty(options.mapsToShow)
    options.mapsToShow = {'linear','black','white','contrast'};
end
if ~isfield(options,'mapMethod') || isempty(options.mapMethod)
    options.mapMethod = 'fitlm'; % fitlm mean
end
if ~isfield(options,'mapSampleRate') || isempty(options.mapSampleRate)
    options.mapSampleRate = 60; % Hz
end
if ~isfield(options,'framesToShow') || isempty(options.framesToShow)
    options.framesToShow = [1 4 6 8 10 12];  % at 60 Hz would be 8 ms, 40, 72 etc
end
if ~isfield(options,'plotflag') || isempty(options.plotflag)
    options.plotflag = 1;
end
if ~isfield(options,'figName') || isempty(options.figName)
    options.figName = 'Dummy';
end

%%%%%
% Cast some variables
figName = options.figName;
plotflag = options.plotflag;
grid_size = options.grid_size;
framesToShow = options.framesToShow;
mapsToShow = options.mapsToShow;
mapMethod = options.mapMethod;
mapSampleRate = options.mapSampleRate;

% Get SVD over responses
[stimulusWeights,overallTimeCourse] = getResponseSVD(responses);

%%%%%
% Now make relevant maps
initMap = NaN(size(stim_matrix,1),size(stim_matrix,2),length(framesToShow)+1,length(mapsToShow));
initPMap = NaN(size(stim_matrix,1),size(stim_matrix,2),length(framesToShow)+1,length(mapsToShow));

for tMap = 1:length(mapsToShow)
    thisMap = mapsToShow{tMap};
    thisMap_stim_matrix = stim_matrix;
    switch(thisMap)
        case 'linear'
            % do nothing
        case 'black'
            thisMap_stim_matrix(thisMap_stim_matrix == 1) = 0;
        case 'white'
            thisMap_stim_matrix(thisMap_stim_matrix == -1) = 0;
        case 'contrast'
            thisMap_stim_matrix(thisMap_stim_matrix == -1) = 1;
    end
    
    %%%%%%
    % NB If we had been logging data about the stimulus on the same time course as the response that we have nw
    % accuumulated (every frame) we could simply take the response vector and the stimulus matrix, above, and compare the two.
    % But we have been changing the stimulus rarely (every 0.5 s usually) and
    % sparsely logging this stimulus so we now need to do the stimulus at each timepoint
    for tFrame = 1:length(framesToShow)
        thisFrame = framesToShow(tFrame);
        [ initMap(:,:,tFrame,tMap),initPMap(:,:,tFrame,tMap)] = makeRFmap(thisMap_stim_matrix,responses(:,thisFrame),0,mapMethod);
    end
    [initMap(:,:,tFrame+1,tMap), initPMap(:,:,tFrame+1,tMap)] = makeRFmap(thisMap_stim_matrix,stimulusWeights,0,mapMethod);
end

%%%%%%%%%
if plotflag
    % Display the maps
    cm = RedWhiteBlue;
    scal_f = 10; % scale image by this before...
    sigma = 3; % ...filtering by this
    frameCentreTimeMs = 1000*framesToShow/mapSampleRate - 0.5/mapSampleRate;
    responsesframeCentreTimeMs = 1000*(1:size(responses,2))/mapSampleRate - 0.5/mapSampleRate;
    
    % First the VEP
    figure('Name',[figName,' : RF responses'])
    
    subplot(221)
    % Plot average responses for each grid location
    on_responses = zeros(grid_size(1),grid_size(2),size(responses,2));
    off_responses = zeros(grid_size(1),grid_size(2),size(responses,2));
    %figure;
    for i=1:grid_size(1)
        for j=1:grid_size(2)
            on_responses(i,j,:) = nanmean(responses(stim_matrix(i,j,:)==1,:),1);
            off_responses(i,j,:) = nanmean(responses(stim_matrix(i,j,:)==-1,:),1);
        end
    end
    on_responses_reshaped = reshape(on_responses,[grid_size(1)*grid_size(2) size(responses,2)]);
    off_responses_reshaped = reshape(off_responses,[grid_size(1)*grid_size(2) size(responses,2)]);
    min_responses = min(min(on_responses_reshaped(:)),min(off_responses_reshaped(:)));
    max_responses = max(max(on_responses_reshaped(:)),max(off_responses_reshaped(:)));
    on_responses_norm = (on_responses_reshaped-min_responses)./(max_responses-min_responses)+0.5;
    on_responses_norm = reshape(on_responses_norm,[grid_size(1),grid_size(2), size(responses,2)]);
    off_responses_norm = (off_responses_reshaped-min_responses)./(max_responses-min_responses)+0.5;
    off_responses_norm = reshape(off_responses_norm,[grid_size(1),grid_size(2), size(responses,2)]);
    
    step = 0.5:1/size(responses,2):1.5;
    step = step(1:end-1);
    for i=1:grid_size(1)
        for j=1:grid_size(2)
            plot(step+j-1, squeeze(on_responses_norm(i,j,:)+i-1),'color','r');
            hold on;
            plot(step+j-1, squeeze(off_responses_norm(i,j,:)+i-1),'color','b');
        end
    end
    title([figName,' ON (r)+OFF(b) Responses'])
    set(gca,'YDir','normal')
    set(gca,'XLim',[0 grid_size(1)+1],'YLim',[0 grid_size(2)+1])
    axis equal
    
    subplot(4,2,5)
    tm = mean(responses,1);
    plot(responsesframeCentreTimeMs,tm,'b-'); hold on
    plot(frameCentreTimeMs,tm(framesToShow),'ro','MarkerFaceColor','w')
    set(gca,'XLim',[0 max(responsesframeCentreTimeMs)],'TickDir','out','Box','off')
    ylabel('mV')
    title('Mean')
    
    subplot(4,2,7)
    plot(responsesframeCentreTimeMs,overallTimeCourse,'b-')
    set(gca,'XLim',[0 max(responsesframeCentreTimeMs)],'TickDir','out','Box','off')
    xlabel('Time (ms)')
    ylabel('Units')
    title('SVD')
    
    subplot(4,2,2)
    plot(1:length(stimulusWeights),max(abs(responses),[],2),'r-')
    set(gca,'XLim',[0 length(stimulusWeights)+1],'TickDir','out','Box','off')
    title('Max |Mean|')
    
    subplot(4,2,4)
    plot(1:length(stimulusWeights),stimulusWeights,'r-'); hold on
    set(gca,'XLim',[0 length(stimulusWeights)+1],'TickDir','out','Box','off')
    title('SVD')
    
    if ~isempty(eye_data)
        subplot(4,2,6)
        plot(1:length(eye_data),eye_data,'b-')
        set(gca,'XLim',[0 length(eye_data)+1],'TickDir','out','Box','off')
        title('Eye area')
    end
    
    if ~isempty(wheel_data)
        subplot(4,2,8)
        plot(1:length(wheel_data),wheel_data,'b-')
        set(gca,'XLim',[0 length(wheel_data)+1],'TickDir','out','Box','off')
        xlabel('Trial')
        title('Speed')
    end
    % drawnow;
    
    % Now the maps
    figure('Name',[figName,' : RF spatial maps'])
    %drawnow;
    
    for tMap = 1:length(mapsToShow)
        tm = initMap(:,:,1:end-1,tMap);
        maxR = max(abs(tm(:)));
        for tFrame = 1:length(framesToShow)
            thisMap = squeeze(initMap(:,:,tFrame,tMap));
            thisMap_s = imresize(thisMap,scal_f);
            thisMap_s = imgaussfilt(thisMap_s,sigma);
            
            subplot(6,length(framesToShow)+2,tFrame+(tMap-1)*(length(framesToShow)+2))
            try
            imagesc(thisMap_s,[-maxR maxR])
            axis equal
            set(gca','XTick',[],'Ytick',[],'Box','off')
            set(gca,'YDir','normal')
            colormap(cm)

            % Annotate
            if tFrame == 3 && tMap == 1
                title(figName)
            end
            if tFrame == 1
                title(sprintf('%s',mapsToShow{tMap}))
            end
            text(0.05,0.95,sprintf('%3.1fms', frameCentreTimeMs(tFrame)),'Units','normalized')
            catch
            end
        end
        
        % And SVD
        thisMap = squeeze(initMap(:,:,end,tMap));
        thisMap_s = imresize(thisMap,scal_f);
        thisMap_s = imgaussfilt(thisMap_s,sigma);
        
        subplot(6,length(framesToShow)+2,tMap*(length(framesToShow)+2)-1)
        imagesc(thisMap_s)
        axis equal
        set(gca','XTick',[],'Ytick',[],'Box','off')
        set(gca,'YDir','normal')
        colormap(cm)
        
        if tMap == 1
            title('SVD')
        end
        
        % And min P value
        thisMap =  initPMap(:,:,tFrame+1,tMap);
        %     thisMap_s = imresize(thisMap,scal_f);
        %     thisMap_s = imgaussfilt(thisMap_s,sigma);
        
        subplot(6,length(framesToShow)+2,tMap*(length(framesToShow)+2))
        imagesc(log10(thisMap),[-3 3])
        axis equal
        set(gca','XTick',[],'Ytick',[],'Box','off')
        set(gca,'YDir','normal')
        colormap(cm)
        
        if tMap == 1
            title('SVD p')
        end
        
        %drawnow;
    end
end
