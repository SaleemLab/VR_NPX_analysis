% function outdat = analyseSparseNoiseNew(ms,pd,options,fileNames,filesPath)
%
% Inputs:
% fileNames: this is a list of filenames in a certain day of interest
% folder. This function will look for a parameters file and for a file with
% the stimulus on each frame info in.
% filesPath: this is the upstream path of those filenames.
% options.figs2plot: this is a 1d vector containing figure numbers, and determines
% which figures are plotted

% Outputs:
% allTrialsMat
% miniscopeMat
% miniscopeMatNorm
% occupancyMap
% stimOnsetsMat
% blackRateMap
% blackRateMapNorm
function outdat = analyseSparseNoiseNew(ms,pd,options,fileNames,filesPath)

outdat = [];

% Some defaults
PSTHWindow = [-1000 5000];
baselinePSTHWindow = [-250 0];
respPSTHWindow = [-0 500];
defaultParams.gridCentreX = -44;
defaultParams.gridWidth = 80;
defaultParams.squareSize = 5; % width of each square in degrees.
defaultParams.gridCentreY = 0;
defaultParams.gridHeight = 60;
contrasts = [0,255];
upsamplingFactor = 10;
gaussSmoothingSigma = 10;


% Get stimulus configs from the .bin file: this uses:
% X:\ibn-vision\Archive - saleemlab\Code\BonVision\Bonsai_analysis\Bonsai_sparsenoise_binread.m
configsDataIndex = contains(fileNames,'stimulusData');
configsFileNames = fileNames(configsDataIndex);

% Also for files from 03/8/19 onwards, get the csv saving a bunch of
% stimulus params.
% Define the centre points of each of the grid of squares in each axis.
% This is not currently used in this function, but is used in the
% plotting function.
% If these centre positions aren't defined, prior to 30/7/19, the
% defaults were; x: to have a gridCentreX of -44 and a gridWidth of 80.
% Therefore the limits of the grid would be -84 and -4. These would be
% the edges of the grid, and with 5 degree diameter squares, this would
% yield grid centres at -84+2.5 = -81.5 and -4-2.5 = -6.5. y:
% gridCentreY = 0, gridHeight = 60. The limits of the grid would
% therefore be -30 to 30. Converted to centres; -27.5 and 27.5.
%     xPositions = -81.5:5:-6.5;
%     yPositions = -27.5:5:27.5;
% This would yield 16 columns and 12 rows. At some point there were 15
% columns.

paramsInfoIndex = contains(fileNames,'sparseNoiseParamsLog');
if any(paramsInfoIndex) % if there is one of these files
    paramsFileNames = fileNames(paramsInfoIndex);
    paramsData = readcell(fullfile(filesPath,paramsFileNames{1}));
    % extract the numbers
    for a = 1 : length(paramsData) % for every column headers
        tempColonInd = strfind(paramsData{a},':'); % find the index of the colon
        tempVal = str2num(paramsData{a}(tempColonInd+1:end)); % take all the characters after that
        switch paramsData{a}(1:tempColonInd-1) % if the characters before it spell...
            case 'gridCentreX' % ... 'gridCentreX' ...
                gridCentreX = tempVal; % ... then log the value appropriately etc
            case 'gridCentreY'
                gridCentreY = tempVal;
            case 'gridWidth'
                gridWidth = tempVal;
            case 'gridHeight'
                gridHeight = tempVal;
            case 'interRepInterval'
                interRepInterval = tempVal;
            case 'squareSize'
                squareSize = tempVal;
        end
    end
else
    gridCentreX = defaultParams.gridCentreX;
    gridWidth = defaultParams.gridWidth;
    squareSize = defaultParams.squareSize; % width of each square in degrees.
    gridCentreY = defaultParams.gridCentreY;
    gridHeight = defaultParams.gridHeight;
end

centreLimsX = gridCentreX + (((gridWidth/2)-squareSize/2)*[-1 1]);
xCentres = (centreLimsX(1):squareSize:centreLimsX(2));
centreLimsY = gridCentreY + (((gridHeight/2)-squareSize/2)*[-1 1]);
yCentres = (centreLimsY(1):squareSize:centreLimsY(2));
numStimRows = size(yCentres,2);
numStimCols = size(xCentres,2);

configsFHandle=fopen(fullfile(filesPath,configsFileNames{1}));
configsData = fread(configsFHandle,inf,'uint8');
configsData = reshape(configsData,...
    [numStimRows,numStimCols,length(configsData)/(numStimCols*numStimRows)]);

if strcmp(options.rigID,'SOLOMON04')
    configsData = flip(configsData,1);
elseif strcmp(options.rigID,'SOLOMON03')
    % then don't do anything
end
configsData(:,:,end) = []; % last frame not shown
fclose(configsFHandle);

% update fig
title([num2str(length(pd.pdUpChangePoints)),' up stimuli found. ',...
    num2str(length(pd.pdDownChangePoints)),' down stimuli found. ',...
    num2str(size(configsData,3)),' total stimuli expected'])

% Align with PD data
% constantISI just takes the first and the last transition to white and
% uses the number of stimuli presented to infer the times of each
% stimulus. It assumes that the times of the onset of each stimulus are
% evenly spaced. For reference, the signal starts off black, then turns
% white with the onset of the first config. It then turns black again
% for the onset of the 2nd config etc. It is important to note that
% when the condition is triggered to stop presenting stimuli, the pd is
% switched to black. This means that if 1 stimulus is presented, the pd
% goes W when it starts and then goes black when the end condition is
% triggered. Similarly, if 2 stimuli are presented, the pd still only
% undergoes 2 transitions. From W to B for the first stimulus, from B
% to W for the second stimulus, and from B to B for the end condition.
% It is important to realise also that the final stimulus may or may
% not have been presented.
timingMethod = 'constantISI';
if mod(size(configsData,3),2) == 1 % if the configs data suggests
    % there was an odd # of stimuli...
    
    % IF YOU HAVE AN ODD NUMBER OF STIMULI, THEN BE WARY OF THIS PIECE OF
    % CODE. ONLY THE EVEN NUMBER OF STIMULI PART IS CONFIRMED.
    
    warning('CODE NOT TESTED, skipping')
    return
    
    % e.g. say we had 4 stimuli (UDUD), that occurred at
    % 20.5, 21, 21.5 and 22s. Then the first and last supThresh
    % values would occur at 20.5 and ~22s (note that the last
    % upChangePoint would be 21.5, and not ~22s), which would
    % indicate the onsets of the first and 4th stimuli
    % respectively. linspace(20.5,22,4) gives us the right onsets
    pdChangesTimesCons = linspace(...
        pd.pdTimebase(find(pd.supThreshPdData,1,'first')),...
        pd.pdTimebase(find(pd.supThreshPdData,1,'last')),...
        size(configsData{1},3));
elseif mod(size(configsData,3),2) == 0 % if the configs data suggests
    % there was an even # of stimuli...
    
    % Take the example where there are 4 stimuli; 0 (W), 0.5 (B), 1
    % (W), 1.5 (B), the upchanges are at 0 and 1. Here, note that
    % linspace(0,1,4) gives incorrect onsets. In contrast, linspace(0,1,3) 
    % gives the right onsets, you just need to create the final one.
    pdChangesTimesCons = linspace(...
        pd.pdTimebase(pd.pdUpChangePoints(1)),...
        pd.pdTimebase(pd.pdUpChangePoints(end)),...
        size(configsData,3)-1);
    pdChangesTimesCons(end+1) = pdChangesTimesCons(end) + mean(diff(pdChangesTimesCons));
end

pdChangesTimesPrec = NaN(size(pdChangesTimesCons));
combPoints = sort(union(pd.pdUpChangePoints,pd.pdDownChangePoints));
pdChangesTimesPrec(1:length(combPoints)) = pd.pdTimebase(combPoints)';
if length(pdChangesTimesCons)~=length(pdChangesTimesPrec)
    warning('pd change detection settings need to be adjusted')
    return
end
nansum(abs(pdChangesTimesCons-pdChangesTimesPrec))

if strcmp(timingMethod,'constantISI')
    pdChangesTimes = pdChangesTimesCons;
elseif strcmp(timingMethod,'precise')
    pdChangesTimes = pdChangesTimesPrec;
end

% Now we have an estimate of the time of each stimulus
PSTHtimebase = PSTHWindow(1) : options.frameDurMs : PSTHWindow(2)-options.frameDurMs;
stimOnsetsMat = NaN(size(configsData,1),size(configsData,2),length(contrasts),1);
respTimebase = respPSTHWindow(1) : options.frameDurMs : respPSTHWindow(2)-options.frameDurMs;


allPSTHmini = NaN(size(configsData,1),size(configsData,2),...
    length(contrasts),1,length(PSTHtimebase));
allPSTHminiNorm = allPSTHmini; % also prep. a normalised version
occupancyMat = NaN(size(configsData,1),size(configsData,2),length(contrasts));
respMat = NaN(size(configsData,1),size(configsData,2),...
    length(contrasts),1,length(respTimebase));
respMatNorm = respMat;


% MOVE THIS TO IN LOOP
allTrialsMat = NaN(length(pdChangesTimes),length(PSTHtimebase));
allTrialsMatNorm = allTrialsMat;
for a = 1:length(pdChangesTimes)
    tempMiniFrames = find((pd.pdTimebase>(pdChangesTimes(a)+PSTHWindow(1))) & ...
        (pd.pdTimebase<=(pdChangesTimes(a)+PSTHWindow(2))));
    tempMiniFrames(tempMiniFrames>size(allTrialsMat,2)) = [];
    allTrialsMat(a,1:length(tempMiniFrames)) = ms.fMiniscopeData(tempMiniFrames,2);
end

% for every row in the grid
for rowNum = 1:size(configsData,1)
    % for every column in the grid
    for colNum = 1:size(configsData,2)
        % i.e. a,b specifies a certain square position.
        % for every contrast at this square position
        for contrastNum = 1:length(contrasts)
            contValTemp = contrasts(contrastNum);
            clear configsIndex tempNumOcc tempMiniFrames
            % find every config where there was the right contrast square here
            % make a list of ones and zeros, where 1s indicate elements with the right contrast
            configsIndex = configsData(rowNum,colNum,:) == contValTemp;
            % find how many of these elements there are
            tempNumOcc = sum(configsIndex);
            % log these in the occupancy map (dims = row, column, contrast)
            occupancyMat(rowNum,colNum,contrastNum) = tempNumOcc;
            % also log all the onset times of these stimuli
            stimOnsetsMat(rowNum,colNum,contrastNum,1:tempNumOcc) = pdChangesTimes(configsIndex);
            % extract the miniscope signal in a perievent window
            for trialNum = 1:tempNumOcc
                tempMiniFrames = find((pd.pdTimebase>(stimOnsetsMat(rowNum,colNum,contrastNum,trialNum)+PSTHWindow(1))) & ...
                    (pd.pdTimebase<=(stimOnsetsMat(rowNum,colNum,contrastNum,trialNum)+PSTHWindow(2))));
                
                if (length(tempMiniFrames) > length(PSTHtimebase))
                    tempMiniFrames = tempMiniFrames(1:length(PSTHtimebase));
                end
                tempMiniFrames(tempMiniFrames>size(ms.fMiniscopeData,1)) = [];
                currentPSTH = ms.fMiniscopeData(tempMiniFrames,2);
                allPSTHmini(rowNum,colNum,contrastNum,trialNum,[1:length(tempMiniFrames)]) = currentPSTH;
                
                % also compile a normalised version
                tempBaseFrames = find((pd.pdTimebase>(stimOnsetsMat(rowNum,colNum,contrastNum,trialNum)+baselinePSTHWindow(1))) & ...
                    (pd.pdTimebase<=(stimOnsetsMat(rowNum,colNum,contrastNum,trialNum)+baselinePSTHWindow(2))));
                if (length(tempBaseFrames) > length(PSTHtimebase))
                    tempBaseFrames = tempBaseFrames(1:length(PSTHtimebase));
                end
                tempBaseline = ms.fMiniscopeData(tempBaseFrames,2);
                
                allPSTHminiNorm(rowNum,colNum,contrastNum,trialNum,[1:length(tempMiniFrames)]) = currentPSTH - nanmean(tempBaseline);
                
                tempRespFrames = find((pd.pdTimebase>(stimOnsetsMat(rowNum,colNum,contrastNum,trialNum)+respPSTHWindow(1))) & ...
                    (pd.pdTimebase<=(stimOnsetsMat(rowNum,colNum,contrastNum,trialNum)+respPSTHWindow(2))));
                if (length(tempRespFrames) > length(PSTHtimebase))
                    tempRespFrames = tempRespFrames(1:length(PSTHtimebase));
                end
                tempRespFrames(tempRespFrames>size(ms.fMiniscopeData,1)) = [];
                tempResp = ms.fMiniscopeData(tempRespFrames,2);

                respMat(rowNum,colNum,contrastNum,trialNum,[1:length(tempRespFrames)]) = tempResp';
                respMatNorm(rowNum,colNum,contrastNum,trialNum,[1:length(tempRespFrames)]) = tempResp'-nanmean(tempBaseline);
            end
        end
    end
end

% respMat = respMat(:,:,:,1:8,:);
% respMatNorm = respMatNorm(:,:,:,1:8,:);
blackRateMap = nanmean(nanmean(respMat(:,:,1,:,:),5),4);
blackRateMapNorm = nanmean(nanmean(respMatNorm(:,:,1,:,:),5),4);

upSampledBlackRateMapNorm = imresize(blackRateMapNorm,[size(blackRateMapNorm,1)*upsamplingFactor,NaN]);
upsampledXPoints = linspace(1,length(xCentres),length(xCentres)*upsamplingFactor);
upsampledYPoints = linspace(1,length(yCentres),length(yCentres)*upsamplingFactor);
upsampledXCentres = interp1(1:length(xCentres),xCentres,upsampledXPoints);
upsampledYCentres = interp1(1:length(yCentres),yCentres,upsampledYPoints);

smoothedBlackRateMapNorm = imgaussfilt(upSampledBlackRateMapNorm,gaussSmoothingSigma);
if length(contrasts) == 2 % if you even presented white squares
    whiteRateMap = nanmean(nanmean(respMatNorm(:,:,2,:,:),5),4);
    whiteRateMapNorm = nanmean(nanmean(respMatNorm(:,:,2,:,:),5),4);
end

ratemapsToAnalyse = {blackRateMapNorm,whiteRateMapNorm,smoothedBlackRateMapNorm};
xCentresToAnalyse = {xCentres,xCentres,upsampledXCentres};
yCentresToAnalyse = {yCentres,yCentres,upsampledYCentres};
for currRatemap = 1:length(ratemapsToAnalyse)
    [~,linInd] = max(ratemapsToAnalyse{currRatemap}(:)); % output the linear index of the max val
    [maxRow,maxColumn] = ind2sub(size(ratemapsToAnalyse{currRatemap}),linInd); % convert to subscripts
    ratemapXpeaks(currRatemap) = xCentresToAnalyse{currRatemap}(maxColumn);
    ratemapYpeaks(currRatemap) = yCentresToAnalyse{currRatemap}(maxRow);
end

outdat.blackRateMapNorm = blackRateMapNorm;
outdat.xCentres = xCentres;
outdat.yCentres = yCentres;

% PLOTTING

% FIGURE 1 ----------------------------------------------------------------
% Plot responses to all configs. Each row is a different config, the x axis
% is time and the color gives the response.
if any(options.figs2plot==1)
    figure('name','Responses_to_all_squares','position',options.figurePosVec);
    imagesc(PSTHtimebase/1000,size(allTrialsMat,2),allTrialsMat);
    xlabel('time relative to stimulus onset (secs)');
    ylabel('trial number');
    box off
    set(gca,'TickDir','out')
    title('Responses to all squares');
end

% FIGURE 2 ----------------------------------------------------------------
% Plot the black ratemap and look at a gaussian fit
% figure('name',['Sparse noise responses for ',doiFiles]);
if any(options.figs2plot==2)
%     figure()
%     title('blackRateMap')
%     imagesc(blackRateMapNorm);
    
    % x = 1 : size(blackRateMap, 2); % Columns.
    % y = 1 : size(blackRateMap, 1); % Rows.
    % [X, Y] = meshgrid(x, y);
    % meanA = mean(blackRateMap(:));
    % centerOfMassX = mean(blackRateMap(:) .* X(:)) / meanA;
    % centerOfMassY = mean(blackRateMap(:) .* Y(:)) / meanA;
    %
    % [fitparams, modelout] = fit2DGauss(blackRateMap,Y,X);
    % theta = 0 : 0.01 : 2*pi;
    % x = fitparams(4)/2 * cos(theta) + fitparams(1);
    % y = fitparams(4)/2 * sin(theta) + fitparams(2);
    % plot(x, y,'w','LineWidth',5);
    % axis equal;
    % axis tight;
    % subplot(4,4,4)
    % imagesc(modelout)
    % set(gca,'YDir','normal')
    % blackRateMapNorm = blackRateMap;
end

% FIGURE 3 ----------------------------------------------------------------
% Plot the black and white ratemaps
if any(options.figs2plot==3)
    sigma = [];
    ratemapLabels = {'black_square_responses','white_square_responses','smoothed_black_square responses'};

    for currRatemap = 1:length(ratemapsToAnalyse)
        figure('name',ratemapLabels{currRatemap},'position',options.figurePosVec);
        if options.fitGaussesFlag && currRatemap<2
            subplot(121) 
        end
%             colormap gray

        imagesc(xCentres,yCentres,ratemapsToAnalyse{currRatemap})
        hold on
        
        set(gca,'YDir','normal')
        cbar2 = colorbar;
        cbar2.Label.String = 'Average miniscope signal';
        set(gca,'TickDir','out')
        xlabel('Position in azimuth, or x, (deg)')
        if ~options.fitGaussesFlag
            set(gca,'XTick',xCentres)
            set(gca,'YTick',yCentres)
        end

        ylabel('Position in elevation, or y, (deg)')
        axis equal
        axis tight
        box off
        
        scatter(ratemapXpeaks(currRatemap),ratemapYpeaks(currRatemap),'kx');

        
        if options.fitGaussesFlag && currRatemap<2
            [X, Y] = meshgrid(xCentres, yCentres);
            fitoptions.sdInit = 2;
            fitoptions.rangeBg = 2; % eg = 2; ie from half to twice the minimum
            fitoptions.interpModelOutMultiplier = 10;
            [fitparams, modelout, fitoptions] = ibn_fit2DGaussLsq(ratemapsToAnalyse{currRatemap},X,Y,fitoptions);
            centreX = fitparams(1); 
            centreY = fitparams(2);
            sigmaDeg = fitparams(4);
            viscircles([centreX,centreY],sigmaDeg);
            subplot(122);
            imagesc(fitoptions.interpx(1,:),fitoptions.interpy(:,2),fitoptions.interpmodelout);
            viscircles([centreX,centreY],sigmaDeg);
            set(gca,'YDir','normal')
            axis equal; axis tight;
        end
        
        sgtitle({filesPath,['Peak @ (x,y): ',num2str(ratemapXpeaks(currRatemap)),...
            ', ',num2str(ratemapYpeaks(currRatemap)),' sigma = ',num2str(sigmaDeg)]});
    end
    
end

% FIGURE 7 ----------------------------------------------------------------
% Plot the the single trials for the ratemap peaks
if any(options.figs2plot==7)
    figure('name','Sparse_noise_traces','position',options.figurePosVec);
    subplot(2,2,1)

    allPSTHminiNorm(allPSTHminiNorm==0) = NaN; % temporarily convert zeros to nans
    for a = 1:size(allPSTHminiNorm,4)
        tempTrace = squeeze(allPSTHminiNorm(maxRow,maxColumn,1,a,:));
        plot(PSTHtimebase,tempTrace)
        hold on
    end
    yLimsPlot1 = get(gca,'YLim');
    plot([0 0],yLimsPlot1,'k--');
    plot([500 500],yLimsPlot1,'k--');
    
    meanBlackTracePeak = nanmean(squeeze(allPSTHminiNorm(maxRow,maxColumn,1,:,:)));
    scatter(PSTHtimebase,meanBlackTracePeak)
    title(['Responses_at_ratemap_peak: x = ',num2str(blackXPeak),', y = ',num2str(blackYPeak)])
    xlabel('time relative to stimulus onset (ms)')
    ylabel('miniscope signal')
    
    if length(contrasts) == 2 % if you even presented white squares
        subplot(2,2,3)
        [~,linIndWhite] = max(whiteRateMapNorm(:)); % output the linear index of the max val
        [maxWhiteRow,maxWhiteColumn] = ind2sub(size(whiteRateMapNorm),linIndWhite); % convert to subscripts
        for a = 1:size(allPSTHminiNorm,4)
            tempTrace = squeeze(allPSTHminiNorm(maxWhiteRow,maxWhiteColumn,1,a,:));
            plot(PSTHtimebase,tempTrace)
            hold on
        end
        yLimsPlot1 = get(gca,'YLim');
        meanWhiteTracePeak = nanmean(squeeze(allPSTHminiNorm(maxWhiteRow,maxWhiteColumn,1,:,:)));
        scatter(PSTHtimebase,meanWhiteTracePeak)
        title(['MiniscopeData for ratemap peak: ',num2str(maxWhiteRow),',',num2str(maxWhiteColumn)])
        xlabel('time relative to stimulus onset (ms)')
        ylabel('miniscope signal')
        plot([0 0],yLimsPlot1,'k--');
        plot([500 500],yLimsPlot1,'k--');
    end
end

% FIGURE 4 ----------------------------------------------------------------
% Plot a figure where each row is a different trial from stimulus
% presentations in and around the receptive field, and where different
% columns are different timepoints, and different color values indicate
% different responses.
if any(options.figs2plot==4)
    figure('name','Sparse_noise_responses','position',options.figurePosVec);
    span = 1;
    rowRange = (maxRow-span:maxRow+span);
    rowRange(rowRange>size(allPSTHminiNorm,1)) = size(allPSTHminiNorm,1);
    rowRange(rowRange<=1) = 1;
    colRange = (maxColumn-span:maxColumn+span);
    colRange(colRange>size(allPSTHminiNorm,2)) = size(allPSTHminiNorm,2);
    colRange(colRange<=1) = 1;
    matOI = squeeze(allPSTHminiNorm(rowRange,colRange,1,:,:));
    imagesc(squeeze(nanmean(nanmean(matOI,1),2)));
    title('Responses to ~intraRF squares')
end

% FIGURE 5 ----------------------------------------------------------------
if any(options.figs2plot==5)
    figure('name','Occupancy_maps','position',options.figurePosVec);
    subplot(1,2,1)
    tempBlackOccMap = occupancyMat(:,:,1,1);
    imagesc(xCentres,yCentres,tempBlackOccMap);
    set(gca,'YDir','normal')
    colorbar
    title(['Occupancy map for black squares: range of reps is ',num2str(range(tempBlackOccMap(:)))])
    if length(contrasts) == 2 % if you even presented white squares
        subplot(1,2,2)
        tempWhiteOccMap = occupancyMat(:,:,2,1);
        imagesc(xCentres,yCentres,tempWhiteOccMap);
        set(gca,'YDir','normal')
        colorbar
        title(['Occupancy map for white squares: range of reps is ',num2str(range(tempWhiteOccMap(:)))])
    end
end

% FIGURE 6 ----------------------------------------------------------------
% Plot the first 4 frames in the same way as the rf maps are plotted to
% show that they're being plotted in the right way
if any(options.figs2plot==6)
    figure('name','First_few_frames','position',options.figurePosVec);
    colormap gray
    for frameNum = 1:4
        subplot(2,2,frameNum)
        imagesc(xCentres,yCentres,configsData(:,:,frameNum));
        set(gca,'YDir','normal')
        caxis([0 255])
        xlabel('position in azimuth (deg)')
        xlabel('position in elevation (deg)')
        title(['frame number ',num2str(frameNum)]);
    end
    sgtitle('First four frames')
end

figure(3)
end
