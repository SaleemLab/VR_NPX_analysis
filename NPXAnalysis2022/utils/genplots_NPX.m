% Functional evaluation of NPX Recordings
% FR wrote this based on analysis code from the lab on 06/22
% Dependencies: NPXAnalysis2022;
%               Visual-response-analysis/Stimuli/sparsenoise

%function genplots_NPX(StimFreq)

% change animals and session list accordingly
animals = {'M22008','M22009','M22021','M22023','M22020'}; % 
sessions = {'20220408','20220412','20220426','20220427','20220505'}; % 

if ismac
    ROOTPATH = '/Users/s.solomon/Filestore/Research2/ibn-vision';
else
    ROOTPATH = 'X:\ibn-vision';
end
DATAPATH = fullfile(ROOTPATH,'DATA','SUBJECTS');

if ~exist('StimFreq')
    StimFreq = 2;
end

if ~exist('StimType')
    StimType = 'Gratings';
end

for i = 1:length(animals)
    SUBJECT = animals{i};
    SESSION = sessions{i};

    EPHYS_DATAPATH = fullfile(DATAPATH,SUBJECT,'ephys',SESSION);
    BONSAI_DATAPATH = fullfile(DATAPATH,SUBJECT,'bonsai',SESSION);

    fileTable = getBoundGAndBonsaiFileNames(SUBJECT,SESSION,DATAPATH);
   
     SNIdx = strcmp(fileTable.StimulusName,'SparseNoise');
     gs = fileTable.gNumber(SNIdx);
     

    glxDir = fullfile(EPHYS_DATAPATH,[SUBJECT, '_', SESSION, '_g', num2str(gs)]);
    glxDir = fullfile(glxDir,[SUBJECT, '_', SESSION, '_g', num2str(gs) '_imec0']);

    [ChannelMapData] = SGLXMetaToCoords_ChannelMap_FR(glxDir);
    channelRangeTH = sort(find(ChannelMapData.ycoords<1000),'ascend');
    channelRangeCx = sort(find(ChannelMapData.ycoords>1000),'ascend');
    nChannelsToBin = 24;

    % Split the data by different contiguous electrodes to allow SparseNoise script to run %
    ShIdxTh = ChannelMapData.shankind(channelRangeTH);    
    ShIdxCx = ChannelMapData.shankind(channelRangeCx);

    % Get Spiking Data for Drifting and Static Gratings
    [Static_tresps,~,Static_tstimData,Static_ttimeVector,Static_KSoptions] = getKSData(SUBJECT,SESSION,'StaticTF');
    StatictrialstoUse = Static_tstimData.StimIndex==StimFreq; %Currently Using the 2Hz Data

    if ismember('DriftingTF',fileTable.StimulusName)
        [Drifting_tresps,~,Drifting_tstimData,Drifting_ttimeVector,Drifting_KSoptions] = getKSData(SUBJECT,SESSION,'DriftingTF');
        DriftingtrialstoUse = Drifting_tstimData.StimIndex==StimFreq; %Currently Using the 2Hz Data
    else
        Drifting_tresps = [];
        Drifting_tstimData = [];
        Drifting_KSoptions = []; 
        DriftingtrialstoUse = [];
    end

    % Get data for Sparse Noise %
    optionsSN.importMode = 'KS'; % LF or MUA or KS
    optionsSN.BinWidth = 1/60; % resolution (in s) of output resps (e.g. 1/60)
    optionsSN.AnalysisTimeWindow = [0 0.5];% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])
    optionsSN.ks_unitType = ''; % 'mua', 'good' or ''

    % Step 2: Find csv files associated with desired stimulus
    [BehaviourDataFiles,EyeDataFiles,TrialParamsFiles] = getBonsaiFileNames('SparseNoise',BONSAI_DATAPATH);

    % Set bonsai data paths
    PERIPHERALS_DATAPATH = fullfile(BONSAI_DATAPATH, BehaviourDataFiles{1});
    EYEDATA_DATAPATH = fullfile(BONSAI_DATAPATH, EyeDataFiles{1});
    TRIALDATA_DATAPATH = fullfile(BONSAI_DATAPATH,TrialParamsFiles{1});
    optionsSN.PERIPHERALS_DATAPATH = PERIPHERALS_DATAPATH; % full path to BehaviourData bonsai logfile
    optionsSN.EYEDATA_DATAPATH = EYEDATA_DATAPATH; % full path to EyeTrackingData bonsai logfile
    optionsSN.TRIALDATA_DATAPATH = TRIALDATA_DATAPATH; % full path to TrialData bonsai logfile

    % Set ephys data path
    optionsSN.gFileNum = gs;
    folderName = findGFolder(EPHYS_DATAPATH,optionsSN.gFileNum);
    T_EPHYS_DATAPATH = fullfile(EPHYS_DATAPATH,folderName,[folderName,'_imec0']);
    optionsSN.EPHYS_DATAPATH = T_EPHYS_DATAPATH; % full path to ephys recording

    % Extract data
    [SNresps,SNotherData,SNstimData,~,~,~,~,optionsSN] = extractAndCollateNPData(optionsSN);  

    figure
    sgtitle([SUBJECT ' ' SESSION ' Thalamus'])
    ThShank = unique(ShIdxTh);
    CxShank = unique(ShIdxCx);
    ShankNo = length(CxShank);

    for jj = 1:length(ThShank)
        ThisShank = ThShank(jj);
        thisRange = channelRangeTH(ShIdxTh==ThisShank);
%         [initMap,tMap, scal_f, sigma] = getSparseNoise_NPX(SNresps,SNotherData,SNstimData,nChannelsToBin,thisRange,optionsSN);
%         thisMap = squeeze(initMap(:,:,end,tMap));
%         thisMap_s = imresize(thisMap,scal_f);
%         thisMap_s = imgaussfilt(thisMap_s,sigma);
% 
%         % Now plot the RF %
%         subplot(3,ShankNo,jj)
%         imagesc(thisMap_s)
%         axis equal
%         set(gca','XTick',[],'Ytick',[],'Box','off')
%         set(gca,'YDir','normal')
%         colormap(cm)

        UnitsToPlot_Static = ismember(Static_KSoptions.peakChannel,thisRange);
        MeanPSTH_Static = mean(mean(Static_tresps(UnitsToPlot_Static,:,StatictrialstoUse),1),3);        

        % Now Plot the PSTHs %

        subplot(2,ShankNo,jj)
        plot(Static_ttimeVector,sum(mean(Static_tresps(UnitsToPlot_Static,:,StatictrialstoUse),3),1))
        xlim([-0.25 8.75])
        set(gca,'TickDir','out')
        box off
        xline(0,'r-')
        xline(3,'r-')

        if ~isempty(Drifting_tresps)
            UnitsToPlot_Drifting = ismember(Drifting_KSoptions.peakChannel,thisRange);
            MeanPSTH_Drifting = mean(Drifting_tresps(UnitsToPlot_Drifting,:,DriftingtrialstoUse),1);
            
            subplot(2,ShankNo,jj+ShankNo)
            plot(Drifting_ttimeVector,sum(mean(Drifting_tresps(UnitsToPlot_Drifting,:,DriftingtrialstoUse),3),1))
            xlim([-0.25 8.75])
            set(gca,'TickDir','out')
            box off
            xline(0,'r-')
            xline(3,'r-')
        else
            MeanPSTH_Drifting = [];
            Drifting_ttimeVector = [];

            subplot(2,ShankNo,jj+ShankNo)
            plot(Drifting_ttimeVector,MeanPSTH_Drifting)
            xlim([-0.25 8.75])
            set(gca,'TickDir','out')
            box off
            
        end

        

%         clear UnitsToPlot_Static UnitsToPlot_Drifting MeanPSTH_Static MeanPSTH_Drifting
    end
%     clear thisRange ThisShank

    figure
    sgtitle([SUBJECT ' ' SESSION ' Cortex'])

    for jj = 1:length(CxShank)
        ThisShank = CxShank(jj);
        thisRange = channelRangeCx(ShIdxCx==ThisShank);
%         [initMap,tMap, scal_f, sigma] = getSparseNoise_NPX(SUBJECT,SESSION,nChannelsToBin,thisRange,gs);
%         thisMap = squeeze(initMap(:,:,end,tMap));
%         thisMap_s = imresize(thisMap,scal_f);
%         thisMap_s = imgaussfilt(thisMap_s,sigma);
% 
%         % Now plot the RF %
%         subplot(3,ShankNo,jj)
%         imagesc(thisMap_s)
%         axis equal
%         set(gca','XTick',[],'Ytick',[],'Box','off')
%         set(gca,'YDir','normal')
%         colormap(cm)       
        
        UnitsToPlot_Static = ismember(Static_KSoptions.peakChannel,thisRange);
        MeanPSTH_Static = mean(Static_tresps(UnitsToPlot_Static,:,StatictrialstoUse),1);

        subplot(2,ShankNo,jj)
        plot(Static_ttimeVector,sum(mean(Static_tresps(UnitsToPlot_Static,:,StatictrialstoUse),3),1))
        xlim([-0.25 8.75])
        set(gca,'TickDir','out')
        box off
        xline(0,'r-')
        xline(3,'r-')

        if ~isempty(Drifting_tresps)
            UnitsToPlot_Drifting = ismember(Drifting_KSoptions.peakChannel,thisRange);
            MeanPSTH_Drifting = mean(Drifting_tresps(UnitsToPlot_Drifting,:,DriftingtrialstoUse),1);
            
            subplot(2,ShankNo,jj+ShankNo)
            plot(Drifting_ttimeVector,sum(mean(Drifting_tresps(UnitsToPlot_Drifting,:,DriftingtrialstoUse),3),1))
            xlim([-0.25 8.75])
            set(gca,'TickDir','out')
            box off
            xline(0,'r-')
            xline(3,'r-')
        else
            MeanPSTH_Drifting = [];
            Drifting_ttimeVector = [];

            subplot(2,ShankNo,jj+ShankNo)
            plot(Drifting_ttimeVector,MeanPSTH_Drifting)
            xlim([-0.25 8.75])
            set(gca,'TickDir','out')
            box off
            
        end

        

%         clear UnitsToPlot_Static UnitsToPlot_Drifting MeanPSTH_Static MeanPSTH_Drifting
    end
%     clear thisRange SNresps SNstimData SNotherData optionsSN Drifting_tresps Drifting_totherData Drifting_tstimData Drifting_ttimeVector ... 
%         Drifting_KSoptions Static_tresps Static_totherData Static_tstimData Static_ttimeVector Static_KSoptions
end

%end