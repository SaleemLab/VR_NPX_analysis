function [es, VRDetails, chans] = ...
    getVRspikes_SL(animal,iseries,expt_list,igroup,iaddinfo,sync_option ...
    ,discreet_steps, flag_load, flag_spont, flag_spkrate,loadTimes,SmthTimeWindow,samplerate,FloadVRfile)

% modified to take multiple iexp -- 2019 MM
%% initialization
global DIRS
DIRS=SetDefaults('animal',animal);
DIRS.ball= 'X:\Archive - saleemlab\Data\Behav'; % need to fix this later 2019-10 MM

if nargin~=5
    flag_load = 1;
end

if nargin<5
    iaddinfo=[];
end

if nargin<6
    discreet_steps = 100;
    sync_option= 'syncpulse'; % current alternative is 'speedsignal' 2019-05 MM
end

if nargin<6
    flag_spont = 1;
end

if nargin<7
    flag_spkrate = 0;
end

if nargin<10
    loadTimes = 0;
end

if nargin<11
    SmthTimeWindow = 150;
end

if nargin<12
    samplerate = 60;
end

if nargin<13
    FloadVRfile = true;
end

for expt=1:length(expt_list)
    
iexp=expt_list(expt)

% fname = [animal '_' num2str(iseries) '_' num2str(iexp)]; probably old format
fname = [animal '_' num2str(iexp)];
dDIRname = [DIRS.Processed filesep num2str(iseries)]; % temp.fix 2019-01 MM
run_rest = 0;

if exist([dDIRname filesep fname '_es' '.mat'],'file') == 2
    try
        load([dDIRname filesep fname '_es','.mat'],'es');
        load([dDIRname filesep fname '_spikeInfo.mat'],'TempChans'); chans = TempChans; clear TempChans;
        load([dDIRname filesep fname '_VRData.mat'],'VRDetails');
        run_rest = 0;
    catch
        run_rest = 1;
    end
else
    disp('No existing file');
    run_rest = 1;
end

if run_rest == 1
    %% get analog channels 
    
    switch sync_option
        case  'syncpulse'
            AC_Info = getSyncPulseSignal_VR(animal, iseries, iexp);
        case 'speedsignal'
            DataFolder=[DIRS.multichanspikes filesep animal];
            ACInfo = OE_get_AC_Signals(DataFolder);
            ACInfo = OE_Compute_Speed_Signal(ACInfo);
    end
    
    %% get all behav data from the day 
    
     [VRdata, VRdata_o, es] = VRWheelLoad_SLv2(animal, iseries, iexp);
    iexp_list = dir([DIRS.ball filesep num2str(iseries) filesep animal '_' num2str(iseries) '*.mat']);
    if length(iexp_list)>1
        for gg = 2:length(iexp_list)
            [VRdata, VRdata_o, es_x] = VRWheelLoad_SLv2(animal, iseries, str2double(iexp_list(gg).name(end-15:end-13)));
            if gg == 2
                VR = combineTwoVRexpts(es, es_x);
            else
                VR = combineTwoVRexpts(VR, es_x);
            end
        end
        es = VR; clear VR
    end
    
    
    %% find behav session recorded during ePhys recording (compare sync pulse signals)
    [VR_SessionOI, ind_M] = SessionFinderForSyncPulseOEVR(es, AC_Info);
    
    [VRdata, VRdata_o, ~] = VRWheelLoad_SLv2(animal, iseries, VR_SessionOI);
    
    %% Setup trimming parameters for either 
    OE_TimeStamps = AC_Info.TimestampsDownsampled;
    OE_SamplingRate = AC_Info.SamplingRateOE/AC_Info.SamplingFactor;
    VR_SessionDuration = es.sampleTimes(max(find(es.iexp==VR_SessionOI))) - es.sampleTimes(min(find(es.iexp==VR_SessionOI))); % seconds
    OE_SessionDuration = max(OE_TimeStamps)-min(OE_TimeStamps);
    % find distance of max correlation point from the beginning of the VR
    % code
    VR_SessionStartDistance = abs(es.sampleTimes(ind_M) - es.sampleTimes(min(find(es.iexp==VR_SessionOI))));    % seconds
    VR_SessionEndDistance = abs(es.sampleTimes(ind_M) - es.sampleTimes(max(find(es.iexp==VR_SessionOI))));      % seconds
    OE_SessionStartDistance = OE_TimeStamps(round(length(OE_TimeStamps)*0.5))-min(OE_TimeStamps);               % seconds
    % get the trimming indexes for the ephys and the VR 
    if VR_SessionStartDistance < OE_SessionStartDistance % distance from end is the same as the comparing chunk is in the middle
        % if VR was started after the ephys recording start
        [v,OE_start_ind]    =  min(abs(OE_TimeStamps-(OE_TimeStamps(round(length(OE_TimeStamps)*0.5))-VR_SessionStartDistance)));
        VR_start_ind        =  min(find(es.iexp==VR_SessionOI));
        if VR_SessionEndDistance < OE_SessionStartDistance
            % if VR was stopped before the ephys recording end
            [v,OE_end_ind]  =  min(abs(OE_TimeStamps-(OE_TimeStamps(round(length(OE_TimeStamps)*0.5))+VR_SessionEndDistance)));
            VR_end_ind      =  max(find(es.iexp==VR_SessionOI));
        else
            % if VR was stopped after the ephys recording end
            [v,VR_end_ind]  =  min(abs(es.sampleTimes-(es.sampleTimes(ind_M)+OE_SessionStartDistance)));
            OE_end_ind      =  length(OE_TimeStamps);
        end
    else
        % if VR was started before the ephys recoring start
        OE_start_ind        = 1;
        [v,VR_start_ind]    = min(abs(es.sampleTimes-(es.sampleTimes(ind_M)-OE_SessionStartDistance)));
        if VR_SessionEndDistance < OE_SessionStartDistance
            % if VR was stopped before the ephys recording end
            [v,OE_end_ind]  =  min(abs(OE_TimeStamps-(OE_TimeStamps(round(length(OE_TimeStamps)*0.5))+VR_SessionEndDistance)));
            VR_end_ind      =  max(find(es.iexp==VR_SessionOI));
        else                                                                
            % if VR was stopped after the ephys recording end
            [v,VR_end_ind]  =  min(abs(es.sampleTimes-(es.sampleTimes(ind_M)+OE_SessionStartDistance)));
            OE_end_ind      =  length(OE_TimeStamps);
        end
    end
    OE_Indexes = [OE_start_ind OE_end_ind];
    VR_Indexes = [VR_start_ind VR_end_ind];
    
    figure
    title('Please check that the 2 sync pulses overlap nicely')
    plot(AC_Info.TimestampsDownsampled(OE_Indexes(1):(OE_Indexes(2)))-min(AC_Info.TimestampsDownsampled(OE_Indexes(1):(OE_Indexes(2)))),...
        AC_Info.SyncPulse(OE_Indexes(1):(OE_Indexes(2)))/max(AC_Info.SyncPulse(OE_Indexes(1):(OE_Indexes(2)))),'k')
    hold on
    plot(es.sampleTimes(VR_Indexes(1):(VR_Indexes(2)))-min(es.sampleTimes(VR_Indexes(1):(VR_Indexes(2)))),...
        es.SyncPulse(VR_Indexes(1):(VR_Indexes(2))),'r.')
    legend({'O.E. sync pulse','VR sync pulse'})
    xlabel('seconds'); ylabel('ON/OFF')
    
    pause(2)
    % 1 select shank from which taking the spike times
%     if nargin<8
%         inp = inputdlg('Please enter space-separated numbers of the group/shank to be loaded:','Sample', [1 20]);
%         igroup = str2num(inp{:});
%     end       
%     
    % 2 save all spike times
    for igroup_idx = 1:length(igroup)
        if exist('iaddinfo')
%             temp = getKwikSpikes_SL(animal, iseries, iexp, igroup(igroup_idx),iaddinfo);
            temp = getKiloSortSpikes_MM(animal, iseries, iexp, igroup(igroup_idx),iaddinfo);
        else
%             temp = getKwikSpikes_SL(animal, iseries, iexp, igroup(igroup_idx));
            temp = getKiloSortSpikes_MM(animal, iseries, iexp, igroup(igroup_idx));
        end
        tempNumCells = length(temp);
        if igroup_idx==1
            chans=temp;
            numCells = tempNumCells;
        else
            for icell = 1:tempNumCells
                chans(numCells+1) = temp(icell);
                numCells = numCells + 1;
            end
        end
    end

    % 3 select only the spiketimes happened during the OE-VR period (kwik spikes include all the spikes of the day!)   
    if ~isempty(iaddinfo)
        FileNameMetaDat = dir([DIRS.ePhys filesep num2str(iseries) filesep iaddinfo filesep '*.dat_meta.mat']);
    else
        FileNameMetaDat = dir([DIRS.ePhys filesep num2str(iseries) filesep '*.dat_meta.mat']);
    end
   
    %FileNameMetaDat = dir([DIRS.multichanspikes filesep num2str(iseries) filesep '*.dat_meta.mat']);
    metaInfo = load([FileNameMetaDat.folder filesep FileNameMetaDat.name]);
    RecordingNumber = find(metaInfo.lims==length(AC_Info.Data));
    EphysExtrems = [sum(metaInfo.lims(1:RecordingNumber-1))+1 sum(metaInfo.lims(1:RecordingNumber))];
    TimeFict = 1/AC_Info.SamplingRateOE:1/AC_Info.SamplingRateOE:sum(metaInfo.lims)/AC_Info.SamplingRateOE;
    TimeStart = TimeFict(EphysExtrems(1));
    TimeEnd   = TimeFict(EphysExtrems(2));
    gc = 1; clear TempChans
%     for nc = 1:length(chans)
%         if strcmp(chans(nc).id(5),'Good')|| strcmp(chans(nc).id(5),'MUA')
%            TempChans(gc) = chans(nc); 
%            gc = gc + 1;
%         end
%     end  
    % making it 'Not Noise' for looking at prelim. data when manual
    % clustering was not completed --2019-03 MM
    for nc = 1:length(chans)
        if ~strcmp(chans(nc).id(5),'Noise')
           TempChans(gc) = chans(nc); 
           gc = gc + 1;
        end
    end
    for gc = 1:length(TempChans)
        TempChans(gc).spiketimes = TempChans(gc).spiketimes((TempChans(gc).spiketimes>(TimeStart+OE_TimeStamps(OE_Indexes(1)))) & (TempChans(gc).spiketimes<(TimeStart+OE_TimeStamps(OE_Indexes(2)))));
        TempChans(gc).spiketimes = TempChans(gc).spiketimes - (TimeStart+OE_TimeStamps(OE_Indexes(1)));
    end
    
    % 4 trim the VR data 
    clear TempES
    esFieldNames = fieldnames(es);
    for fn = 1:length(esFieldNames)
        if length(es.(esFieldNames{fn}))== length(es.sampleTimes)
            TempES.(esFieldNames{fn}) = es.(esFieldNames{fn})(VR_Indexes(1):VR_Indexes(2));
        else
            TempES.(esFieldNames{fn}) = es.(esFieldNames{fn});
        end
    end
    TempES.sampleTimes = TempES.sampleTimes-min(TempES.sampleTimes);
    
    % 5 compute the firing rate in the same bins as the VR  
% in case we want to smooth the firing rate already now
%     sigma = 0.05/(max(TempES.sampleTimes)/length(TempES.sampleTimes)); % standard deviation in number of samples (converted from time in seconds)
%     Width = 0.05*3/(max(TempES.sampleTimes)/length(TempES.sampleTimes)); % convert size from seconds to number of samples
%     x_g = linspace(-Width/2, Width/2, Width);
%     gaussFilter = exp(-x_g.^2/(2*sigma^2));
%     gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize

    sampleTimesedges = [TempES.sampleTimes;TempES.sampleTimes(end)+TempES.sampleTimes(end)-TempES.sampleTimes(end-1)];
    for gc = 1:length(TempChans)
        TimeStamps = TempChans(gc).spiketimes;
        ifr = histcounts(TimeStamps,sampleTimesedges)/(max(TempES.sampleTimes)/length(TempES.sampleTimes));
        TempChans(gc).FiringRate = ifr';
        TempES.SpikeTrain(gc,:) = ifr';
        TempES.spikeIDs{gc} = chans(gc).id;
        TempES.chanIDs{gc} = chans(gc).ichan;
%         TempES.bestchan{gc} = chans(gc).bestchan;
%         TempES.waveform{gc} = chans(gc).waveform;
        try 
            TempES.ProbeIDs{gc} = chans(gc).ProbeIDs;
        catch
            TempES.ProbeIDs{gc} = 1;
        end
        TempES.spikeTrain(:,gc) = ifr(:);
        TempES.spikeIDs{gc} = chans(gc).id;
        TempES.chanIDs{gc} = chans(gc).ichan;
        TempES.ProbeIDs{gc} = 1;
%         TempES.bestchan{gc} = chans(gc).bestchan;
%         TempES.waveform{gc} = chans(gc).waveform;
        clear ifr
        % convolve firing rate with gaussian filter
        % TempChans(gc).FiringRate = conv(ifr', gaussFilter_, 'same');
    end
    es = TempES;

    ACInfo.OE_TimeStampsDS = AC_Info.TimestampsDownsampled(OE_Indexes(1):OE_Indexes(2))-min(AC_Info.TimestampsDownsampled(OE_Indexes(1):OE_Indexes(2)));
    ACInfo.OE_SyncPulse = AC_Info.SyncPulse(OE_Indexes(1):OE_Indexes(2));
    OE_PhDi = interp1qr(AC_Info.Timestamps-min(AC_Info.Timestamps),...
                                            AC_Info.Data(:,2),...
                                            AC_Info.TimestampsDownsampled);
    OE_speed(:,1) = interp1qr(AC_Info.Timestamps-min(AC_Info.Timestamps),...
                                            AC_Info.Data(:,3),...
                                            AC_Info.TimestampsDownsampled);
    OE_speed(:,2) = interp1qr(AC_Info.Timestamps-min(AC_Info.Timestamps),...
                                            AC_Info.Data(:,4),...
                                            AC_Info.TimestampsDownsampled);                                        
    ACInfo.OE_speed = OE_speed(OE_Indexes(1):OE_Indexes(2),:);
    ACInfo.OE_PhotoDiode = OE_PhDi(OE_Indexes(1):OE_Indexes(2));
    % 6 save all the data into processed file
    % new server format (from 2019)
     if exist([DIRS.Processed filesep num2str(iseries)],'dir') == 7
        save([DIRS.Processed filesep num2str(iseries) filesep animal '_' num2str(iexp) '_' iaddinfo '_spikeInfo.mat'],'TempChans');
        save([DIRS.Processed filesep num2str(iseries) filesep animal '_' num2str(iexp) '_' iaddinfo '_es.mat'],'es');
        save([DIRS.Processed filesep num2str(iseries) filesep animal '_' num2str(iexp) '_' iaddinfo '_OE.mat'],'ACInfo');
        VRDetails = VRdata.EXP;
        save([DIRS.Processed filesep num2str(iseries) filesep animal '_' num2str(iexp) '_' iaddinfo '_VRData.mat'],'VRDetails');
    else
        mkdir([DIRS.Processed filesep num2str(iseries)])
        save([DIRS.Processed filesep num2str(iseries) filesep animal '_' num2str(iexp) '_' iaddinfo '_spikeInfo.mat'],'TempChans');
        save([DIRS.Processed filesep num2str(iseries) filesep animal '_' num2str(iexp) '_' iaddinfo '_es.mat'],'es'); 
        save([DIRS.Processed filesep num2str(iseries) filesep animal '_' num2str(iexp) '_' iaddinfo '_OE.mat'],'ACInfo');
        VRDetails = VRdata.EXP;
        save([DIRS.Processed filesep num2str(iseries) filesep animal '_' num2str(iexp) '_' iaddinfo '_VRData.mat'],'VRDetails');
    end
    
    disp ('*****Saving done!*****')

%     % old server format (until 2018)
%     if exist([DIRS.data filesep animal filesep num2str(iseries)],'dir') == 7
%         save([DIRS.data filesep animal filesep num2str(iseries) filesep animal '_' num2str(iexp) '_spikeInfo.mat'],'TempChans');
%         save([DIRS.data filesep animal filesep num2str(iseries) filesep animal '_' num2str(iexp) '_es.mat'],'es');
%         save([DIRS.data filesep animal filesep num2str(iseries) filesep animal '_' num2str(iexp) '_OE.mat'],'ACInfo');
%         VRDetails = VRdata.EXP;
%         save([DIRS.data filesep animal filesep num2str(iseries) filesep animal '_' num2str(iexp) '_VRData.mat'],'VRDetails');
%     else
%         mkdir([DIRS.data filesep animal filesep num2str(iseries)])
%         save([DIRS.data filesep animal filesep num2str(iseries) filesep animal '_' num2str(iexp) '_spikeInfo.mat'],'TempChans');
%         save([DIRS.data filesep animal filesep num2str(iseries) filesep animal '_' num2str(iexp) '_es.mat'],'es'); 
%         save([DIRS.data filesep animal filesep num2str(iseries) filesep animal '_' num2str(iexp) '_OE.mat'],'ACInfo');
%         VRDetails = VRdata.EXP;
%         save([DIRS.data filesep animal filesep num2str(iseries) filesep animal '_' num2str(iexp) '_VRData.mat'],'VRDetails');
%     end
end 
end
end

