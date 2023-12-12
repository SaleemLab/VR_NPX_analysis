% Program to provide visualisation of log / peripherals files for Kilosort
% recordings to ensure that all dat been saved as appropriate
%
% SGS 03/04/2022
% Program to provide visualisation of log / peripherals files for
% neuropixel recordings to ensure that all data been acquired as appropriate
%
% SGS 03/04/2022

%%%%%%
SUBJECT = 'M22009';
SESSION = '20220412';
StimulusName = 'StaticTF'; % name of stimulus here or leave empty if all

%%%%%%
% Set up paths
if ismac
    ROOTPATH = '/Users/s.solomon/Filestore/Research2/ibn-vision';
else
    ROOTPATH = 'X:\ibn-vision\';
end

% Get  paths
DATAPATH = fullfile(ROOTPATH,'DATA','SUBJECTS');
BONSAI_DATAPATH = fullfile(DATAPATH,SUBJECT,'bonsai',SESSION);

% After 7th April 2022 we changed the way in shwich stimuli were
% encoded, and saved the photdiode signal, but before that we are totally reliant on QuadState
if etime([str2num(SESSION(1:4)) str2num(SESSION(5:6)) str2num(SESSION(7:8)) 00 00 00],[2022 04 07 00 00 00]) < 0 % Happened before April 7th
    loadflag = 'old';
else
    loadflag = 'new';
end

% If not specified stimulus file, find all with behavioural data
if isempty(StimulusName)
    switch(loadflag)
        case 'new'
            td = dir(fullfile(BONSAI_DATAPATH,'*WheelLog*.csv')); % >7th April 2022
        case 'old'
            td = dir(fullfile(BONSAI_DATAPATH,'*BehaviourData*.csv')); % For < 7th April 2022
    end
    tx = {}; [tx{1:length(td)}] = deal(td.name);
    StimulusName = strtok(tx,'_');
    %StimulusName = setdiff(StimulusName,{'Grey','SparseNoise'}); % Not include grey or sparsenoise yet
else
    if ischar(StimulusName)
        StimulusName = {StimulusName};
    end
end

for thisStimulus = 1:length(StimulusName)
    % Find csv files associated with desired stimulus
    fileTable = getBoundGAndBonsaiFileNames(SUBJECT,SESSION,DATAPATH);
    theseFiles = find(contains(fileTable.StimulusName,StimulusName{thisStimulus}));
    BehaviourDataFiles = fileTable.WheelLog(theseFiles);
    EyeDataFiles = fileTable.EyeLog(theseFiles);
    TrialParamsFiles = fileTable.Trial(theseFiles);
    PDFiles = fileTable.PD(theseFiles);
    gNums = fileTable.gNumber(theseFiles);
    if ~isempty(PDFiles)
        options.PD_FLAG = 1;    % Have saved the Photodiode output
    else
        options.PD_FLAG = 0;
    end

    % Set bonsai data paths
    PERIPHERALS_DATAPATH = fullfile(BONSAI_DATAPATH, BehaviourDataFiles{1});
    EYEDATA_DATAPATH = fullfile(BONSAI_DATAPATH, EyeDataFiles{1});
    TRIALDATA_DATAPATH = fullfile(BONSAI_DATAPATH,TrialParamsFiles{1});
    if ~isempty(PDFiles)
        PHOTODIODE_DATAPATH = fullfile(BONSAI_DATAPATH,PDFiles{1});
        pd_flag = 1;
    else
        pd_flag = 0;
    end
    
    %%%%%%
    % Import various files
    stimData = import_BonVisionParams(TRIALDATA_DATAPATH);
    peripherals = import_bonsai_peripherals(PERIPHERALS_DATAPATH);
    eyeData = import_bonsai_eyedata(EYEDATA_DATAPATH);
    % For when photodiode is encoded in the peripherals file
    if ~pd_flag
        photodiode.Time = peripherals.Time;
        photodiode.Photodiode = peripherals.Photodiode;
    else
        % ...or encoded separately
        [photodiode, photodiode_sync] = import_bonsai_photodiode(PHOTODIODE_DATAPATH);
%         [photodiode, photodiode_sync] = import_bonsai_photodiode(PHOTODIODE_ANALOG_DATAPATH,SYNCPULSE_ANALOG_DATAPATH);
%         photodiode_sync = alignSyncPulses(photodiode_sync,peripherals);
%         photodiode.Time = interp1(photodiode_sync.Time, photodiode_sync.realignedTime, photodiode.Time,'linear','extrap');
    end
    
    switch(loadflag)
        case 'new'
            % Convert all times to reference to first stimulus onset
            tt = peripherals.FrameComputerDateVec(1,:); tt(4:end)=0; % Get the date, make time = 0 at midnight
            wheelData.FrameTime = etime(peripherals.FrameComputerDateVec,tt)*1000; % Make frame time in ms since start of day
            eyeData.FrameTime = etime(eyeData.FrameComputerDateVec,tt)*1000;
            photodiode.FrameTime = etime(photodiode.FrameComputerDateVec,tt)*1000;
            
            switch(StimulusName{thisStimulus})
                case 'Grey'
                    ; % do nothing
                case 'SparseNoise'
                    % We would like to move towards encoding the quadstate change
                    % from the photodiode, but we are still limited to finding it
                    % in the stimData if present, or in the peripherals
                    idx_trial_start = find(diff(peripherals.QuadState)>0.5) + 1;
                    if length(idx_trial_start) ~= size(stimData,1)
                        % try to see if can retrieve by looking at down states too
                        try
                            idx_trial_start = sort(union( ...
                                find(diff(peripherals.QuadState)>0.5) + 1, ...
                                find(diff(peripherals.QuadState)<-0.5) + 1));
                        catch
                            error(sprintf('Expected %01d stimuli from TrialParams log, found %01d in Peripherals QuadState',size(stimData,1),length(idx_trial_start)))
                        end
                    end
                    % For sparsenoise, prior to 12th April 2022, we did not save the time
                    % of the stimulus in a Trial log. We can recover from the quadstate transitions recorded in wheel data though
                    stimData.FrameTime = etime(peripherals.FrameComputerDateVec(idx_trial_start,:),tt)*1000;
                    stimData.Time = stimData.FrameTime;
                otherwise
                    % DELETED FOLLOWING LINE
                    % stimData.Time = peripherals.Time(idx_trial_start); % in milliseconds
                    stimData.FrameTime = etime(stimData.FrameComputerDateVec,tt)*1000;
                    stimData.Time = stimData.FrameTime;
            end
            
        case 'old'
            %%%%%%
            % Get stimulus onset from QuadState in Peripherals (this also helps
            % check if the right number of stimuli)
            idx_trial_start = find(diff(peripherals.QuadState)>0.5) + 1;
            if length(idx_trial_start) ~= size(stimData,1)
                % try to see if can retrieve by looking at down states too
                try
                    idx_trial_start = sort(union( ...
                        find(diff(peripherals.QuadState)>0.5) + 1, ...
                        find(diff(peripherals.QuadState)<-0.5) + 1));
                catch
                    error(sprintf('Expected %01d stimuli from TrialParams log, found %01d in Peripherals QuadState',size(stimData,1),length(idx_trial_start)))
                end
            end
            
            %stimData.Time = stimData.Timestamp;
            switch(StimulusName{thisStimulus})
                case 'SparseNoise'
                    % For sparsenoise, prior to 12th April 2022, we did not save the time
                    % of the stimulus in a Trial log. We can recover from the quadstate transitions recorded in wheel data though
                    stimData.FrameTime = peripherals.Time(idx_trial_start);
                    stimData.Time = stimData.FrameTime;
            end
    end
    
    % Get wheel speed
    wheelData.pos = peripherals.Wheel;
    wheelData.Time = peripherals.Time;
    wheelData = import_processWheelEphys(wheelData,'gaussian',12);
    
    %%%%%%%%
    % Summarise the data
    %%%%%%%%
    
    % Step 1: plot the pupil area and wheel speed over the
    % whole session, and indicate when stimulus started on the same plot
    f1 = figure('Name',sprintf('%s :: %s :: %s :: peripheral data, time series',SUBJECT,SESSION,StimulusName{thisStimulus}));
    
    % Wheel speed output
    h1 = subplot(3,1,1);
    plot((wheelData.Time-wheelData.Time(1))/1000,wheelData.smthSpeed,'g-');
    set(gca,'TickDir','out','Box','off')
    ylabel('Speed')
    
    % Pupil area output
    h2 = subplot(3,1,2);
    plot((eyeData.Time-eyeData.Time(1))/1000,eyeData.Area,'g-');
    set(gca,'TickDir','out','Box','off')
    ylabel('Pupil')
    
    % Photodiode output
    h3 = subplot(3,1,3);
    plot((photodiode.Time-photodiode.Time(1))/1000,photodiode.Photodiode,'g-');
    set(gca,'TickDir','out','Box','off')
    ylabel('Photodiode')
    xlabel('Time (s)')
    
    switch(StimulusName{thisStimulus})
        case 'Gray'
            ; % do nothing
        otherwise
            % Plot stimulus onsets if there are a manageable number of them
            hold on;
            if length(stimData.Time) <= 100
                for kk = 1:length(stimData.Time)
                    plot(((stimData.Time(kk)-photodiode.Time(1))/1000)*ones(1,2),[min(photodiode.Photodiode) max(photodiode.Photodiode)],'r-');
                end
            else
                text(0.05,0.95,'Too many stimuli to plot - only plotting first and last','Units','normalized','Color','r')
                plot(((stimData.Time(1)-photodiode.Time(1))/1000)*ones(1,2),[min(photodiode.Photodiode) max(photodiode.Photodiode)],'r-');
                plot(((stimData.Time(end)-photodiode.Time(1))/1000)*ones(1,2),[min(photodiode.Photodiode) max(photodiode.Photodiode)],'r-');
            end
            hold off;
    end
    linkaxes([h1 h2 h3],'x');
    drawnow
    
    switch(StimulusName{thisStimulus})
        case 'Grey'
            ; % do nothing - not trial based
        otherwise    % Now plot average
            f2 = figure('Name',sprintf('%s :: %s :: %s :: peripheral data, triggered averages',SUBJECT,SESSION,StimulusName{thisStimulus}));
            obsWindow = [-0.2 0.5];    % s to view around stimulus
            ts_t = stimData.Time; % Time at which Bonsai recorded the output of the QuadState flip
            
            % Wheel speed
            h1 = subplot(2,3,1);
            tw = []; tt = [];
            for thisStim = 1:size(stimData,1)
                [~,idx0] = min(abs(wheelData.Time-(ts_t(thisStim) + obsWindow(1)*1000)));
                if thisStim == 1
                    [~,idx1] = min(abs(wheelData.Time-(ts_t(thisStim) + obsWindow(2)*1000)));
                    lIndx = idx1-idx0;
                end
                if idx0+lIndx <= length(wheelData.smthSpeed)
                    tw(thisStim,:) = wheelData.smthSpeed(idx0:idx0+lIndx);
                    tt(thisStim,:) = (wheelData.Time(idx0:idx0+lIndx)-ts_t(thisStim))/1000;
                end
            end
            imagesc(tt(1,:),1:size(tw,1),tw)
            set(gca,'TickDir','out','Box','off')
            title('Speed')
            
            h2 = subplot(2,3,4);
            trace = nanmean(tw,1);
            plot(tt(1,:),trace,'b-'); hold on
            plot([0 0],[min(trace) max(trace)],'r-'); hold off
            set(gca,'TickDir','out','Box','off')
            linkaxes([h1 h2],'x')
            
            % Pupil area output
            h1 = subplot(2,3,2);
            tw = []; tt = [];
            for thisStim = 1:size(stimData,1)
                [~,idx0] = min(abs(eyeData.Time-(ts_t(thisStim) + obsWindow(1)*1000)));
                if thisStim == 1
                    [~,idx1] = min(abs(eyeData.Time-(ts_t(thisStim) + obsWindow(2)*1000)));
                    lIndx = idx1-idx0;
                end
                if idx0+lIndx <= length(eyeData.Area)
                    tw(thisStim,:) = eyeData.Area(idx0:idx0+lIndx);
                    tt(thisStim,:) = (eyeData.Time(idx0:idx0+lIndx)-ts_t(thisStim))/1000;
                end
            end
            imagesc(tt(1,:),1:size(tw,1),tw)
            set(gca,'TickDir','out','Box','off')
            title('Pupil')
            
            h2 = subplot(2,3,5);
            trace = nanmean(tw,1);
            plot(tt(1,:),trace,'b-'); hold on
            plot([0 0],[min(trace) max(trace)],'r-'); hold off
            set(gca,'TickDir','out','Box','off')
            linkaxes([h1 h2],'x')
            xlabel('Time (ms)')
            
            % Photodiode output
            h1 = subplot(2,3,3);
            tw = []; tt = [];
            for thisStim = 1:size(stimData,1)
                [~,idx0] = min(abs(photodiode.Time-(ts_t(thisStim) + obsWindow(1)*1000)));
                if thisStim == 1
                    [~,idx1] = min(abs(photodiode.Time-(ts_t(thisStim) + obsWindow(2)*1000)));
                    lIndx = idx1-idx0;
                end
                if idx0+lIndx <= length(photodiode.Photodiode)
                    tw(thisStim,:) = photodiode.Photodiode(idx0:idx0+lIndx);
                    tt(thisStim,:) = (photodiode.Time(idx0:idx0+lIndx)-ts_t(thisStim))/1000;
                end
            end
            imagesc(tt(1,:),1:size(tw,1),tw)
            set(gca,'TickDir','out','Box','off')
            title('Photodiode')
            
            h2 = subplot(2,3,6);
            trace = nanmean(tw,1);
            plot(tt(1,:),tw,'c-'); hold on
            plot(tt(1,:),trace,'b-','Linewidth',2);
            plot([0 0],[min(trace) max(trace)],'r-'); hold off
            set(gca,'TickDir','out','Box','off')
            linkaxes([h1 h2],'x')
            drawnow
    end
end