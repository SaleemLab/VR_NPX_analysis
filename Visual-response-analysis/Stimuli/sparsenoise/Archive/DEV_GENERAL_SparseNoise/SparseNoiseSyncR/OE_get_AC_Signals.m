%% Tomaso Muzzu - UCL - 12/12/2018

%% version for Bonvision stimulus

% Analog channel reader. Loads all the channels the photodiode, speed, and sync pulse
% signals

function ACInfo = OE_get_AC_Signals(directory)
% Select the folders where the OpenEphys data are stored
% Gets anologue channels at a rate of 1kHz. 
% ACInfo = getEPhysAnalogSignals(expObject, downsampled_rate); (in kHz)

    FoldersList = uigetfile_n_dir(directory,'Select the folder containing the recording of interest');
    ACInfo.ExpDate = FoldersList{1,1}(end-18:end);
    FileList = dir([FoldersList{1} filesep '101_ADC*.continuous']);
    if isempty(FileList)
        FileList = dir([FoldersList{1} filesep '100_ADC*.continuous']);
        if isempty(FileList)
            FileList = dir([FoldersList{1} filesep '103_ADC*.continuous']);
        end
    end
    % check how many channels were used and load them onto memory
    for j = 1:length(FileList)
        % Actually loads the analog channel
        [data(:,j), timestamps(:,j), info] = load_open_ephys_data_faster([FileList(j).folder filesep FileList(j).name]);
        fprintf(['\nLoading channel ' num2str(j) ' of ' num2str(length(FileList)) '.\n']);
        ACInfo.AnChannelsOE(j) = str2num(info.header.channel(end));  % physical channel used
    end
    ACInfo.Data = data; 
    ACInfo.Timestamps = timestamps(:,j)-min(timestamps(:,j));
    ACInfo.SamplingRateOE = info.header.sampleRate;    % sampling rate in Hz
    % check if there is more than one recording in the file (appended)
    SepIndexes = find(diff(ACInfo.Timestamps)>1);
    if ~isempty(SepIndexes)
        ACInfo.SessionStarts = SepIndexes;
    end
        
    % %%%%%%%% Keep sampling rate and 'channel' as metadata recorded into
    % the ACInfo
    
    % data(:,1) = sync pulse signal
    % data(:,2) = photodiode signal
    % data(:,3) = signal A from rotary encoder
    % data(:,4) = signal B from rotary encoder
    
    clear FileList info
end
   

    
    
