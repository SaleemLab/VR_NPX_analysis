function [thisOEData, ttl, wheelData, fileinfo, timestamps] = readCONandPDandWheel(thisOEfileName,options)


% Function to load the open ephys data, photodiode pulses and wheel data.

% AP 10/07/19

if ~isfield(options,'chanNumbersToInclude') || isempty(options.chanNumbersToInclude)
    channels = -1;    % Channels number associated with the stimulus, to include in the export (-1 = all) or [1 2 3 4]
else
    channels = options.chanNumbersToInclude;
end
if length(channels) == 1 && channels == -1
    channels = 1:32;
end

%%%%%%%%%%%%%%%%---------------------------------%%%%%%%%%%%%%%%
% Fix sample rates
mainSampleRate = 30000;  % Main sample rate for OE files is 30 kHz

%%%%%%%%%%%%%%%%---------------------------------%%%%%%%%%%%%%%%
if options.multiplierDat ~= 1
    if rem(options.multiplierDat,1)
        error('Must use integers (-N...-1,1,...N) for multiplier')
    else
        fprintf('\nMultiplying Ephys signal by %01d    ',options.multiplierDat)
    end
end
if options.referenceDat
    fprintf('\nReferencing to average of other channels')
    if ~isempty(options.ignoreElectrodes)
        fprintf(', while ignoring Channel(s) ')
        fprintf('%01d,',options.ignoreElectrodes)
    end
    fprintf('...')
end

%% Set paths and memory
%folder_dir = dir(thisOEfileName);
%folder_name = folder_dir.name;
%oe_path = fullfile(conToppath,folder_name);
oe_path = thisOEfileName;
oe_dir = dir(oe_path);
oe_ch = cellfun(@any, regexp({oe_dir.name},'CH\d.*?.continuous'));

oe_ch_filenames = natsortfiles(cellfun(@(x) [oe_path filesep x], ...
    {oe_dir(oe_ch).name},'uni',false));

% dmem = memory; 
% memToLeaveFree = 4 * 2^30; % num of GB to keep free
% memToAllocate = dmem.MemAvailableAllArrays - memToLeaveFree;
% memToAllocate = max(0, memToAllocate);
% nint16s = memToAllocate/2;

%% Load the data
chs = oe_ch_filenames(channels);
nCH = length(chs);
% chunkSizeSamps = nint16s/nCH;
% nChunks = ceil(totalSamps/chunkSizeSamps);

% Load in first channel to set memory restrictions
fprintf(1, 'Reading channel 1 for initialization\n');
[data, timestamps, data_info] = load_open_ephys_data_faster(chs{1});
%[data, timestamps, data_info] = load_open_ephys_data(chs{1});
totalSamps = length(data);

% OEDataAllocate = [];


OEDataAllocate = zeros(nCH, totalSamps, 'int16');

for n = 1:nCH
    fprintf(1, 'reading CH %d/%d\n', n, nCH);
    
    if n>1 % can skip this for the first channel, it was loaded above
        [data, timestamps, info] = load_open_ephys_data_faster(chs{n});
    end
    
    OEDataAllocate(n,:) = int16(data);
    
end
    
% Invert if necessary
if options.multiplierDat ~= 1
    meanOEvals = int16(mean(OEDataAllocate,2));
    
    tmean = repmat(meanOEvals,1,size(OEDataAllocate,2));
    OEDataAllocate = OEDataAllocate-tmean;
    OEDataAllocate = OEDataAllocate.*options.multiplierDat;
    OEDataAllocate = OEDataAllocate+tmean;
end

% Process data if requested (subtract average of other electrodes)
if options.referenceDat
    thisDat = zeros(length(channels),size(OEDataAllocate,2),'int16');
    for thisChannel = 1:length(channels)
        
        % And not include 'dead' or otherwise corrupted channels
        [~, ignElectrodes] = setdiff(channels, options.ignoreElectrodes);
        
        % Comparison electrodes should only be those outside this
        % tetrode...
        refElectrodes = setdiff(1:length(channels),thisChannel);
        refElectrodes = intersect(refElectrodes, ignElectrodes);
        % Debug report
        fprintf('\nCh %01d: referenced to Chs',thisChannel)
        fprintf(' %01d',refElectrodes)
        
        % Create average of z-scores across included electrodes
        zscored= zscore(double(OEDataAllocate(refElectrodes,:))')';
        meanZscore = mean(zscored,1);
        mChDat = mean(OEDataAllocate(thisChannel,:));
        thisChDat = double(OEDataAllocate(thisChannel,:))-mChDat;
        % Regress against target channel
        weight = meanZscore'\thisChDat';
        % Subtract weighted version of mean
        thisDat(thisChannel,:) = int16((thisChDat - weight*meanZscore)+mChDat);
    end
    OEDataAllocate = thisDat;
end



thisOEData =  OEDataAllocate;
fprintf(1, 'done\n')


%% Load photodiote or TTL data
ttl.ttlStartTimes = [];

switch options.sync_input
    
    % Not sure if this code works but I keep it here in case it is ever needed
    % in the future
    case 'ttl'
        
        % Get events filename
        ttl_file = cellfun(@any, regexp({oe_dir.name},'all_channels.*?.events'));
        ttl_filename = [oe_path filesep oe_dir(ttl_file).name];
        [data, timestamps, info] = load_open_ephys_data_faster(ttl_filename);
        
        % Structure of file: data = channel
        ttl_timestamps = timestamps(data == (sync_channel - 1));
        
        % First timestamp is arbitrary with regards to when open ephys was
        % started, so subtract the first timestamp (should be the start of
        % recording timestamp, and be equal to the continuous file
        % timestamps(1)/sample rate)
        ttl_timestamps = ttl_timestamps - timestamps(1);
        
    case 'adc'
        
        % Get ADC sync channel filename
        oe_adc = cellfun(@any, strfind({oe_dir.name},['ADC' num2str(options.sync_channel)]));
        oe_adc_filename = [oe_path filesep oe_dir(oe_adc).name];
        nADC = 1;
        
        %         chunkSizeSamps = nint16s;
        %         nChunks = ceil(totalSamps/chunkSizeSamps);
        
        
        ttlData = zeros(nADC, totalSamps, 'int16');
        
        fprintf(1, 'reading PD ADC %d/%d\n', options.sync_channel, nADC);
        
        [data, timestamps, sync_info] = load_open_ephys_data_faster(oe_adc_filename);
        %[data, timestamps, sync_info] = load_open_ephys_data(oe_adc_filename);
        ttlData = int16(data);
        
        %ttlind = find(ttlData>1.5 | ttlData<0.5);
        %ttlind = find(ttlData>1.5);
        %ttl_diff = diff(ttlind);
        %ttlind2 = find(ttl_diff > 10000);
        %ttl_timestamps = timestamps([ttlind(1); ttlind(ttlind2+1)]);
        %ttl_ts = [ttlind(1); ttlind(ttlind2+1)]./mainSampleRate;
        %ttl_ts = ttl_ts(1:end-1);
        if options.photo_br ==1
            % White squares
            ttlind = find(ttlData>0.5);
            ttl_diff = diff(ttlind);
            ttlind2 = find(ttl_diff > 2500); % ~80ms difference
            upPhases = [ttlind(1); ttlind(ttlind2+1)];
            
            %Black squares
            ttlind = find(ttlData<0.5);
            ttl_diff = diff(ttlind);
            ttlind2 = find(ttl_diff > 2500);
            downPhases = [ttlind(1); ttlind(ttlind2(1:end-1)+1)];
            
        else
            
            % White squares
            ttlind = find(ttlData>0.5);
            ttl_diff = diff(ttlind);
            firstCrossing = ttlind(1);
            sampleRate = sync_info.header.sampleRate;
            ttl_diff_ts = ttl_diff./sampleRate;
            upPhases = [ttlind(find(ttl_diff_ts > 0.08 &  ttl_diff_ts < 1.02)+1)]; % find sample points in range of phase reversals
            
            
            
%             upPhases = [firstCrossing;upPhases];
            
            %Some times I get an extra tll at the end and others I don't, not
            %sure why but this is a quick (not good) fix ftm
            %          if mod(length(upPhases(find(upPhases==blockStarts(end)):end)),2)==1
            %              warning('One detected ttl was removed from the end');
            %              upPhases = upPhases(1:end-1);
            %          end
            
            % Black squares
            ttlind2 = find(ttlData<0.5);
            c = find(diff(ttlind2)>2);
            downPhases = ttlind2(c(diff(c)>2000)+1);
            downPhases = downPhases(2:end);
%            ttl.blockStarts = blockStarts;
        end
         
       %ttl.ttlStartTimes = ttl_ts;
        ttl.upPhases = upPhases;
        ttl.downPhases = downPhases;
        
end

%% Load wheel data

if options.includeWheel==1
        
        nADC = length(options.wheel_channels);
        wheelData = zeros(nADC, totalSamps, 'int16');
        for i=1:nADC
            
            % Get ADC wheel channel filename
            oe_adc = cellfun(@any, strfind({oe_dir.name},['ADC' num2str(options.wheel_channels(i))]));
            oe_adc_filename = [oe_path filesep oe_dir(oe_adc).name];
            
            fprintf(1, 'reading Wheel ADC %d/%d\n', i, nADC);
            [data, timestamps, wheel_info] = load_open_ephys_data_faster(oe_adc_filename);
            %[data, timestamps, wheel_info] = load_open_ephys_data(oe_adc_filename);
            wheelData(i,:) = int16(data);
            
        end
 
else
    wheelData = [];
end

% Save some info for later
fileinfo.data_info = data_info;
fileinfo.sync_info = sync_info;
if options.includeWheel==1
fileinfo.wheel_info = wheel_info;
end

end
