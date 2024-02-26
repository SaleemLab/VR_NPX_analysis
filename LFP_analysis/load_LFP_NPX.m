function [raw_LFP tvec new_SR chan_config sorted_config] = load_LFP_NPX(options,column,varargin)


options.importMode = 'KS';
[file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column);% Since it is LF

if ~contains(imecMeta.acqApLfSy,'384,0') % NPX2 only has AP but NPX1 has AP and LF
    options.importMode = 'LF';
    [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column);% Since it is LF
    probe_type = 1;
else
    probe_type = 2;
end

% Default values
p = inputParser;
addParameter(p,'selected_channels',[chan_config.Channel],@isnumeric) % Select channels to load for LFP
addParameter(p,'start_sec',[1],@isnumeric) % Start timepoint in seconds (default is [1])
addParameter(p,'duration_sec',[imecMeta.fileTimeSecs-2],@isnumeric) % Duration of recording to load (default is whole session - 2 seconds)
addParameter(p,'desired_SR',[1250],@isnumeric) % Desired SR for downsampling. By default, it is 1250 Hz
addParameter(p,'probe_no',[1],@isnumeric) % Probe 1 or 2
addParameter(p,'probe_1_total_sample',[],@isnumeric) % total number of LF band sample

% addParameter(p,'stdev',[],@isnumeric)
% addParameter(p,'show','on',@isstr)
% addParameter(p,'noise',[],@ismatrix)
% addParameter(p,'passband',[125 300],@isnumeric)
% addParameter(p,'EMGThresh',.9,@isnumeric);
% addParameter(p,'saveMat',false,@islogical);
% addParameter(p,'minDuration',20,@isnumeric)
% addParameter(p,'plotType',1,@isnumeric)
% addParameter(p,'savepath',[],@isstr)

% assign parameters (either defaults or given)
parse(p,varargin{:});
selected_channels = p.Results.selected_channels;
start_sec = p.Results.start_sec;
duration_sec = p.Results.duration_sec;
desired_SR = p.Results.desired_SR;
probe_no = p.Results.probe_no;
probe_1_total_sample = p.Results.probe_1_total_sample;

% Load Data
BinWidth = 1/desired_SR;
% start_sec = 1;
start_samp = round(start_sec*imecMeta.imSampRate);
% duration_sec = imecMeta.fileTimeSecs-2;
nSamp     = round(duration_sec*imecMeta.imSampRate);
chanTypes = str2double(strsplit(imecMeta.acqApLfSy, ','));
nEPhysChan = chanTypes(1);
downSampleRate = fix(imecMeta.imSampRate*BinWidth); % Donwsample rate

% Read the data
% cd(options.ANALYSIS_DATAPATH)
clipDur = 10; % seconds
nClipSamps = round(imecMeta.imSampRate*clipDur);
nClips = floor(nSamp/nClipSamps);
samples_to_pass = 0;

% imecMeta.selected_channels = selected_channels;


% % 
% raw_LFP = [];
% for clip = 1:10
%     tic
%     disp(sprintf('loading clip %i',clip))
%     temp_LFP = ReadNPXBin(start_samp+samples_to_pass, nClipSamps, imecMeta, file_to_use, options.EPHYS_DATAPATH);
%     temp_LFP(nEPhysChan+1:end,:) = []; % Get rid of sync channel
%     temp_LFP = GainCorrectIM(temp_LFP, 1:nEPhysChan, imecMeta);
% 
%     % % Downsample the data
%     raw_LFP = [raw_LFP downsample(temp_LFP(selected_channels,:)',downSampleRate)'];
%     samples_to_pass = samples_to_pass + nClipSamps;
% %     samples_to_pass = samples_to_pass + nClipSamps + 1;
%     toc
% end


% Loading LFP using neuropixel utilities function (load slightly faster)
DIR = dir(fullfile(options.EPHYS_DATAPATH,'*ChanMap.mat'));
if isempty(DIR) % If empty load ephys meta 
    metafile = dir(fullfile(options.EPHYS_DATAPATH,'*.ap.meta'));
    SGLXMetaToCoords_ChannelMap_masa(metafile);
    DIR = dir(fullfile(options.EPHYS_DATAPATH,'*ChanMap.mat'));
end

channel_map_filename = fullfile(DIR.folder,DIR.name);

imec = Neuropixel.ImecDataset(options.EPHYS_DATAPATH,'ChannelMap',channel_map_filename);
samples_to_pass = 0;
raw_LFP = [];

disp(sprintf('%i number of %i clip to process',nClips,clipDur));
for clip = 1:nClips
    tic
    disp(sprintf('loading clip %i',clip))

    timeWindow = [1+start_samp+samples_to_pass:start_samp+samples_to_pass+nClipSamps]; % in seconds
    %     [temp_LFP, ~] = imec.readAP_timeWindow(timeWindow);
    if probe_type == 1
        [temp_LFP] = double(imec.readLF_idx(timeWindow)); % if lf.bin load LF data, if ap.bin load AP data (for NPX2, only ap)
    else
        [temp_LFP] = double(imec.readAP_idx(timeWindow)); % if lf.bin load LF data, if ap.bin load AP data (for NPX2, only ap)
        %     temp_LFP = GainCorrectIM(temp_LFP, 1:nEPhysChan, imecMeta);
    end
    % % Downsample the data
    raw_LFP = [raw_LFP downsample(temp_LFP(selected_channels,:)',downSampleRate)'];
    %     time_to_pass = time_to_pass + clipDur;
    samples_to_pass = samples_to_pass + nClipSamps;
    toc
end

new_SR = imecMeta.imSampRate/downSampleRate;

if probe_no == 2
    probe_1_tvec = (start_samp:nClipSamps*nClips+start_samp-1)/imecMeta.imSampRate;

    %     If downsized to 1250 Hz, then 1-2ms drift is not an issue (for now record the probe 2 tvec in terms of probe 1 times)
    [~,fname] = fileparts(options.EPHYS_DATAPATH);
    load(fullfile(options.EPHYS_DATAPATH,[fname,'_aligned_AP_sample_number.mat']));
    aligned_LFP_sample_number = round(aligned_AP_sample_number(1:12:end)/12); % LFP is every 12 sample points of the AP band data
    aligned_LFP_tvec = aligned_LFP_sample_number(start_samp:nClipSamps*nClips+start_samp-1)/imecMeta.imSampRate;
    
    for nchannel = 1:length(selected_channels)
        raw_LFP(nchannel,:) = interp1(aligned_LFP_tvec(1:2:end),raw_LFP(nchannel,:),probe_1_tvec(1:2:end),'linear'); %probe 1 LFP time to match probe 2 LFP time
        if sum(isnan(raw_LFP(nchannel,:)))>0 % if there are nan points due to interpolation. replace it with the previous datapoint
            nan_index = find(isnan(raw_LFP(nchannel,:)));
            diff_nan_index = diff(nan_index); non_continuous_index = find(diff_nan_index>1);
            for n = 1:length(nan_index)
                if nan_index(n) == 1 && isempty(non_continuous_index)
                    raw_LFP(nchannel,nan_index(n)) = raw_LFP(nchannel,nan_index(end)+1); %temp solution to nans in the beginning
                else
                raw_LFP(nchannel,nan_index(n)) = raw_LFP(nchannel,nan_index(n)-1);
                end
            end
        end
    end

end

tvec = start_samp/imecMeta.imSampRate:1/new_SR:start_sec+(length(raw_LFP(1,:))-1)/new_SR;
end

