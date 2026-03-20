function [resps] = perievent_LFP_response(stimTimes,AnalysisTimeWindow,options,varargin)

p = inputParser;
% addParameter(p,'doDetrend',false,@islogical);
addParameter(p,'filter_freq',[],@isnumeric);
% addParameter(p,'filter_freq',[0.5 3;4 12;9 17;30 60;60 100;125 300; 300 600],@isnumeric) % Frequency range: [freq1(1) freq1;freq2 freq2;......]
% addParameter(p,'CSD_V1_CA1_normalisation',0,@isnumeric) % Normalised CSD within region or not
% addParameter(p,'x_col',1,@isnumeric) % which x_col to use
addParameter(p,'selected_channels',[],@isnumeric) % which x_col to use
% addParameter(p,'LFP',[],@isnumeric) % which x_col to use
% addParameter(p,'best_channels',[],@isstruct) % which x_col to use


parse(p,varargin{:})
% filter_type = p.Results.filter_type;
filter_freq = p.Results.filter_freq;
selected_channels = p.Results.selected_channels;
extractedTimeWindow = AnalysisTimeWindow;
% extractedTimeWindow(1) = extractedTimeWindow(1)-5;
% extractedTimeWindow(2)= extractedTimeWindow(2)+5;


% 
% if isfield(clusters,'merged_spike_id')
%     clusters.cluster_id = unique(clusters.merged_cluster_id);
%     for ncell = 1:length(unique(clusters.merged_cluster_id))
%         tempt_peak_channel = clusters.peak_channel(clusters.merged_cluster_id == clusters.cluster_id(ncell));
%         tempt_peak_depth = clusters.peak_depth(clusters.merged_cluster_id == clusters.cluster_id(ncell));
%         tempt_peak_waveform = clusters.peak_channel_waveforms(clusters.merged_cluster_id == clusters.cluster_id(ncell),:);
%         tempt_cell_type = clusters.cell_type(clusters.merged_cluster_id == clusters.cluster_id(ncell));
% 
%         if length(tempt_peak_depth)== 2
%             tempt_peak_channel = tempt_peak_channel(1);
%             tempt_peak_depth = tempt_peak_depth(1);
%             tempt_peak_waveform = tempt_peak_waveform(1,:);
%             tempt_cell_type = tempt_cell_type(1);
%         else % find median peak depth assign that value to the unit
%             [~,index]= min(tempt_peak_depth - median(tempt_peak_depth));
%             tempt_peak_channel = tempt_peak_channel(index);
%             tempt_peak_depth = tempt_peak_depth(index);
%             tempt_peak_waveform = tempt_peak_waveform(index,:);
%             tempt_cell_type = tempt_cell_type(index);
%         end
% 
%         merged_peak_channel(ncell) = tempt_peak_channel;
%         merged_peak_depth(ncell) = tempt_peak_depth;
%         merged_peak_waveform(ncell,:) = tempt_peak_waveform;
%         merged_cell_type(ncell,:) = tempt_cell_type;
%     end
% 
%     clusters.peak_channel = merged_peak_channel;
%     clusters.peak_depth = merged_peak_depth;
%     clusters.peak_channel_waveforms = merged_peak_waveform;
%     clusters.cell_type = merged_cell_type;
% 
%     clusters.spike_id = clusters.merged_spike_id;
%     clusters.cluster_id = unique(clusters.merged_cluster_id);
% end

BinWidth = 1/1250;
% options.importMode = 'LF'; % LF or MUA or KS
% AnalysisTimeWindow two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])

%             options.ks_unitType = 'good'; % 'mua', 'good' or ''
nprobe = options.probe_id+1; % probe_no is [1,2] based on options.probe_id (0 and 1)

column = 1;
options.importMode = 'KS';
[file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column);% Since it is LF

if ~contains(imecMeta.acqApLfSy,'384,0') % NPX2 only has AP but NPX1 has AP and LF
    options.importMode = 'LF';
    [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column);% Since it is LF
    probe_type = 1;
    disp('Using lf.bin for NP1')
else

    DIR = dir(fullfile(options.EPHYS_DATAPATH,'*lf*'));
    if ~isempty(DIR)
        options.importMode = 'LF';
        [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column);% Since it is LF
        probe_type = 1;
        disp('Using preprocessed lf.bin for NP2')
    else
        probe_type = 2;
        disp('Using preprocessed ap.bin for NP2. Will take longer time to process')
    end
end




% Extract LFP data relative to the EVENT onset time
% Preallocate matrix for speed - so use the first stimulus to work
% out how much space we will need for the resampled data
first_stim = stimTimes(1);

%%%% OLD - use Neuropixels library
%         downSampleRate = fix(imec.fsLF*BinWidth); % the rate at which we downsample depends on the acquisition rate and target binwidth
%         tresps = imec.readLF_timeWindow([first_stim+AnalysisTimeWindow(1) first_stim+AnalysisTimeWindow(2)]);
%         % Design low pass filter (with corner frequency determined by the
%         % desired output binwidth
%         d1 = designfilt('lowpassiir','FilterOrder',12, ...
%             'HalfPowerFrequency',(1/BinWidth)/2,'DesignMethod','butter','SampleRate', imec.fsLF);

%%%% NEW - use SpikeGLX_Datafile_Tools
imecMeta = ReadMeta(fullfile(options.EPHYS_DATAPATH,file_to_use));
nSamp = fix(SampRate(imecMeta)*range(extractedTimeWindow));
chanTypes = str2double(strsplit(imecMeta.acqApLfSy, ','));
nEPhysChan = chanTypes(1);
downSampleRate = round(SampRate(imecMeta)*BinWidth); % the rate at which we downsample depends on the acquisition rate and target binwidth
% Design low pass filter (with corner frequency determined by the
% desired output binwidth)
d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',(1/BinWidth)/2,'DesignMethod','butter','SampleRate',SampRate(imecMeta));

% Additional filtering
if ~isempty(filter_freq)
    filter_type  = 'bandpass';
    filter_width = filter_freq;                 % range of frequencies in Hz you want to filter between
    filter_order = round(6*SampRate(imecMeta)/(max(filter_width)-min(filter_width)));  % creates filter for ripple
    norm_freq_range = filter_width/(SampRate(imecMeta)/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
    b_filter = fir1(filter_order, norm_freq_range,filter_type);
end

if ~isempty(selected_channels)
    % Read the data
    tresps = ReadBin(fix(SampRate(imecMeta)*first_stim), nSamp, imecMeta, file_to_use, options.EPHYS_DATAPATH);

    if strcmp(imecMeta.typeThis, 'imec')
        tresps = GainCorrectIM(tresps, 1:nEPhysChan, imecMeta);
    else
        tresps = GainCorrectNI(tresps, 1:nEPhysChan, imecMeta);
    end
    tresps= tresps(selected_channels,:); % Get rid of sync channel
else
    % Read the data
    tresps = ReadBin(fix(SampRate(imecMeta)*first_stim), nSamp, imecMeta, file_to_use, options.EPHYS_DATAPATH);
    tresps(nEPhysChan+1:end,:) = []; % Get rid of sync channel

    if strcmp(imecMeta.typeThis, 'imec')
        tresps = GainCorrectIM(tresps, 1:nEPhysChan, imecMeta);
    else
        tresps = GainCorrectNI(tresps, 1:nEPhysChan, imecMeta);
    end

end

% Downsample
tresps = downsample(tresps',downSampleRate)';
% Create vector for PSTH times at the downsampled rate (NB I am not
% sure this is completely correct - it depends on how the imec
% class does this [in the imec class, it does:
%   idxWindow = imec.closestSampleLFForTime(timeWindowSec);
%   sampleIdx = idxWindow(1):idxWindow(2);
timeVector = linspace(extractedTimeWindow(1),extractedTimeWindow(2),size(tresps,2));

% Initialise matrix for ephys data
resps = zeros([size(tresps,1) size(tresps,2) length(stimTimes)]);
% For each trial
H = waitbar(0,'Importing trials');
for thisTrial = 1:length(stimTimes)
    waitbar(thisTrial/length(stimTimes),H)
    % Establish what the stimulus time would be in ephys
    stim_on = stimTimes(thisTrial);

    % Get the LFP from IMEC file
    %%%% OLD - use Neuropixels library
    %             tresps = imec.readLF_timeWindow([stim_on+AnalysisTimeWindow(1) stim_on+AnalysisTimeWindow(2)]);
    %%%% NEW - use SpikeGLX_Datafile_Tools
    if ~isempty(selected_channels)
        % Read the data
        tresps = ReadBin(fix(SampRate(imecMeta)*stim_on), nSamp, imecMeta, file_to_use, options.EPHYS_DATAPATH);
        tresps(nEPhysChan+1:end,:) = []; % Get rid of sync channel

        if strcmp(imecMeta.typeThis, 'imec')
            tresps = GainCorrectIM(tresps, 1:nEPhysChan, imecMeta);
        else
            tresps = GainCorrectNI(tresps, 1:nEPhysChan, imecMeta);
        end
        
        tresps= tresps(selected_channels,:); % Get rid of sync channel
    else
        % Read the data
        tresps = ReadBin(fix(SampRate(imecMeta)*stim_on), nSamp, imecMeta, file_to_use, options.EPHYS_DATAPATH);
        tresps(nEPhysChan+1:end,:) = []; % Get rid of sync channel

        if strcmp(imecMeta.typeThis, 'imec')
            tresps = GainCorrectIM(tresps, 1:nEPhysChan, imecMeta);
        else
            tresps = GainCorrectNI(tresps, 1:nEPhysChan, imecMeta);
        end
    end

    % Zero-mean the data ...
    tresps = tresps-repmat(mean(tresps,2),1,size(tresps,2));
    % ... and lowpass filter ...
    if ~(contains(imecMeta.acqApLfSy,'384,0') & probe_type == 1) % if NP2 data but contains lf.bin, it is a preprocessed lf with low pass filter already applied
        tresps = filtfilt(d1,double(tresps)')';
    end

    if ~isempty(filter_freq)
        tresps = filtfilt(b_filter,1,tresps')';
    end
    % ... then downsample that data
    ttresps = downsample(tresps',downSampleRate)';
    % ... and save
    if size(ttresps,2) == size(resps,2)
        resps(:,:,thisTrial) = single(ttresps);
    else
        resps(:,:,thisTrial) = nan(size(resps,1),size(resps,2));
    end
end
close(H)

end
