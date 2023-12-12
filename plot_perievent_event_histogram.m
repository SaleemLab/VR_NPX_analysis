function plot_perievent_event_histogram(event1,event2,varargin)

%INPUTS
%   all_spike_data      spike
%               
%   (options)
%   'lfp'               -A buzcode-style lfp structure... if you would
%                        rather just input the lfp instead of loading from
%                        basepath
%                           Default: load from basePath with bz_GetLFP
%   'spikes'            -A buzcode-style spike structure 
%                           Default: load from basePath with bz_GetSpikes
%   'NREMInts'          -Interval of times for NREM (seconds) 
%                        (Default: loaded from SleepState.states.mat, 
%                                   run SleepScoreMaster if not exist)
%                        use [0 Inf] to detect over all time points
%   'DetectionChannel'  -Channel with the most robust Slow Waves. (0-Indexing a la neuroscope). 
%                        (Default: 'autoselect')
%                        'useold' to use channel from existing SlowWaves.events.mat 
%                        If providing lfp via the 'lfp' input, make sure to
%                        give a detection channel here.
%   'noSpikes'          -true/false - set to true to not use spike information
%                        (default: false)
%   'MUAspikes'       -true/false - use MUA peaks (500-5000Hz) extracted 
%                        from the .dat file instead of spikes
%   'CTXChans'          -LFP channels that are in the cortex...  
%                        default: region 'CTX' from baseName.sessionInfo.mat or xml
%   'sensitivity'       -sensititivity (0-1) for determining LFP thresholds
%                        sensitivity for setting gamma/delta thresholds.
%                        lower sensitivity will result in fewer False Positives,
%                        but more Missed slow waves. (default 0.6)
%   'filterparms'       -filtering parameters structure with fields:
%           .deltafilter    [low high] bounds (default: [0.5 8]Hz)
%           .gammafilter    [low high] bounds (default: [100 400]Hz)
%           .gammasmoothwin  window for smoothing gamma power (default: 0.08 s)
%           .gammanormwin    window for normalizing gamma power (default: 20s)
%   'showFig'           -true/false show a quality control figure (default: true)
%   'saveMat'           -logical (default=true) to save in buzcode format
%   'forceReload'       -logical (default: false) to re-detect
%   'noPrompts'         -true/false disable any user prompts (default: false)
%
%
%OUTPUTS
%   SlowWaves    a buzcode structure
%   VerboseOut   extra output stuff for detection quality checks/figures
%
%
%
%DLevenstein 2016/2017
%If used, please cite: Levenstein et al 2018, currently on bioRxiv
%TO DO
%-incorporate multiple channels for detection of slow wave, which is robust
%on all (deep) lfp channels in the local cortical population
%-update input parameter list
%NOTE: requires 2017a or higher. sad. (functions: movingmad and movingmedian)
% Modified by Masahiro Takigawa 2023

% Default values
p = inputParser;
% addParameter(p,'group','by region',@isstr) % 
addParameter(p,'mode',1,@isnumeric) % Hanning window for pwelch analysis in seconds
% addParameter(p,'group_name',{},@iscell) % Frequency range: [freq1(1) freq1;freq2 freq2;......]
addParameter(p,'event_name','event',@isstr) % Frequency range: [freq1(1) freq1;freq2 freq2;......]
addParameter(p,'twin',[-1.2 1.2],@isnumeric) % time window around the event for analysis
addParameter(p,'plot_option',1,@isnumeric) % Powers Selected frequency for plotting


% assign parameters (either defaults or given)
parse(p,varargin{:});
twin = p.Results.twin;

MUA_filter_length = 50;
SD_alpha = 5; %2 std width
MUA_filter_alpha = (MUA_filter_length-1)/SD_alpha;
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,4
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width

bin_size = 0.01;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(event1, event2, twin, bin_size);
psth = filtfilt(w,1,mean(binnedArray./bin_size)); % normalize to Hz
psth_se = filtfilt(w,1,std(binnedArray./bin_size)./sqrt(length(event2)));

hold on
plot(bins,psth,'b')
hold on
% patch([time fliplr(time)], [Ymax fliplr(Ymin)], 'g')
patch([bins fliplr(bins)],[psth+psth_se fliplr(psth-psth_se)],'b','FaceAlpha','0.3','LineStyle','none');
hold on
plot([0 0],get(gca,'ylim'),'r-')


end

