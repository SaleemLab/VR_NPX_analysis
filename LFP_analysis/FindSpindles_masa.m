function [spindles] = FindSpindles_masa(lfp,timevec,varargin)
%FindRipples - Find cortical spindles (9-17Hz oscillations).
%
% USAGE
%    [spindles] = FindSpindles_masa(lfp.data,lfp.timestamps,<options>)
%    OR
%    [ripples] = bz_FindRipples(basepath,channel,<options>)
%
%    spindles are detected using the normalized squared signal (NSS) by
%    thresholding the baseline, merging neighboring events, thresholding
%    the peaks, and discarding events with excessive duration.
%    Thresholds are computed as multiples of the standard deviation of
%    the NSS. Alternatively, one can use explicit values, typically obtained
%    from a previous call.  The estimated EMG can be used as an additional
%    exclusion criteria
%
% INPUTS
%    lfp            unfiltered LFP (one channel) to use (nsample x 1)
%	 timestamps	    timestamps to match filtered variable
%    <options>      optional list of property-value pairs (see table below)
%
%    OR
%
%    basepath       path to a single session to run findRipples on
%    channel      	Ripple channel to use for detection
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'thresholds'  thresholds for ripple beginning/end and peak, in multiples
%                   of the stdev (default = [2 5]); must be integer values
%     'durations'   min inter-ripple interval and max ripple duration, in ms
%                   (default = [30 100]). 
%     'minDuration' min ripple duration. Keeping this input nomenclature for backwards
%                   compatibility
%     'restrict'    interval used to compute normalization (default = all)
%     'frequency'   sampling rate (in Hz) (default = 1250Hz)
%     'stdev'       reuse previously computed stdev
%     'show'        plot results (default = 'off')
%     'noise'       noisy unfiltered channel used to exclude ripple-
%                   like noise (events also present on this channel are
%                   discarded)
%     'passband'    N x 2 matrix of frequencies to filter for ripple detection 
%                   (default = [130 200])
%     'EMGThresh'   0-1 threshold of EMG to exclude noise
%     'saveMat'     logical (default=false) to save in buzcode format
%     'plotType'   1=original version (several plots); 2=only raw lfp
%    =========================================================================
%
% OUTPUT
%
%    ripples        buzcode format .event. struct with the following fields
%                   .timestamps        Nx2 matrix of start/stop times for
%                                      each ripple
%                   .detectorName      string ID for detector function used
%                   .peaks             Nx1 matrix of peak power timestamps 
%                   .stdev             standard dev used as threshold
%                   .noise             candidate ripples that were
%                                      identified as noise and removed
%                   .peakNormedPower   Nx1 matrix of peak power values
%                   .detectorParams    struct with input parameters given
%                                      to the detector
% SEE ALSO
%
%    See also bz_Filter, bz_RippleStats, bz_SaveRippleEvents, bz_PlotRippleStats.

% Copyright (C) 2004-2011 by Michaël Zugaro, initial algorithm by Hajime Hirase
% edited by David Tingley, 2017
% modified by Masahiro Takigawa, 2023

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Reference for spindle detection
% https://www.pnas.org/doi/10.1073/pnas.1805517115
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7363445/
% https://www.pnas.org/doi/10.1073/pnas.1103612108

% Based on Tingley and Buzsáki 2020

% Default values
p = inputParser;
addParameter(p,'thresholds',[1 3],@isnumeric)
addParameter(p,'durations',[500 3000],@isnumeric)
addParameter(p,'restrict',[],@isnumeric)
addParameter(p,'frequency',1000,@isnumeric)
addParameter(p,'stdev',[],@isnumeric)
addParameter(p,'show','on',@isstr)
addParameter(p,'noise',[],@ismatrix)
addParameter(p,'passband',[9 17],@isnumeric)
addParameter(p,'EMGThresh',.9,@isnumeric);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'plotType',1,@isnumeric)
addParameter(p,'savepath',[],@isstr)
addParameter(p,'behaviour',[],@isstruct)

lfp; 
timevec;

% assign parameters (either defaults or given)
parse(p,varargin{:});
frequency = p.Results.frequency;
show = p.Results.show;
restrict = p.Results.restrict;
sd = p.Results.stdev;
noise = p.Results.noise;
% ThresholdFactor = p.Results.thresholds;
lowThresholdFactor = p.Results.thresholds(1);
highThresholdFactor = p.Results.thresholds(2);

minSpindleDuration = p.Results.durations(1);
maxSpindleDuration = p.Results.durations(2);
plotType = p.Results.plotType;
passband = p.Results.passband;
EMGThresh = p.Results.EMGThresh;
savepath = p.Results.savepath;
behaviour = p.Results.behaviour;

filter_type  = 'bandpass';
filter_order = round(6*frequency/(max(passband)-min(passband)));  % creates filter for ripple
norm_freq_range = passband/(frequency/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_spindle = fir1(filter_order, norm_freq_range,filter_type);
signal = filtfilt(b_spindle,1,lfp);

%% Detect spindle and calculate noise
% zscored_spindle = zscore(abs(hilbert(signal)));
zscored_spindle = zscore(envelope(signal,round(frequency/5),'rms'));
% zscored_spindle = zscore(rms(signal,2));
zscored_spindle = smoothdata(zscored_spindle,'gaussian',round(frequency/5));
% lowThresholdFactor = 1.5;

if isfield(behaviour,'mobility_zscore')
    speed_threshold = 3; % zscore of three in terms of movement
    behaviour.speed=behaviour.mobility_zscore;
else
    speed_threshold = 5; %cm
end

speed = interp1(behaviour.sglxTime,behaviour.speed,timevec,'nearest');
speed_thresholded = speed<speed_threshold; % find events with speed below 5

if size(speed_thresholded,1)<size(speed_thresholded,2)
    speed_thresholded=speed_thresholded';
end

if size(zscored_spindle,1)<size(zscored_spindle,2)
    zscored_spindle=zscored_spindle';
end

% speed_thresholded
thresholded = zscored_spindle > lowThresholdFactor;
thresholded = thresholded + speed_thresholded;
thresholded = thresholded-1;
thresholded(thresholded ~= 1) = 0;

% Detect spindle periods by thresholding zscored ripple power
start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);
% Exclude last spindle if it is incomplete
if length(stop) == length(start)-1
	start = start(1:end-1);
end
% Exclude first spindle if it is incomplete
if length(stop)-1 == length(start)
    stop = stop(2:end);
end
% Correct special case when both first and last spindles are incomplete
if start(1) > stop(1)
	stop(1) = [];
	start(end) = [];
end
firstPass = [start,stop];
if isempty(firstPass)
    spindles.onset =[];
    spindles.offset =[];
    spindles.peak_zscore =[];
    spindles.peaktimes = [];
	disp('Detection by thresholding failed');
	return
else
	disp(['Detecting: ' num2str(length(firstPass)) ' spindle events after thresholding.']);
end



% Discard spindles with a peak power < highThresholdFactor
secondPass = [];
peakNormalizedPower = [];
for i = 1:size(firstPass,1)
    [maxValue,maxIndex] = max(zscored_spindle([firstPass(i,1):firstPass(i,2)]));
    if maxValue > highThresholdFactor
        secondPass = [secondPass ; firstPass(i,:)];
        peakNormalizedPower = [peakNormalizedPower ; maxValue];
    end
end

if isempty(secondPass)
    spindles.onset =[];
    spindles.offset =[];
    spindles.peak_zscore =[];
    spindles.peaktimes = [];
    disp('Peak thresholding failed.');
    return
else
    disp(['Detecting: ' num2str(length(secondPass)) ' spindle events after spindle peak thresholding.']);
end

% Calculate peak power
peakNormalizedPower = [];
for i = 1:size(secondPass,1)
    [maxValue,maxIndex] = max(zscored_spindle([secondPass(i,1):secondPass(i,2)]));
    peakNormalizedPower(i,1) = maxValue;
end

% Detect negative peak position for each spindle
peakPosition = zeros(size(secondPass,1),1);
for i=1:size(secondPass,1)
	[minValue,minIndex] = min(signal(secondPass(i,1):secondPass(i,2)));
	peakPosition(i) = minIndex + secondPass(i,1) - 1;
end

% Discard spindles that are way too long
spindles = [timevec(secondPass(:,1)) timevec(peakPosition) ...
           timevec(secondPass(:,2)) peakNormalizedPower];
duration = spindles(:,3)-spindles(:,1);
spindles(duration>maxSpindleDuration/1000,:) = NaN;
% disp(['Detecting: ' num2str(size(spindles,1)) ' spindle events after discarding short and long events.']);

% Discard spindles that are too short
spindles(duration<minSpindleDuration/1000,:) = NaN;
spindles = spindles((all((~isnan(spindles)),2)),:);
disp(['Detecting: ' num2str(size(spindles,1)) ' spindle events after discarding short and long events.']);

% If a noise channel was provided, find spindle-like events and exclude them
bad = [];
if ~isempty(noise)

    noiselfp = filtfilt(b_spindle,1,noise);
    zscored_noise = zscore(abs(hilbert(noise)));

	excluded = logical(zeros(size(spindles,1),1));
	% Exclude spindles when concomittent noise crosses high detection threshold
	for i = 1:size(spindles,1)

        time_index = find(timevec > spindles(i,1) & timevec < spindles(i,3));

		if any(zscored_noise(time_index(1):time_index(2))>5)
			excluded(i) = 1;
        end
	end
	bad = spindles(excluded,:);
	spindles = spindles(~excluded,:);
	disp(['Detecting: ' num2str(size(spindles,1)) ' spindles events after removing spindle-band noise.']);
end


% Optionally, plot results
if strcmp(show,'on')
    if plotType == 1
    	figure;
    	if ~isempty(noise)
            subplot(4,1,1)
            plot(timevec,signal)
            title('Spindle filtered signal LFP (9-17Hz)')
            subplot(4,1,2)
            plot(timevec,zscored_spindle)
            title('z-scored Spindle power (9-17Hz)')
            subplot(4,1,3)
            plot(timevec,noise)
            title('Spindle filtered noise LFP (9-17Hz)')
            subplot(4,1,4)
            plot(timevec,zscored_noise)
            title('z-scored Spindle power noise (9-17Hz)')
    	else
            subplot(4,1,1)
            plot(timevec,signal)
            subplot(4,1,2)
            plot(timevec,zscored_spindle)
        end

        figure
        spindle_first = round(size(spindles,1)/2);
        spindle_last = round(size(spindles,1)/2)+5;

        subplot(3,1,1)
        time_index = find(timevec>=spindles(spindle_first,1)-1 &timevec<=spindles(spindle_last,1)+1);
        plot(timevec(time_index),signal(time_index))
        title('spindle filtered LFP (9 - 17Hz)')

        subplot(3,1,2)
        plot(timevec(time_index),zscored_spindle(time_index))
        hold on
%         plot([timevec(time_index(1)) timevec(time_index(end))],[ThresholdFactor ThresholdFactor],'k','linestyle','--');
        plot([timevec(time_index(1)) timevec(time_index(end))],[highThresholdFactor highThresholdFactor],'k','linestyle','--');
        plot([timevec(time_index(1)) timevec(time_index(end))],[lowThresholdFactor lowThresholdFactor],'k');
        title('spindle filtered RMS zscore (9 - 17Hz)')

        for j=spindle_first:spindle_last
            hold on
            plot([spindles(j,1) spindles(j,1)],[min(zscored_spindle(time_index)) spindles(j,4)],'g-');
            plot([spindles(j,2) spindles(j,2)],[spindles(j,4) spindles(j,4)],'k-');
            plot([spindles(j,1) spindles(j,3)],[spindles(j,4) spindles(j,4)],'k-');
            plot([spindles(j,3) spindles(j,3)],[min(zscored_spindle(time_index)) spindles(j,4)],'r-');
        end
        
        subplot(3,1,3)
        plot(timevec(time_index),speed(time_index))
        title('Movement (pixel change zscored)')
        hold on
%         plot([timevec(time_index(1)) timevec(time_index(end))],[highThresholdFactor highThresholdFactor],'k','linestyle','--');
        for j=spindle_first:spindle_last
            hold on
            plot([spindles(j,1) spindles(j,1)],[min(speed(time_index)) lowThresholdFactor],'g-');
            plot([spindles(j,2) spindles(j,2)],[lowThresholdFactor lowThresholdFactor],'k-');
            plot([spindles(j,1) spindles(j,3)],[lowThresholdFactor lowThresholdFactor],'k-');
            plot([spindles(j,3) spindles(j,3)],[min(speed(time_index)) lowThresholdFactor],'r-');
        end

    end
end


%% BUZCODE Struct Output
temp = spindles; clear spindles

spindles.onset = temp(:,1);
spindles.offset = temp(:,3);
spindles.peaktimes = temp(:,2);            %peaktimes? could also do these as timestamps and then ripples.ints for start/stops?
spindles.peak_zscore = temp(:,4);  %amplitudes?
% ripples.stdev = sd;
% if ~isempty(bad)
%     ripples.noise.times = bad(:,[1 3]);
%     ripples.noise.peaks = bad(:,[2]);
%     ripples.noise.peakNormedPower = bad(:,[4]);
% else
%     ripples.noise.times = [];
%     ripples.noise.peaks = [];
%     ripples.noise.peakNormedPower = [];
% end

%The detectorinto substructure
detectorinfo.detectorname = 'Findspindles_masa';
detectorinfo.detectiondate = datetime('today');
detectorinfo.detectionintervals = restrict;
detectorinfo.detectionparms = p.Results;
detectorinfo.detectionparms = rmfield(detectorinfo.detectionparms,'noise');

%Put it into the spindle structure
spindles.detectorinfo = detectorinfo;

%Save
if p.Results.saveMat
    save(fullfile(savepath,'extracted_spindles_events.mat'),'spindles')
end


function y = Filter0(b,x)

if size(x,1) == 1
	x = x(:);
end

if mod(length(b),2)~=1
	error('filter order should be odd');
end

shift = (length(b)-1)/2;

[y0 z] = filter(b,1,x);

y = [y0(shift+1:end,:) ; z(1:shift,:)];

function [U,stdA] = unity(A,sd,restrict)

if ~isempty(restrict),
	meanA = mean(A(restrict));
	stdA = std(A(restrict));
else
	meanA = mean(A);
	stdA = std(A);
end
if ~isempty(sd),
	stdA = sd;
end

U = (A - meanA)/stdA;

