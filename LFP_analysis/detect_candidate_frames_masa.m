function [frames] = detect_candidate_events_masa(tvec,spike_times,time_epoches,bin_size,options)
%   based on Davidson, Kloosterman, and Wilson(2009,Neuron) criteria.  
%   - events detected by MUA only
%   - zscore of 3 minimum, and edges detected by zscore of 0.  
%     modificiation- If event is too long, then edge criteria shifted to
%     zscore of 0.25, then 0.5
%   - events within 50 ms are combined together
%   - replay events less than and equal 100ms are removed. 
% thus bayesian decoded candidate replay events have at least 5
% consecutive 20 ms bins, and 5 different neurons.

% things to decide:
%   - durations: shortest 40/50ms, above 100ms use for decoding. stored as
%   'reactivations' right now these are not merged or otherwise processed
%   - speed filter: 5cm/s
%   - theta filter (may have some bleeding from sharp wave band filtering) added 
%     for ranking: with csc high theta low ripple power
%   - adaptive thresholding: tag those that have been truncated and have a
%   look: start and end may have different min zscore as a result
options;

parameters= list_of_parameters;
SR = 1/mean(diff(tvec));

% Define Gaussian window for smoothing
windowWidth = 0.03; % seconds
gaussianWindow = gausswin(windowWidth/bin_size);
% Normalize to have an area of 1 (i.e., to be a probability distribution)
gaussianWindow = gaussianWindow / sum(gaussianWindow);

tvec = interp1(tvec,tvec,tvec(1):bin_size:tvec(end));
timevec_edge = (tvec(1)-(tvec(2)-tvec(1))/2....
    :bin_size:...
    tvec(end)+(tvec(end)-tvec(end-1))/2)';


mua= filtfilt(gaussianWindow,1,histcounts(spike_times,timevec_edge)')';

if ~isempty(time_epoches)
    [tvec_epoches,idx] = RestrictInts(tvec',time_epoches); %Replace with InInterval
    mua = mua(idx)';
end

end


