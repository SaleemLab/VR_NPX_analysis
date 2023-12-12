function ttl = convertPDtoTimeStamps_FR(ttlData,thresh_pos,thresh_neg,sampleRate,options)
if ~isfield(options,'stim_dur')
    options.stim_dur = 3;
    warning('Using default options.stim_dur of 0.25s')
end
if options.paradigm==1
    %%%%%
    % White squares
    ttlind = find(ttlData>thresh_pos);
    ttl_diff = diff(ttlind);
    ttl_diff_ts = ttl_diff./sampleRate;
    upPhases = [ttlind(find(ttl_diff_ts > 1.5)+1)]; % find sample points in range of phase reversals

    %%%%
    % Black squares
    ttlind2 = find(ttlData<thresh_neg);
    ttl_diff2 = diff(ttlind2);
    c = find(ttl_diff2>10);
    %                 downPhases = ttlind2(c(diff(c)>100)+1); % SGS replaced with below on 19/03/2020
    downPhases = ttlind2(c(diff(c)>0.1*sampleRate)+1);
    
elseif strcmp(options.paradigm,'sn')
    % White squares
    ttlind = find(ttlData>thresh_pos);
    ttl_diff = diff(ttlind);
    ttl_diff_ts = ttl_diff./sampleRate;
    upPhases = [ttlind(find(ttl_diff_ts > 0.05)+1)]; % find sample points in range of phase reversals
    
    %%%%
    % Black squares
    ttlind2 = find(ttlData<thresh_neg);
    ttl_diff2 = diff(ttlind2);
    c = find(ttl_diff2>10);
    %                 downPhases = ttlind2(c(diff(c)>100)+1); % SGS replaced with below on 19/03/2020
    downPhases = ttlind2(c(diff(c)>0.1*sampleRate)+1);
    
elseif strcmp(options.paradigm,'masa')
    %%%%%
    % White squares
   
    ttlind = find(ttlData>thresh_pos);
    ttl_diff = diff(ttlind);
    ttl_diff_ts = ttl_diff./sampleRate;
    upPhases = [ttlind([1;find(ttl_diff_ts > 0.5*options.stim_dur)+1])]; % find sample points in range of phase reversals
    
    %%%%
    % Black squares
    ttlind2 = find(ttlData<thresh_neg);
    ttl_diff2 = diff(ttlind2);
    ttl_diff_ts2 = ttl_diff2./sampleRate;
    downPhases = [ttlind2(find(ttl_diff_ts2 > 0.4*options.stim_dur)+1)]; % find sample points in range of phase reversals
    
% Checking
%     start = 2030000;
%     finish = 2050000;
%     plot([start:finish],ttlData(start:finish),'k')
%     hold on
%     scatter(upPhases(find(upPhases>start&upPhases<finish)),400,'r')
%     scatter(downPhases(find(downPhases>start&downPhases<finish)),100,'b')
%    
end
%%%%%
if isfield(options,'paradigm')
    switch options.paradigm
        case 'srp'
            firstCrossing = ttlind(1);
            blockStarts = [firstCrossing; ttlind(find(ttl_diff_ts > 5)+1)]; % find sample points in range of block changes
            if length(blockStarts)>10
                warning('More than 10 blocks detected. Please check');
                %blockStarts(diff(blockStarts)<(200*sampleRate)) = [];
            end
            upPhases = sort(union(blockStarts,upPhases));
    end
end


%%%%
% Now try to find the rise in the peak from the lowest
% point and replace values
% Note that occassionally it can fail - so iterate the threshold in 5 steps between the
% the difference between the original theshold crossing (thresh_neg) and
% the minimum values
absMin = min(ttlData(upPhases(end-3):upPhases(end-2)))+0.02;
thresholds = linspace(absMin, thresh_neg, 5);
for thisThreshold = 1:length(thresholds)
    ttlind = find(ttlData<=thresholds(thisThreshold));
    IND = nearestpoint(upPhases, ttlind,'previous');
    tmpPhases = nan(length(upPhases),1);
    tmpPhases(~isnan(IND)) = ttlind(IND(~isnan(IND)));

    % Check to see if we are OK (ie the difference between the original
    % threshold crossing and the new one is not more than 50 ms, or 3 frames)
    if length(find (tmpPhases-upPhases < -1*0.05*sampleRate)) == 0
        break
    end
end

% If this is an SRP paradigm (presumabl though the same issue may
% apply in other paradigms
% The problem with the above is that the attempt to find the minimum
% can yield points that are NaNs (for the first TTL) or are
% displaced by 30s (for the subsequent blockstarts).
%
% Over write the first tmpPhase, and each subsequent
% tmpPhase that exceeds 5s, with the original upPhase
% point

% switch options.paradigm
%     case 'srp'
tp = find(diff(upPhases)>5*sampleRate);
tmpPhases([1; tp+1]) = upPhases([1; tp+1]);
% end

% Reassign
upPhases = tmpPhases;
upPhases = unique(upPhases);
upPhases(isnan(upPhases)) = [];

%%%
% Now do the same for downPhases
absMax = max(ttlData(upPhases(end-2):upPhases(end-1)))-0.02;
thresholds = linspace(absMax, thresh_pos, 5);
for thisThreshold = 1:length(thresholds)
    ttlind = find(ttlData>=thresholds(thisThreshold));
    IND = nearestpoint(downPhases, ttlind,'previous');
    tmpPhases = nan(length(downPhases),1);
    tmpPhases(~isnan(IND)) = ttlind(IND(~isnan(IND)));

    % Check to see if we are OK (ie the difference between the original
    % threshold crossing and the new one is not more than 50 ms, or 3 frames)
    if length(find (tmpPhases-downPhases < -1*0.05*sampleRate)) == 0
        break
    end
end

tp = find(diff(downPhases)>5*sampleRate);
tmpPhases([1; tp+1]) = downPhases([1; tp+1]);

% Reassign
downPhases = tmpPhases;
downPhases = unique(downPhases);
downPhases(isnan(downPhases)) = [];
%%%
% AP commented this out. Will address this issue at a later stage in the
% analysis.
% switch options.paradigm
%     case 'srp'
%         % Note we may have to remove the last upPhase to account for issues with returning to
%         % white photodiode at end of stimulus presentation - this should
%         % only be done if the last phase is an upPhase, and there is
%         % exactly one more upPhase than downPhases
%         if (length(upPhases) == 1001) || (length(upPhases) == 2001)
%             upPhases = upPhases(1:end-1);
%         end
%         % Also remove any downphases that occur before upPhases
%         downPhases(downPhases < min(upPhases)) = [];
% end

% Remove any downphases that occur before upPhases (that should probably be
% generic and not specific to the srp experiment (AP it's not)
switch options.paradigm
    case 'srp'
        downPhases(downPhases < min(upPhases)) = [];
end

%%%%
% Output
ttl.upPhases = upPhases;
ttl.downPhases = downPhases;


