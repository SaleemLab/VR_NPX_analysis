function dataStruct = alignToSglxTime(dataStruct, sglx_sync_times, plotFlag)

%% generic function to align bonsai csv file with SGLX times usign async pulses

%% input vars
% 1) dataStruct is a struct which has fields 
% 'Sync' 0 and 1 of async pulses
% 'Time' which has a timebase, assumed to be in ms
% 2) sglx_sync_times is a vector of sync pulse times from sglx, assumed to
% be in s

if nargin<3
    plotFlag = false;
end

%% function

% find idx and times where sync pulse goes from 0 to 1 in bonsai data
idx = find(diff(dataStruct.Sync)==1);
syncTimes_bonsai = dataStruct.Time(idx+1)./1000; % convert from ms to s

% find peak in corrlation to align sync pulses from bonsai and sglx
[r, lags] = xcorr(diff(syncTimes_bonsai), diff(sglx_sync_times));
[val, idx] = max(r);
best_lag = lags(idx);


% next we use the sign of best lag in combination with the numbers of sync pulses
% available from sglx and bonsai to determine which programme started first
% and which programme ended first. We keep sync pulses from bonsai and sglx
% only when both are available in order to align time bases.

if best_lag < 0 
    % if sglx started first, we skip sglx sync pulses until bonsai sync
    % pulses are also available to align with.
    nSyncOffset = -best_lag+1;
    t_npix = sglx_sync_times(nSyncOffset:end); 
    
    if numel(syncTimes_bonsai)>=numel(t_npix)
        % if bonsai finished last, we keep bonsai sync pulses up to the
        % last sglx sync pulse
        t_bonsai = syncTimes_bonsai(1:numel(t_npix));
    else
        % if sglx finished last, we keep sglx sync pulses up to the last
        % bonsai sync pulse
        t_npix = t_npix(1:numel(syncTimes_bonsai));
        t_bonsai = syncTimes_bonsai;
    end
        
else % if bonsai started first we keep all the sglx sync pulses for now
    nSyncOffset = best_lag+1;
    t_npix = sglx_sync_times(1:end);
    
    if numel(syncTimes_bonsai(nSyncOffset:end))>=numel(t_npix) 
    % if sglx finished first, use bonsai sync pulses up to the last sglx time
        t_bonsai = syncTimes_bonsai(nSyncOffset:numel(t_npix)+nSyncOffset-1);
        
    else 
    % if bonsai finished first, we keep sglx sync pulses up to the last
    % bonsai sync time
        t_bonsai = syncTimes_bonsai(nSyncOffset:numel(syncTimes_bonsai));
        t_npix = t_npix(1:numel(t_bonsai));
    end
    
end

if plotFlag
% optionally ploting result of syncing
figure
plot(diff(t_npix),'k')
hold on
plot(diff(t_bonsai),'r:')
end

% now we do the actual interpolation of bonsai time to sglx time
% this asumes that data struct time is in ms and sglx sync times are in s
dataStruct.sglxTime = interp1(t_bonsai, t_npix, dataStruct.Time./1000,'linear', 'extrap'); 

end
