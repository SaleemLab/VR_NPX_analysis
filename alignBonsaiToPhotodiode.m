function [bonsai_data] = alignBonsaiToPhotodiode(bonsai_data,syncTimes_Photodiode,option)

if strcmp(option,'replay') % If visual stimuli during sleep
    [pks,locs]= findpeaks(abs(bonsai_data.pos-200));

    stimulus_time = bonsai_data.sglxTime(locs);

    bonsai_idx = sort([find(diff(bonsai_data.quad)==-1); find(diff(bonsai_data.quad)==1)]);
    syncTimes_Quad = bonsai_data.sglxTime(bonsai_idx);

    for n = 1:length(stimulus_time)
        this_stimulus_time = stimulus_time(n);
        temp = find(syncTimes_Quad>this_stimulus_time & syncTimes_Quad<this_stimulus_time + 4);
        quad_after(n) = syncTimes_Quad(temp(1));

        temp = find(syncTimes_Photodiode>this_stimulus_time & syncTimes_Photodiode<this_stimulus_time + 4);
        photodiode_after(n) = syncTimes_Photodiode(temp(1));

        temp = find(syncTimes_Quad<this_stimulus_time & syncTimes_Quad>this_stimulus_time - 4);
        quad_before(n) = syncTimes_Quad(temp(end));

        temp = find(syncTimes_Photodiode<this_stimulus_time & syncTimes_Photodiode>this_stimulus_time - 4);
        photodiode_before(n) = syncTimes_Photodiode(temp(end));
    end

    delays_after = quad_after-photodiode_after;
    delays_before = quad_before - photodiode_before;
    mean_delays = mean([delays_after; delays_before],1);

    bonsai_data.stimuli_onset = bonsai_data.sglxTime(locs)-mean_delays';
    bonsai_data.stimuli_onset_raw = bonsai_data.sglxTime(locs);
    bonsai_data.stimuli_id = bonsai_data.pos(locs);
    bonsai_data.stimuli_track = bonsai_data.pos(locs);
    bonsai_data.stimuli_track(bonsai_data.pos(locs)<200) = 2;
    bonsai_data.stimuli_track(bonsai_data.pos(locs)>200) = 1;
elseif strcmp(option,'mean delay') % If visual stimuli during sleep
    bonsai_idx = sort([find(diff(bonsai_data.QuadState)==-1)+1; find(diff(bonsai_data.QuadState)==1)+1]);
    syncTimes_Quad = bonsai_data.sglxTime(bonsai_idx);
    mean_delay = nanmean(syncTimes_Quad - syncTimes_Photodiode);
    bonsai_data.corrected_sglxTime = bonsai_data.sglxTime - mean_delay;

    %     for n = 1:length(bonsai_data.sglxTime)
    %         this_stimulus_time = bonsai_data.sglxTime(n);
    %         temp = find(syncTimes_Quad>this_stimulus_time & syncTimes_Quad<this_stimulus_time + 0.5);
    %         quad_after(n) = syncTimes_Quad(temp(1));
    %
    %         temp = find(syncTimes_Photodiode>this_stimulus_time & syncTimes_Photodiode<this_stimulus_time + 4);
    %         photodiode_after(n) = syncTimes_Photodiode(temp(1));
    %
    %         temp = find(syncTimes_Quad<this_stimulus_time & syncTimes_Quad>this_stimulus_time - 4);
    %         quad_before(n) = syncTimes_Quad(temp(end));
    %
    %         temp = find(syncTimes_Photodiode<this_stimulus_time & syncTimes_Photodiode>this_stimulus_time - 4);
    %         photodiode_before(n) = syncTimes_Photodiode(temp(end));
    %     end
    %
    %     delays_after = quad_after-photodiode_after;
    %     delays_before = quad_before - photodiode_before;
    %     mean_delays = mean([delays_after; delays_before],1);

else
    bonsai_idx = sort([find(diff(bonsai_data.QuadState)==-1)+1; find(diff(bonsai_data.QuadState)==1)+1]);
    syncTimes_Quad = bonsai_data.sglxTime(bonsai_idx);

    
    % Especially in early recordings in a session there may be significant
    % delay between the start of Bonsai measurements and the start of ephys
    % measurements. This can make xcorr less reliable as there are a large
    % number of potential correspondences
    % We therefore first attempt to align the overall sync pulse traces ...
    [r, lags] = xcorr(diff(syncTimes_Quad), diff(syncTimes_Photodiode));
    [~, joint_idx] = max(r);
    best_lag = lags(joint_idx);
    % ...then remove the syncpulses from bonsai trace that occur before the
    % optimum lag
    if best_lag > 0
        syncTimes_Quad = syncTimes_Quad(syncTimes_Quad>syncTimes_Quad(best_lag));
    end
    % ...and rerun the xcorr
    [r, lags] = xcorr(diff(syncTimes_Quad), diff(syncTimes_Photodiode));
    [~, joint_idx] = max(r);
    best_lag = lags(joint_idx);
    % check which async pulse sequence is available first (ie. which was turned
    % on first)
    if best_lag < 0
        nSyncOffset = -best_lag+1;
        t_npix = syncTimes_Photodiode(nSyncOffset:end); % sync sglx times
        if length(t_npix) > length(syncTimes_Quad)
            t_npix = t_npix(1:length(syncTimes_Quad));
        end
        t_bonsai = syncTimes_Quad(1:numel(t_npix)); %sync bonsai times
    else
        nSyncOffset = best_lag+1;
        t_npix = syncTimes_Photodiode(1:end);
        if length(t_npix) > length(syncTimes_Quad)
            t_npix = t_npix(1:length(syncTimes_Quad));
        end
        t_bonsai = syncTimes_Quad(nSyncOffset:numel(t_npix)+nSyncOffset-1);
    end
     tt = interp1(t_bonsai, t_npix, bonsai_data.sglxTime,'linear','extrap');
%     tt = interp1(unique(bonsai_data.corrected_sglxTime),unique(bonsai_data.corrected_sglxTime),linspace(bonsai_data.corrected_sglxTime(1),bonsai_data.corrected_sglxTime(end),length(bonsai_data.corrected_sglxTime)),'linear')';
    [unique_t,index,~] = unique(tt);
     bonsai_data.corrected_sglxTime = interp1(index,unique_t,1:length(tt));
end