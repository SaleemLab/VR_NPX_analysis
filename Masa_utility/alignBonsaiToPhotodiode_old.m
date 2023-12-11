function [bonsai_data] = alignBonsaiToPhotodiode_old(bonsai_data,syncTimes_Photodiode,option)

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
else
    bonsai_idx = sort([find(diff(bonsai_data.quad)==-1)+1; find(diff(bonsai_data.quad)==1)+1])
    syncTimes_Quad = bonsai_data.sglxTime(bonsai_idx);
    mean_delay = mean(syncTimes_Quad - syncTimes_Photodiode);
    bonsai_data.corrected_sglxTime = bonsai_data.sglxTime - mean_delay;

%     for n = 1:length(bonsai_data.sglxTime)
%         this_stimulus_time = stimulus_time(n);
%         temp = find(syncTimes_Quad>this_stimulus_time & syncTimes_Quad<this_stimulus_time + 4);
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
end