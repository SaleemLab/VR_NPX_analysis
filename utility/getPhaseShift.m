function [phaseShiftM,phaseShiftSe] = getPhaseShift(signal,freq1,freq2,fs,refChan)
%GETPHASESHIFT gets the phase shift down a neuropixel probe at a filtered
%signal between freq1 and 2.
% Give a column matrix of the phase difference between two rows of signal.
% Window - want to get the middle bunch of peaks to compare
window = 10; % Seconds
bin = ceil(window*fs);
binStart = floor((size(signal,1)/2)-(bin/2));
pbin = window*((freq2+freq1)/2); % Number of peaks I should be looking for 
% Had issues with using the butter function to create a filter (spitting
% out NaNs) so using this instead which seems to work well enough
d = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',3,'HalfPowerFrequency2',7, ...
    'SampleRate',fs);

for i = 1:size(signal,2)
    filtSignal(i,:) = filtfilt(d,signal(10000:end,i));
    phaseSignal(i,:) = rad2deg(angle(hilbert(filtSignal(i,:))));
%     [pks,locs] = findpeaks(phaseSignal(i,:),'MinPeakDistance',(1/freq2)*fs); % Has to have cycle of smallest freq 
%     idx = min(find(locs > binStart));
%     storeLocs(i,:) = locs(idx:idx+pbin-1);
end

% Then randomly choose 1000 time points to look at the phase across and
% look at the diff at each time point to get a spread
iterations = 1000;
for i = 1:iterations
    selectedPhases(:,i) = phaseSignal(:,randi(size(phaseSignal,2)));
    selectedPhases(:,i) = selectedPhases(:,i) - selectedPhases(refChan,i);
end

phaseShiftM = mean(selectedPhases');
phaseShiftSe = std(selectedPhases')/sqrt(size(selectedPhases,2));
end

