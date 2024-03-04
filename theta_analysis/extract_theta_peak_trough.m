% Extract theta cycles peaks and troughs
% SD, 2020
% Select from two methods:
    % 1st - MATLAB FUNCTION
    % 2nd - Manual extraction
% OUTPUT: a matrix for peaks and troughs information. 1st column: indices of peaks/troughs; 2nd column: local max or local min; 3rd col: index for
% peak (1) or trough (-1); 4th coulmn: timestamps of each peak/trough (based on CSC time)

function extract_theta_peak_trough(method)

load extracted_CSC.mat
load extracted_position.mat

best_theta_lower_ripple_CSC =  find(strcmp({CSC.channel_label},'best_theta_low_ripple')==1);
lfp   = CSC(best_theta_lower_ripple_CSC).CSCraw;
theta = CSC(best_theta_lower_ripple_CSC).theta;
t     = CSC(best_theta_lower_ripple_CSC).CSCtime;
clear CSC

if method == 1
    % method 1 - using MATLAB findpeaks function
    [pks(:,2), pks(:,1)] = findpeaks(theta);
    [trs(:,2), trs(:,1)] = findpeaks(-theta);
    trs(:,2) = - trs(:,2);
    % clean both with 2 conditions:
    % - peaks must be positive, troughs negative
    % - it has to go peak-trough-peak etc
    extrema1 = clean_peaks_and_troughs(theta,pks,trs);
    
    % separate back into peaks and troughs
    theta_peaks   = extrema1(extrema1(:,3) == 1,:);
    theta_troughs = extrema1(extrema1(:,3) == -1,:);
    theta_peaks(:,4)  = t(theta_peaks(:,1));
    theta_troughs(:,4) = t(theta_troughs(:,1));

    
elseif method ==2
    % method 2 - using local extrema detection over sliding window
    [pks2,trs2,minChange] = findMinMax(theta, 0.25);
    
    % clean both with 2 conditions:
    % - peaks must be positive, troughs negative
    % - it has to go peak-trough-peak etc
    extrema2 = clean_peaks_and_troughs(theta,pks2,trs2);
    
    % separate back into peaks and troughs
    theta_peaks  = extrema2(extrema2(:,3) == 1,:);
    theta_troughs = extrema2(extrema2(:,3) == -1,:);
    theta_peaks(:,4)  = t(theta_peaks(:,1));
    theta_troughs(:,4) = t(theta_troughs(:,1));

end

save('theta_peak_trough.mat','theta_peaks','theta_troughs','-v7.3')

end