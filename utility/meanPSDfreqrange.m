function [P, f] = meanPSDfreqrange(pxx,fxx,F)
% finds mean psd power in frequency ranges
% F is size [n,2] with n = number of ranges
%
for n = 1:size(F,1)
    
    [~,f1idx] = min(abs(fxx-F(n,1)));
    [~,f2idx] = min(abs(fxx-F(n,2)));
    
    P(:,n) = mean(pxx(f1idx:f2idx,:));
    
    f(n,:) = [fxx(f1idx),fxx(f2idx)];
    
end