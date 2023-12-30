function [time,Data,indI,indF] = segmentData(t,eeg,tWindow,fs,n0)
L = length(eeg);
lWindow = fs*tWindow;
nWindows = floor(L/lWindow);
for n = 1:nWindows
    winInd = (n-1)*lWindow + (1:lWindow);
    time{n} = t(winInd);
    Data{n} = eeg(winInd);
    indI{n} = winInd(1)+n0-1;
    indF{n} = winInd(end)+n0-1;
end