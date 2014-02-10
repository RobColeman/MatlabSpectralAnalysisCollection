function out = fBandpassButterworth(x,sr,LowPassElb,LowPassStop,HighPassElb,HighPassStop)
%
%
%   out = fBandpassButterworth(x,sr,LowPassElb,LowPassStop,HighPassElb,HighPassStop)
%
%
%   x: single channel of time series
%   sr: sampling rate

nyq = sr/2;
Rp = 1;
Rs = 48;

%% high pass filter
[nHi,WnHi]      = buttord(HighPassElb/nyq,HighPassStop/nyq,Rp,Rs);
[zHi,pHi,kHi]   = butter(nHi,WnHi,'high'); % design filter
[sosHi,gHi]     = zp2sos(zHi,pHi,kHi); % convert to SOS form
TdHi            = dfilt.df2tsos(sosHi,gHi); % make filter object
tbufHi          = filter(TdHi,x)'; % apply filter

%% lowpass filter
[nLo,WnLo]      = buttord(LowPassElb/nyq,LowPassStop/nyq,Rp,Rs);
[zLo,pLo,kLo]   = butter(nLo,WnLo,'low'); % design filter
[sosLo,gLo]     = zp2sos(zLo,pLo,kLo); % convert to SOS form
TdLo            = dfilt.df2tsos(sosLo,gLo); % make filter object
out             = filter(TdLo,tbufHi)'; % apply filter

end