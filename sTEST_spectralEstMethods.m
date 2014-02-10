%% script - test fWaveletALAdmz 
clear all;
close all;
clc;

sr          = 256;
time        = 20; 
sigfreqs    = linspace(10,20,11); %[12 18];%
sl          = 1/3;
meanamp     = 100;
varamp      = 35; 

%% deterministic constants
n   = sr*time;
d   = length(sigfreqs);
x   = linspace(0,2*pi*time,n);
Y   = zeros(1,n);
blanks   = zeros(d,n); 
%% deterministic RVs
siglen  = randi(round(n*sl),d,1);
SoC     = randi(2,1,n)-1;
amp     = abs(rand(d,1)*varamp+meanamp);

for i = 1:d
    imax(i)     = (n-siglen(i))-1;    
    locs(i,1)   = randi(imax(i));
    locs(i,2)   = locs(i,1)+siglen(i);
end % over componants

count = 0;
for h = sigfreqs
    count = count+1;
    ihz     = h;
    isig    = x(locs(count,1):locs(count,2));

    if SoC(count) % sin or cos
        blanks(count,locs(count,1):locs(count,2))  = ...
            amp(count)*sin(ihz.*(isig));
    else
        blanks(count,locs(count,1):locs(count,2))  = ...
            amp(count)*sin(ihz.*(isig));
    end
end % over frequencies used

X = sum(blanks,1);

%% spectral estimation methods





% Morlet wavelet
scales = linspace(100,700,50);
envcwt = 1;
[cwt,cfreqs,WAHM]= fRCWT(X','morlet',sr,scales,envcwt);
% STFT 
F      = linspace(min(sigfreqs)-2,max(sigfreqs)+2,50);
[S,F,T,P] = spectrogram(X,[],[],F,sr);
% STFT for 12 and 18 Hz only
[S2,F2,T2,P2] = spectrogram(X,[],[],[12 18],sr);
% FFT for integrated power
[xPower,freqs] = fGetAmpandPower(X,sr,'fft',[],[],0);
xPower = xPower((min(sigfreqs)-2)<freqs & freqs < (max(sigfreqs)+2));
freqs = freqs((min(sigfreqs)-2)<freqs & freqs < (max(sigfreqs)+2));


% Bandpass butterworth filters
env = 1;
Filt1Hz      = [10 11; 13 14];
Filt2Hz      = [16 17; 19 20];
LowHiPass1      = [Filt1Hz(1,2) Filt1Hz(2,1)];
LowHiStop1      = [Filt1Hz(1,1) Filt1Hz(2,2)];
LowHiPass2      = [Filt2Hz(1,2) Filt2Hz(2,1)];
LowHiStop2      = [Filt2Hz(1,1) Filt2Hz(2,2)];

BFilters(2,:)    = fSSVEP_BandPass(X,sr,LowHiPass1,LowHiStop1,env,[],4,'butter');
BFilters(1,:)    = fSSVEP_BandPass(X,sr,LowHiPass2,LowHiStop2,env,[],4,'butter');

EFilters(2,:)    = fSSVEP_BandPass(X,sr,LowHiPass1,LowHiStop1,env,[],4);
EFilters(1,:)    = fSSVEP_BandPass(X,sr,LowHiPass2,LowHiStop2,env,[],4);

%% plot stiff
figure;
subplot(8,1,1);imagesc(flipud(blanks));title('Individual signals');
subplot(8,1,2);plot(X);title('Combined signals');
subplot(8,1,3);plot(freqs,xPower);title('Power via FFT');
subplot(8,1,4);imagesc(flipud(cwt));title('Morlet wavelet');
subplot(8,1,5);imagesc(flipud(P));title('Spectrogram via STFT');
subplot(8,1,6);imagesc(flipud(P2));title('Spectrogram via STFT - 12 and 18 hz only');
subplot(8,1,7);imagesc(BFilters);title('Bandpass Butterworth filters - 12 and 18 hz only');
subplot(8,1,8);imagesc(EFilters);title('Bandpass Elliptical filters - 12 and 18 hz only');

disp(sigfreqs)
disp(cfreqs);
%%
[BPower,freqs] = fGetAmpandPower(BFilters,sr,'fft',[],[],0);
[FPower,freqs] = fGetAmpandPower(BFilters,sr,'fft',[],[],0);
BPower  = BPower(freqs<64);
FPower  = FPower(freqs<64);
freqs   = freqs(freqs<64);
figure;
subplot(4,1,1);plot(BFilters');title('Bandpass Butterworth filters - 12 and 18 hz only');
subplot(4,1,2);plot(freqs,BPower);title('Power spectra of Butterworth envelope');
subplot(4,1,3);plot(EFilters');title('Bandpass elliptical filters - 12 and 18 hz only');
subplot(4,1,4);plot(freqs,FPower);title('Power spectra of elliptical envelope');
%% keyboard;
% %% generate filters
% 
% Filt1Hz      = [10 11; 13 14];
% LowHiPass1      = [Filt1Hz(1,2) Filt1Hz(2,1)];
% LowHiStop1      = [Filt1Hz(1,1) Filt1Hz(2,2)];
% nyq = sr/2;Rp = 1;Rs = 48;
% [n,Wn]      = buttord(LowHiStop1/nyq,LowHiPass1/nyq,Rp,Rs);
% [z,p,k]     = butter(n,Wn,'bandpass'); % design filter
% [sos,g]     = zp2sos(z,p,k); % convert to SOS form
% Filter12          = dfilt.df2tsos(sos,g); % make filter object
% fvtool(Filter12)
% 
% 
% 
% Filt2Hz      = [16 17; 19 20];
% LowHiPass2      = [Filt2Hz(1,2) Filt2Hz(2,1)];
% LowHiStop2      = [Filt2Hz(1,1) Filt2Hz(2,2)];
% nyq = sr/2;Rp = 1;Rs = 48;
% [n,Wn]      = buttord(LowHiStop2/nyq,LowHiPass2/nyq,Rp,Rs);
% [z,p,k]     = butter(n,Wn,'bandpass'); % design filter
% [sos,g]     = zp2sos(z,p,k); % convert to SOS form
% Filter18          = dfilt.df2tsos(sos,g); % make filter object
% fvtool(Filter18)
% 
% keyboard;
% 
% 
% 
%% elliptical filter design

% Filt1Hz         = [10 11; 13 14];
% LowHiPass1      = [Filt1Hz(1,2) Filt1Hz(2,1)];
% LowHiStop1      = [Filt1Hz(1,1) Filt1Hz(2,2)];
% nyq = sr/2;Rp = 1;Rs = 48;
% [n,Wp]          =ellipord(LowHiStop1/nyq,LowHiPass1/nyq,Rp,Rs);
% [z,p,k]         = ellip(n,Rp,Rs,Wp,'bandpass');
% [sos,g]         = zp2sos(z,p,k); % convert to SOS form
% EFilter12       = dfilt.df2tsos(sos,g);
% fvtool(EFilter12)
% 
% 
% Filt2Hz         = [16 17; 19 20];
% LowHiPass2      = [Filt2Hz(1,2) Filt2Hz(2,1)];
% LowHiStop2      = [Filt2Hz(1,1) Filt2Hz(2,2)];
% nyq = sr/2;Rp = 1;Rs = 48;
% [n,Wp]          =ellipord(LowHiStop2/nyq,LowHiPass2/nyq,Rp,Rs);
% [z,p,k]         = ellip(n,Rp,Rs,Wp,'bandpass');
% [sos,g]         = zp2sos(z,p,k); % convert to SOS form
% EFilter18       = dfilt.df2tsos(sos,g);
% fvtool(EFilter18)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% examples with EEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('SSVEPAT_Su12_17-Sep-2012_Data2.mat');
X = squeeze(data.EEG(31,257:end-256));
X = fFilterEEGdata(X,[],sr,[]);
clear data;

%% spectral estimation methods


% FFT for integrated power
[xPower,freqs] = fGetAmpandPower(X,sr,'fft',[],[],0);
xPower = xPower((min(sigfreqs)-2)<freqs & freqs < (max(sigfreqs)+2));
FFTfreqs = freqs((min(sigfreqs)-2)<freqs & freqs < (max(sigfreqs)+2));

% Bandpass butterworth filters
env = 1;
Filt1Hz      = [10 11; 13 14];
Filt2Hz      = [16 17; 19 20];
LowHiPass1      = [Filt1Hz(1,2) Filt1Hz(2,1)];
LowHiStop1      = [Filt1Hz(1,1) Filt1Hz(2,2)];
LowHiPass2      = [Filt2Hz(1,2) Filt2Hz(2,1)];
LowHiStop2      = [Filt2Hz(1,1) Filt2Hz(2,2)];

EEGBFilters(2,:)    = fSSVEP_BandPass(X,sr,LowHiPass1,LowHiStop1,env,[],4,'butter');
EEGBFilters(1,:)    = fSSVEP_BandPass(X,sr,LowHiPass2,LowHiStop2,env,[],4,'butter');

EEGEFilters(2,:)    = fSSVEP_BandPass(X,sr,LowHiPass1,LowHiStop1,env,[],4);
EEGEFilters(1,:)    = fSSVEP_BandPass(X,sr,LowHiPass2,LowHiStop2,env,[],4);

%%
[BPower,freqs] = fGetAmpandPower(EEGBFilters,sr,'fft',[],[],0);
[FPower,freqs] = fGetAmpandPower(EEGEFilters,sr,'fft',[],[],0);
BPower  = BPower(freqs<64);
FPower  = FPower(freqs<64);
freqs   = freqs(freqs<64);
figure;
subplot(6,1,1);plot(X);title('EEG');
subplot(6,1,2);plot(FFTfreqs,xPower);title('Power vis FFT');
subplot(6,1,3);plot(EEGBFilters');title('Bandpass Butterworth filters - 12 and 18 hz only');
subplot(6,1,4);plot(freqs,BPower);title('Power spectra of Butterworth envelope');
subplot(6,1,5);plot(EEGEFilters');title('Bandpass elliptical filters - 12 and 18 hz only');
subplot(6,1,6);plot(freqs,FPower);title('Power spectra of elliptical envelope');
