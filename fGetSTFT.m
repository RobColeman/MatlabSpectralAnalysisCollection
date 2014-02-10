function [Pxx,Txx,Fxx,Sxx] = fGetSTFT(X,sr,F,Window,noverlap)


% Usage: [Pxx,Txx,Fxx] = fGetSTFT(X,sr,F,Window,noverlap)
%
%   Inputs: 
%       X: Matrix of signals, channels along rows
%       sr: sampling rate of signal
%       F: Frequencies to sample
%       Window: The window size of the STFT
%       noverlap: How much to overlap each window when doing STFT, must be
%                   strictly smaller than Window
%
%   Outputs:
%       Pxx: Power Spectral Estimate of Signals, i.e. spectrograms
%       Txx: Central timepoints of each Spectrogram sample
%       Fxx: Frequencies of Spectrograms
%       Sxx: Spetrograms
%
%

%% defaults for EEG

if size(X,1) > size(X,2)
    X = X'; % channels along rows
end

if ~exist('sr','var')
    sr = 1024;
elseif isempty(sr)
    sr = 1024;
end
if ~exist('F','var')
    F = 1:45; % integers from 1 to 45
elseif isempty(F)
    F = 1:45;
end

if ~exist('Window','var')
    Window = 128;
elseif isempty(Window)
    Window = 128;
end

if ~exist('noverlap','var')
    noverlap = 100;
elseif isempty(noverlap)
    noverlap = 100;
end

dsamp = 0;
nodetrend = 1;
%% main



for ch = 1:size(X,1)
    chX = double(squeeze(X(ch,:)));
    chX = fFilterEEGdata(chX,dsamp,sr,nodetrend); % filter data - don't downsample
    [Sx(ch,:,:),Fx,Tx,Px(ch,:,:)]=spectrogram(chX,Window,noverlap,F,sr);    
end % over singals

Pxx = squeeze(Px); % if we're only running one signal
Sxx = squeeze(Sx);
Txx = Tx;
Fxx = Fx;
end % function