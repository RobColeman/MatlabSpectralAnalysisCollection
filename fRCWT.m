function [cwt,cfreqs,WAHM]= fRCWT(x,wavelet,sr,scales,EV)
%
%   Usage:  [cwt cfreqs WAHM]= fRCWT(x,wavelet,sr,scales,EV)
%
%   Perfoms Continuous Wavelet Transform on 1D input signal 'x'
%   Inputs:
%       x:        Signal, dyadic length n=2^J, real-valued, ensure input is even
%       wavelet:  String 'Gauss', 'DerGauss','Sombrero', 'Morlet'
%       EV:       Calculate CWT envelope for outputs
%       scales:   Scales at which to perform CWT
%  Outputs:
%       cwt:      Wavelet transformed data, size length(scales) by length(x)
%       cfreqs:   Pseudo-Frequencies of each scale
%       WAHM:     Pseudo-Frequency Width of Wavelet Scale at Half Maximum
%                 response

%% defaults

if ~exist('scales','var');
    scales = linspace(9,150,32);
elseif isempty(scales)
    scales = linspace(9,150,32);
end

if ~exist('EV','var');
    EV = 1;
elseif isempty(EV)
    EV = 1;
end

if ~exist('intType','var');
    intType = 'Hermite';
elseif isempty(intType)
    intType = 'Hermite';
end

if ~exist('sr','var');
    sr = 256;
elseif isempty(sr)
    sr = 256;
end

if ~exist('wavelet','var');
    wavelet = 'Morlet';
elseif isempty(wavelet)
    wavelet = 'Morlet';
end

%% preparation
nscales = length(scales);
[nc n] = size(x);
if n < nc
    x = x';
end
[nc n] = size(x);
if mod(n,2) > 0
    x = x(1:end-1);
end
[nc n] = size(x);
omega0 = 5;

xhat = fft(x);
xi   = [ (0: (n/2)) (((-n/2)+1):-1) ] .* (2*pi/n);

cwt = nan(n,nscales);
cfreqs = nan(1,nscales);
WAHM = nan(2,nscales);
% scale to hz calc
ihzfactor = sr/n;
%% main loop

for s = 1:length(scales)
    S = scales(s);
    omega =  n .* xi ./ S ;
    
    % generate mother wavelet
    if strcmpi(wavelet,'Gauss'),
        motherwavelet = exp(-omega.^2 ./2);
    elseif strcmpi(wavelet,'DerGauss'),
        motherwavelet = i.*omega.*exp(-omega.^2 ./2);
    elseif strcmpi(wavelet,'Sombrero'),
        motherwavelet = (omega.^2) .* exp(-omega.^2 ./2);
    elseif strcmpi(wavelet,'Morlet'),
        motherwavelet = exp(-(omega - omega0).^2 ./2) - exp(-(omega.^2 + omega0.^2)/2);
    end % gen mother wavelet
    
    % Renormalization
    motherwavelet = motherwavelet ./ sqrt(S);
    what = motherwavelet .* xhat;
    w    = ifft(what);
    iCWT = real(w)';
    
    % Envelope interpolation
    if EV
        H = hilbert(iCWT);
        iCWT = real(H).^2+imag(H).^2;
    end
    cwt(:,s)  =  iCWT;
    % record peak freq of motherwavelet for Pseudo-Frequencies
    [Y,mind] = max(motherwavelet);
    %WAHM(1,s) = find(motherwavelet>Y/2,1,'first');
    %WAHM(2,s) = find(motherwavelet>Y/2,1,'last');
    cfreqs(1,s) = mind;
end % nvoice

cfreqs = cfreqs.*ihzfactor;
%WAHM = WAHM.*ihzfactor;
cwt = cwt';

end % function