function [xPower,freqs] = fGetAmpandPower(X,sr,method,p,freqs,showplot)
%   [xPower,xabsAmp,xFFT,freqs] = fGetAmpandPower(X,sr,method,p,freqs,showplot)
%
%   Inputs:
%       X       = signal, single or multidimensional
%                       signals along rows
%       sr      = sampling rate
%       method  = 'fft' or 'ar'
%       p       = order of the Autoregressive Model, a power of 2, default = 64
%       showplot= plot the amplitude and power
%
%
%
%  Outputs:
%       xPower      = Power Spectral Density
%       xabsAmps    = Absoluts Amps
%       xFFT        = raw FFT of signal
%       freqs       = index of frequencies of output
%

%% defaults
nyq = sr/2;
if ~exist('showplot','var');
    showplot = 0;
elseif isempty(showplot)
    showplot = 0;
end
if ~exist('method','var')
    method = 'burg';
elseif isempty(method)
    method = 'burg';
end
L = length(X);
if ~exist('p','var')
    powsof2 = 2.^(1:5);
    p = min(64,powsof2(find(powsof2<L,1,'last')));
elseif isempty(p)
    powsof2 = 2.^(1:5);
    p = min(64,powsof2(find(powsof2<L,1,'last')));
end
if ~exist('freqs','var')
    freqs = linspace(1,nyq,nyq*5);
elseif isempty(freqs)
    freqs = linspace(1,nyq,nyq*5);
end
if strcmpi(method,'fft')
   NFFT = 2^nextpow2(L);            % up to next power of two, what is my
   freqs = sr/2*linspace(0,1,NFFT/2+1); % index of frequencies    
end
if size(X,1) > size(X,2)
    X = X';
end
X = detrend(X')';
%% main
parfor i = 1:size(X,1)
    x = double(X(i,:));
    switch lower(method)
        case 'fft'
            L = length(x);
            NFFT = 2^nextpow2(L);            % up to next power of two, what is my 
            xfft = fft(x,NFFT)./L;            
            xPower(i,:) = real(xfft(1:NFFT/2+1)).^2;
        case 'burg'
            [xPower(i,:)] = pburg(x,p,freqs,sr);
        case 'ar'
            [xPower(i,:)] = pburg(x,p,freqs,sr);
        case 'cov'
            [xPower(i,:)] = pcov(x,p,freqs,sr);
    end % use Burg Autoregressive method
end % over signals
if showplot
   figure;plot(freqs,xPower);title('Power Spectrum');   
end % plot
end % function

