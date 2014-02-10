function [cwt cfreqs] = fPlotCWT(x,type,HzLow,HzHigh,sr,EV,intType,nvoice,oct,scale)
%
% fPlotCWT(X,HzLow,HzHigh,sr,EV,intType,nvoice,oct,scale)
%
%  Plots the CWT of a signal
%
% types:
%   imagesc
%   countourf
%   surf

if ~exist('type','var');
    type = 'imagesc';
elseif isempty(type)
    type = 'imagesc';
end

if ~exist('HzLow','var');
    HzLow = 4;
elseif isempty(HzLow)
    HzLow = 4;
end

if ~exist('HzHigh','var');
    HzHigh = 60;
elseif isempty(HzHigh)
    HzHigh = 60;
end

if ~exist('scale','var');
    scale = 4;
elseif isempty(scale)
    scale = 4;
end

if ~exist('oct','var');
    oct = 1;
elseif isempty(oct)
    oct = 1;
end

n = length(x);omega0 = 5;
noctave = floor(log2(n))-oct;

if ~exist('nvoice','var');
    nvoice = 16;
elseif isempty(nvoice)
    nvoice = 16;
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
    sr = 512;
elseif isempty(sr)
    sr = 512;
end

if ~exist('wavelet','var');
    wavelet = 'Morlet';
elseif isempty(wavelet)
    wavelet = 'Morlet';
end


%% Wavelet

[cwt cfreqs]= fRCWT(x,wavelet,sr,EV,intType,nvoice,oct,scale);

%% figure
H1 = figure;
switch lower(type)
    case lower('imagesc')
        H1 = imagesc(flipud(cwt(cfreqs > HzLow & cfreqs < HzHigh,:)));
    case lower('countourf')
        H1 = contourf(flipud(cwt(cfreqs > HzLow & cfreqs < HzHigh,:)));
    case lower('surf')
        H1 = surf(flipud(cwt(cfreqs > HzLow & cfreqs < HzHigh,:)));
end % type of plot


end % function