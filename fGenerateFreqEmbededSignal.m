function [out blanks] = fGenerateFreqEmbededSignal(t,sr,sigfreqs,siglength,muamp,varamp,noisevar,plotstuff)
%
%  Usage: [out blanks] = fGenerateFreqEmbededSignal(t,sr,sigfreqs,siglength,muamp,varamp,noisevar)
%
%
%
%
%
%
%
if ~exist('t','var');
    t = 20;
elseif isempty(t)
    t = 20;
end % sampling rate
if ~exist('sr','var');
    sr = 256;
elseif isempty(sr)
    sr = 256;
end % sampling rate
if ~exist('sigfreqs','var');
    sigfreqs = linspace(4,35,32);
elseif isempty(sigfreqs)
    sigfreqs = linspace(4,35,32);
end % sampling rate
if ~exist('siglength','var');
    siglength = 1/4;
elseif isempty(siglength)
    siglength = 1/4;
end % sampling rate
if ~exist('muamp','var');
    muamp = 100;
elseif isempty(muamp)
    muamp = 100;
end % sampling rate
if ~exist('varamp','var');
    varamp = 35;
elseif isempty(varamp)
    varamp = 35;
end % sampling rate
if ~exist('noisevar','var');
    noisevar = 10;
elseif isempty(noisevar)
    noisevar = 10;
end % sampling rate
if ~exist('plotstuff','var');
    plotstuff = 0;
elseif isempty(plotstuff)
    plotstuff = 0;
end % sampling rate


%% deterministic constants
n   = sr*t;
d   = length(sigfreqs);
x   = linspace(0,2*pi*t,n);
Y   = zeros(1,n);
blanks   = zeros(d,n);
%% deterministic RVs
siglen  = randi(round(n*siglength),d,1);
SoC     = randi(2,1,n)-1;
amp     = abs(rand(d,1)*varamp+muamp);

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

out = sum(blanks,1);
out = (out+randn(size(out))*noisevar^2);
out = fCenterSphereData(out);

if plotstuff
    figure;
    subplot(3,1,1);plot(blanks);title('Componants Signals');
    subplot(3,1,2);plot(sum(blanks,1));title('Combined Signals');
    subplot(3,1,2);plot(out);title('Combined Signals with added noise');
end


end % function
