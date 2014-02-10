function [out] = fEnvelopeCWT(X,intType)
%
%   
%
%
%   intType: Type of interpolation to use, string
%           'Hermite': Hermite cubic interpolaiton
%           'CubSplEndCond': Cubic spline with end conditions
%           'CubSpl': Cubic spline
%           'CubSmoothSpl': Cubic Smoothing spline
%           
%

if ~exist('intType','var');
    intType = 'Hermite';
end

%% find peaks
t = 1:length(X);
[~,pt] = findpeaks(X);py   = X(pt);
[~,vt] = findpeaks(-X);vy  = X(vt);

%% interpolate

switch lower(intType)
    case lower('Hermite');
        out = pchip(pt,py,t)-pchip(vt,vy,t); % difference between upper and lower envelope of CWT input
        
    case lower('CubSplEndCond');
        out = fnval(csape(pt,py),t)-fnval(csape(vt,vy),t);
        
    case lower('CubSpl');
        out = csapi(pt,py,t)-csapi(vt,vy,t);
end % interptype

end % function