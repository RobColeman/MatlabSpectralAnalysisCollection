function X = fGenerateFakeEEG(length,numchannels,freq,noisemax)
%
%   X = fGenerateFakeEEG(length,numchannels,freq,noisemax)
%
%
%
%
%
%
%

X = 0:((2*pi)/length):(2*pi);
X = repmat(X,numchannels,1);

for i = 1:numchannels
   gnoise = normrnd(0,rand*noisemax,1,size(X,2));
   t = randi(3)-1;
   if t == 1
      X(i,:) = sin(freq*X(i,:))+gnoise;       
   elseif t == 2
      X(i,:) = cos(freq*X(i,:))+gnoise;
   else
      X(i,:) = X(i,:)+gnoise;
   end  
end

