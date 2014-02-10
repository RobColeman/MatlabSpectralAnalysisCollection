function out = fmdzCausalCWT(data,sr,frequencies,cwtlength)
%
%   Usage: out = fmdzCausalCWT(data,sr,frequencies,cwtlength)
%
%
%
%
%
%
%

if size(data,1) > size(data,2)
    data = data';
end
[nc,ns] = size(data);
cwtsamps = cwtlength*sr;



for e = 1:nc
    edata = squeeze(data(e,:));
    for f = frequencies
        fftcwtS = sin(f.*linspace(0,2*pi*cwtlength,cwtsamp));
        fftcwtC = Cos(f.*linspace(0,2*pi*cwtlength,cwtsamp));
        for t = 1:ns            
            tlast  = t;
            tfirst = tlast-cwtsamps+1;
            if tfirst < 1
                tfX = [zeros(1,abs(tfirst)+1)  edata(1:tlast)];
            else
                tfX = edata(tfirst:tlast);
            end
            
            efrespS = conv(tfX,motherS,'same'); 
            efrespC = conv(tfX,motherS,'same'); 
            

            out(e,f,t) = tr;
        end
    end % frequencies
end




end % function