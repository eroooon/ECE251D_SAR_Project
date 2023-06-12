function [data, rnggrid] = doRangeCompression(rxsig, c, fs)

speedOfLight = physconst('LightSpeed'); % speed of light

fftRes = size(rxsig,1) + length(c) + 1;

% fft of matched filter
w =  zeros(fftRes,1);
w(1:length(c)) = c;
C = fft(w,fftRes);

% output
cdata = zeros(size(rxsig));
for ii=1:size(rxsig,2)
    
    % fft of range swath
    swath = zeros(fftRes,1);
    swath(1:size(rxsig,1)) = rxsig(:,ii);
    S = fft(swath,fftRes);
    
    Y = S .* C;
    y = ifft(Y,fftRes);
    cdata(:,ii) = y(1:size(rxsig,1));
end

data = zeros(size(rxsig));
data(1:(size(rxsig,1)-length(c)+1),:) = cdata(length(c):size(rxsig,1),:);
rnggrid = (0:1:size(rxsig,1)-1)*speedOfLight/fs/2;