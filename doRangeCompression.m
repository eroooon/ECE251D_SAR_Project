function [data, rnggrid] = doRangeCompression(rxsig, tpd, bw, fs)

% speed of light
speedOfLight = physconst('LightSpeed');

% create time steps to generate reference waveform coefficients
t = 0:1/fs:(tpd-1/fs);

% chirp rate, LFM signal characteristic
tau = bw/tpd;

% generate reference lfm waveform
ref = exp(1i*pi*tau*(t.*t));

% generate range compressed image
data = zeros(size(rxsig));
for ii=1:size(rxsig,2)

    % correlate each range line with the complex conjugate of the
    % transmitted chirp
    convdata = conv(rxsig(:,ii),fliplr(conj(ref)));
    data(:,ii) = convdata(end-size(rxsig,1)+1:end);
end

%rnggrid = (0:1:size(rxsig,1)-1)*speedOfLight/(2*fs);
rnggrid = (0:1:size(rxsig,1)-1)*speedOfLight/(2*tau*2*tpd);