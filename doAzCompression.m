function [data] = doAzCompression(pcdata, fc, fs, prf, nPulses, R, maxRange, platSpeed)

speedOfLight = physconst('LightSpeed'); % speed of light

% Set up fast time grid
nRangeSamples = ceil((2*maxRange/speedOfLight)*fs);
fastTime = (0:1/fs:(nRangeSamples-1)/fs);

% Range frequency span and wavenumber
fRange = linspace(fc-fs/2,fc+fs/2,length(fastTime));
kRange = 2*(2*pi*fRange)/speedOfLight;

% convert sar to frequency domain
X = fftshift(fft2(pcdata));

% Cross-range wavenumber
kAz = 2*pi*linspace(-prf/2,prf/2,nPulses)./platSpeed;

% Matrix of the cross-range wavenumbers
kazimuth = kAz.';

% Build up matrix of normalized azimuth reference functions
kx = zeros(length(kazimuth),length(kRange));
ref = zeros(length(kazimuth),length(kRange));
for ii=1:length(kazimuth)
    for jj=1:length(kRange)
        kx(ii,jj) = sqrt(kRange(jj)^2-kazimuth(ii)^2);
        ref(ii,jj) = exp(-1i*R*kx(ii,jj));
        ref(ii,jj) = abs(kRange(jj))/sqrt(kRange(jj)^2+kazimuth(ii)^2)*ref(ii,jj);
    end
end

% Apply 2D matched filter
Y = (ref').*X;

% Stolt interpolation
out = Y.';
for ii=1:length(kazimuth)
        out(ii,:) = interp1(kx(ii,:),out(ii,:),kRange);
end
out(isnan(out)) = 1e-30;

% Applying last phase offset after interpolation to avoid interpolating
% zeros
out = out.*exp(-1i*R*kRange);

% ifft & return compressed image
data = ifft2(out);
