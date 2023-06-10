function azcompresseddata = helperRangeMigration(sigData,fastTime,fc,fs,prf,speed,numPulses,c,Rc)
% This function demonstrates the range migration algorithm for imaging the
% side-looking synthetic aperture radar. The pulse compressed synthetic
% aperture data is considered in this algorithm.

% Set the range frequency span.
frequencyRange = linspace(fc-fs/2,fc+fs/2,length(fastTime));
krange = 2*(2*pi*frequencyRange)/c;

% Set the cross-range wavenumber.
kaz = 2*pi*linspace(-prf/2,prf/2,numPulses)./speed;

% Generate a matrix of the cross-range wavenumbers to match the size of
% the received two-dimensional SAR signal.
kazimuth = kaz.';
kx = krange.^2-kazimuth.^2;

% Set the final wavenumber to achieve azimuth focusing.
kx = sqrt(kx.*(kx > 0));
kFinal = exp(1i*kx.*Rc);

% Perform a two-dimensional FFT on the range compressed signal.
sdata =fftshift(fft(fftshift(fft(sigData,[],1),1),[],2),2);

% Perform bulk compression to get the azimuth compression at the reference
% range. Perform filtering of the 2-D FFT signal with the new cross-range
% wavenumber to achieve complete focusing at the reference range and as a
% by-product, partial focusing of targets not lying at the reference range.
fsmPol = (sdata.').*kFinal;

% Perform Stolt interpolation to achieve focusing for targets that are not
% lying at the reference range.
stoltPol = fsmPol;
for i = 1:size((fsmPol),1)
    stoltPol(i,:) = interp1(kx(i,:),fsmPol(i,:),krange(1,:));
end
stoltPol(isnan(stoltPol)) = 1e-30;
stoltPol = stoltPol.*exp(-1i*krange.*Rc);
azcompresseddata = ifft2(stoltPol);

end
