% Source:
% https://www.mathworks.com/help/radar/ug/stripmap-synthetic-aperture-radar-sar-image-formation.html
clear all; close all;

% Constants
c = physconst('LightSpeed'); % speed of light

%% Radar Configuration

% Input parameters
fc        = 4e9;  % SAR center frequency
rangeRes  = 3;    % desired RR, m
xRangeRes = 3;    % desired CRR, m
prf       = 1000; % pulse repitition frequency
                  % PRF determines max unambiguous range and is the cross-
                  % range sampling frequency. PRF <= 2*fd
tpd       = 3*10^-6;  % pulse width
fs        = 120*10^6; % sample rate
aperture  = 4;  % sq aperterture, m^2
maxRange  = 2500; % m
Rf = 1000; % reference range, m

% Platform parameters
speed   = 100;  % m/s
tflight = 4; % flight duration, s

% Derived parameters
bw = c/(2*rangeRes); % signal bandwidth
slowTime = 1/prf; % Cross-range time information.
                  % Defines the time instances at which the pulses are
                  % transmitted along the flight path.
numpulses = tflight/slowTime +1;
truncrangesamples = ceil((2*maxRange/c)*fs);
fastTime = (0:1/fs:(truncrangesamples-1)/fs); % time duration for the operation of each pulse.
lambda = c/fc;

% radar LFM signal source
plat     = phased.Platform('InitialPosition', [0;-200;500],...
                           'Velocity', [0; speed; 0]);
waveform = phased.LinearFMWaveform('SampleRate',fs,...
                                   'PulseWidth', tpd,...
                                   'PRF', prf,...
                                   'SweepBandwidth', bw);
                               
% SAR transmitter, antenna looks orthogonal to flight track
N = 4;
antenna = phased.CosineAntennaElement('FrequencyRange', [1e9 6e9]); % fc = 4e9
%antenna = phased.ULA('Element',element,'NumElements',N); % fc = 4e9
antennaGain = aperture2gain(aperture,c/fc); 

transmitter = phased.Transmitter('PeakPower', 50e3,...
                                 'Gain', antennaGain);
radiator = phased.Radiator('Sensor', antenna,...
                           'OperatingFrequency', fc,...
                           'PropagationSpeed', c);

% SAR receiver
collector = phased.Collector('Sensor', antenna,...
                             'PropagationSpeed', c,...
                             'OperatingFrequency', fc);
receiver = phased.ReceiverPreamp('SampleRate', fs,...
                                 'NoiseFigure', 30);

% propagation channel
channel = phased.FreeSpace('PropagationSpeed', c,...
                           'OperatingFrequency', fc,...
                           'SampleRate', fs,...
                           'TwoWayPropagation', true);
                       
%% Target Configuration

targetpos = [800,0,0; 1000,0,0; 1300,0,0]'; 
targetvel = [0,0,0;   0,0,0;    0,0,0]'; % stationary

target = phased.RadarTarget('OperatingFrequency',fc,...
                            'MeanRCS', [1,1,1]); % m^2
pointTargets = phased.Platform('InitialPosition',targetpos,...
                               'Velocity',targetvel);
                           
% The figure below describes the ground truth based on the target
% locations.
figure(1); h = axes;
plot(targetpos(2,1),targetpos(1,1),'*g'); grid on; hold on;
plot(targetpos(2,2),targetpos(1,2),'*r');
plot(targetpos(2,3),targetpos(1,3),'*b'); hold off;
set(h,'Ydir','reverse');
xlim([-10 10]); ylim([700 1500]);
title('Ground Truth Target Locations');
ylabel('Range'); xlabel('Cross-Range');
legend('Target 1','Target 2','Target 3','location','best');


%% SAR Simulation

% Define the broadside angle
refangle = zeros(1,size(targetpos,2));
rxsig = zeros(truncrangesamples,numpulses);
for ii = 1:numpulses
    % Update radar platform and target position
    [radarpos, radarvel] = plat(slowTime);
    [targetpos,targetvel] = pointTargets(slowTime);
    
    % Get the range and angle to the point targets
    [targetRange, targetAngle] = rangeangle(targetpos, radarpos);
    
    % Generate the LFM pulse
    sig = waveform();
    % Use only the pulse length that will cover the targets.
    sig = sig(1:truncrangesamples);
    
    % Transmit the pulse
    sig = transmitter(sig);
    
    % Define no tilting of beam in azimuth direction
    targetAngle(1,:) = refangle;
    
    % Radiate the pulse towards the targets
    sig = radiator(sig, targetAngle);
    
    % Propagate the pulse to the point targets in free space
    sig = channel(sig, radarpos, targetpos, radarvel, targetvel);
    
    % Reflect the pulse off the targets
    sig = target(sig);
    
    % Collect the reflected pulses at the antenna
    sig = collector(sig, targetAngle);
    
    % Receive the signal  
    sig = receiver(sig);
    
    % store
    rxsig(:,ii) = sig;
end

figure(2); imagesc(abs(rxsig)); title('SAR Raw Data, Stripmap Scan')
xlabel('Cross-Range Samples');
ylabel('Range Samples');

%% Post-Processing

% Perform range compression. Each row of the received signal, which
% contains all the information from each pulse, can be matched filtered
% to get the dechirped or range compressed signal.
[cdata, rnggrid] = ...
    doRangeCompression(rxsig, waveform.PulseWidth, waveform.SweepBandwidth, fs);

figure(3); imagesc(abs(cdata)); title('SAR Range Compressed Data')
xlabel('Cross-Range Samples')
ylabel('Range Samples')

% Perform azimuth Compression.
rma_processed = doAzCompression(cdata,fc,fs,prf,numpulses,Rf,maxRange,speed);

ranges = (0:1:truncrangesamples-1).*c/(2*fs);
az = ((0:1:numpulses-1) - (numpulses-1)/2)*lambda/2;

figure(4); %imagesc((abs((rma_processed(1700:2300,600:1400).'))));
imagesc(abs(rma_processed.'));
title('SAR Data focused using Range Migration algorithm ');
xlabel('Cross-Range Samples');
ylabel('Range Samples');
