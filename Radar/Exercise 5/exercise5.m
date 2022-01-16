% Exercise 5 - Phase-migration & Array-alignment
%
% Luiz Felipe da Silveira Coelho - luizfelipe.coelho@smt.ufrj.br
% Nov 24, 2021
%


clc
clear
close all


% DEFINITIONS

% Signal
carrierFrequency = 24*1e9;
speedOfLight = 299792458;
waveLength = carrierFrequency/speedOfLight;
bandWidth = 500*1e8;
pulseDuration = 40*1e-6;
pulseRepetitionFrequency = 1/pulseDuration;
chirpDuration = .5*pulseDuration;
chirpSlope = bandWidth/chirpDuration;
samplingRate = 2*bandWidth;
rangeSwath = 10;
rangeResolution = speedOfLight/(2*bandWidth);
fastTimeLength = round(rangeSwath/rangeResolution);

% Target
xAxisDistance = 10;
yAxisDistance = 0;
radialDistance = sqrt(xAxisDistance^2 + yAxisDistance^2);
maximumRadialVelocity = waveLength/(4*pulseDuration);
targetRadialVelocity1 = -.5*maximumRadialVelocity;
targetRadialVelocity2 = -1*maximumRadialVelocity;


% SIMULATION

% Definitions
azimuthAxis = linspace(-1, 1, 500);

% Transmitted signal
expectedDelay = 2*radialDistance/speedOfLight;
delayInSamples = round(expectedDelay*samplingRate);
chirpTimeAxis = linspace(0, chirpDuration, chirpDuration*samplingRate);
transmittedSignal = exp(1j*pi*0.*chirpTimeAxis).*exp(1j*pi*chirpSlope.* ...
    chirpTimeAxis.^2);

% First scenario
numberOfAntennasInTx = 2;
numberOfAntennasInRx = 7;
distanceBetweenTxAntennas = waveLength/5;  % d<<wavelength
distanceBetweenRxAntennas = numberOfAntennasInTx*distanceBetweenTxAntennas;
firstScenarioSpecs = struct('numberOfTxAntennas', numberOfAntennasInTx, ...
    'numberOfRxAntennas', numberOfAntennasInRx, 'TxAntennaDistance', ...
    distanceBetweenTxAntennas, 'RxAntennaDistance', ...
    distanceBetweenRxAntennas);

% Second scenario

% Third scenario

% Fourth scenario


% Run - Scenario generation
channelFirstScenario1 = zeros(numberOfAntennasInTx, length(azimuthAxis));
channelFirstScenario2 = zeros(numberOfAntennasInTx, length(azimuthAxis));

for azimuthIndex = 1:length(azimuthAxis)
    % Velocity 1
    channelFirstScenario1(:, azimuthIndex) = generate_channel_matrix( ...
        radialDistance, targetRadialVelocity1, maximumRadialVelocity, ...
        waveLength, pi*azimuthAxis(azimuthIndex), firstScenarioSpecs) * ...
        ones(numberOfAntennasInRx, 1);
    % Velocity 2
    channelFirstScenario2(:, azimuthIndex) = generate_channel_matrix( ...
        radialDistance, targetRadialVelocity2, maximumRadialVelocity, ...
        waveLength, pi*azimuthAxis(azimuthIndex), firstScenarioSpecs) * ...
        ones(numberOfAntennasInRx, 1);
end


% Run - Signal through channel
% Scenario 1
dechirpedFirstScenario1 = zeros(fastTimeLength, length(azimuthAxis), ...
    numberOfAntennasInTx);
dechirpedFirstScenario2 = zeros(fastTimeLength, length(azimuthAxis), ...
    numberOfAntennasInTx);
for txIndex = 1:numberOfAntennasInTx
    % Velocity 1
    receivedFirstScenario1 = transmittedSignal(1:fastTimeLength)' * ...
        channelFirstScenario1(txIndex, :);
    % Velocity2
    receivedFirstScenario2 = transmittedSignal(1:fastTimeLength)' * ...
        channelFirstScenario2(txIndex, :);

    % Run - Dechirping
    dechirpMatrix = diag( ...
        transmittedSignal(delayInSamples+1:delayInSamples+fastTimeLength));
    % Velocity 1
    dechirpedFirstScenario1(:, :, txIndex) = dechirpMatrix*conj(...
        receivedFirstScenario1);
    % Velocity 2
    dechirpedFirstScenario2(:, :, txIndex) = dechirpMatrix*conj(...
        receivedFirstScenario2);
end

% sampleChoice = 1;
% dechirpedFirstScenario1 = [dechirpedFirstScenario1(:, :, 1); ...
%     dechirpedFirstScenario1(:, :, 2)];
% dechirpedFirstScenario2 = [dechirpedFirstScenario2(:, :, 1); ...
%     dechirpedFirstScenario2(:, :, 2)];

% Test
dechirpedFirstScenario1 = .5*dechirpedFirstScenario1(:, :, 1) + ...
    .5*dechirpedFirstScenario1(:, :, 2);
dechirpedFirstScenario2 = .5*dechirpedFirstScenario2(:, :, 1) + ...
    .5*dechirpedFirstScenario2(:, :, 2);


% FFT
fftSize = 2^10;
% freqAxis = fftshift([linspace(0, fftSize/2-1, fftSize/2) ...
%     linspace(-fftSize/2, -1, fftSize/2)]*pulseRepetitionFrequency/fftSize);
freqAxis = fftshift([linspace(0, fftSize/2-1, fftSize/2) ...
    linspace(-fftSize/2, -1, fftSize/2)]*samplingRate/fftSize);
rangeAxis = freqAxis*speedOfLight/(2*chirpSlope);

rangeAzimuthFirstScenario1 = abs(fftshift(fft(dechirpedFirstScenario1, ...
    fftSize, 1)));
rangeAzimuthFirstScenario2 = abs(fftshift(fft(dechirpedFirstScenario2, ...
    fftSize, 1)));


% PLOTTING
figure,
mesh(azimuthAxis, rangeAxis, rangeAzimuthFirstScenario1, 'FaceColor', ...
    'flat')
ylim([0 115])
view(2)

figure,
mesh(azimuthAxis, rangeAxis, rangeAzimuthFirstScenario2, 'FaceColor', ...
    'flat')
ylim([0 115])
view(2)


% FUNCTIONS
function channelMatrix = generate_channel_matrix(targetDistance, ...
    targetSpeed, maximumTargetSpeed, waveLength, azimuthAngle, SystemSpecs)

transmitterNo = SystemSpecs.numberOfTxAntennas;
transmitterDistance = SystemSpecs.TxAntennaDistance;
receiverNo = SystemSpecs.numberOfRxAntennas;
receiverDistance = SystemSpecs.RxAntennaDistance;

channelMatrix = zeros(transmitterNo, receiverNo);

for transmitterIndex = 1:transmitterNo
    targetMotion = 4*targetDistance*pi/waveLength + ...
        (transmitterIndex-1)*pi*targetSpeed/maximumTargetSpeed;
    for receiverIndex = 1:receiverNo
        phaseSum = 2*pi*sin(azimuthAngle) * ...
            ((transmitterIndex-1)*transmitterDistance + ...
            (receiverIndex-1)*receiverDistance)/waveLength;
        channelMatrix(transmitterIndex, receiverIndex) = exp(1j*(...
            phaseSum + targetMotion));
    end
end

end



% EoF

