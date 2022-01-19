% Fast-time x Slow-time FFT
% 
% Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
% Jonathan S. Cuellar - j.silvacuellar@poli.ufrj.br
% Jan. 2022
%

clc
clear
close all

%% Definitions
pulseDuration = 40*1e-6;  % Seconds
dutyCycle = .5;
carrierFrequency = 30*1e9;  % Hz
bandWidth = 500*1e6;  % Hz
kmph2mps = @(kmphSpeed) kmphSpeed/3.6;
mps2kmph = @(mpsSpeed) mpsSpeed*3.6;
radialVelocity = kmph2mps(200);  % m/s
expectedRadialDistance = 500;  % m
targetRadialDistance = 300;  % m
speedOfLight = 299792458;  % m/s
rangeSwath = 20;  % m
numberOfPulses = 128;
fastTimeFFTLength = 2^10;
slowTimeFFTLength = 2^9;
timeFFTLength = 2^6;
% signalToNoiseRatio = -13.5;  % dB
signalToNoiseRatio = 2^20;  % dB

%% Pre Calculations:
pulseRepetitionFrequency = 1/pulseDuration;
chirpDuration = dutyCycle*pulseDuration;
chirpSlope = bandWidth/chirpDuration;
rangeResolution = speedOfLight/(2*bandWidth);
samplesInFastTime = round(rangeSwath/rangeResolution);
% samplingRate = 2*bandWidth*1.5;
% samplingRate = 9*1e11;
samplingRate = 9*1e8;
dopplerEffect = carrierFrequency*2*radialVelocity/speedOfLight;
dopplerRatio = dopplerEffect/carrierFrequency;
expectedDelayTime = 2*expectedRadialDistance/speedOfLight;
expectedDelaySamples = round(expectedDelayTime*samplingRate);

%% Transmitted Signal
timeAxis = linspace(0, pulseDuration, pulseDuration*samplingRate);

range2freq = @(distance) 2*chirpSlope*distance/speedOfLight;

signalPulse = exp(1j*2*pi*(carrierFrequency).*(timeAxis)).*...
        exp(1j*pi*chirpSlope.*(timeAxis).^2);

transmittedChirp = signalPulse(1:round(chirpDuration*samplingRate));

transmittedPulse = [transmittedChirp zeros(1, round((pulseDuration- ...
    chirpDuration)*samplingRate))];

referenceSignal = transmittedPulse(expectedDelaySamples+1: ...
    samplesInFastTime+expectedDelaySamples);

%% Generating received pulses and dechirp process 
% (in fast-time x slow-time matrix)

fastSlowMatrix = zeros(fastTimeFFTLength/2, numberOfPulses);
radialDistanceCorrected = zeros(1, numberOfPulses);
for pulse = 1:numberOfPulses
    radialDistanceCorrected(pulse) = targetRadialDistance + ((pulse-1)* ...
        radialVelocity*pulseDuration);
    receiverDelaySample = round(2*samplingRate*radialDistanceCorrected(pulse)/ ...
        speedOfLight);
    
    tau = (2*radialDistanceCorrected(pulse))/speedOfLight;
    
    tauTimeAxis = linspace(0,chirpDuration+tau,round((chirpDuration+tau)*samplingRate));
    
    % Radar Channel
    channelTau = exp(1j.*pi.*(-2.*carrierFrequency.*tau-2.*dopplerEffect.*tauTimeAxis+chirpSlope.*((tauTimeAxis.^2).*((1-dopplerRatio).^2-1)+(tau.^2)-2.*tau.*tauTimeAxis.*(1-dopplerRatio))));

    % Received Signal
    receivedPulse = channelTau.*signalPulse(1:round((chirpDuration+tau)*samplingRate));

    receivedPulse = [zeros(1, receiverDelaySample) ...
        receivedPulse(receiverDelaySample:end) ...
        zeros(1, round((pulseDuration-chirpDuration-tau)*samplingRate)-...
        receiverDelaySample)];
    
    receivedPulse = awgn(receivedPulse,signalToNoiseRatio);
    
    % Plots begin
    if (pulse==numberOfPulses)
        figure
%         plot(timeAxis(1:length(transmittedPulse)), real(transmittedPulse)), grid on
        plot(timeAxis, real(transmittedPulse(1:length(timeAxis)))), grid on
        xlabel('Time [s]', 'interpreter', 'latex')
        ylabel('Amplitude', 'interpreter', 'latex')
        title('Reference Signal', 'interpreter', 'latex')
        
        figure
        plot(timeAxis(1:length(receivedPulse)), real(receivedPulse)), grid on
        xlabel('Time [s]', 'interpreter', 'latex')
        ylabel('Amplitude', 'interpreter', 'latex')
        title('Received Signal', 'interpreter', 'latex')
        
        referenceSpectrum = fftshift(fft(transmittedPulse, timeFFTLength));
        referenceFFT = referenceSpectrum(timeFFTLength/2+1:end);
        receivedSpectrum = fftshift(fft(receivedPulse, timeFFTLength));
        receivedFFT = receivedSpectrum(timeFFTLength/2+1:end);
        
        if rem(timeFFTLength, 2) == 0
            freqAxis = fftshift([linspace(0, timeFFTLength/2-1, ...
                timeFFTLength/2) linspace(-timeFFTLength/2, -1, ...
                timeFFTLength/2)] * samplingRate/timeFFTLength);
            freqAxis = freqAxis(timeFFTLength/2+1:end);
        else
            freqAxis = fftshift([linspace(0, (timeFFTLength-1)/2, ...
                timeFFTLength/2+1) linspace(-(timeFFTLength-1)/2, -1, ...
                timeFFTLength/2)] * samplingRate/timeFFTLength);
        end
        
        figure
        plot(freqAxis, abs(referenceFFT)), grid on
        xlabel('Freq [Hz]', 'interpreter', 'latex')
%         xlim([300*1e8 (300+5)*1e8])
        ylabel('Amplitude', 'interpreter', 'latex')
        title('Reference Spectrum', 'interpreter', 'latex')
        
        figure
        plot(freqAxis, abs(receivedFFT)), grid on
        xlabel('Freq [Hz]', 'interpreter', 'latex')
%         xlim([296.6*1e8 (296.6+5)*1e8])
        ylabel('Amplitude', 'interpreter', 'latex')
        title('Received Spectrum', 'interpreter', 'latex')
    end
    % Plots end
        
    dechirpedSignal = referenceSignal.*conj(receivedPulse(...
        expectedDelaySamples+1:samplesInFastTime+expectedDelaySamples));
%     dechirpedSignal = referenceSignal.*conj(receivedPulse(1:round((expectedDelayTime+chirpDuration)*samplingRate)));
    dechirpedSpectrum = fftshift(fft(dechirpedSignal, fastTimeFFTLength));
    fastSlowMatrix(:, pulse) = dechirpedSpectrum(fastTimeFFTLength/2+1:end);
end

%% Fast-time processing

integratedPulses = sqrt(sum(abs(fastSlowMatrix).^2, 2));
if rem(fastTimeFFTLength, 2) == 0
    fastFreqAxis = fftshift([linspace(0, fastTimeFFTLength/2-1, ...
        fastTimeFFTLength/2) linspace(-fastTimeFFTLength/2, -1, ...
        fastTimeFFTLength/2)] * samplingRate/fastTimeFFTLength);
    fastFreqAxis = fastFreqAxis(fastTimeFFTLength/2+1:end);
else
    fastFreqAxis = fftshift([linspace(0, (fastTimeFFTLength-1)/2, ...
        fastTimeFFTLength/2+1) linspace(-(fastTimeFFTLength-1)/2, -1, ...
        fastTimeFFTLength/2)] * samplingRate/fastTimeFFTLength);
end
freq2range = @(freq) freq*speedOfLight/(2*chirpSlope);
rangeAxis = freq2range(fastFreqAxis);
[maxValueFastTime, indexFastTime] = max(integratedPulses);

%% Slow-time processing

fastSlowMatrix2 = zeros(fastTimeFFTLength/2,slowTimeFFTLength);
for k = 1:1:length(fastSlowMatrix)
    fastSlowMatrix2(k,:) = fftshift(fft(fastSlowMatrix(k,:), slowTimeFFTLength));
end

lambda = speedOfLight/carrierFrequency;
Vmax = lambda/(4*pulseDuration);
% Vmax = lambda/(4*chirpDuration);
deltaV = lambda/(2*pulseDuration*numberOfPulses);
velocityAxisX = linspace(-Vmax,Vmax,numberOfPulses);

dopplerAnalysis = fastSlowMatrix(indexFastTime, :).';
slowTimeSpectrum = periodogram(dopplerAnalysis, [], ...
    slowTimeFFTLength, 'centered');
if rem(slowTimeFFTLength, 2) == 0
    slowFreqAxis = fftshift([linspace(0, slowTimeFFTLength/2-1, ...
        slowTimeFFTLength/2) linspace(-slowTimeFFTLength/2, -1, ...
        slowTimeFFTLength/2)] * pulseRepetitionFrequency/slowTimeFFTLength);
else
    slowFreqAxis = fftshift([linspace(0, (slowTimeFFTLength-1)/2, ...
        slowTimeFFTLength/2+1) linspace(-(slowTimeFFTLength-1)/2, -1, ...
        slowTimeFFTLength/2)] * pulseRepetitionFrequency/slowTimeFFTLength);
end
freq2velocity = @(freq) (freq*speedOfLight/(2*carrierFrequency)).*(1 - ...
    freq/(2*carrierFrequency));
velocityAxis = mps2kmph(freq2velocity(slowFreqAxis));

%% Plotting

figure
mesh(0:pulseDuration:(numberOfPulses-1)*pulseDuration, rangeAxis, ...
    abs(fastSlowMatrix), 'FaceColor', 'flat')
view(2)
xlabel('Time [$s$]', 'interpreter', 'latex')
xlim([0 (numberOfPulses-1)*pulseDuration])
ylabel('Range [$m$]', 'interpreter', 'latex')
title('Fast-time vs. slow-time matrix', 'interpreter', 'latex')

figure
mesh(velocityAxis, rangeAxis, ...
    abs(fastSlowMatrix2), 'FaceColor', 'flat')
view(2)
xlabel('Velocity [$km/h$]', 'interpreter', 'latex')
xlim([min(velocityAxis) max(velocityAxis)])
ylabel('Range [$m$]', 'interpreter', 'latex')
ylim([0 expectedRadialDistance+1])
title('Fast-time vs. slow-time matrix', 'interpreter', 'latex')

figure
plot(rangeAxis, 10*log10(integratedPulses)), hold on
plot(rangeAxis(indexFastTime), 10*log10(integratedPulses(indexFastTime) ...
    ), 'x'), hold off, grid on
xlabel('Range [$m$]', 'interpreter', 'latex')
ylabel('Power [dB]', 'interpreter', 'latex')
title(['Fast-time using noncoherent integration of $', ...
    num2str(numberOfPulses), '$ pulses.'], 'interpreter', 'latex')

figure
plot(velocityAxis, 10*log10(slowTimeSpectrum)), grid on
xlabel('Velocity [km/h]', 'interpreter', 'latex')
ylabel('Power [dB]', 'interpreter', 'latex')
title('Slow-time DFT', 'interpreter', 'latex')

%% Functions

function noisySignal = add_noise(inputSignal, signalToNoiseRatio)
    % Adds white Gaussian noise to the signal, adjuted by the SNR (Signal
    % to Noise Ratio).
    signalPower = mean(abs(inputSignal).^2);
    misadjustedNoise = randn(size(inputSignal));
    adjustedNoise = sqrt(signalPower*10^(-.1*signalToNoiseRatio) / ...
        mean(abs(misadjustedNoise).^2))*misadjustedNoise;
    noisySignal = inputSignal + adjustedNoise;
end


% EoF
