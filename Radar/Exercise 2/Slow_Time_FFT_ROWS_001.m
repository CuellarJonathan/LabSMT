% Slow Time FFT
% 
% Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
% Jonathan S. Cuellar - j.silvacuellar@poli.ufrj.br
% Dec. 2021
%

clc
clear
close all

% Definitions
pulseDuration = 40*1e-6;  % Seconds
dutyCycle = .5;
carrierFrequency = 30*1e9;  % Hz
bandWidth = 500*1e6;  % Hz
kmph2mps = @(kmphSpeed) kmphSpeed/3.6;
mps2kmph = @(mpsSpeed) mpsSpeed*3.6;
radialVelocity = kmph2mps(200);  % m/s
% radialVelocity = kmph2mps(10*1e8/275000);  % m/s
% radialVelocity = kmph2mps(100);  % m/s
expectedRadialDistance = 2500;  % m
targetRadialDistance = 300;  % m
speedOfLight = 299792458;  % m/s
rangeSwath = 20;  % m
numberOfPulses = 128;
fastTimeFFTLength = 2^10;
% slowTimeFFTLength = numberOfPulses;
slowTimeFFTLength = 2^9;
timeFFTLength = 2^6;
% timeFFTLength = 2^25;
% signalToNoiseRatio = -13.5;  % dB
signalToNoiseRatio = 2000000;  % dB

% scaleError = 1e8/(2^10*2^8)
% scaleError = 1e8/(fastTimeFFTLength*numberOfPulses*4)
% radialVelocity = radialVelocity*scaleError; 

% Pre Calculations:
pulseRepetitionFrequency = 1/pulseDuration;
chirpDuration = dutyCycle*pulseDuration;
chirpSlope = bandWidth/chirpDuration;
rangeResolution = speedOfLight/(2*bandWidth);
samplesInFastTime = round(rangeSwath/rangeResolution);
% samplingRate = 2*bandWidth*1.5;
samplingRate = 9*1e8;
% samplingRate = 100*1e8;
dopplerEffect = carrierFrequency*2*radialVelocity/speedOfLight;
dopplerRatio = dopplerEffect/carrierFrequency;
expectedDelayTime = 2*expectedRadialDistance/speedOfLight;
expectedDelaySamples = round(expectedDelayTime*samplingRate);

chirpTimeAxis = linspace(0, chirpDuration, chirpDuration*samplingRate);
timeAxis = linspace(0, pulseDuration, pulseDuration*samplingRate);
% 0 (zero) here is because we are working in base band.
% transmittedChirp = exp(1j*2*pi*0*chirpTimeAxis).*exp(1j*pi*chirpSlope* ...
%     chirpTimeAxis.^2);

% receivedChirp = exp(-1j*pi*2*dopplerEffect*chirpTimeAxis).*...
%     exp(1j*pi*chirpSlope*(1-dopplerRatio)*chirpTimeAxis.^2);  % Eq. 9.23

% Received signal
receivedChirp = exp(1j*2*pi*(carrierFrequency+dopplerEffect).*chirpTimeAxis).*exp(1j*pi*chirpSlope*(1+dopplerRatio).*chirpTimeAxis.^2);

range2freq = @(distance) 2*chirpSlope*distance/speedOfLight;

% referenceSignal = transmittedChirp(expectedDelaySamples+1: ...
%     samplesInFastTime+expectedDelaySamples);

transmittedChirp = [exp(1j*2*pi*(carrierFrequency).*(chirpTimeAxis)).*...
        exp(1j*pi*chirpSlope*(1).*(chirpTimeAxis).^2) ...
        zeros(1, round((pulseDuration-chirpDuration)*samplingRate))];
    
transmittedPulse = [transmittedChirp zeros(1, round((pulseDuration- ...
    chirpDuration)*samplingRate))];
referenceSignal = transmittedChirp(expectedDelaySamples+1: ...
    samplesInFastTime+expectedDelaySamples);

% Generating received pulses and dechirp process (in fast-time x slow-time
% matrix)

fastSlowMatrix = zeros(fastTimeFFTLength/2, numberOfPulses);
radialDistanceCorrected = zeros(1, numberOfPulses);
for pulse = 1:numberOfPulses
    radialDistanceCorrected(pulse) = targetRadialDistance + ((pulse-1)* ...
        radialVelocity*pulseDuration);
    receiverDelaySample = round(2*samplingRate*radialDistanceCorrected(pulse)/ ...
        speedOfLight);
%     receivedPulse = [zeros(1, receiverDelaySample) ...
%         add_noise(receivedChirp, signalToNoiseRatio) ...
%         zeros(1, round((pulseDuration-chirpDuration)*samplingRate)- ...
%         receiverDelaySample)];
    
    tau = (2*radialDistanceCorrected(pulse))/speedOfLight;
%     receivedTimeAxis = linspace(0,chirpDuration+tau,(chirpDuration+tau)*samplingRate);
%     receivedPulse = [exp(1j*2*pi*(carrierFrequency+dopplerEffect).*(receivedTimeAxis-tau)).*...
%         exp(1j*pi*chirpSlope*(1+dopplerRatio).*(receivedTimeAxis-tau).^2) ...
%         zeros(1, round((pulseDuration-chirpDuration-tau)*samplingRate)-...
%         receiverDelaySample)];
    
    tauTimeAxis = linspace(chirpDuration,chirpDuration+tau,tau*samplingRate);
    
%     receivedPulseTau = exp(1j*2*pi*(carrierFrequency+dopplerEffect).*(tauTimeAxis-tau)).*...
%         exp(1j*pi*chirpSlope*(1+dopplerRatio).*(tauTimeAxis-tau).^2);
% 
%     receivedPulse = [exp(1j*2*pi*(carrierFrequency+dopplerEffect).*(chirpTimeAxis-tau)).*...
%         exp(1j*pi*chirpSlope*(1+dopplerRatio).*(chirpTimeAxis-tau).^2) ...
%         zeros(1, round((pulseDuration-chirpDuration)*samplingRate)-...
%         receiverDelaySample)];

    receivedPulseTau = exp(1j*2*pi*(carrierFrequency.*(tauTimeAxis-tau)-dopplerEffect.*dopplerRatio.*tauTimeAxis)).*...
        exp(1j*pi*chirpSlope.*(tauTimeAxis-tau-(dopplerRatio/2).*tauTimeAxis).^2);

    receivedPulse = [exp(1j*2*pi*(carrierFrequency.*(chirpTimeAxis-tau)-dopplerEffect.*dopplerRatio.*chirpTimeAxis)).*...
        exp(1j*pi*chirpSlope.*(chirpTimeAxis-tau-(dopplerRatio/2).*chirpTimeAxis).^2) ...
        zeros(1, round((pulseDuration-chirpDuration)*samplingRate)-...
        receiverDelaySample)];
    
    receivedPulse = [zeros(1, receiverDelaySample) ...
        receivedPulse(receiverDelaySample:(round(chirpDuration*samplingRate))) ...
        receivedPulseTau ...
        zeros(1, round((pulseDuration-chirpDuration-tau)*samplingRate)-...
        receiverDelaySample)];
    
    receivedPulse = awgn(receivedPulse,signalToNoiseRatio);
    
    if (pulse==numberOfPulses)
        figure
%         plot(timeAxis(1:length(referenceSignal)), real(referenceSignal)), grid on
        plot(timeAxis(1:length(transmittedChirp)), real(transmittedChirp)), grid on
        xlabel('Time [s]', 'interpreter', 'latex')
        ylabel('Amplitude', 'interpreter', 'latex')
        title('Reference Signal', 'interpreter', 'latex')
        
        figure
%         plot(timeAxis(1:length(receivedPulse(...
%         expectedDelaySamples+1:samplesInFastTime+expectedDelaySamples))), real(receivedPulse(...
%         expectedDelaySamples+1:samplesInFastTime+expectedDelaySamples))), grid on
        plot(timeAxis(1:length(receivedPulse)), real(receivedPulse)), grid on
        xlabel('Time [s]', 'interpreter', 'latex')
        ylabel('Amplitude', 'interpreter', 'latex')
        title('Received Signal', 'interpreter', 'latex')
        
        referenceSpectrum = fftshift(fft(transmittedChirp, timeFFTLength));
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
        
    
    dechirpedSignal = referenceSignal.*conj(receivedPulse(...
        expectedDelaySamples+1:samplesInFastTime+expectedDelaySamples));
%     dechirpedSignal = referenceSignal.*conj(receivedPulse(...
%         expectedDelaySamples+1:samplesInFastTime+expectedDelaySamples));
    dechirpedSpectrum = fftshift(fft(dechirpedSignal, fastTimeFFTLength));
    fastSlowMatrix(:, pulse) = dechirpedSpectrum(fastTimeFFTLength/2+1:end);
end

fastSlowMatrix2 = zeros(fastTimeFFTLength/2,slowTimeFFTLength);
for k = 1:1:length(fastSlowMatrix)
    fastSlowMatrix2(k,:) = fftshift(fft(fastSlowMatrix(k,:), slowTimeFFTLength));
end

% Noncoherent integration
integratedPulses = sqrt(sum(abs(fastSlowMatrix).^2, 2));
% fastSpectrum = fftshift(fft(integratedPulses, samplesInFastTime));
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

% Slow-time processing

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

% Plotting

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
