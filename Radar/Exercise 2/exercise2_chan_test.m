% Exercise 2:
% 
% In automotive communications, the carrier frequency range lies between 76
% and 81 GHz. Assume a bandwidth of BW = 500 MHz, a carrier of f_c = 80
% GHz, a pulse duration of T = 40 \mus, and a range of detection d = 300 m.
% Assume the target is moving at v = 300 km/h. Consider that we need to
% cover a range swath of 20 m, so that you can choose L.
% 
% In this case, the Doppler effect originates a frequency shift of f_D =
% f_c*(2*v/c) = 44444.44 Hz.
% 
% Use a pulse burst FMCW with one and four pulses and estimate the range
% and the Doppler shift in both cases. Use a DFT in the fast-time axis to
% estimate the range, and then calculate the Doppler shift from the data
% delay for M = 1 and by exploring the slow-time information for M = 4.
%
% You can automatically detect a single target's delay by comparing the
% outcome of the matched filter with some prescribed threshold. Compare the
% results with those you know precisely the reflection delay, with the
% estimated delay, and with the exact delay +/- 5% of its nominal value.
%
%
% Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
% Jonathan S. Cuellar - j.silvacuellar@poli.ufrj.br
% Jan. 2022
%

clc
clear
close all


% Definitions
pulseDuration = 40*1e-6;  % Seconds
dutyCycle = .5;
carrierFrequency = 80*1e9;  % Hz
bandWidth = 500*1e6;  % Hz
kmph2mps = @(kmphSpeed) kmphSpeed/3.6;
mps2kmph = @(mpsSpeed) mpsSpeed*3.6;
radialVelocity = kmph2mps(60);  % m/s
expectedRadialDistance = 300;  % m
targetRadialDistance = 300;  % m
speedOfLight = 299792458;  % m/s
rangeSwath = 20;  % m
numberOfPulses = 128;
fastTimeFFTLength = 2^10;
slowTimeFFTLength = 2^10;
signalToNoiseRatio = -13.5;  % dB


% Pre Calculations:
pulseRepetitionFrequency = 1/pulseDuration;
waveLength = speedOfLight/carrierFrequency;
chirpDuration = dutyCycle*pulseDuration;
chirpSlope = bandWidth/chirpDuration;
rangeResolution = speedOfLight/(2*bandWidth);
samplesInFastTime = round(rangeSwath/rangeResolution);
samplingRate = 2*bandWidth*1.5;
dopplerEffect = carrierFrequency*2*radialVelocity/speedOfLight;
dopplerRatio = dopplerEffect/carrierFrequency;
expectedDelayTime = 2*expectedRadialDistance/speedOfLight;
expectedDelaySamples = round(expectedDelayTime*samplingRate);
maximumRadialVelocity = waveLength/chirpDuration;

assert(maximumRadialVelocity > abs(radialVelocity), ['The chosen ', ...
    'radial velocity is greater in magnitude than the maximum radial', ...
    'velocity. Decrease the radial velocity in magnitude.'])

if pulseRepetitionFrequency/2 < dopplerEffect
    warning(['The Doppler frequency is out of bounds. ', ...
        'Reduce pulse duration or radial velocity.'])
end

% Generating transmitted pulse:
chirpTimeAxis = linspace(0, chirpDuration, chirpDuration*samplingRate);
pulseTimeAxis = linspace(0, pulseDuration, pulseDuration*samplingRate);

transmittedChirp = exp(1j*2*pi*(carrierFrequency).*chirpTimeAxis).*...
    exp(1j*pi*chirpSlope*chirpTimeAxis.^2);

assert(length(transmittedChirp) == length(chirpTimeAxis), ['The len', ...
    'gth of the transmitted chirp is not correct.'])
    
transmittedPulse = [transmittedChirp zeros(1, round((pulseDuration- ...
    chirpDuration)*samplingRate))];

assert(length(transmittedPulse) == length(pulseTimeAxis), ['The len', ...
    'gth of the transmitted pulse is not correct.'])

referenceSignal = transmittedPulse(1:samplesInFastTime);

assert(length(referenceSignal) == samplesInFastTime, ['Length of the ', ...
    'referece signal is not correct.'])


% Generating received pulses and dechirp process (in fast-time x slow-time
% matrix)
fastSlowMatrix = zeros(fastTimeFFTLength/2, numberOfPulses);
radialDistanceCorrected = zeros(1, numberOfPulses);

for pulse = 1:numberOfPulses
    radialDistanceCorrected(pulse) = targetRadialDistance + ((pulse-1)* ...
        radialVelocity*pulseDuration);
    receivedDelayTime = 2*radialDistanceCorrected(pulse)/speedOfLight;
    receivedDelaySamples = round(samplingRate*receivedDelayTime);
    
    channel_m = exp(-1j*2*pi*(carrierFrequency*receivedDelayTime+ ...
        (dopplerEffect*chirpTimeAxis))).*exp(1j*pi*chirpSlope*(( ...
        receivedDelayTime+dopplerRatio*chirpTimeAxis).^2) - ...
        2*receivedDelayTime*chirpTimeAxis - 2*dopplerRatio*chirpTimeAxis.^2);
    
    receivedChirp = channel_m.*transmittedChirp;
    
%     receivedChirp = exp(1j*2*pi*(carrierFrequency*(chirpTimeAxis- ...
%         receivedDelayTime)-dopplerEffect*chirpTimeAxis)).*exp(1j* ...
%         chirpSlope*pi*(chirpTimeAxis-receivedDelayTime - dopplerRatio* ...
%         chirpTimeAxis).^2);
    
    assert(length(receivedChirp) == length(chirpTimeAxis), ['Incor', ...
        'rect length for received chirp.'])

    receivedPulse = add_noise([zeros(1, receivedDelaySamples) ...
        receivedChirp ...
        zeros(1, round((pulseDuration-chirpDuration)*samplingRate) - ...
        receivedDelaySamples)], signalToNoiseRatio);
    
    assert(length(receivedPulse) == length(pulseTimeAxis), ['Incor', ...
        'rect length for received pulse.'])
            
    dechirpedSignal = referenceSignal.*conj(receivedPulse(...
        expectedDelaySamples+1:samplesInFastTime+expectedDelaySamples));
    dechirpedSpectrum = fftshift(fft(dechirpedSignal, fastTimeFFTLength));
    fastSlowMatrix(:, pulse) = dechirpedSpectrum(fastTimeFFTLength/2+1:end);
end


% Noncoherent integration
integratedPulses = sqrt(sum(abs(fastSlowMatrix).^2, 2));
fastFreqAxis = fftshift([linspace(0, fastTimeFFTLength/2-1, ...
    fastTimeFFTLength/2) linspace(-fastTimeFFTLength/2, -1, ...
    fastTimeFFTLength/2)] * samplingRate/fastTimeFFTLength);
fastFreqAxis = fastFreqAxis(fastTimeFFTLength/2+1:end);

freq2range = @(freq) freq*speedOfLight/(2*chirpSlope);
rangeAxis = freq2range(fastFreqAxis);
[~, indexFastTime] = max(integratedPulses);

fprintf("Range result: %.4f m.\n", rangeAxis(indexFastTime))
textDistance = {['Distance: ' num2str(rangeAxis(indexFastTime)) ' m.']};

% Slow-time processing
dopplerAnalysis = fastSlowMatrix(indexFastTime, :).';
slowTimeSpectrum = periodogram(dopplerAnalysis, [], ...
    slowTimeFFTLength, 'centered');
slowFreqAxis = fftshift([linspace(0, slowTimeFFTLength/2-1, ...
    slowTimeFFTLength/2) linspace(-slowTimeFFTLength/2, -1, ...
    slowTimeFFTLength/2)] * pulseRepetitionFrequency/slowTimeFFTLength);
freq2velocity = @(freq) (freq*speedOfLight/(2*carrierFrequency)).*(1 - ...
    freq/(2*carrierFrequency));
velocityAxis = mps2kmph(freq2velocity(slowFreqAxis));
[~, indexSlowTime] = max(abs(slowTimeSpectrum));

fprintf("Velocity result: %.4f km/h.\n", velocityAxis(indexSlowTime))
textVelocity = {['Velocity: ' num2str(velocityAxis(indexSlowTime)) 'km/h.']};


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
plot(rangeAxis, 10*log10(integratedPulses)), hold on
plot(rangeAxis(indexFastTime), 10*log10(integratedPulses(indexFastTime) ...
    ), 'x', 'linewidth', 2, 'markersize', 10), hold off, grid on
text(rangeAxis(indexFastTime)*2, ...
    10*log10(integratedPulses(indexFastTime)), textDistance, ...
    'interpreter', 'latex', 'fontsize', 12)
xlabel('Range [$m$]', 'interpreter', 'latex')
ylabel('Power [dB]', 'interpreter', 'latex')
title(['Fast-time using noncoherent integration of $', ...
    num2str(numberOfPulses), '$ pulses.'], 'interpreter', 'latex')

figure
plot(velocityAxis, 10*log10(slowTimeSpectrum)), hold on
plot(velocityAxis(indexSlowTime), ...
    10*log10(slowTimeSpectrum(indexSlowTime)), 'x', 'linewidth', 2, ...
    'markersize', 10), hold off, grid on
text(-20, 10*log10(slowTimeSpectrum(indexSlowTime))-2.5, textVelocity, ...
    'interpreter', 'latex', 'fontsize', 12)
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
