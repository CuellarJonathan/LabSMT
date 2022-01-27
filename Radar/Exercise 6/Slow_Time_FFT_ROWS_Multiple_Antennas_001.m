% Multiple Antennas - Fast-time x Slow-time FFT
%
% Jonathan S. Cuellar - j.silvacuellar@poli.ufrj.br
% Jan. 2022
%

clc
clear
close all

%% Definitions

Mt = 1; % Transmitter number
Mr = 5; % Receiver number
d_antennas_t = 0.1; % m % Transmitter antennas distance
d_antennas_r = 0.1; % m % Receiver antennas distance

kmph2mps = @(kmphSpeed) kmphSpeed/3.6;
mps2kmph = @(mpsSpeed) mpsSpeed*3.6;

targetRadialDistance_y = 100; % m % Vehicle towards distance
targetRadialDistance_x = 200; % m % Vehicle lateral distance
% targetRadialDistance_x = 57.735; % m % Vehicle lateral distance

angle_target_rad = atan2(targetRadialDistance_y,targetRadialDistance_x);
angle_target = rad2deg(angle_target_rad);

radialVelocity = kmph2mps(10); % km/h Vehicle velocity

radialVelocity_y = radialVelocity*sin(angle_target_rad); % m/s % Vehicle towards velocity
radialVelocity_x = radialVelocity*cos(angle_target_rad); % m/s % Vehicle lateral velocity

% angle_velocity_rad = atan2(radialVelocity_y,radialVelocity_x);

% signalV = radialVelocity_x/abs(radialVelocity_x);

% radialVelocity = signalV*hypot(radialVelocity_x,radialVelocity_y);

pulseDuration = 40*1e-6;  % Seconds
dutyCycle = .5;
carrierFrequency = 30*1e9;  % Hz
bandWidth = 500*1e6;  % Hz
% radialVelocity = kmph2mps(200);  % m/s
expectedRadialDistance = 1000;  % m
% targetRadialDistance = 300;  % m
speedOfLight = 299792458;  % m/s
rangeSwath = 20;  % m
numberOfPulses = 128;
fastTimeFFTLength = 2^10;
% slowTimeFFTLength = 2^9;
slowTimeFFTLength = numberOfPulses;
timeFFTLength = 2^6;
% timeFFTLength = 2^25;
% signalToNoiseRatio = -27;  % dB
signalToNoiseRatio = -13.5;  % dB
% signalToNoiseRatio = 2^20;  % dB

%% Pre Calculations:

pos_t = positionAntennas(Mt,d_antennas_t);
pos_r = positionAntennas(Mr,d_antennas_r);

rad2deg = @(x) x*180/pi;
deg2rad = @(x) x*pi/180;

% target_angle_rad = atan2(targetRadialDistance_x,targetRadialDistance_y); % Azimuth in radians
% target_angle_deg = rad2deg(target_angle_rad); % Azimuth in degrees

pulseRepetitionFrequency = 1/pulseDuration;
waveLength = speedOfLight/carrierFrequency;
chirpDuration = dutyCycle*pulseDuration;
chirpSlope = bandWidth/chirpDuration;
rangeResolution = speedOfLight/(2*bandWidth);
samplesInFastTime = round(rangeSwath/rangeResolution);
% samplingRate = 2*bandWidth*1.5;
samplingRate = 9*1e8;
% samplingRate = 1*1e11;
dopplerEffect = carrierFrequency*2*radialVelocity/speedOfLight;
dopplerRatio = dopplerEffect/carrierFrequency;
expectedDelayTime = 2*expectedRadialDistance/speedOfLight;
expectedDelaySamples = round(expectedDelayTime*samplingRate);
maximumRadialVelocity = waveLength/(4*pulseDuration*Mt);

assert(maximumRadialVelocity > abs(radialVelocity), ['The chosen ', ...
    'radial velocity is greater in magnitude than the maximum radial', ...
    'velocity. Decrease the radial velocity in magnitude or number of transmitter antennas.'])

if pulseRepetitionFrequency/2 < dopplerEffect
    warning(['The Doppler frequency is out of bounds. ', ...
        'Reduce pulse duration or radial velocity.'])
end

%% Transmitted Signal

timeAxis = linspace(0, pulseDuration, pulseDuration*samplingRate);
range2freq = @(distance) 2*chirpSlope*distance/speedOfLight;
signalPulse = exp(1j*2*pi*(carrierFrequency).*(timeAxis)).*...
        exp(1j*pi*chirpSlope.*(timeAxis).^2);
    
assert(length(signalPulse) == length(timeAxis), ['The len', ...
    'gth of the transmitted pulse is not correct.'])

transmittedChirp = signalPulse(1:round(chirpDuration*samplingRate));

transmittedPulse = [transmittedChirp zeros(1, round((pulseDuration- ...
    chirpDuration)*samplingRate))];
referenceSignal = transmittedPulse(expectedDelaySamples+1: ...
    samplesInFastTime+expectedDelaySamples);

assert(length(referenceSignal) == samplesInFastTime, ['Length of the ', ...
    'referece signal is not correct.'])

%% Generating received pulses and dechirp process 
% (in fast-time x slow-time matrix)

fastSlowMatrix = zeros(Mt,Mr,fastTimeFFTLength/2, numberOfPulses);

fprintf("Range Receiver:\n")

for pulse = 1:numberOfPulses
    for mt = 1:1:Mt
        fprintf("\nPulse %d:\n\n", pulse)
        for mr = 1:1:Mr
            radialDistanceCorrected_x = targetRadialDistance_x + (((pulse-1)*Mt+mt-1)* ...
                radialVelocity_x*pulseDuration);
            radialDistanceCorrected_t_y = targetRadialDistance_y - pos_t(mt) + (((pulse-1)*Mt+mt-1)* ...
                radialVelocity_y*pulseDuration);
            radialDistanceCorrected_r_y = targetRadialDistance_y - pos_r(mr) + (((pulse-1)*Mt+mt-1)* ...
                radialVelocity_y*pulseDuration);

            radialDistanceCorrected_t = hypot(radialDistanceCorrected_t_y,radialDistanceCorrected_x);
            radialDistanceCorrected_r = hypot(radialDistanceCorrected_r_y,radialDistanceCorrected_x);

            radialDistanceCorrected = radialDistanceCorrected_t+radialDistanceCorrected_r;
            
            fprintf("Range mt = %d e mr = %d: %.45f m.\n", mt, mr, radialDistanceCorrected/2)
                
            tau = radialDistanceCorrected/speedOfLight;
            receiverDelaySample = round(tau*samplingRate);

            tauTimeAxis = linspace(0,chirpDuration+tau,round((chirpDuration+tau)*samplingRate));

            % Radar Channel
            channelTau = exp(1j.*pi.*(-2.*carrierFrequency.*tau-2.*dopplerEffect.*tauTimeAxis+chirpSlope.*((tauTimeAxis.^2).*((1-dopplerRatio).^2-1)+(tau.^2)-2.*tau.*tauTimeAxis.*(1-dopplerRatio))));

            assert(length(channelTau) == length(tauTimeAxis), ['Incor', ...
                'rect length for channel pulse.'])

            % Received Signal
            receivedPulse = channelTau.*signalPulse(1:round((chirpDuration+tau)*samplingRate));
            receivedPulse = [zeros(1, receiverDelaySample) ...
                receivedPulse(receiverDelaySample:end) ...
                zeros(1, round((pulseDuration-chirpDuration-tau)*samplingRate)-...
                receiverDelaySample)];

            % Noise addition
            receivedPulse = awgn(receivedPulse,signalToNoiseRatio);

            dechirpedSignal = referenceSignal.*conj(receivedPulse(...
                expectedDelaySamples+1:samplesInFastTime+expectedDelaySamples));
        %     dechirpedSignal = referenceSignal.*conj(receivedPulse(1:round((expectedDelayTime+chirpDuration)*samplingRate)));
            dechirpedSpectrum = fftshift(fft(dechirpedSignal, fastTimeFFTLength));
            fastSlowMatrix(mt,mr,:, pulse) = dechirpedSpectrum(fastTimeFFTLength/2+1:end);
        end
    end
end

%% Fast-time processing

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

fprintf("\nEstimated Range Result: \n \n")

% lambda = speedOfLight/carrierFrequency;
% Vmax = lambda/(4*pulseDuration);
lambda = waveLength;
Vmax = maximumRadialVelocity;
deltaV = lambda/(2*pulseDuration*numberOfPulses);
freqMax = 2*Vmax*carrierFrequency/speedOfLight;
% velocityAxisX = linspace(-Vmax,Vmax,numberOfPulses);

fastSlowMatrix2 = zeros(Mt,Mr,fastTimeFFTLength/2,numberOfPulses);
fastSlowMatrix3 = zeros(Mt,Mr,fastTimeFFTLength/2,numberOfPulses);
fastSlowMatrix4 = zeros(Mt,fastTimeFFTLength/2,numberOfPulses);

for mt = 1:1:Mt
    for mr = 1:1:Mr
        fastSlowMatrixT = num2cell(fastSlowMatrix(mt,mr,:,:));
        fastSlowMatrixTemp = cell2mat(reshape(fastSlowMatrixT,fastTimeFFTLength/2,numberOfPulses));
        integratedPulses = sqrt(sum(abs(fastSlowMatrixTemp).^2, 2));
        [maxValueFastTime, indexFastTime] = max(integratedPulses);
        [maxValueFastTime2, indexFastTime2] = max(maxValueFastTime);

        fprintf("Range result mt = %d e mr = %d: %.45f m.\n", mt, mr, rangeAxis(indexFastTime(indexFastTime2)))
%         textDistance = {['Distance: ' num2str(rangeAxis(indexFastTime)) ' m.']};


%% Slow-time processing

        for k = 1:1:length(fastSlowMatrixTemp)
            fastSlowMatrix2(mt,mr,k,:) = fftshift(fft(fastSlowMatrixTemp(k,:), slowTimeFFTLength));
        end
        
        fastSlowMatrixT2 = num2cell(fastSlowMatrix2(mt,mr,:,:));
        fastSlowMatrixTemp2 = cell2mat(reshape(fastSlowMatrixT2,fastTimeFFTLength/2,numberOfPulses));

        dopplerAnalysis = fastSlowMatrixTemp(indexFastTime, :).';
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
        freq2velocity = @(freq) (freq*speedOfLight/(2*Mt*carrierFrequency)).*(1 - ...
            freq/(2*carrierFrequency));
        velocityAxis = mps2kmph(freq2velocity(slowFreqAxis));
        [maxSlowTime, indexSlowTime] = max(abs(slowTimeSpectrum));
        [maxSlowTime2, indexSlowTime2] = max(maxSlowTime);

        fprintf("Velocity result: %.4f km/h.\n", velocityAxis(indexSlowTime(indexSlowTime2)))
%         textVelocity = {['Velocity: ' num2str(velocityAxis(indexSlowTime)) 'km/h.']};

        % freq2velocity = @(freq) freq;
        % velocityAxis = freq2velocity(slowFreqAxis);

        % v = freq*c/(2*fc)
        % freq = 2*v*fc/c
        
   
    end
    
    %% Angle FFT processing
    
    dopplerAnalysis2 = fastSlowMatrix2(mt,:,indexSlowTime, indexSlowTime2).';
    slowTimeSpectrum2(mt,:,:) = periodogram(dopplerAnalysis2, [], ...
        Mr, 'centered');
    
    if rem(Mr, 2) == 0
        angleFreqAxis = fftshift([linspace(0, Mr/2-1, ...
            Mr/2) linspace(-Mr/2, -1, ...
            Mr/2)] * pulseRepetitionFrequency/Mr);
    else
        angleFreqAxis = fftshift([linspace(0, (Mr-1)/2, ...
            Mr/2+1) linspace(-(Mr-1)/2, -1, ...
            Mr/2)] * pulseRepetitionFrequency/Mr);
    end        
    
    for k = 1:1:fastTimeFFTLength/2
        for i = 1:1:slowTimeFFTLength
            fastSlowMatrix3(mt,:,k,i) = fftshift(fft(fastSlowMatrix2(mt,:,k,i), Mr));
%             fastSlowMatrix4(mt,k,i) = sqrt(sum(fastSlowMatrix3(mt,:,k,i).^2))/Mr;
            fastSlowMatrix4(mt,k,i) = (sum(fastSlowMatrix3(mt,:,k,i).^2));
        end
    end
    
%     for k = 1:1:fastTimeFFTLength/2
%         for i = 1:1:slowTimeFFTLength
%              matrixTemp = fftshift(fft(fastSlowMatrix2(mt,:,k,i), Mr));
%              [maxAngleR, indexAngleR] = max(abs(matrixTemp));
%              fastSlowMatrix3(k,i) = matrixTemp(indexAngleR);
%         end
%     end

%     fastSlowMatrix3(mt,:,indexSlowTime,indexSlowTime2) = fftshift(fft(fastSlowMatrix2(mt,:,indexSlowTime,indexSlowTime2), Mr));
%         
    fastSlowMatrixT3 = num2cell(fastSlowMatrix4(mt,:,:));
    fastSlowMatrixTemp3 = cell2mat(reshape(fastSlowMatrixT3,fastTimeFFTLength/2,numberOfPulses));

    phase2angle = @(phase) asin(phase*lambda/(2*pi*d_antennas_r));
    angleAxis = rad2deg(angle(phase2angle(angleFreqAxis)));
    [maxAngle, indexAngle] = max(angle(dopplerAnalysis2));
%     [maxAngle2, indexAngle2] = max(maxAngle);
%     [maxAngle3, indexAngle3] = max(maxAngle2);

%     fprintf("Angle result: %.2fº\n", angleAxis(round(maxAngle)))
    
end

[maxRange, indexRange] = max(abs(fastSlowMatrixTemp3));
[maxRange2, indexRange2] = max(maxRange);

fprintf('\n')
fprintf('Estimated Range Result: %.4f m \n', rangeAxis(indexRange(indexRange2)))
fprintf('Estimated Velocity Result:  %.4f km/h \n', velocityAxis(indexRange2))

%% Plotting

% figure
% mesh(0:pulseDuration:(numberOfPulses-1)*pulseDuration, rangeAxis, ...
%     abs(fastSlowMatrix), 'FaceColor', 'flat')
% view(2)
% xlabel('Time [$s$]', 'interpreter', 'latex')
% xlim([0 (numberOfPulses-1)*pulseDuration])
% ylabel('Range [$m$]', 'interpreter', 'latex')
% title('Fast-time vs. slow-time matrix', 'interpreter', 'latex')
% set(groot,'DefaultAxesTickLabelInterpreter','Tex');
% set(gca,'TickLabelInterpreter','latex')

figure
mesh(velocityAxis, rangeAxis, ...
    (abs(fastSlowMatrixTemp2)), 'FaceColor', 'flat')
view(2)
c = colorbar('TickLabelInterpreter','latex');
c.Label.String = 'Axis';
c.Label.Interpreter = 'latex';
xlabel('Velocity [$km/h$]', 'interpreter', 'latex')
xlim([min(velocityAxis) max(velocityAxis)])
ylabel('Range [$m$]', 'interpreter', 'latex')
ylim([0 expectedRadialDistance+1])
title('Fast-time vs. slow-time matrix', 'interpreter', 'latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(gca,'TickLabelInterpreter','latex')

figure
mesh(velocityAxis, rangeAxis, ...
    abs(fastSlowMatrixTemp3), 'FaceColor', 'flat')
view(2)
c = colorbar('TickLabelInterpreter','latex');
c.Label.String = 'Axis';
% c.Label.Interpreter = 'latex';
xlabel('Velocity [$km/h$]', 'interpreter', 'latex')
xlim([min(velocityAxis) max(velocityAxis)])
ylabel('Range [$m$]', 'interpreter', 'latex')
ylim([0 expectedRadialDistance+1])
title(['Fast-time vs. slow-time - Multiple Antennas (Mr = ' num2str(Mr) ')'], 'interpreter', 'latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(gca,'TickLabelInterpreter','latex')

% figure
% plot(rangeAxis, 10*log10(integratedPulses)), hold on
% plot(rangeAxis(indexFastTime), 10*log10(integratedPulses(indexFastTime) ...
%     ), 'x'), hold off, grid on
% xlabel('Range [$m$]', 'interpreter', 'latex')
% ylabel('Power [dB]', 'interpreter', 'latex')
% title(['Fast-time using noncoherent integration of $', ...
%     num2str(numberOfPulses), '$ pulses.'], 'interpreter', 'latex')
% set(groot,'DefaultAxesTickLabelInterpreter','Tex');
% set(gca,'TickLabelInterpreter','latex')

% figure
% plot(velocityAxis, 10*log10(slowTimeSpectrum)), grid on
% xlabel('Velocity [km/h]', 'interpreter', 'latex')
% ylabel('Power [dB]', 'interpreter', 'latex')
% title('Slow-time DFT', 'interpreter', 'latex')
% set(groot,'DefaultAxesTickLabelInterpreter','Tex');
% set(gca,'TickLabelInterpreter','latex')

% figure
% % plot(real(angleAxis), 10*log10(slowTimeSpectrum2(1,:,:))), grid on
% plot(real(angleAxis), (slowTimeSpectrum2(1,:,:))), grid on
% xlabel('Angle [Degrees]', 'interpreter', 'latex')
% ylabel('Power [dB]', 'interpreter', 'latex')
% title('Angle Slow-time DFT', 'interpreter', 'latex')
% set(groot,'DefaultAxesTickLabelInterpreter','Tex');
% set(gca,'TickLabelInterpreter','latex')

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

function d = distanceAntenna(pos,dx,dy,Mtr,M)
    d = zeros(Mtr,M);
    for n=1:Mtr
        d(n,1) = hypot(pos(n)-dx,dy); %sqrt((d_T(i-1,1)^2)+(d_antenas^2)-(2*d_T(i-1,1)*d_antenas*cos(angle)));  % Distance in meters, target
    end
end

function  pos = positionAntennas(Mtr,d_antenas)
    pos = zeros(1,Mtr);
    pos(1) = (Mtr-1)*d_antenas/2;
    for n = 2:Mtr
        pos(n) = pos(n-1)-d_antenas;
    end
end


% EoF
