% Exercise 7
% In automotive communications, a vehicle is located at the front view and using a virtual array
% of radar with 20 elements is in place, implementing a beamformer aimed at five degrees from the
% source vehicle movement direction corresponding to the target position. Assume the target vehicle
% changes its position by three degrees with the beamforming unchanged in the source vehicle. What
% will be the change in the range and Doppler shift estimation?
% Assume a bandwidth of BW = 500 MHz, a carrier of fc = 80 GHz, a pulse duration of
% T = 40µs, and a range of detection of d = 100 m. Assume the target is moving at v = 60 km/h.
% Consider that we need to cover a range swath of 10m, so that you can choose L and M = 4.

%Mateus Goncalves Guina
%DRE: 115061923
clc
clear
close all

% -------------------------------------------------------------------------
%                               Definitions
% -------------------------------------------------------------------------
T = 40*1e-6;  % Pulse duration PRI
PRF = 1/T;  % Pulse repetition frequency
dc = 1;  % Duty-cycle
T_ch = dc*T;  % Chirp duration
BW = 500*1e6;  % Signal bandwidth
f_c = 80*1e9;  % Carrier frequency
v = 60;  % Relative target speed, km/h
ac = 10;  % Aceleration m/s^2
d = 100;  % Distance in meters, target
c = 299792458;  % Speed of light, m/s
d_w = 10;  % Range swath
delta_d = c/(2*BW);  % Range resolution
L = round(d_w/delta_d); % Total number of samples in fast-time
kmh2ms = @(x) 1000*x/(60*60);  % Function to convert km/h to m/s
elements=20;
% Doppler effects
f_D = ((2*kmh2ms(v))/(c)) * f_c;  % Doppler effect
nu = f_D/f_c;  % Doppler ratio
f_c2 = 0;  % Reduce computational burden.
a = BW/T_ch;  % Chirp slope

% Sampling rate estimation:
% To be able to detect a velocity of 300 km/h, one must have a minimal
% sampling rate. The sampling period is given by the following equation.
Fs = 2*BW*(1.5);
M = 10;  % Number of pulses
tau = 2*d/c;  % How much time to m=0 peak
n0 = round(tau*Fs - 4*L/5);  % I want n0 as small as possible
fft_size = 2^10;  % High FFT size is a problem in range estimation




% -------------------------------------------------------------------------
%                    Beamforming
% -------------------------------------------------------------------------

angle_of_arrival = [5;0];
lambda = physconst('LightSpeed')/f_c;
array = phased.ULA('NumElements',elements,'ElementSpacing',lambda/2);
array.Element.FrequencyRange = [f_c-BW/2 f_c+BW/2];

 beamformer = phased.PhaseShiftBeamformer('SensorArray',array,...
    'OperatingFrequency',f_c,'Direction',angle_of_arrival,...
    'WeightsOutputPort',true);


% -------------------------------------------------------------------------
%                    Signal generation and processing
% -------------------------------------------------------------------------
t_ch = linspace(0, T_ch, T_ch*Fs);

% Transmitted signal (reference signal for xcorr)
tx = exp(1j*pi*(f_c2).*t_ch).*exp(1j*pi*a.*t_ch.^2);
tx = [tx zeros(1, round((T-T_ch)*Fs))];

% Received signal
rx = exp(1j*pi*(f_c2+f_D).*t_ch).*exp(1j*pi*a*(1+nu).*t_ch.^2);

% Memory allocation
A = zeros(2*L + 1, M);
N = L+n0;
snr = 150000;
for m = 1:M
    % Target
    d_m = d + (m-1)*kmh2ms(v)*(1+ac)*T;
    t_samp_m = round(2*(d_m/c)*Fs);
    rx_m = awgn([zeros(1, t_samp_m) rx zeros(1, round((T-T_ch)*Fs) - t_samp_m)], snr, 'measured');

    x = collectPlaneWave(array,rx_m',angle_of_arrival,f_c);%sinal com 5 graus
    x2 = collectPlaneWave(array,rx_m',[3;0],f_c);%sinal com 3 graus


    [y,w] = beamformer(x);
    [y2,w2] = beamformer(x2);
    A(:, m) = xcorr(y(n0+1:N,1)', tx(1, n0+1:N), L);%sinal recebido com 5 graus
    A2(:, m) = xcorr(y2(n0+1:N,1)', tx(1, n0+1:N), L); %sinal recebido com 3 graus 
end


% -------------------------------------------------------------------------
%                           Fast-time DFT
% -------------------------------------------------------------------------
ms2kmh = @(x) (60*60)*x/1000;  % Function to convert m/s to km/h
fprintf('Distance estimation, (DFT):\n')
freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
                 linspace(-fft_size/2, -1, fft_size/2)])*(Fs/fft_size);

% d_fft = (n0*c)/(2*Fs);  % My range estimation a priori

% f_tt = (d_fft*2*a)/c;
A = A(L+2:end, :);
A2 = A2(L+2:end, :);
for m = 1:M
%     d_fft_past = d_fft;
    AA = fftshift(fft((A(:, m)), fft_size));
    AA2 = fftshift(fft((A2(:, m)), fft_size));
    [value, f_idx] = max(abs(AA));
    [value2, f_idx2] = max(abs(AA2));
    f_t = freq(f_idx);
    f_t2 = freq(f_idx2);
    d_fft = (f_t*c)/(2*a);

    txt = {['Estimated distance: ' num2str(d_fft) ' m'],...
        ['$\hat{f}_t = ' num2str(f_t) '$ Hz']};
    txt2 = {['Estimated distance: ' num2str(d_fft) ' m'],...
        ['$\hat{f}_t = ' num2str(f_t) '$ Hz']};
    
    if (m == 1)||(m==M)
        figure
        plot(freq, 10*log10(abs(AA)), 'linewidth', 2), hold on,
        plot(f_t, 10*log10(value), 'o', 'linewidth', 2, 'markersize', 10)
        xline(f_t, '--', 'linewidth', 2), hold off, grid on
        text(-8*1e8, 22, txt, 'interpreter', 'latex', 'fontsize', 12)
        xlabel('Frequency, $f$ Hz', 'interpreter', 'latex')
        ylabel('Maginitude, dB', 'interpreter', 'latex')
        title(['Fast-time FFT, angulo 5' ], 'interpreter', 'latex')
        fprintf('valor %f m',d_fft)
        figure
        plot(freq, 10*log10(abs(AA2)), 'linewidth', 2), hold on,
        plot(f_t2, 10*log10(value2), 'o', 'linewidth', 2, 'markersize', 10)
        xline(f_t, '--', 'linewidth', 2), hold off, grid on
        text(-8*1e8, 20, txt2, 'interpreter', 'latex', 'fontsize', 12)
        xlabel('Frequency, $f$ Hz', 'interpreter', 'latex')
        ylabel('Maginitude, dB', 'interpreter', 'latex')
        title(['Fast-time FFT, angulo 3'], 'interpreter', 'latex')
    end
end

% -------------------------------------------------------------------------
%                          Plot Beamforming
% -------------------------------------------------------------------------
figure
pattern(array,f_c,[-90:90],0,'PropagationSpeed',physconst('LightSpeed'),'Type',...
    'powerdb','CoordinateSystem','polar','Weights',w)

azang = -90:5:90;
figure
subplot(211)
pattern(array,f_c,[-90:90],0,'CoordinateSystem','rectangular',...
    'Type','powerdb','PropagationSpeed',physconst('LightSpeed'))
set(gca,'xtick',azang);
title('Recepção sem Beamforming Weights')
subplot(212)
pattern(array,f_c,[-90:90],0,'CoordinateSystem','rectangular',...
    'Type','powerdb','PropagationSpeed',physconst('LightSpeed'),...
    'Weights',w)

set(gca,'xtick',azang);
title('Recepção com Beamforming Weights')

