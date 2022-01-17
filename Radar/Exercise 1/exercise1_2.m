% Show the real part of a radar chirp waveform with aT = 20, T = 2Tch and
% M = 5 for illustration.

clc
clear
close all

% Definitions
T = 1;  % Time duration
a = 20;  % Chrip rate
d_c = .5;  % Duty-cycle
f_c = 0;  % Carrier frequency
Fs = 1e3;  % Sampling frequency
M = 5;  % Number of pulses

t = linspace(0, T, d_c*T*Fs);
chirp = exp(1j*pi*f_c.*t).*exp(1j*(pi*a.*t.^2));
chirp_pulse = [chirp zeros(1, (1-d_c)*T*Fs)];

t_M = linspace(0, M*T, M*T*Fs);
x = zeros(1, M*T*Fs);
for m = 1:5
    x((m-1)*T*Fs+1:m*T*Fs) = chirp_pulse;
end

figure(1),
plot(t_M, real(x)), grid on
xlabel('Time, $t$ [sec]', 'interpreter', 'latex')
ylabel('Amplitude, $x(t)$', 'interpreter', 'latex')
title('Radar Chirp for $M=5$', 'interpreter', 'latex')