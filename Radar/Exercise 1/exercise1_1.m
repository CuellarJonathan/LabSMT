% • For a chirp supported for a normalized period of time Tch = 1 = T plot
% the real part of the radar chirp for aT = 20.
% • Plot the ambiguity function.
% • Plot the spectrum of the radar signal.
% • Repeat the above plots for sinusoid pulse and for a square wave pulse.

clc
clear
close all

% Global Definitions
T = 1;  % Pulse duration
Fs = 1e3;  % Sampling rate
t = linspace(0, T, T*Fs);  % Time axis

% -------------------------------------------------------------------------
%                               Chirp Signal
% -------------------------------------------------------------------------
% Definitions
Tch = 1;
a = 20;
f_c = 0;

% Chrip Signal
chirp = exp(1j*pi*(f_c).*t).*exp(1j*(pi*a.*t.^2));

% Ambiguity function (auto-correlation)
[r_chirp, lags] = xcorr(chirp, 10, 'normalized');

% Fast Fourier Transform
fft_size = 2^ceil(log2(length(chirp)));
CHIRP = fftshift(fft(chirp, fft_size))/(fft_size);
freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
                 linspace(-fft_size/2, -1, fft_size/2)])*(Fs/fft_size);

% Plotting
figure(1),
plot(t, real(chirp), 'linewidth', 2), grid on
xlabel('Time, $t$ [sec]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Amplitude, $x(t)$', 'interpreter', 'latex', 'fontsize', 14)
title('In-phase component of chirp signal', 'interpreter', 'latex',...
    'fontsize', 14)
a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', 14)

figure(2),
plot(lags, abs(r_chirp), 'linewidth', 2), grid on
xlabel('Lags [samples]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$A(t,f_{\rm D})$', 'interpreter', 'latex', 'fontsize', 14)
title('Ambiguity Function', 'interpreter', 'latex', 'fontsize', 14)
a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', 14)

figure(3),
plot(freq, abs(CHIRP), 'linewidth', 2), grid on
xlabel('Frequency, $f$ [Hz]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Magnitude $|X(f)|$', 'interpreter', 'latex', 'fontsize', 14)
title('Spectrum of the Radar Signal', 'interpreter', 'latex',...
    'fontsize', 14)
xlim([-150 150]);
a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', 14)

% -------------------------------------------------------------------------
%                           Sinusoid Pulse
% -------------------------------------------------------------------------
% Definitions
f_c = 10;  % Sinusoid frequency
phi = 0;  % Phase

% Sinusoidal signal
sinusoid = exp(1j*(2*pi*f_c*t + phi));

% Ambiguity function
[r_sin, lags] = xcorr(sinusoid, 10, 'normalized');

% FFT
fft_size = 2^ceil(log2(length(sinusoid)));
SINUSOID = fftshift(fft(sinusoid, fft_size))/(fft_size);
freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
                 linspace(-fft_size/2, -1, fft_size/2)])*(Fs/fft_size);

% Plotting
figure(4),
plot(t, real(sinusoid), 'linewidth', 2), grid on
xlabel('Time, $t$ [sec]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Amplitude, $x(t)$', 'interpreter', 'latex', 'fontsize', 14)
title('Real part of complex sinusoidal signal', 'interpreter', 'latex',...
     'fontsize', 14)
a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', 14)

figure(5),
plot(lags, abs(r_sin), 'linewidth', 2), grid on
xlabel('Lags [samples]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$A(t,f_{\rm D})$', 'interpreter', 'latex', 'fontsize', 14)
title('Ambiguity Function', 'interpreter', 'latex', 'fontsize', 14)
a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', 14)

figure(6),
plot(freq, abs(SINUSOID), 'linewidth', 2), grid on
xlabel('Frequency, $f$ [Hz]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Magnitude $|X(f)|$', 'interpreter', 'latex', 'fontsize', 14)
title('Spectrum of the Radar Signal', 'interpreter', 'latex',...
    'fontsize', 14)
xlim([-50 50]);
a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', 14)

% -------------------------------------------------------------------------
%                             Square Wave
% -------------------------------------------------------------------------
% Definitions
d_c = .4;  % Duty-cycle

% Squared wave
sqr = [zeros(1, ((1-d_c)*T*Fs/2))...
       ones(1, (d_c*T*Fs))...
       zeros(1, ((1-d_c)*T*Fs)/2)];

% Ambiguity function
[r_sqr, lags] = xcorr(sqr, 10, 'normalized');

% FFT
fft_size = 2^ceil(log2(length(sqr)));
SQR = fftshift(fft(sqr, fft_size))/(fft_size);
freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
                 linspace(-fft_size/2, -1, fft_size/2)])*(Fs/fft_size);

% Plotting
figure(7),
plot(t, real(sqr), 'linewidth', 2), grid on
xlabel('Time, $t$ [sec]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Amplitude, $x(t)$', 'interpreter', 'latex', 'fontsize', 14)
title('Square wave signal', 'interpreter', 'latex',...
    'fontsize', 14)
ylim([-0.5 1.5])
a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', 14)
xlim([0 1])

figure(8),
plot(lags, abs(r_sqr), 'linewidth', 2), grid on
xlabel('Lags [samples]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$A(t,f_{\rm D})$', 'interpreter', 'latex', 'fontsize', 14)
title('Ambiguity Function', 'interpreter', 'latex', 'fontsize', 14)
a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', 14)

figure(9),
plot(freq, abs(SQR), 'linewidth', 2), grid on
xlabel('Frequency, $f$ [Hz]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Magnitude $|X(f)|$', 'interpreter', 'latex', 'fontsize', 14)
title('Spectrum of the Radar Signal', 'interpreter', 'latex',...
    'fontsize', 14)
xlim([-150 150])
a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', 14)


