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
% may 2021
%

clc
clear
close all

%% -------------------------------------------------------------------------
%                               Definitions
%% -------------------------------------------------------------------------

% Global definitions
c = 299792458;  % Speed of light, m/s

% Chirp definitions
T = 80*1e-6;  % Pulse duration PRI
PRF = 1/T;  % Pulse repetition frequency
dc = 0.5;  % Duty-cycle
T_ch = dc*T;  % Chirp duration
BW = 500*1e6;  % Signal bandwidth
f_c = 80*1e9;  % Carrier frequency
a = BW/T_ch;  % Chirp slope

% Target definitions
v = 100;  % Relative target speed, km/h
d = 300;  % Distance in meters, target

ac = 0;  % Aceleration m/s^2

% Clutter definitions
v_clutter = 200;  % Relative clutter speed, km/h
d_clutter = 300;  % Distance in meters, clutter

d_w = 20;  % Range swath
delta_d = c/(2*BW);  % Range resolution
L = round(d_w/delta_d); % Total number of samples in fast-time
kmh2ms = @(x)x/(3.6);  % Function to convert km/h to m/s

% Doppler effects
f_D = ((2*kmh2ms(v))/(c)) * f_c;  % Doppler effect
nu = f_D/f_c;  % Doppler ratio
f_D_clutter = ((2*kmh2ms(v_clutter))/(c)) * f_c;
nu_clutter = f_D_clutter/f_c;

f_c2 =0;  % Reduce computational burden.

%
% Sampling rate estimation:
% To be able to detect a velocity of 300 km/h, one must have a minimal
% sampling rate. The sampling period is given by the following equation.
% Ts = v*(360/1000)*T*(1/c).
% To alleviate the high sampling rate, one can multiply the this time by M
% and use the 1st and the M-th pulse to estimate the velocity.
%
M = 128;  % Number of pulses
Ts = M*kmh2ms(abs(v))*T*(1/c);
Fs = 1/Ts;  % Sampling rate
% Fs = 2*BW;

tau = 2*d/c;  % How much time to m=0 peak
n0 = zeros(M);
n0(1) = round(tau*Fs-(4*L/5));  % I want n0 as small as possible
fft_size = 2^12;  % High FFT size is a problem in range estimation


%% -------------------------------------------------------------------------
%                    Signal generation and processing
% -------------------------------------------------------------------------
% Time axis
t_ch = linspace(0, T_ch, T_ch*Fs);

% Transmitted signal (reference signal for xcorr)
tx = exp(1j*pi*(f_c2).*t_ch).*exp(1j*pi*a.*t_ch.^2);
tx = [tx zeros(1, round((T-T_ch)*Fs))];

% Memory allocation
A = zeros(2*L + 1, M);
N = L+n0(1);
snr = 150000;
for m = 1:M
    % Target
    d_m = d + (m-1)*kmh2ms(v)*(1+ac)*T;
    t_samp_m = round(2*(d_m/c)*Fs);
%     rx_m = awgn([zeros(1, t_samp_m) rx zeros(1, round((T-T_ch)*Fs) - t_samp_m)] +...
%         rx_clutter, snr, 'measured');

    % Doppler effects
    f_D = ((2*kmh2ms(v))/(c)) * f_c;  % Doppler effect
    nu = f_D/f_c;  % Doppler ratio
    f_D_clutter = ((2*kmh2ms(v_clutter))/(c)) * f_c;
    nu_clutter = f_D_clutter/f_c;

    tau = 2*d_m/c;  % How much time to m=0 peak
    n0(m) = round(tau*Fs-(4*L/5));  % I want n0 as small as possible
    N = L+n0(m);

    % Received signal
    xx = exp(1j*pi*(f_c2+f_D).*t_ch).*exp(1j*pi*a*(1+f_D/f_c2).*t_ch.^2);
    t_clutter = round(2*(d_m/c)*Fs);
    rx_clutter = [zeros(1, t_clutter)...
              exp(1j*pi*(f_c2+f_D_clutter).*t_ch).*exp(1j*pi*a*(1+nu_clutter).*t_ch.^2)...
              zeros(1, round((T-T_ch)*Fs) - t_clutter)];

    rx_m = awgn(rx_clutter, snr, 'measured');
    % Dechirping process
    A(:, m) = xcorr(rx_m(1, n0(m)+1:N), tx(1, n0(m)+1:N), L);
end
%% -------------------------------------------------------------------------
%                   Slow-time vs. Fast-time matrix
%% -------------------------------------------------------------------------

% A = A(L+2:end, :);
Afft=zeros(fft_size,M);
for m = 1:M
    Afft(:,m)=fftshift(fft(A(:,m),fft_size));
end
% LL = linspace(0, L-1, L);
LL = linspace(0, L-1, fft_size);
MM = linspace(0, M-1, M);

figure(9),
surf(MM, LL,abs(Afft))
xlabel('Slow-time, $m$', 'interpreter', 'latex')
ylabel('Fast-time, $l$', 'interpreter', 'latex')
zlabel('Normalized dechirped signal', 'interpreter', 'latex')
title('Slow-time vs. Fast-time matrix', 'interpreter', 'latex')
% saveas(figure(9), "2_slowTimeVsFastTimeMatrix.png");


%% -------------------------------------------------------------------------
%                           Fast-time DFT
%% -------------------------------------------------------------------------
ms2kmh = @(x) 3.6*x;  % Function to convert m/s to km/h
fprintf('Distance estimation, (DFT):\n')
freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
                 linspace(-fft_size/2, -1, fft_size/2)])*(Fs/fft_size);

d_fft = zeros(M);
d_fft(1) = (n0(1)*c)/(2*Fs);  % My range estimation a priori
% Where should it be?
f_tt = (d_fft(1)*2*a)/c;
for m = 1:M
%     d_fft_past = d_fft(m);
    AA = fftshift(fft((A(:, m)), fft_size));
    [value, f_idx] = max(abs(AA));
    f_t = freq(f_idx);
    d_fft(m) = (f_t*c)/(2*a);
%     fprintf('Estimated range, for m = %d: %.4f m.\n', m, d_fft)
    txt = {['distance: ' num2str(d_fft(m)) ' m'],...
        ['$\hat{f}_t = ' num2str(f_t) '$ Hz']};
    if (m == 1) || (m == 2) || (m == M)
        figure(m),
        plot(freq, 10*log10(abs(AA)), 'linewidth', 2), hold on,
        plot(f_t, 10*log10(value), 'o', 'linewidth', 2, 'markersize', 10)
        xline(f_tt, '--', 'linewidth', 2), hold off, grid on
        text(1e8, 30, txt, 'interpreter', 'latex', 'fontsize', 12)
        xlabel('Frequency, $f$ Hz', 'interpreter', 'latex')
        ylabel('Maginitude, dB', 'interpreter', 'latex')
        xlim([-4e8 4e8]);
        ylim([0 40]);
        
        title(['Fast-time FFT, m = ' num2str(m)], 'interpreter', 'latex')
%         saveas(figure(m), strcat("2_spectrum",num2str(m),".png"));
        
    end
end

v_hat = (d_fft(M) - d_fft(3))/(M*T);
fprintf('Estimated velocity M = %d: %.4f kmh.\n', M, ms2kmh(v_hat))


fprintf('\n')


%% -------------------------------------------------------------------------
%                           Distance Estimation
%% -------------------------------------------------------------------------
fprintf('Distance estimation (maximum value in matrix):\n')
for m = 1:M
    [~, idx] = max(abs(A(:, m)));
    d_hat1 = ((LL(idx)+n0(m))/Fs)*c/2;
    
    if (m == 1) || (m == 2) || (m == M)
    fprintf('Estimated distance m = %d: %.4f m.\n', m, d_hat1)
    end
end




%% -------------------------------------------------------------------------
%                       DFT from slow-time
%% -------------------------------------------------------------------------
m = round(M);
[~, idx] = max(abs(A(:, m)));
A_slw = A(idx, :);

% My simple high-pass filters:
b1 = [1 -1];
b2 = [1 -2 1];

y1 = filter(b1, 1, A_slw);
y2 = filter(b2, 1, A_slw);

fft_size = 2^ceil(log2(length(A_slw)));
freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
                 linspace(-fft_size/2, -1, fft_size/2)])*(PRF/fft_size);
Y1 = fftshift(fft(y1, fft_size));
Y2 = fftshift(fft(y2, fft_size));


% figure,
% plot(freq, abs(Y1), 'linewidth', 2), grid on
% xline(f_D, '--', 'linewidth', 2)
% xlabel('Frequency, $f$ HZ', 'interpreter', 'latex')
% ylabel('Magnitude, dB', 'interpreter', 'latex')
% title(['Slow-time FFT using two-pulse canceller for max value in m = ' num2str(m)], 'interpreter', 'latex')


% EoF

