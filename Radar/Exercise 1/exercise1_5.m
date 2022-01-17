%  Try different values of M ranging from 1 to 10

clc
clear
close all

V = 1247;  % Total amount of target velocity samples
c = 299792458;  % Speed of light, m/s
% K band = 18:26.5 GHz
% K_a band = 26.5:40 GHz
f_c = 27*1e9;  % Carrier frequency
v = 360.*linspace(-360, 360, V)./1e3;  % Reference velocity
f_D = ((2*v)./(c - v))*f_c;  % Reference f_D
nu = f_D/f_c;
T_ch = 10/max(f_D);  % Chirp duration, with normalization
T = 2*T_ch;  % pulse duration in seconds
a = 20/T;
Fs = (1/T_ch)*1e3;  % Sampling rate (1e3 is my PCs limit for ok graph)
maxlag = round((T_ch*Fs));  % max number of lags
t_ch = linspace(0, T_ch, T_ch*Fs);
M = 10;

% Targets info:
v_targ = 360*120*1e-3;  % Target speed in m/s
d_targ = 15000;  % Target's distance in m
f_D_targ = ((2*v_targ)/(c-v_targ))*f_c;
nu_targ = f_D_targ/f_c;
t_d = 2*d_targ/c;  % Time interval between tx and rx
t_samp = round(t_d*Fs);
snr = -15;  % Signal to noise ratio, dB

rx_ch = exp(1j*pi*(0+f_D_targ).*t_ch).*exp(1j*(pi*a*(1+nu_targ).*t_ch.^2));
rx = [zeros(1, t_samp) rx_ch zeros(1, round((T-T_ch)*Fs) - t_samp)];

fft_size = 2^14;
RX = fftshift(fft(rx, fft_size))/(fft_size);
freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
                 linspace(-fft_size/2, -1, fft_size/2)])*(Fs/fft_size);

tic;
parfor idx = 1:V
    % Reference signal for sweep
    ref = exp(1j*pi*(0+f_D(idx)).*t_ch).*exp(1j*(pi*a*(1+nu(idx)).*t_ch.^2));
    for m = 1:M
        % Received signal
        d_m = d_targ + (m-1)*T*v_targ;
        t_samp_m = round(2*(d_m/c)*Fs);
        rx = [zeros(1, t_samp_m) rx_ch zeros(1, round((T-T_ch)*Fs) - t_samp_m)];
        % Reference signal
        d_ref_m = (m-1)*T*v(idx);
        t_ref_samp_m = round(2*(d_ref_m/c)*Fs);
        ref_pulse = [zeros(1, t_ref_samp_m) ref zeros(1, round((T-T_ch)*Fs) - t_ref_samp_m)];
        % Ambiguity function
        A_m = xcorr(awgn(rx, snr, 0), ref_pulse, maxlag, 'normalized');
        if m > 1
            A_sum = A_sum + A_m;
        else
            A_sum = A_m;
        end
        A(: , idx, m) = A_sum/m;
    end
end
toc()

% Distance and velocity estimation
[val, idx1] = max(A(:, :, 3));
[~, f_D_hat_idx] = max(val);  % Velocity
t_samp_hat = idx1(f_D_hat_idx) - (maxlag+1);
d_hat = ((t_samp_hat/Fs)*c)/2;  % estimated distance
f_D_hat = f_D(f_D_hat_idx);
v_hat = (f_D_hat*c)/(2*f_c + f_D_hat);  % Estimated velocity

fprintf('\n')
fprintf('RESULTS for M = 3:\n')
fprintf('---------------------------------------\n')
fprintf('Estimated velocity: %.4f km/h.\n', (v_hat/360)*1e3)
fprintf('Estimated distance: %.4f m.\n', d_hat)

[val, idx1] = max(A(:, :, 6));
[~, f_D_hat_idx] = max(val);  % Velocity
t_samp_hat = idx1(f_D_hat_idx) - (maxlag+1);
d_hat = ((t_samp_hat/Fs)*c)/2;  % estimated distance
f_D_hat = f_D(f_D_hat_idx);
v_hat = (f_D_hat*c)/(2*f_c + f_D_hat);  % Estimated velocity

fprintf('\n')
fprintf('RESULTS for M = 6:\n')
fprintf('---------------------------------------\n')
fprintf('Estimated velocity: %.4f km/h.\n', (v_hat/360)*1e3)
fprintf('Estimated distance: %.4f m.\n', d_hat)

[val, idx1] = max(A(:, :, 10));
[~, f_D_hat_idx] = max(val);  % Velocity
t_samp_hat = idx1(f_D_hat_idx) - (maxlag+1);
d_hat = ((t_samp_hat/Fs)*c)/2;  % estimated distance
f_D_hat = f_D(f_D_hat_idx);
v_hat = (f_D_hat*c)/(2*f_c + f_D_hat);  % Estimated velocity

fprintf('\n')
fprintf('RESULTS for M = 10:\n')
fprintf('---------------------------------------\n')
fprintf('Estimated velocity: %.4f km/h.\n', (v_hat/360)*1e3)
fprintf('Estimated distance: %.4f m.\n', d_hat)

% Plot
lags = linspace(-1, 1, 2*maxlag+1);
fd_norm = linspace(-10, 10, V);
[XX, YY] = meshgrid(fd_norm, lags);

figure(1),
s = surf(XX, YY, abs(A(:, :, 3)));
set(s, 'linestyle', 'none')
xlabel('$F_D \tau$', 'interpreter', 'latex')
ylabel('$t/\tau$', 'interpreter', 'latex')
zlabel('$A(F_D, t)$', 'interpreter', 'latex')
title('Ambiguity Function for $M=3$', 'interpreter', 'latex')

figure(2),
s = surf(XX, YY, abs(A(:, :, 6)));
set(s, 'linestyle', 'none')
xlabel('$F_D \tau$', 'interpreter', 'latex')
ylabel('$t/\tau$', 'interpreter', 'latex')
zlabel('$A(F_D, t)$', 'interpreter', 'latex')
title('Ambiguity Function for $M=6$', 'interpreter', 'latex')

figure(3),
s = surf(XX, YY, abs(A(:, :, 10)));
set(s, 'linestyle', 'none')
xlabel('$F_D \tau$', 'interpreter', 'latex')
ylabel('$t/\tau$', 'interpreter', 'latex')
zlabel('$A(F_D, t)$', 'interpreter', 'latex')
title('Ambiguity Function for $M=10$', 'interpreter', 'latex')

figure(4),
plot(freq, abs(RX), 'linewidth', 2), grid on
xlabel('Frequency, $f$ [Hz]', 'interpreter', 'latex')
ylabel('Magnitude, $|X(f)|$', 'interpreter', 'latex')
title('Spectrum of the Radar Signal', 'interpreter', 'latex')
xlim([-2e4 2e4])


% EoF