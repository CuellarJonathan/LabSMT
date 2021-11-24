% Exercise 2:
%
% Jonathan S. Cuellar - j.silvacuellar@poli.ufrj.br
% November 2021
%

clc
clear
close all

%% -------------------------------------------------------------------------
%                               Definitions
%% -------------------------------------------------------------------------
M = 10;  % Number of pulses
Mt = 1; % Número de antenas
T = 40*1e-6;  % Pulse duration PRI
PRF = 1/T;  % Pulse repetition frequency
dc = 0.5;  % Duty-cycle
T_ch = dc*T;  % Chirp duration
BW = 500*1e6;  % Signal bandwidth
f_c = 80*1e9;  % Carrier frequency

v = 300;  % Relative target speed, km/h
d = 300;  % Em metros

v_clutter = 0;  % Relative clutter speed, km/h
ac = 0;  % Aceleration m/s^2

d_clutter = 0;  % Distance in meters, clutter
c = 299792458;  % Speed of light, m/s
lamb_c = c/f_c; % Wave-Length
d_w = 20;  % Range swath
delta_d = c/(2*BW);  % Range resolution
L = round(d_w/delta_d); % Total number of samples in fast-time
kmh2ms = @(x)x/(3.6);  % Function to convert km/h to m/s

f_c2 =0;  % Reduce computational burden.
a = BW/T_ch;  % Chirp slope
%
% Sampling rate estimation:
% To be able to detect a velocity of 300 km/h, one must have a minimal
% sampling rate. The sampling period is given by the following equation.
% Ts = v*(360/1000)*T*(1/c).
% To alleviate the high sampling rate, one can multiply the this time by M
% and use the 1st and the M-th pulse to estimate the velocity.
%

Ts = M*kmh2ms(abs(v))*T*(1/c);
Fs = 1/Ts;  % Sampling rate
% Fs = 2*BW;

tau = 2*d/c;  % How much time to m=0 peak
Tr = 0.5*1e-6; %  Properly choosing the sampling period Tr
n0 = zeros(M);
% n0 = round(tau*Fs-(4*L/5));  % I want n0 as small as possible
% n0 = round(tau/Tr);
fft_size = 2^12;  % High FFT size is a problem in range estimation

%% -------------------------------------------------------------------------
%                    Signal generation and processing
% -------------------------------------------------------------------------
% Time axis
t_ch = linspace(0, T_ch, T_ch*Fs);

% Transmitted signal (reference signal for xcorr)
tx1 = exp(1j*pi*(f_c2).*t_ch).*exp(1j*pi*a.*t_ch.^2);
tx = [tx1 zeros(1, round((T-T_ch)*Fs))];

% Chrip Signal
chirp = tx;

% Ambiguity function (auto-correlation)
[r_chirp, lags] = xcorr(chirp, 150, 'normalized');

% Memory allocation
A = zeros(2*L + 1, M);
snr = 150000;

d_m = zeros(M);

rad2deg = @(x) x*180/pi;

fprintf('Car Movement:\n\n');

for m = 1:M
    
    d_m(m) =  d + (m-1)*kmh2ms(v)*(1+ac)*T;
    
    fprintf('New Distance: t = %f, d_m = %.4f m, v = %.4f kmh. \n', (m-1)*T, d_m(m), v);
    
    % Doppler effects
    f_D = ((2*kmh2ms(v))/(c)) * f_c;  % Doppler effect
    nu = f_D/f_c;  % Doppler ratio
    f_D_clutter = ((2*kmh2ms(v_clutter))/(c)) * f_c;
    nu_clutter = f_D_clutter/f_c;

    tau = 2*d_m(m)/c;  % How much time to m=0 peak
    n0(m) = round(tau*Fs-(4*L/5));

    N = L+n0(m);

    rx = exp(1j*(2*pi*(f_c2*(t_ch-m*T-tau)-f_D*nu*(t_ch-m*T-tau)+(2*a*pi*f_c2*((t_ch-m*T-tau-nu*(t_ch-m*T-tau)).^2))/2)));

    rx_m = awgn(rx, snr, 'measured');
    A(:, m) = xcorr(rx_m(1, n0(m)+1:N), tx(1, n0(m)+1:N), L);
end
%% -------------------------------------------------------------------------
%                   Slow-time vs. Fast-time matrix
%% -------------------------------------------------------------------------

A = A(L+2:end, :);
Afft=zeros(fft_size,M);
for m = 1:M
    Afft(:,m)=fftshift(fft(A(:,m),fft_size));
end
% LL = linspace(0, L-1, L);
LL = linspace(0, L-1, fft_size);
MM = linspace(0, M-1, M);

CHIRP = fftshift(fft(chirp, fft_size))/(fft_size);
freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
                 linspace(-fft_size/2, -1, fft_size/2)])*(Fs/fft_size);

t = t_ch;

fig3 = figure;
% Plotting
subplot(3,1,1);
plot(t, real(tx1), 'linewidth', 1), grid on
xlim([0 T]);
xlabel('Time, $t$ [sec]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Amplitude, $x(t)$', 'interpreter', 'latex', 'fontsize', 14)
title('In-phase component of chirp signal', 'interpreter', 'latex',...
    'fontsize', 14)
% ab = get(gca, 'XTickLabel');
% set(gca, 'XTickLabel', ab, 'fontsize', 14)

% Plotting
subplot(3,1,2);
plot(lags, abs(r_chirp), 'linewidth', 1), grid on
xlabel('Lags [samples]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$A(t,f_{\rm D})$', 'interpreter', 'latex', 'fontsize', 14)
title('Ambiguity Function', 'interpreter', 'latex', 'fontsize', 14)
xlim([-150 150]);
% ab = get(gca, 'XTickLabel');
% set(gca, 'XTickLabel', ab, 'fontsize', 14)

% Plotting
subplot(3,1,3);
plot(freq, abs(CHIRP), 'linewidth', 1), grid on
xlabel('Frequency, $f$ [Hz]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Magnitude $|X(f)|$', 'interpreter', 'latex', 'fontsize', 14)
title('Spectrum of the Radar Signal', 'interpreter', 'latex',...
    'fontsize', 14)
xlim([-1e9 1e9]);
% ab = get(gca, 'XTickLabel');
% set(gca, 'XTickLabel', ab, 'fontsize', 14)

figure;
surf(MM, LL,abs(Afft(:,:)))
xlabel('Slow-time, $m$', 'interpreter', 'latex')
ylabel('Fast-time, $l$', 'interpreter', 'latex')
zlabel('Normalized dechirped signal', 'interpreter', 'latex')
%     ylim([0 LL]);
%     xlim([0 MM]);
title(['Slow-time vs. Fast-time matrix'], 'interpreter', 'latex')


% Frequência x Tempo

txA = 1e8;
[fsh1,t1] = fshift(T,T_ch,BW,0,f_c,M*Mt,txA,tau(1));
[fsh2,t2] = fshift(T,T_ch,BW,tau(1),f_c,M*Mt,txA,tau(1));
fsh3 = fsh1+fsh2;

fig1 = figure;

% subplot(2,1,1);
plot(t1,fsh1);
hold on
plot(t2,fsh2);
hold off
title('Frequency shift');
xlim([0 M*Mt*T+tau(1)]);
xlabel('Tempo (s)');ylabel('Frequência (Hz)');
legend('Transmitido','Refletido')

% subplot(2,1,2);
% plot(t2,fsh3);
% title('Frequency shift SUM');
% xlim([0 M*Mt*T+tau(1)]);
% xlabel('Tempo (s)');ylabel('Frequência (Hz)');


%% -------------------------------------------------------------------------
%                           Fast-time DFT
%% -------------------------------------------------------------------------
ms2kmh = @(x) 3.6*x;  % Function to convert m/s to km/h
freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
                 linspace(-fft_size/2, -1, fft_size/2)])*(Fs/fft_size);

v_all = zeros(M);
d_hat_fft = zeros(M);

for m = 1:M
    d_fft_past = (n0(m)*c)/(2*Fs);  % My range estimation a priori
    % Where should it be?
    f_tt = (d_fft_past*2*a)/c;
    AA = fftshift(fft((A(:, m)), fft_size));
    [value, f_idx] = max(abs(AA));
    d_hat_fft(m) = ((LL(f_idx)+n0(m))/Fs)*c/2;
    f_t = freq(f_idx);
%         v_t = c*abs(f_t - f_tt)/2/f_c;
    d_fft = (f_t*c)/(2*a);
%         v_t = (c*f_t/2/f_c)-((2*d_fft*a)/2/f_c);
%         v_all(n,m) = (d_fft - d_fft_past)/T;
%         v_t = abs(d_fft(n) - d_fft_past(n))/2/T;
%     fprintf('Estimated range, for m = %d: %.4f m.\n', m, d_fft)
    txt = {['Distance Antena: ' num2str(d_hat_fft(m)) ' m'],...
        ['$f_t = ' num2str(f_t) '$ Hz']};
    if (m == 1) || (m == 2) || (m == round(M/2)) || (m == M)
        figure,
        plot(freq, 10*log10(abs(AA)), 'linewidth', 2); 
        hold on
        plot(f_t, 10*log10(value), 'o', 'linewidth', 2, 'markersize', 10);
        xline(f_tt, '--', 'linewidth', 2); 
        hold off
        grid on
        text(7e6, 50, txt, 'interpreter', 'latex', 'fontsize', 12);
        xlabel('Frequency, $f$ Hz', 'interpreter', 'latex');
        ylabel('Maginitude, dB', 'interpreter', 'latex');
        xlim([-4e8 4e8]);
        ylim([0 60]);

        title(['Fast-time FFT, m = ' num2str(m)], 'interpreter', 'latex');

    end
end


fprintf('\n')


%% -------------------------------------------------------------------------
%                           Distance Estimation
%% -------------------------------------------------------------------------
fprintf('Distance estimation (maximum value in matrix):\n\n')
d_hat = zeros(M);
A_max = zeros(M);
v_hat = zeros(M);
d_past = d;
for m = 1:M
    [~, idx] = max(abs(A(:, m)));
    d_hat(m) = ((LL(idx)+n0(m))/Fs)*c/2;
    if(d_hat(m) < 0) 
        d_hat(m) = 0; 
    end
    A_max(m) = max(abs(A(:, m)));

    if (((m >= 1)&&(m <= 4)) || (m == round(M/2)) || (m == M))
        fprintf('Estimated Distance Antena com m = %d: %.4f m.\n', m, d_hat(m))
    end
    if(m==1)
        d_past = d_hat(m);
    end
    v_hat(m) = (d_hat(m)-d_past)/(m*T);
    d_past = d_hat(1);
    if (((m >= 1)&&(m <= 4)) || (m == round(M/2)) || (m == M))
        fprintf('Estimated Distance m = %d: %.4f m.\n', m, d_hat(m))
        fprintf('Estimated Velocity m = %d: %.4f kmh.\n', m, ms2kmh(v_hat(m)));
    end
end

figure,
plot(linspace(0,M*T,M), d_hat(:,1), 'linewidth', 2), grid on
text(0.25*max(M*T), max(d_hat(:,1))-0.8*(max(d_hat(:,1))-min(d_hat(:,1))), ['Estimated Velocity M = ' num2str(M) ': ' num2str(ms2kmh(v_hat(M))) ' km/h'],...
    'interpreter', 'latex', 'fontsize', 12);
xlabel('Time, Seconds', 'interpreter', 'latex')
ylabel('Range, Meters', 'interpreter', 'latex')
title('Slow-Time Configuration', 'interpreter', 'latex')


%% -------------------------------------------------------------------------
%                       DFT from slow-time
%% -------------------------------------------------------------------------

v_f = ms2kmh((d_hat_fft(M)-d_hat_fft(1))/(M*T));
fprintf('\nFFT Estimated Velocity M = %d: %.4f kmh.\n', M, v_f);

% m = round(M);
% [~, idx] = max(abs(A(1,:, m)));
% % A_slw = A(1,idx, :);
% 
% Acell1 = num2cell(A(1,idx,:));
% Acell2 = reshape(Acell1,1, M);
% A_slw = cell2mat(Acell2);
% 
% % My simple high-pass filters:
% b1 = [1 -1];
% b2 = [1 -2 1];
% 
% y1 = filter(b1, 1, A_slw);
% y2 = filter(b2, 1, A_slw);
% 
% fft_size = 2^ceil(log2(length(A_slw)));
% freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
%                  linspace(-fft_size/2, -1, fft_size/2)])*(PRF/fft_size);
% Y1 = fftshift(fft(y1, fft_size));
% Y2 = fftshift(fft(y2, fft_size));
% 
% 
% figure,
% plot(freq, abs(Y2), 'linewidth', 2), grid on
% xline(f_D, '--', 'linewidth', 2)
% xlabel('Frequency, $f$ HZ', 'interpreter', 'latex')
% ylabel('Magnitude, dB', 'interpreter', 'latex')
% title(['Slow-time FFT using two-pulse canceller for max value in m = ' num2str(m)], 'interpreter', 'latex')


% EoF

function [fsh,t] = fshift(T,T_ch,BW,tau,fc,M,txA,tend)
    t = 0:1/txA:(M*T+tend)+1/txA;
    fsh = fc;
    for m = 0:1:(M-1)
        w = @(tw) ((tw > (m*T+tau)) & (tw <= (m*T+T_ch+tau)));
        fsh = fsh + (((BW)/(T)).*(t-m*T-tau)).*w(t);
    end
end

