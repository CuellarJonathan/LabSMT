clc
clear
close all

%% -------------------------------------------------------------------------
%                               Definitions
%  -------------------------------------------------------------------------

T = 40*1e-6;  % Pulse duration PRI
BW = 500*1e6;  % Signal bandwidth
f_c = 80*1e9;  % Carrier frequency
v = 150;  % Relative target speed, km/h
d = 300;  % Distance in meters, target
d_w = 10;  % Range swath
c = 299792458;  % Speed of light, m/s
delta_d = c/(2*BW);  % Range resolution
L = round(d_w/delta_d); % Total number of samples in fast-time
a = BW/T;

sigma_rcs=linspace(0.3,9,50); %Area de reflexao(m^2). Considerando 0.3 como uma placa de transito e 9 com um caminhao
r=linspace(100,500,50); %distancia entre os radares

%% -------------------------------------------------------------------------
%                               Max Radar
%  -------------------------------------------------------------------------

disp('Tempo execução max radar:');
tic
K = zeros(length(sigma_rcs),length(r));
for m = 1:length(r)
    K(:,m)=T*BW*sqrt(sigma_rcs(m))./r;
end
toc

K_min = min(min(K));
K_max = max(max(K));

%% -------------------------------------------------------------------------
%                               Leakage
%  -------------------------------------------------------------------------

leakage_area = 1./((pi.*BW.*T./K).^2);

delta_tau = linspace(T/K_min,T/K_max,length(K));

leakage = 1./((pi.*BW.*delta_tau).^2);

syms t;

leakage_int = zeros(1,length(delta_tau));

disp('Tempo execução leakage:');
tic
for k = 1:1:length(delta_tau)
    leakage_int(k) = int((1/T)*exp(1j*2*pi*(f_c*t+(1/2)*a*t^2))*exp(-1j*2*pi*(f_c*(t-delta_tau(k))+(1/2)*a*(t-delta_tau(k))^2)),t,0,T)^2;
end
toc

%% -------------------------------------------------------------------------
%                               Plots
%  -------------------------------------------------------------------------

fig1 = figure;
subplot(2,2,1);
surf(sigma_rcs,r,K)
colorbar
xlabel('\sigma, m^{2}')
ylabel('r, m')
zlabel('K','interpreter', 'latex')
xlim([0 9]);
title('Quantidade maxima de radares', 'interpreter', 'latex')

subplot(2,2,2);
surf(sigma_rcs,r,leakage_area)
colorbar
xlabel('\sigma, m^{2}')
ylabel('r, m')
zlabel('Leakage Max','interpreter', 'latex')
xlim([0 9]);
title('Vazamento por Delta tau', 'interpreter', 'latex')

subplot(2,2,3);
plot(delta_tau,leakage);
% xlabel('\Delta\tau, s')
% ylabel('Leakage Max')
% xlim([0 9]);
% title('Vazamento por Delta tau', 'interpreter', 'latex')

% subplot(2,2,4);
hold on
% surf(sigma_rcs,r,abs(real(leakage_int)))
plot(delta_tau,abs(real(leakage_int)))
% colorbar
% zlabel('Leakage','interpreter', 'latex')
hold off

xlim([0 1e-6]);
% ylim([0 400]);
xlabel('\Delta\tau, s')
ylabel('Leakage')
title('Vazamento por Delta tau', 'interpreter', 'latex')
legend('Leakage Max','Leakage \Delta\tau');

set(fig1, 'Position',  [200, 50, 1024, 720])
