%   In automotive communications, assume a bandwidth of total_BW = 600 MHz, a 
% carrier of f c = 80 GHz, a pulse duration of T = 40μs, and a range of 
% detection of d = 300 m. Assume the target is moving at v = 150 km/h.
% Consider that we need to cover a range swath of 10m, so that you can
% choose L and M = 4.
%   Divide the bandwidth among three vehicles each using its own OFDM 
% waveform to estimate their relative distances and Doppler shifts.
% Simulate the situation and discuss the results, all with an SNR of 50 dBs.

% For more information on the properties of Stepped Frequency OFDM radar solutions,
% please check:
% B. Schweizer, C. Knill, D. Schindler and C. Waldschmidt,
% "Stepped-Carrier OFDM-Radar Processing Scheme to Retrieve High-Resolution
% Range-Velocity Profile at Low Sampling Rate," in IEEE Transactions on Microwave
% Theory and Techniques, vol. 66, no. 3, pp. 1610-1618, March 2018,
% doi: 10.1109/TMTT.2017.2751463.
% on: https://ieeexplore.ieee.org/document/8053909

clc
clear
close all


%% Global Properties
c = 299792458; % m/s

%% Radio Resources
total_BW = 600*1e6; % 600 MHz
fc = 80*1e9; % 80 GHz
%T = 40*1e-6; % 40 micro s was not possible because of vua

%% Radar Properties
d = 300; % range of detection = 300 m
% range_swath = 10; % 10 m

%% OFDM Properties
BW = total_BW/3;

ofdm_N = 800; % Number of OFDM subcarriers
n_subcarriers = ofdm_N;

delta_f = BW/n_subcarriers
Ts = 1/delta_f;

numberOfDataBitsPerSymbol = n_subcarriers;

% Setup cycle prefix
Tcp = ((2*d)/c)*(2); % Double the delay relative to max distance
Tcp = min(Tcp, Ts);

% Symbol duration
T = Ts + Tcp

% Setup OFDM Stepped Frequency
ofdm_L = 64; % Number of transmitted sumbols
n_symbols = ofdm_L;

ofdm_M = 2; % Number of subsymbols in a block

subsymbol_L = ofdm_M*ofdm_L; % number of symbols per frame
subsymbol_N = round(ofdm_N/ofdm_M); % Number of carriers per frame

subsymbol_BW = BW/ofdm_M;

%% Transmission Properties
SNR = -10;

%% Target Properties
numberOfVehicles = 3;
v = 150*(1/3.6) % 150 km/h = 41.67 m/s;

%% Check constraints
vua = c/(4*fc*T*ofdm_M) % unambiguosly measurable velocity

rangeResolution = c/(2*BW)

velocityResolution = c/(2*fc*T*(ofdm_M*(ofdm_L-1)+1))

Rua = c/(2*delta_f) % > d

Rmax = (Tcp*c/2) % > d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sampling data
dataInput = randi([0 3], numberOfDataBitsPerSymbol*n_symbols, 1);

%% OFDM Transmission
% QPKS Modulation
n = 0:pi/4:2*pi-pi/4;
in_phase = cos(n+pi/4);
quadrature = sin(n+pi/4);
symbol_map = (in_phase + quadrature*1i);

qpsk_modulatedData = (symbol_map(dataInput+1))';

%Series/Paralel
padding = mod(n_subcarriers-mod(length(qpsk_modulatedData),n_subcarriers),n_subcarriers);
X_padded = [qpsk_modulatedData;zeros(padding,1)];
numberOfSymbols = length(X_padded)/n_subcarriers;

OFDM_blocks = reshape(X_padded, n_subcarriers, numberOfSymbols);

% Sent Subsymbols
sent_subsymbols = zeros(ofdm_M, subsymbol_N, numberOfSymbols);
for m = 1:1:ofdm_M
    sent_subsymbols(m,:,:) = OFDM_blocks((m-1)*subsymbol_N+1:m*subsymbol_N,:);
end

% Received Subsymbols
tau0 = (2*d)/c;

received_subsymbols = zeros(); % IFFT data
for m = 1:1:ofdm_M
    for n = 1:1:subsymbol_N
        for l = 1:1:numberOfSymbols
            fm = (fc+(m-1)*(subsymbol_N)*delta_f);
            received_subsymbols(m,n,l) = sent_subsymbols(m,n,l)...
                *exp(-1j*2*pi*tau0*(fm+(n-1)*delta_f))...
                *exp(1j*2*pi*((2*v*((m-1)+(l-1)*ofdm_M)*T)/c)*(fm+(n-1)*delta_f));
        end
    end
end


%% Add Reception Noise
power = sum(abs(sent_subsymbols(:)).^2)/length(sent_subsymbols(:))
sigma = sqrt((power/(10^(SNR/10)/2))) % sqrt(n0/2)
noise = (randn(size(sent_subsymbols)))*sigma; % Additive noise
received_subsymbols = received_subsymbols + noise;

%% Radar Signal Analysis
D = received_subsymbols./sent_subsymbols;

% Recovered OFDM environment
OFDM_stepped_blocks = zeros(subsymbol_N*ofdm_M, numberOfSymbols);
for m=1:1:ofdm_M
    for n=1:1:subsymbol_N
        for l=1:1:numberOfSymbols
            OFDM_stepped_blocks((n)+(m-1)*subsymbol_N,(m)+(l-1)*ofdm_M) = D(m,n,l);
        end
    end
end

% DFT of length M*B is applied in direction of time
% to obtain the velocity profile
velocityProfile = fft(OFDM_stepped_blocks, ofdm_M*ofdm_L, 2);

% An IDFT is applied on the velocity profile in frequency direction 
% to obtain the desired two-dimensional high resolution range-velocity profile.
rangeVelocityProfile = ifft(velocityProfile, ofdm_M*subsymbol_N, 1);

range = linspace(0, c/(2*delta_f), subsymbol_N*ofdm_M);
velocity = linspace(0, c/(2*fc*T), ofdm_M*ofdm_L);

% Find peak (get estimate)
peak = 0;
estimative_range = 0;
estimative_velocity = 0;
rangeVelocityProfile_length = size(rangeVelocityProfile);
for i = 1:1:rangeVelocityProfile_length(1)
    for j =1:1:rangeVelocityProfile_length(2)
        value = abs(rangeVelocityProfile(i,j));
        if value > peak
            peak = value;
            estimative_range = i;
            estimative_velocity = j;
        end
    end
end

disp(' ');
disp('Range estimative: ');
disp([mat2str(range(estimative_range), 6), ' m']);
disp(' ');
disp('Velocity estimative: ');
disp([mat2str(velocity(estimative_velocity), 6), ' m/s']);

%% Plotting

figure(1),
surf(velocity, range, abs(rangeVelocityProfile))
xlabel('Velocity $m/s$', 'interpreter', 'latex')
ylabel('Range, $m$', 'interpreter', 'latex')
title('Range-velocity profile', 'interpreter', 'latex')
%saveas(figure(1), "4_rangeVelocityProfile.png")

% Eof