%   In automotive communications, assume a bandwidth of BW = 600 MHz, a 
% carrier of f c = 80 GHz, a pulse duration of T = 40μs, and a range of 
% detection of d = 300 m. Assume the target is moving at v = 150 km/h.
% Consider that we need to cover a range swath of 10m, so that you can
% choose L and M = 4.
%   Divide the bandwidth among three vehicles each using its own OFDM 
% waveform to estimate their relative distances and Doppler shifts.
% Simulate the situation and discuss the results, all with an SNR of 10 dBs.

clc
clear
close all


%% Global Properties
c = 299792458; % m/s

%% Radio Resources
BW = 600*1e6; % 600 MHz
fc = 80*1e9; % 80 GHz
T_symbol = 40*1e-6; % 40 μs

%% OFDM Properties
user_BW = BW/3;

ofdm_N = 64; % Number of subcarriers per user

delta_f = user_BW/ofdm_N;

% Relate T = 40micro s to the number of transmited data symbols
fft_samplePeriod = (1/(2*user_BW)); % time between samples
fft_numberOfSamplesInASymbol = round(T_symbol/fft_samplePeriod);
fft_numberOfTransmittedDataPoints = floor(fft_numberOfSamplesInASymbol/2);
fft_numberOfCyclePrefixSamples = fft_numberOfSamplesInASymbol - fft_numberOfTransmittedDataPoints;

% Setup OFDM Stepped Frequency
ofdm_M = 4; % Number of OFDM frames

block_L = fft_numberOfTransmittedDataPoints/ofdm_N; % number of symbols per frame
block_N = round(ofdm_N/ofdm_M); % Number of carriers per frame

block_BW = user_BW/ofdm_M;
block_samplingRate = 2*block_BW;

% Symbols duration
Tdata = fft_numberOfTransmittedDataPoints*fft_samplePeriod;
Tcp = fft_numberOfCyclePrefixSamples*fft_samplePeriod; % Half the channel width (OFDM recover)

fft_sizeOfCyclePrefix = 0; %ofdm_N*(Tcp/Tdata); % Number of symbols allocated to cyclic prefix

%% Radar Properties
d = 300; % range of detection = 300 m
range_swath = 10; % 10 m
minimumRangeResolution = c/(2*user_BW); % delta_D

%% Transmission Properties
SNR = 50; % 10 dBs

%% Target Properties
numberOfVehicles = 3;
v = 150*(1/3.6); % 150 km/h = 41.67 m/s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sampling chirp function samples
qpsk_dataInput = randi(3, fft_numberOfTransmittedDataPoints, 1);

%% QPKS Modulation
qpsk_Modulator = comm.QPSKModulator;

qpsk_modulatedData = qpsk_Modulator(qpsk_dataInput);

%% OFDM Transmission
n_subcarriers = ofdm_N;

fft_padding = mod(n_subcarriers-mod(length(qpsk_modulatedData),n_subcarriers),n_subcarriers);
X_padded = [qpsk_modulatedData;zeros(fft_padding,1)];
padding = length(X_padded)/n_subcarriers;
X_blocks = reshape(X_padded, n_subcarriers, padding);
x = ifft(X_blocks);

%% Build transmitted data tensor (OFDM frame)
data = zeros(ofdm_M, block_N, block_L); % IFFT data
for m = 1:1:ofdm_M
    for n = 1:1:block_N
        for l = 1:1:block_L
            data(m,:,l) = x((m-1)*block_N+1:m*block_N,l);
        end
    end
end

%% Build received data tensor (OFDM frame)
T = Tdata;%+Tcp;
tau0 = (2*d)/c;

dataReceived = zeros(ofdm_M, block_N, block_L); % IFFT data
for m = 1:1:ofdm_M
    for n = 1:1:block_N
        for l = 1:1:block_L
            fm = (fc+(m-1)*(block_N)*delta_f);
            dataReceived(m,n,l) = data(m,n,l)*exp(-1j*2*pi*tau0*(fm+(n-1)*delta_f))...
                *exp(1j*2*pi*((2*v*((m-1)+(l-1)*ofdm_M)*T)/c)*(fm+(n-1)*delta_f));
        end
    end
end

%% Add Reception Noise
noise = awgn(zeros(size(dataReceived)), SNR);
dataReceived = dataReceived+noise;

%% Radar Signal Analysis

D = dataReceived./data;

% Recovered OFDM environment (freq/time) (64 carriers)
OFDM_matrix = zeros(block_N*ofdm_M, block_L*ofdm_M);
for m=1:1:ofdm_M
    for n=1:1:block_N
        for l=1:1:block_L
            OFDM_matrix((n)+(m-1)*block_N,(m)+(l-1)*ofdm_M) = D(m,n,l);
        end
    end
end

velocityProfile = fft(OFDM_matrix, ofdm_M*block_L, 1);

rangeVelocityProfile = ifft(velocityProfile, ofdm_M*block_N, 2);

% Find peak (get estimate)
peak = 0;
estimative_delay = 0;
estimative_frequency = 0;
for i = 1:1:ofdm_M*block_L
    for j =1:1:ofdm_M*block_N
        value = abs(rangeVelocityProfile(i,j));
        if value > peak
            peak = value;
            estimative_delay = i;
            estimative_frequency = j;
        end
    end
end

disp('Delay estimative: ');
disp(estimative_delay);
disp(' ');
disp('Doppler shift estimative: ');
disp(estimative_frequency);

%% Plotting

figure(1),
x_plot = 1:1:(ofdm_M*block_N);
y_plot = 1:1:(block_L*ofdm_M);
surf(x_plot, y_plot, abs(rangeVelocityProfile))
xlabel('Frequency $Hz$', 'interpreter', 'latex')
ylabel('Time, $s$', 'interpreter', 'latex')
zlabel('Normalized dechirped signal', 'interpreter', 'latex')
title('Range-velocity profile', 'interpreter', 'latex')
%saveas(figure(1), "4_rangeVelocityProfile.png")

%% Functions
function retval = u(t)
  if t < 0
      retval = 0;
  else
      retval = 1;
  end
end

% Eof