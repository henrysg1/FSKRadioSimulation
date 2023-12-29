%==========================================================================
%|      Project:  University of Nottingham CubeSat Radio Simulation
%|
%|       Author:  Henry Sands-Grant
%|     Language:  MATLAB
%|   To Compile:  -
%|
%|       To Run:  Open MATLAB script in MATLAB and click run button
%|
%+-----------------------------------------------------------------------------
%|
%|  Description:  Simulation of FSK radio modulation and demodulation of
%|                the AstroDev Lithium-2 radio used in a CubeSat
%|                application
%|
%|      License:  GNU General Public License v3.0
%|
%==========================================================================

clc;
clear all;
close all;

% ====== Configurations of radio ======

% Carrier center frequency
center_freq = 433.5e6; % 433.5MHz
% Amplitude of carrier signal
A = 5;
% Amplitude of intefering signal
Ai = 5;
% Frequency deviation
f_deviation = 48e3; % 48kHz
% Frequency deviations
f1 = center_freq + f_deviation; % Frequency for binary '1'
f2 = center_freq - f_deviation; % Frequency for binary '0'
% Binary information (for now, currently random 1s and 0s but will be
% updated in the future to represent packet data)
n = 8; % Number of bits to generate.
bit_data = randi([0, 1], [1, n]);
% Bit rate
bit_rate = 9600;
% Bit period
bit_period = 1 / bit_rate;         
% Sample rate
sample_rate = 4 * bit_rate;


% ====== Creating initial digital signal ======


digital_signal = [];

for n = 1 : 1 : length(bit_data)
    if bit_data(n) == 1
       signal_element = ones(1,sample_rate);
    else
        signal_element = zeros(1,sample_rate);
    end
     digital_signal = [digital_signal signal_element];
end


% ====== Plotting initial digital signal ======


t1 = bit_period / sample_rate : bit_period / sample_rate : sample_rate * length(bit_data) * (bit_period/sample_rate);

subplot(5, 1, 1);
plot(t1, digital_signal, 'lineWidth', 2.5); grid on;
axis([0 bit_period * length(bit_data) -0.5 1.5]);
ylabel('Amplitude(V)');
xlabel('Time(s)');
title('Initial digital signal');


% ====== Creating noise-free modulated signal ======


t2 = bit_period / 99 : bit_period /99 : bit_period;      

ss = length(t2);

modulated_signal = [];
for (i = 1 : 1 : length(bit_data))
    if (bit_data(i) == 1)
        y = A * cos(2 * pi * f1 * t2);
    else
        y = A * cos(2 * pi * f2 * t2);
    end
    modulated_signal = [modulated_signal y];
end


% ====== Plotting noise-free modulated signal ======


t3= bit_period / 99 : bit_period / 99 : bit_period * length(bit_data);

subplot(5, 1, 2);
plot(t3, modulated_signal);
ylabel('Amplitude(V)');
xlabel('Time(s)');
title('Binary-FSK modulated signal coresponding to initial digital signal');


% ====== Adding noise and interference to modulated signal ======

% AWGN
SNR = 30; % Signal to Noise Ratio
noisy_signal = awgn(modulated_signal, SNR, 'measured');

% Simulate Doppler Shift
doppler_freq_shift = 100; % Doppler frequency shift in Hz
t_doppler = 0 : 1 / sample_rate : length(noisy_signal) / sample_rate - 1 / sample_rate;
doppler_signal = cos(2 * pi * doppler_freq_shift * t_doppler);
modulated_signal_with_doppler = noisy_signal .* doppler_signal;

% Simulate Interference
interference_freq = 434e6; % Interfering signal frequency
interference = Ai * cos(2 * pi * interference_freq * t_doppler);
modulated_signal_with_interference = modulated_signal_with_doppler + interference;


% ====== Plotting noisy modulated signal ======


subplot(5, 1, 3);
plot(t3, modulated_signal_with_interference);
ylabel('Amplitude(V)');
xlabel('Time(s)');
title('Binary-FSK modulated signal including interference');


% ====== Designing FIR Filter ======


% Filter specifications
order = 60;                      % Order of the filter
cutoff_freq = 0.2;               % Normalized cutoff frequency (0.2 corresponds to 20% of the Nyquist rate)
filter_type = 'low';             % Type of filter: 'low' for lowpass

% Designing the FIR filter using the windowing method with a Hamming window
b = fir1(order, cutoff_freq, filter_type, hamming(order+1));


% ====== Applying FIR Filter to the received signal ======


% Apply the FIR filter to the noisy signal
filtered_signal = filter(b, 1, modulated_signal_with_interference);


% ====== Plotting filtered modulated signal ======

subplot(5, 1, 4); 
plot(t3, filtered_signal);
ylabel('Amplitude(V)');
xlabel('Time(s)');
title('Filtered Binary-FSK modulated signal');


% ====== Demodulating Binary-FSK signal ======


demodulated_data = [];
for n = ss : ss : length(modulated_signal_with_interference)
  t = bit_period / 99 : bit_period / 99 : bit_period;

  % Calculating carrier signals for 1s and 0s based on modulated frequency
  y1 = cos(2 * pi * f1 * t);                   
  y2 = cos(2 * pi * f2 * t);   

  product_y1 = y1.*modulated_signal_with_interference((n - (ss - 1)) : n);
  product_y2 = y2.*modulated_signal_with_interference((n - (ss - 1)) : n);

  t4 = bit_period / 99 : bit_period / 99 : bit_period;

  z1 = trapz(t4, product_y1);                     
  z2 = trapz(t4, product_y2);

  integrated_product_y1 = round(2 * z1 / bit_period);
  integrated_product_y2 = round(2 * z2 / bit_period);

  % Determining logic level based on signal being greater or less than A/2
  if(integrated_product_y1 > A/2)      
    decoded_bit = 1;
  else
    decoded_bit = 0;
  end
  demodulated_data = [demodulated_data decoded_bit];
end


% ====== Plotting interpreted digital signal after demodulation ======


digital_signal = [];
for n = 1 : length(demodulated_data)
    if demodulated_data(n) == 1
       signal_element = ones(1, sample_rate);
    else
       signal_element = zeros(1, sample_rate);
    end
     digital_signal = [digital_signal signal_element];
end

t4 = bit_period / sample_rate : bit_period / sample_rate : sample_rate * length(demodulated_data) * (bit_period / sample_rate);

subplot(5, 1, 5)
plot(t4, digital_signal, 'LineWidth', 2.5); grid on;
axis([0 bit_period * length(demodulated_data) -0.5 1.5]);
ylabel('Amplitude(V)');
xlabel('Time(s)');
title('Interpreted digital signal received after binary FSK demodulation');


% ====== Calculating Bit Error Rate (BER) ======


bit_errors = sum(bit_data ~= demodulated_data);
BER = bit_errors/length(bit_data);
fprintf('Bit Error Rate is %f\n', BER);
