clc;
clear all;
close all;

% Carrier center frequency
center_freq = 433.5e6; % 433.5MHz

% Amplitude of carrier signal
A=5;                                          

% Frequency deviation
f_deviation = 48e3; % 48kHz

% Sample rate
sample_rate = 100;

% Frequency deviations
f1 = center_freq + f_deviation; % Frequency for binary '1'
f2 = center_freq - f_deviation; % Frequency for binary '0'

% Binary information
%bit_data=[ 1 0 0 1 1 0 1];

% Binary information (for now, currently random 1s and 0s but will be
% updated in the future to represent packet data

n = 8; % Number of bits to generate.
bit_data = randi([0, 1], [1, n]);

% Bit rate
bit_rate = 9600;

% Bit period
bit_period = 1 / bit_rate;                                          

% Calculating transmission as a binary signal
digital_signal = [];

for n = 1 : 1 : length(bit_data)
    if bit_data(n) == 1
       signal_element = ones(1,sample_rate);
    else
        signal_element = zeros(1,sample_rate);
    end
     digital_signal = [digital_signal signal_element];
end

% Plot digital signal
t1 = bit_period / sample_rate : bit_period / sample_rate : sample_rate * length(bit_data) * (bit_period/sample_rate);

subplot(3, 1, 1);
plot(t1, digital_signal, 'lineWidth', 2.5); grid on;
axis([0 bit_period * length(bit_data) -0.5 1.5]);
ylabel('Amplitude(V)');
xlabel('Time(s)');
title('Initial digital signal');

% Binary-FSK modulation
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

% Plot Binary-FSK modulation
t3= bit_period / 99 : bit_period / 99 : bit_period * length(bit_data);

subplot(3, 1, 2);
plot(t3, modulated_signal);
ylabel('Amplitude(V)');
xlabel('Time(s)');
title('Binary-FSK modulated signal coresponding to initial digital signal');


%Binary-FSK demodulation
demodulated_data = [];
for n = ss : ss : length(modulated_signal)
  t = bit_period / 99 : bit_period / 99 : bit_period;

  % Calculating carrier signals for 1s and 0s based on modulated frequency
  y1 = cos(2 * pi * f1 * t);                   
  y2 = cos(2 * pi * f2 * t);   

  product_y1 = y1.*modulated_signal((n - (ss - 1)) : n);
  product_y2 = y2.*modulated_signal((n - (ss - 1)) : n);

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

% Plot binary information as digital signal received after demodulation
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

subplot(3, 1, 3)
plot(t4, digital_signal, 'LineWidth', 2.5); grid on;
axis([0 bit_period * length(demodulated_data) -0.5 1.5]);
ylabel('Amplitude(V)');
xlabel('Time(s)');
title('Interpreted digital signal received after binary FSK demodulation');
