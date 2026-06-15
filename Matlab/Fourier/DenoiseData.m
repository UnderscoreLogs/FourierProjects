clc, clearvars, close all

startTime = 0;
endTime = 1;
Fs = 1080; %Sample rate
N = (endTime - startTime) * Fs;
freqRes = Fs/N;

freq1 = 120;
freq2 = 50;

t = (0:N-1)/Fs;

clean_function = sin(2 * pi * freq1 * t) + sin(2 * pi * freq2 * t);

% Add noise to the clean function
noise = 0.5 * randn(size(clean_function));
noisy_function = clean_function + noise;

figure(1)
plot(t,noisy_function,'r')
hold on
plot(t,clean_function,'b')
title('Raw Data Vs. Clean Data')
xlabel('Time')
ylabel('Amplitude')
legend('Raw Data', 'Clean Data')

fhat = fft(noisy_function);

frequency_axis = (0:N/2)*(freqRes);

PSD = abs(fhat/N).^2;
PSD = PSD(1:N/2 + 1);

figure(2)
plot(frequency_axis,PSD);
title('PSD')
xlabel('Frequency')
ylabel('Power')

fhat_filtered = fhat;

threshold = 0.1*max(abs(fhat));

fhat_filtered = fhat .* (abs(fhat) > threshold);

denoised_function = ifft(fhat_filtered);

figure(3)
plot(t, real(denoised_function), 'g')
hold on
plot(t, clean_function, 'r')
title('Denoised Signal')
xlabel('Time')
ylabel('Amplitude')
legend('Denoised Signal', 'Expected Signal')


