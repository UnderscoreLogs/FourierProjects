clc, clearvars, close all
%-------------------------Settings--------------------------------------
startTime = 0;
endTime = 2;
Fs = 1000; %Sample rate
window_size = 2^7;

%------------------------------Variables-------------------------------
N = (endTime - startTime) * Fs;
t = (0:N-1)/Fs;

df = Fs / window_size
dt = window_size / Fs

n = (0:window_size-1)';
W = 0.5*(1 - cos(2*pi*n/(window_size-1))); %Hann window

%--------------------------Data Option 1----------------------------------
% Creating a function that changes as time goes on
freq1 = 120;
freq2 = 50;
freq3 = 220;


% Piecewise function using "logical indexing"
clean_function = zeros(size(t));

clean_function(t < 3) = sin(2*pi*freq1*t(t < 3));

clean_function(t >= 3 & t < 6) = sin(2*pi*freq2*t(t >= 3 & t < 6));

clean_function(t >= 6 & t < 7) = sin(2*pi*freq2*t(t >= 6 & t < 7)) + ...
                                 sin(2*pi*freq3*t(t >= 6 & t < 7));

clean_function(t >= 7) = sin(2*pi*freq3*t(t >= 7));

%--------------------------------Noise------------------------------------
% Add noise to the clean function
noise = 0.5 * randn(size(clean_function));

%---------------------------Chirp Data Option 2--------------------------
f0 = 50;
f1 = 250;
T = endTime;

chirp = 2*pi*( f0*t + (f1 - f0)/(3*T^2) * t.^3 );
raw_data = sin(chirp) + noise;

%--------------Plot Option 1-------------------------------------------
% Plot original data
%figure(1)
%plot(t,raw_data,'r')
%hold on
%plot(t,clean_function,'b')
%title('Raw Data Vs. Clean Data')
%xlabel('Time')
%ylabel('Amplitude')
%legend('Raw Data', 'Clean Data')

%-----------------------Plot Option 2-------------------------------
sound(2*raw_data,Fs) %Plays the sound of the raw_data

figure(1)
plot(t,raw_data,'r')
xlabel('Time')
ylabel('Amplitude')
title('Raw Data')

%-----------------------STFT------------------------------------------
% Adding zeros so that a whole number of windows fits within the data
padSize = window_size - rem(length(raw_data), window_size);
raw_data_padded = raw_data;
if padSize ~= window_size
    raw_data_padded = [raw_data, zeros(1, padSize)];
end

windowSampleStep = round(window_size * 0.5); % Must be a whole number
numWindows = floor((length(raw_data_padded) - window_size)/windowSampleStep) + 1;

S = zeros(window_size/2 + 1, numWindows); %Spectrogram data

for k = 1:numWindows
    
    % The start/end index position of each window
    startIndex = (k-1)*windowSampleStep + 1;
    endIndex = startIndex + window_size - 1;    

    % Hann window to the current segment of raw data
    windowData = raw_data_padded(startIndex:endIndex); 
    windowData = windowData(:);
    windowedData = windowData .* W;
    
    % Calculating the FFT of the window
    fftData = fft(windowedData); 

    % Storing magnitude spectrum
    S(:,k) = abs(fftData(1:window_size/2 + 1)); 
end

%--------------------------Visuallizing------------------------------
% Making spectrogram

freqAxis = (0:window_size/2) * Fs / window_size;

timeAxis = (0:numWindows-1) * windowSampleStep  / Fs;

figure(2)
imagesc(timeAxis,freqAxis, S)
axis xy
ylim([0, 300])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Spectrogram')
colorbar

%--------------------Comparing to Regular FFT------------------------
% Compare to regular FFT

fhat_reg = fft(raw_data);

frequency_axis = (0:N/2)*(Fs/N);

PSD = abs(fhat_reg/N).^2;
PSD = PSD(1:N/2 + 1);

figure(3)
plot(frequency_axis,PSD);
title('PSD')
xlabel('Frequency')
ylabel('Power')

%---------------------Comparing to built in functions-----------------
% Comparing to the built in function
figure(4)

% Parameters matching your implementation AI STUFF
window = W;                          % your Hann window
noverlap = window_size - windowSampleStep;  % overlap
nfft = window_size;           % FFT size

spectrogram(raw_data, 128, 64, 128,1000, 'yaxis')
title('Built-in Spectrogram (MATLAB)')
