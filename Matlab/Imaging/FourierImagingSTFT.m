clc, clearvars, close all
%----------------------Settings------------------------------------
win_height = 2^6
win_width = 2^6

overlap = 0.5;

% Step size between windows
windowSampleStep_Y = round(win_height * overlap);
windowSampleStep_X = round(win_width  * overlap);

%-----------------------------------Load Image----------------------------------
f = imread('C:\Users\smegl\OneDrive\Desktop\Code\Matlab\Imaging\standard_test_images\NoiseTest_Test_Image_Crop.jpg'); %zeros(size,size) + noise;
f = rgb2gray(f); %convert to gray scale
f = im2double(f); %convert to double

[rows, cols] = size(f);

% Optional noise
% noise = 0.2 * randn(rows, cols);
% f = f + noise;

figure(1)
imshow(f, [])

%-----------------------2D Window-----------------------------
[X, Y] = meshgrid(1:win_width, 1:win_height);

W = (sin(pi*X/win_width)).^2 .* (sin(pi*Y/win_height)).^2;

%--------------------Padding Data------------------------------------------
% Adding zeros so that a whole number of windows fits within the data
padSize_Y = mod(-rows, windowSampleStep_Y);
padSize_X = mod(-cols, windowSampleStep_X);

raw_data_padded = padarray(f, [padSize_Y padSize_X], 0, 'post');
raw_data_padded = padarray(raw_data_padded, [windowSampleStep_Y windowSampleStep_X], 0, 'pre');

[padded_rows, padded_cols] = size(raw_data_padded);

%-------------------Number of Windows-----------------
numWindows_X = floor((padded_cols - win_width)/windowSampleStep_X) + 1;

numWindows_Y = floor((padded_rows - win_height)/windowSampleStep_Y) + 1;

%--------------------Outputs Arrays----------------------
f_denoised = zeros(padded_rows, padded_cols);
window_weight = zeros(padded_rows,padded_cols);

%------------------Mask for Denoising-----------------------
[x,y] = meshgrid(-win_width/2:win_width/2-1, -win_height/2:win_height/2-1); %creating a grid of all frequencies -N/2 to N/2-1

radius = sqrt(x.^2 + y.^2); % comparison radius

cutoff = 20; %Frequency cutoff
    
mask = radius < cutoff; %Mask, at what points do you cut

%-----------------------2D STFT--------------------------------------
for j = 1:numWindows_Y
    
    % The start/end index position of each window
    startY = (j-1)*windowSampleStep_Y + 1;
    endY   = startY + win_height - 1;    

    for k = 1:numWindows_X
    
        % The start/end index position of each window
        startX = (k-1)*windowSampleStep_X + 1;
        endX   = startX + win_width - 1;  
    
        %--------------------Window of Data--------------------
        windowData = raw_data_padded(startY:endY, startX:endX); 
        
        %-------------Apply Window---------------------------
        windowedData = windowData .* W;

        %------------------------ FFT -----------------------------
       
        % Calculating the FFT of the window
        fhat_data = fft2(windowedData); 

        % Shift FFT to center
        fhat_shifted = fftshift(fhat_data);

        %---------------- Low Pass Filter -----------------------------
        fhat_filtered = fhat_shifted .* mask;
        
        %-----------------Inverse FFT--------------------------------
        f_window_denoised = real(ifft2(ifftshift(fhat_filtered)));
   
        %------------------------- Reconstruction ---------------------    
        f_denoised(startY:endY, startX:endX) = f_denoised(startY:endY, startX:endX) + f_window_denoised .* W;

        % Update the window weight for the current segment
        window_weight(startY:endY, startX:endX) = window_weight(startY:endY, startX:endX) + W.^2;
    end
end

%---------------------- Normalize ---------------------------------------
window_weight(window_weight < 1e-6) = 1;

f_denoised = f_denoised ./ window_weight;

%---------------------- Trim Back To Original Size ----------------------
f_denoised = f_denoised(windowSampleStep_Y:rows+windowSampleStep_Y-1, windowSampleStep_X:cols+windowSampleStep_X-1);

%---------------------- Display Image ----------------------------------
figure(2)
imshow(f_denoised, []);
title('Denoised Image');
