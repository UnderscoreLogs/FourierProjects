clc, clearvars, close all


f = imread('C:\Users\smegl\OneDrive\Desktop\Code\Matlab\Imaging\standard_test_images\lena_gray_512.tif'); %zeros(size,size) + noise;

[rows, cols] = size(f);


%noise = 0.5 * randn(rows,cols);
%f = f + noise;
imshow(f, [])

fhat_image = fft2(f);

figure (2)
absF = abs(fhat_image);
logabsF = log(absF);
imshow(logabsF,[-1 5]);
colormap(jet)
colorbar

%--------------------------------Low Pass Filter----------------------
fhat_filtered_shifted = fftshift(fhat_image); % reordering the fft

[x,y] = meshgrid(-cols/2:cols/2-1, -rows/2:rows/2-1); %creating a grid of all frequencies -N/2 to N/2-1
radius = sqrt(x.^2 + y.^2); % comparison radius
cutoff = 80; %Frequency cutoff

mask = radius < cutoff; %Mask, at what points do you cut

fhat_filtered_shifted = fhat_filtered_shifted .* mask; %Applying mask

%-------------------- Frequency domain image-----------------------
figure (3)
absF_denoised = abs(fhat_filtered_shifted);
logabsF_denoised = log(absF_denoised);
imshow(logabsF_denoised,[-1 5]);
colormap(jet)
colorbar
%------------------------- Inversing back to normal domain--------
f_denoised = ifft2(ifftshift(fhat_filtered_shifted));

%-------------------------- Denoised image---------------------
figure (4)
imshow(f_denoised,[])