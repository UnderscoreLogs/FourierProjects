clc, clearvars

data = h5read('C:\Users\smegl\OneDrive\Desktop\Code\Matlab\MRI\MRI_K_Space_Raw_Datasets\NYUFastMRI\knee_singlecoil_test\knee_singlecoil_test','/kspace');

h5disp('C:\Users\smegl\OneDrive\Desktop\Code\Matlab\MRI\MRI_K_Space_Raw_Datasets\NYUFastMRI\knee_singlecoil_test\knee_singlecoil_test');

[rows, cols, slices] = size(data.r);

image = zeros(rows, cols);

kspace = double(data.r) + 1i * double(data.i);



for slice = 40
    k = kspace(:,:,slice);
    image = image + abs(fftshift(ifft2(ifftshift(k))));

    figure(1)
    imshow(abs(fftshift(ifft2(ifftshift(k)))),[])
    title(['slice = ', num2str(slice)])
    pause(0.1)
end


