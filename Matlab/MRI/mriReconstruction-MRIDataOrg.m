clc
clearvars

filename = 'C:\Users\smegl\OneDrive\Desktop\Code\Matlab\MRI\MRI_K_Space_Raw_Datasets\Hand.h5';

h5disp(filename);

data = h5read(filename, '/dataset/data');

% Number of acquisitions
numAcq = length(data.data);

% Samples in one acquisition
samples = length(data.data{1});

% Complex samples
Nx = samples/2;

% Find maximum ky index
ky_all = data.head.idx.kspace_encode_step_1;
Ny = max(ky_all) + 1;

% Allocate k-space
kspace = zeros(Ny, Nx);

for k = 1:numAcq

    line = data.data{k};

    % Convert interleaved re/im -> complex
    realPart = line(1:2:end);
    imagPart = line(2:2:end);

    complexLine = realPart + 1i*imagPart;

    % Get ky location
    ky = data.head.idx.kspace_encode_step_1(k);

    % Store in correct row
    kspace(ky+1,:) = complexLine;

end

% Reconstruct image
image = fftshift(ifft2(ifftshift(kspace)));

figure
imshow(abs(image),[])