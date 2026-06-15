clc, clearvars, close all

image = imread('C:\Users\smegl\OneDrive\Desktop\Code\Matlab\Imaging\standard_test_images\lena_gray_512.tif');

[rows, cols] = size(image);

rho = double(image);
rho = rho/max(rho(:));

%----------------settings------------------------------------

gamma = 2*pi*42.577e6; % Gyromagnetic ratio in rad/s/T

B0 =3; % Tesla base magnetic field

Gx = 30e-3; % x gradient in T/m
Gy = 25e-3; % y gradient in T/m

xRes = 50; % Resolutions
yRes = 50;

FOVx = 0.5; % meter, Length of the FOV
FOVy = 0.5; % meter, Height of the FOV

%------------------- Calculations -----------------------------------

%Spatial coordinates
x = linspace(-FOVx/2, FOVx/2, cols);
y = linspace(-FOVy/2, FOVy/2, rows);

Fs = (gamma / (2*pi)) * Gx * FOVx %Sampling rate has to be > band width of RF Pulse
pulseDuration = (2*pi) / (gamma * Gy * FOVy) %2pi phase wrap range

T = xRes/Fs % Tread time, Total readout diration 

t = linspace(-T/2, T/2 - T/xRes, xRes);

%-------------------Simulation--------------------------

kspace = zeros(yRes,xRes);

% x gradient
B_Gx = Gx * x;

% Magnetic Field
B = zeros(rows, cols);
B = B0 + repmat(B_Gx, rows, 1);

% Creating a phase step vector
phaseSteps = (-yRes/2):(yRes/2 - 1);

% Quadrature Detection  
TargetLarmor = gamma * B0;
RefCos= cos(TargetLarmor*t);
RefSin= -sin(TargetLarmor*t);

for idx = 1:yRes

    %Phase Encoding Pulse
    phaseStep = phaseSteps(idx);
    B_Gy = Gy * phaseStep * y;
    

    %Shielded 
    sigma = 10e-6; %shielding cosntant NOT USED CURRENTLY
    
    B_eff = B; %effective magnetic field with shielding from electrons
    
    %--------------------------------------------------------
    
    %Larmor  Frequencies
    larmorFreqs = gamma * B_eff; % Larmor frequencies
        
    %Voltage Signal (AI optimized part this because it was taking 10 years 
    % to calculate with a nested loop)
    B_Gy_mat = repmat(B_Gy(:), 1, cols);
    
    rho_vec = rho(:);
    larmor_vec = larmorFreqs(:);
    
    B_Gy_vec = gamma .* B_Gy_mat(:) .* pulseDuration;
    
    t = reshape(t, 1, []);  % row vector
    
    phi = larmor_vec * t + B_Gy_vec;   % (pixels × time)
    
    S = sum(rho_vec .* exp(-1i * phi), 1);
        
    %------------------Signal Processing---------------------------
    
    % Compute the correlation between the voltage signal and reference signals
    correlationSignal = (S .* RefCos) + 1i * (S .* RefSin);
    
    %-----------K-Space----------------------------------------------
        
    kspace(idx,:) = correlationSignal;
end 

%%%-----------------IFFT---------------------------

img_rec = fftshift(ifft2(ifftshift(kspace)));

figure;
imshow(abs(img_rec), []);
title(['Reconstructed Image ',num2str(xRes),'x',num2str(yRes)]);