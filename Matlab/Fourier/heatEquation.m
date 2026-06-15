clc, clearvars, close all

alpha = 0.5; % Thermal diffusivity constant
L = 100; % Length 
N = 1000; % number of samples
dx = L/N;
x = -L/2:dx:L/2-dx; % x domain

%Intital Condition
u0 = 0*x;
u0((L/2 - L/10)/dx:(L/2+L/10)/dx) = 1;
u0hat = fft(u0);

kappa = (2*pi/L) * [0:N/2-1 -N/2:-1];

% Time parameters
dt = 50; % Time step
T = 600; % Total time
t = 0:dt:T; % t domain

uhat = u0hat.' .* exp(-alpha .* kappa.' .^2 .* t);

u = real(ifft(uhat)); %heat equation

figure (2)
imagesc(t, x, u);
xlabel('Time'); 
ylabel('x');
title('Heat Equation Solution');
colorbar;

figure (1)
h=waterfall(x,t(1:1:end),u(:,1:1:end).');
set(h,'LineWidth',5,'FaceAlpha',.5);
colormap(flipud(jet)/1.5)
xlabel('Space');
ylabel('Time');
zlabel('Temperature');
