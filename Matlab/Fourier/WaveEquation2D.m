clc, clearvars, close all

v = 1; % 
L = 100; % Length
W = 100; %Width
N = 100; % number of samples

dx = L/N;
x = -L/2:dx:L/2-dx; % x domain
dy = W/N;
y = -W/2:dy:W/2-dy; % y domain

[X, Y] = meshgrid(x,y);

%Intital Condition
f = zeros(N,N); %u(x,0) %Gaussian Pulse: exp((-(X).^2 - (Y).^2)/3) SIne Wave: sin(4pix/L)
%f((L/2 - L/10)/dx:(L/2+L/10)/dx) = 1;
fhat = fft2(f);

g = exp((-(X).^2 - (Y).^2)/3); %u'(x,0)
ghat = fft2(g);

kappa_x = (2*pi/L) * [0:N/2-1 -N/2:-1];
kappa_y = (2*pi/W) * [0:N/2-1 -N/2:-1];

[KX, KY] = meshgrid(kappa_x, kappa_y);

% Time parameters
dt = 1; % Time step
T = 300; % Total time
t = (0:dt:T).'; % t domain

u = zeros(N, N, length(t));
omega = sqrt((v.^2 .* KX.^2) + (v.^2 .* KY.^2));
omega(omega == 0) = 1; %prevents division by zero
for i = 1:length(t) 
    uhat = (fhat .* cos(omega .* t(i))) + ((ghat ./ (omega)) .* sin(omega .* t(i)));
    uhat(1,1) = fhat(1,1) + ghat(1,1) * t(i);
    u(:,:,i) = real(ifft2(uhat)); %wave equation
end

c_min = min(u, [], 'all') - 0.5;
c_max = max(u, [], 'all') + 0.5;

% figure(1)
%     surf(X,Y,u(:,:,2));
%     shading interp;
%     colormap jet;
%     colorbar;
%     clim([c_min, c_max]);
%     ylabel("Y");
%     xlabel("X");
%     zlabel("Amplitude");
%     xlim([min(x),max(x)]);
%     ylim([min(y), max(y)]);
%     zlim([c_min,c_max]);
%     title(['t=', num2str(t(2))]);
% figure(2)
%     surf(X,Y,u(:,:,41));
%     shading interp;
%     colormap jet;
%     colorbar;
%     clim([c_min, c_max]);
%     ylabel("Y");
%     xlabel("X");
%     zlabel("Amplitude");
%     xlim([min(x),max(x)]);
%     ylim([min(y), max(y)]);
%     zlim([c_min,c_max]);
%     title(['t=', num2str(t(41))]);
% 
% figure(3)
%     surf(X,Y,u(:,:,101));
%     shading interp;
%     colormap jet;
%     colorbar;
%     clim([c_min, c_max]);
%     ylabel("Y");
%     xlabel("X");
%     zlabel("Amplitude");
%     xlim([min(x),max(x)]);
%     ylim([min(y), max(y)]);
%     zlim([c_min,c_max]);
%     title(['t=', num2str(t(101))]);

figure (2)

for i = 1:length(t)
    surf(X,Y,u(:,:,i));
    shading interp;
    colormap jet;
    colorbar;
    clim([c_min, c_max]);
    ylabel("Y");
    xlabel("X");
    zlabel("Amplitude");
    xlim([min(x),max(x)]);
    ylim([min(y), max(y)]);
    zlim([c_min,c_max]);
    title(['t=', num2str(t(i))]);
    drawnow limitrate;
    pause(0.5);
end