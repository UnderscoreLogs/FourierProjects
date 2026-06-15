clc, clearvars, close all
k = 1; % thermal conductivity
rho = 1; % Density
c_p = 1; % Specific heat capacity

alpha_calculated = k/(rho*c_p); % Thermal diffusivity constant alpha = k/(rho*c_p)

% Silver: 1.66 × 10⁻⁴
% Copper: 1.17 × 10⁻⁴
% Aluminum: 9.70 × 10⁻⁵
% Iron: 2.30 × 10⁻⁵
% Steel: 1.20 × 10⁻⁵
% Concrete: 7.00 × 10⁻⁷
% Glass: 3.40 × 10⁻⁷
% Water: 1.43 × 10⁻⁷
% Air: 1.90 × 10⁻⁵

alpha = 1.17e-4;

% Size parameters
L = 0.5; % Length (x) (meters)
W = 0.5; % Width (y) (meters)
N = 100; % number of samples

% Time parameters
T = 100; % Total time (seconds)
dt = 1; % Time step
t = 0:dt:T; % t domain

% Temperature parameters (°C or K)
hotspot = 100; 
ambient = 20;

dx = L/N;
x = -L/2:dx:L/2-dx; % x domain
dy = W/N;
y = -W/2:dy:W/2-dy; % y domain

[X, Y] = meshgrid(x,y);

%Intital Condition
u0 = hotspot - (hotspot - ambient) * (X.^2 + Y.^2) / max(max(X.^2 + Y.^2)) * 1.5;
%u0 = zeros(N,N);
% u0(abs(X) < L/5 & abs(Y) < W/5) = hotspot - ambient;
u0hat = fft2(u0); %2d fft

kappa_x = (2*pi/L) * [0:N/2-1 -N/2:-1];
kappa_y = (2*pi/W) * [0:N/2-1 -N/2:-1];

[KX, KY] = meshgrid(kappa_x, kappa_y);

u = zeros(length(kappa_x),length(kappa_y),length(t));
for i = 1:length(t)
    uhat = u0hat .* exp(-alpha .* (KX.^2 + KY.^2) .* t(i));
    u(:,:,i) = real(ifft2(uhat)); % temperature equation for x,y,t
    %total_energy(i) = sum(sum(u(:,:,i))) % checking that total energy
end


figure (1)
    surf(X, Y, u(:,:,1));
    shading interp;
    colormap jet;
    colorbar;
    clim([ambient, hotspot]); 
    xlim([-L/2 L/2]);
    ylim([-W/2 W/2]);
    zlim([0, hotspot]);
    xlabel('x'); ylabel('y'); 
    zlabel('Temperature (above ambient)');
    title(['Time = ' num2str(t(1))]);

figure (3)
    surf(X, Y, u(:,:,100));
    shading interp;
    colormap jet;
    colorbar;
    clim([ambient, hotspot]); 
    xlim([-L/2 L/2]);
    ylim([-W/2 W/2]);
    zlim([0, hotspot]);
    xlabel('x'); ylabel('y'); 
    zlabel('Temperature (above ambient)');
    title(['Time = ' num2str(t(100))]);    

figure (4)
    surf(X, Y, u(:,:,51));
    shading interp;
    colormap jet;
    colorbar;
    clim([ambient, hotspot]); 
    xlim([-L/2 L/2]);
    ylim([-W/2 W/2]);
    zlim([0, hotspot]);
    xlabel('x'); ylabel('y'); 
    zlabel('Temperature (above ambient)');
    title(['Time = ' num2str(t(51))]); 

% figure (1)
% for i = 1:length(t)
%     surf(X, Y, u(:,:,i));
%     shading interp;
%     colormap jet;
%     colorbar;
%     clim([ambient, hotspot]); 
%     xlim([-L/2 L/2]);
%     ylim([-W/2 W/2]);
%     zlim([0, hotspot]);
%     xlabel('x'); ylabel('y'); 
%     zlabel('Temperature (above ambient)');
%     title(['Time = ' num2str(t(i))]);
%     drawnow limitrate;
%     pause(0.05);
% end
