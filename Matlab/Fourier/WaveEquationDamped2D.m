clc, clearvars, close all

v = 1; % 
gamma = 0.01; % damping coeffiecent 
L = 100; % Length
W = 100; %Width
N = 100; % number of samples

dx = L/N;
x = -L/2:dx:L/2-dx; % x domain
dy = W/N;
y = -W/2:dy:W/2-dy; % y domain

[X, Y] = meshgrid(x,y);

%Intital Condition
f = -exp((-(X).^2 - (Y).^2)/3); %u(x,0) %Gaussian Pulse: exp((-(X).^2 - (Y).^2)/3) SIne Wave: sin(4pix/L)
%f((L/2 - L/10)/dx:(L/2+L/10)/dx) = 1;
fhat = fft2(f);

g = zeros(N,N); %u'(x,0)
ghat = fft2(g);

kappa_x = (2*pi/L) * [0:N/2-1 -N/2:-1];
kappa_y = (2*pi/W) * [0:N/2-1 -N/2:-1];

[KX, KY] = meshgrid(kappa_x, kappa_y);

% Time parameters
dt = 0.01; % Time step
T = 0.01; % Total time
t = (0:dt:T).'; % t domain

%u = zeros(N, N, length(t));
omega = (v.^2 .* KX.^2) + (v.^2 .* KY.^2); %should be named omega^2?
threshold = 1e-10;
under_mask = (gamma.^2 - 4 .* omega) < threshold;
over_mask = (gamma.^2 - 4 .* omega) > threshold;
crit_mask = ~under_mask & ~over_mask;

omega_d_over  = zeros(N,N);
omega_d_under = zeros(N,N);
omega_d_over(over_mask)   = sqrt(gamma.^2 - 4.*omega(over_mask))   ./ 2; %Overdamped omega (regular)
omega_d_under(under_mask) = sqrt(4.*omega(under_mask) - gamma.^2)  ./ 2;%Underdamped omega (switch the sign)

for i = 1:length(t) 
    uhat = zeros(N,N);
    %Underdamped
    uhat(under_mask) = exp((-gamma/2) .* t(i)) .* ((fhat(under_mask) .* cos(omega_d_under(under_mask) .* t(i))) + ((ghat(under_mask) + (gamma/2) .* fhat(under_mask)) ./ omega_d_under(under_mask)) .* sin(omega_d_under(under_mask) .* t(i)));

    %Overdamped
    uhat(over_mask) = ((fhat(over_mask) - ((fhat(over_mask) .* (-gamma./2 + omega_d_over(over_mask)) - ghat(over_mask))./(2.*omega_d_over(over_mask)))) .* exp(((-gamma./2) + omega_d_over(over_mask)) .* t(i))) + (((fhat(over_mask) .* (-gamma./2 + omega_d_over(over_mask)) - ghat(over_mask))./(2.*omega_d_over(over_mask)))) .* exp(((-gamma./2) - omega_d_over(over_mask)) .* t(i));

    %Critically damped
    uhat(crit_mask) = fhat(crit_mask) .* exp((-gamma./2).*t(i)) + (ghat(crit_mask) + (gamma./2) .* fhat(crit_mask)).*t(i).*exp((-gamma./2).*t(i));

    u_frame = real(ifft2(uhat)); %wave equation
    
    figure (1)
    surf(X,Y,u_frame);
    shading interp;
    colormap jet;
    colorbar;
    clim([-0.1, 0.1]);
    ylabel("Y");
    xlabel("X");
    zlabel("Amplitude");
    xlim([min(x),max(x)]);
    ylim([min(y), max(y)]);
    zlim([-1,1]);
    title(['t=', num2str(t(i))]);
    drawnow limitrate;
    pause(0.01);
end


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

c_min = min(u_frame, [], 'all');
c_max = max(u_frame, [], 'all');
