clc, clearvars, close all

v = 1; % 
L = 100; % Length 
N = 1000; % number of samples
dx = L/N;
x = -L/2:dx:L/2-dx; % x domain

%Intital Condition
f = exp(-x.^2 / 10); % sin(4*pi*x/L); %u(x,0) %Gaussian Pulse: exp(-x.^2 / 10)
%f((L/2 - L/10)/dx:(L/2+L/10)/dx) = 1;
fhat = fft(f);

g = zeros(1,N); %u'(x,0)
ghat = fft(g);

kappa = (2*pi/L) * [0:N/2-1 -N/2:-1];

% Time parameters
dt = 1; % Time step
T = 50; % Total time
t = (0:dt:T).'; % t domain

uhat = (fhat .* cos(v .* kappa .* t)) + ((ghat ./ (v .* kappa)) .* sin(v .* kappa .* t));
uhat(:,1) = fhat(1)+ghat(1)*t;
u = real(ifft(uhat, [], 2)); %wave equation

y_min = min(u, [], 'all') - 0.5;
y_max = max(u, [], 'all') + 0.5;

% figure (2)
% 
% for i = 1:length(t)
%     plot(x,u(i,:),'b');
%     ylabel("Amplitude");
%     xlabel("Distance");
%     xlim([min(x),max(x)]);
%     ylim([y_min, y_max]);
%     title(['t=', num2str(t(i))]);
%     drawnow limitrate;
%     pause(0.05);
% end


figure (1)
h=waterfall(x(1:5:end), t(1:10:end), u(1:10:end,1:5:end));
set(h,'LineWidth',5,'FaceAlpha',.5);
colormap(flipud(jet)/1.5)
xlabel('Space');
ylabel('Time');
zlabel('Displacement');