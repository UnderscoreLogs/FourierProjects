clc, clearvars, close all

x = -5:0.001:5;
y = exp(-1*x.^2);

graph = sqrt(pi)*exp(-(2*pi^2)*x.^2); %0.5* ((x==15) + (x==-15));

f= -length(y)/2:length(y)/2-1;

figure(1)
plot(x,graph);
title('Fourier Space: Frequency vs. Magnitude');
ylim([0 max(real(graph)+0.1)]);


figure(2)
plot(x,y);
title('Gaussian: Position vs Magnitude');
ylim([0 max(y+0.1)]);