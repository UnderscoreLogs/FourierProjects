clc, clearvars, close all

n = 2^7;
L = 10; % Define the length of the domain
dx = 10/(n-1); % Define the step size

x = 0:dx:L;
y = 0:dx:L;

[X, Y] = meshgrid(x, y);

f = (sin(pi*X/L))^2 .* (sin(pi*Y/L))^2; %((1 - cos(2*pi*X/L))/2 .* (1 - cos(2*pi*Y/L))/2)/3; %Function to test

correctValues = 0;
shift = n/2;

for i = 1:(n-shift)
    for j = 1:(n-shift)
        val = f(i,j) + f(i+shift, j) + f(i, j+shift);
        
        if abs(val - 1) < 1e-6
            correctValues = correctValues + 1;
        end
    end
end

figure (1)
surf(X, Y, f);
colorbar();
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Function Value');
title(['3D Surface Plot of f(x, y) ' ...
    '# of correct values: ' num2str(correctValues)]);

