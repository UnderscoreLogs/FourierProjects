clc, clearvars, close all

N = 64; % number of samples
L = 30; % Total Length

dx = L/N;
x = -L/2:dx:L/2-dx;

f = cos(x) .* exp((-x.^2)/25); % Original Function
realDerivative = -((sin(x) .* exp((-x.^2)/25)) + ((2/25)*x.*f)); % Actual Derivative

numericalDerivative = zeros(1,length(realDerivative));
for kappa = 1:length(realDerivative) - 1
    numericalDerivative(kappa) = (f(kappa + 1) - f(kappa)) / (dx);
end

fhat = fft(f);
kappa = (2*pi/L) * [0:N/2-1 -N/2:-1];
derivativeFhat = 1i*kappa.*fhat;
fftDerivative = real(ifft(derivativeFhat));

figure (1)

subplot(3,1,1)
plot(x,realDerivative,'b')
hold on
plot(x,numericalDerivative,'--r')
hold on
plot(x,fftDerivative,'--y')
legend('Calculated Derivative', 'Numerical Derivative', 'FFT Derivative')
title('Derivative Methods Comparison')

%% Comparison as n increases

max_n = 20;
min_n = 4;
stop_display = 150;
n_axis = min_n:2:max_n;

errorNumerical = zeros(1,length(n_axis));
errorFFT = zeros(1,length(n_axis));
index = 1;

figure(1)
subplot(3,1,2)
h1 = plot(x, realDerivative, 'b'); hold on
h2 = plot(x, numericalDerivative, '--r');
h3 = plot(x, fftDerivative, '--y');
legend('Exact','Numerical','FFT')

for n = min_n:2:max_n % Only even n(s)
    dx_n = L/n;
    x_n = -L/2:dx_n:L/2-dx_n;
    
    f_n = cos(x_n) .* exp((-x_n.^2)/25); % Original Function
    realDerivative_n = -((sin(x_n) .* exp((-x_n.^2)/25)) + ((2/25)*x_n.*f_n)); % Actual Derivative
    
    numericalDerivative_n = zeros(1,length(realDerivative_n));
    for k = 1:length(realDerivative_n) - 1
        numericalDerivative_n(k) = (f_n(k + 1) - f_n(k)) / (dx_n);
    end
    
    fhat_n = fft(f_n);
    k = (2*pi/L) * [0:n/2-1 -n/2:-1]; % fft in mat lab has the negative values in the second half of the matrix
    derivativeFhat_n = 1i*k.*fhat_n;
    fftDerivative_n = real(ifft(derivativeFhat_n));

    errorNumerical(index) = max(abs(realDerivative_n - numericalDerivative_n)); % Largest difference between the numerical derivative and the real derivative
    errorFFT(index) = max(abs(realDerivative_n - fftDerivative_n)); % Largest difference between the fft derivative and the real derivative
    index = index + 1;
    
    if n <= stop_display
        set(h1,'XData',x_n,'YData',realDerivative_n)
        set(h2,'XData',x_n,'YData',numericalDerivative_n)
        set(h3,'XData',x_n,'YData',fftDerivative_n)
    
        title(['n = ', num2str(n)])
    
        drawnow limitrate;
        pause(0.05);
    end

end



figure(1)
subplot(3,1,3)
plot(n_axis,errorNumerical, '*r')
hold on
plot(n_axis,errorFFT, '*y')
legend('Numerical Error', 'FFT Error')
title('Error Comparison')
