clear; clc; close all;

%% Parameters
N = 8;                 % Number of antennas
c = 3e8;                % Speed of light
f = 10e9;               % Frequency
lambda = c/f;
d = lambda/2;           % Spacing

theta0 = 0;            % Steering angle (deg)
theta = linspace(-90,90,4000);

k = 2*pi*f/c;

theta_rad  = deg2rad(theta);
theta0_rad = deg2rad(theta0);

%% Compute Gain

psi = k*d*(sin(theta_rad) - sin(theta0_rad));

% Avoid division by zero
numerator   = sin(N*psi/2);
denominator = sin(psi/2);
denominator(abs(denominator)<1e-12) = 1e-12;

G = (numerator ./ denominator).^2;

%% Plot
% G_dB = 10*log10(G);
% plot(theta, G_dB);
% ylabel('Gain (dB)');
% xlabel('\theta (deg)');
% ylim([-60 10*log10(N^2)]);
% xlim([-90 90]);



figure('Color','w');
plot(theta, G,'LineWidth',2);
xlabel('\theta (deg)');
ylabel('Gain = |h^H(\theta)\omega|^2');
title(['Beamforming Gain (N = ' num2str(N) ', Peak = N^2)']);
grid on;
% xlim([-90 90]);

