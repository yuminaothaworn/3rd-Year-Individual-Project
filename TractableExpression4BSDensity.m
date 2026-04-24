

% clear all;
% close all;

clc;
hold on;

% Parameters
lambda = linspace(1e-4,1e-2,25);   % BS density
alpha  = 2.5;                        % path loss exponent
mu     = 1;                        % fading parameter

noise_db = 100; %SNR
noise_linear = 10^(noise_db/10);
sigma2 = 1 / (mu * noise_linear);  % noise power

% SINR threshold (fixed)
threshold_dB = -10;
T = 10^(threshold_dB/10);

% Compute coverage probability using arrayfun
pc = arrayfun(@(lam) pc_fun(T, lam, alpha, mu, sigma2), lambda);

% Plot
plot(lambda, pc, 'LineWidth', 1.5);
xlabel('BS Density \lambda');
ylabel('Coverage Probability');
title(['Coverage Probability vs BS Density (T = ', num2str(threshold_dB), ' dB)']);
grid on;

save("ExpressionResults4_BSDensity_7_1","pc","lambda");



function pc = pc_fun(T, lambda, alpha, mu, sigma2)

    beta = beta_fun(T, alpha, mu);

    integrand_v = @(v) exp( ...
        -pi * lambda * v * beta ...
        -mu * T * sigma2 * v.^(alpha/2) );

    pc = pi * lambda * integral(integrand_v, 0, Inf);

end


function beta = beta_fun(T, alpha, mu)

    % Rayleigh fading: g ~ Exp(1)
    fg = @(g) exp(-g);

    a = -2 / alpha;

    integrand_g = @(g) ...
        g.^(2/alpha) .* ...
        ( arrayfun(@(x) upperGamma(a, mu*T*x), g) ...
        - gamma(a) ) .* fg(g);

    Eg = integral(integrand_g, 0, Inf, ...
                  'RelTol',1e-6,'AbsTol',1e-10);

    beta = 2 * (mu*T)^(2/alpha) / alpha * Eg;

end


function G = upperGamma(a, x)

    G = integral(@(t) t.^(a-1).*exp(-t), x, Inf, ...
                 'RelTol',1e-8,'AbsTol',1e-12);

end