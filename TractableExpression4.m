%correct version for expression

clear all;
close all;

hold on;

% Parameters
lambda = 1e-4;      % BS density
alpha  = 2.5;         % path loss exponent
mu     = 1;         % fading parameter

noise_db=10000000;
noise_linear= 10^(noise_db/10);

sigma2 = 1 / (mu*noise_linear); % noise


% T range (SINR threshold)
threshold_dB = linspace(-10,20,50);
Tvec = 10.^(threshold_dB./10);


pc = zeros(size(Tvec));

for i = 1:length(Tvec)
    T = Tvec(i);
    pc(i) = pc_fun(T, lambda, alpha, mu, sigma2);
end

p_c=pc;
plot(threshold_dB, p_c, 'LineWidth', 1.5);
xlabel('Threshold T (dB)');
ylabel('Coverage Probability');
grid on;

save("ExpressionResults2_area","p_c","threshold_dB");


function pc = pc_fun(T, lambda, alpha, mu, sigma2)

    beta = beta_fun(T, alpha, mu);

    integrand_v = @(v) exp( ...
        -pi*lambda*v.*beta ...
        -mu*T*sigma2.*v.^(alpha/2) );

    pc = pi*lambda * integral(integrand_v, 0, Inf);

end

function beta = beta_fun(T, alpha, mu)

    % Rayleigh fading: g ~ Exp(1)
    fg = @(g) exp(-g);

    a = -2/alpha;

    integrand_g = @(g) ...
        g.^(2/alpha) .* ...
        ( arrayfun(@(x) upperGamma(a, mu*T*x), g) ...
        - gamma(a) ) .* fg(g);

    Eg = integral(integrand_g, 0, Inf, ...
                  'RelTol',1e-6,'AbsTol',1e-10);

    beta = 2*(mu*T)^(2/alpha)/alpha * Eg;

end

function G = upperGamma(a, x)
    G = integral(@(t) t.^(a-1).*exp(-t), x, Inf, ...
                 'RelTol',1e-8,'AbsTol',1e-12);
end

