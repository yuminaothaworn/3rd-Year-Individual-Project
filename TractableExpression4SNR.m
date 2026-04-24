clear all;
%hold on;


lambda = 1;      
alpha  = 4;
mu     = 1;


noise_db = linspace(-10, 100, 50);
noise_linear = 10.^(noise_db/10);
sigma2_vec = 1 ./ (mu .* noise_linear);


threshold_dB = -10;
T = 10^(threshold_dB/10);


pc = arrayfun(@(sig) pc_fun(T, lambda, alpha, mu, sig), sigma2_vec);


plot(noise_db, pc, 'LineWidth', 1.5);
xlabel('SNR (dB)');
ylabel('Coverage Probability');
title(['Coverage Probability vs SNR (T = ', num2str(threshold_dB), ' dB)']);
grid on;

save("ExpressionResults4_SNR_4_1","pc","noise_db");



function pc = pc_fun(T, lambda, alpha, mu, sigma2)

    beta = beta_fun(T, alpha, mu);

    integrand_v = @(v) exp( ...
        -pi * lambda * v * beta ...
        -mu * T * sigma2 * v.^(alpha/2) );

    pc = pi * lambda * integral(integrand_v, 0, Inf);

end


function beta = beta_fun(T, alpha, mu)

    a = -2 / alpha;

    integrand_g = @(g) ...
        g.^(2/alpha) .* ...
        ( arrayfun(@(x) upperGamma(a, mu*T*x), g) ...
        - gamma(a) ) .* exp(-g);

    Eg = integral(integrand_g, 0, Inf, ...
                  'RelTol',1e-6,'AbsTol',1e-10);

    beta = 2 * (mu*T)^(2/alpha) / alpha * Eg;

end


function G = upperGamma(a, x)

    G = integral(@(t) t.^(a-1).*exp(-t), x, Inf, ...
                 'RelTol',1e-8,'AbsTol',1e-12);

end