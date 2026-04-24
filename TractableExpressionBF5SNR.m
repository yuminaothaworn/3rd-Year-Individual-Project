
% clear all;
% close all;

hold on;


lambda = 1e-0;
alpha  = 4;
mu     = 1;


f = 1000;
c = 3e8;
lambdaBF = c/f;
k = 2*pi/lambdaBF;

d = lambdaBF/2;
n = 8;

theta0 = 0;
theta = 0;

delta = sind(theta) - sind(theta0);

if abs(delta) < 1e-12
    BFgain = n^2;
else
    BFgain = (sin(n*(k*d*delta)/2) / sin((k*d*delta)/2))^2;
end


nbrPoints = 50;
SNR_dB_vec = linspace(-10,100,nbrPoints);
SNR_linear_vec = 10.^(SNR_dB_vec./10);


threshold_dB = -10;
T = 10^(threshold_dB/10);

pc = zeros(size(SNR_dB_vec));

for i = 1:length(SNR_dB_vec)

    SNR_linear = SNR_linear_vec(i);
    sigma2 = 1 / (mu * SNR_linear);

    pc(i) = pc_fun(T, lambda, alpha, mu, sigma2, BFgain);

    fprintf('Simulating SNR = %.1f dB\n', SNR_dB_vec(i));
end

plot(SNR_dB_vec, pc, '*');

xlabel('SNR (dB)');
ylabel('Coverage Probability');
grid on;

ylim([0 1.1]);

save("ExpressionResultsBF_SNR_4", "pc", "SNR_dB_vec");



function pc = pc_fun(T, lambda, alpha, mu, sigma2, BFgain)

    beta = beta_fun(T, alpha, mu, BFgain);

    integrand_v = @(v) exp( ...
        -pi*lambda*v.*beta ...
        -mu*T/BFgain*sigma2.*v.^(alpha/2) );

    pc = pi*lambda * integral(integrand_v, 0, Inf);

end


function beta = beta_fun(T, alpha, mu, BFgain)

    fg = @(g) exp(-g);

    a = -2/alpha;

    integrand_g = @(g) ...
        g.^(2/alpha) .* ...
        ( arrayfun(@(x) upperGamma(a, mu*T/BFgain*x), g) ...
        - gamma(a) ) .* fg(g);

    Eg = integral(integrand_g, 0, Inf, ...
                  'RelTol',1e-6,'AbsTol',1e-10);

    beta = 2*(mu*T/BFgain)^(2/alpha)/alpha * Eg;

end


function G = upperGamma(a, x)
    G = integral(@(t) t.^(a-1).*exp(-t), x, Inf, ...
                 'RelTol',1e-8,'AbsTol',1e-12);
end