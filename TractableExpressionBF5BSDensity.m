% Coverage probability vs BS density (lambda) - expression model with BF



hold on;

alpha  = 2.5;
mu     = 1;

noise_db = 1000000;
noise_linear = 10^(noise_db/10);
sigma2 = 1 / (mu * noise_linear);


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


nbrPoints = 30;
lambda_vec = linspace(1e-4, 1e-2, nbrPoints);


threshold_dB = -10;
T = 10^(threshold_dB/10);

pc = zeros(size(lambda_vec));

for i = 1:length(lambda_vec)

    lambda = lambda_vec(i);

    pc(i) = pc_fun(T, lambda, alpha, mu, sigma2, BFgain);

    fprintf('Simulating lambda = %.2e\n', lambda);
end

plot(lambda_vec, pc, '*');

xlabel('Base station density \lambda');
ylabel('Coverage Probability');
grid on;

save("ExpressionResultsBF_BS_4", "pc", "lambda_vec");



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