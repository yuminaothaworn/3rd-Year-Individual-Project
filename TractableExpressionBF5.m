%match with MCBF3
% 
% clear all;
% close all;


hold on;

% Parameters
lambda = 1e-4;      % BS density
alpha  = 4;         % path loss exponent
mu     = 1;         % fading parameter

noise_db=100000;
noise_linear= 10^(noise_db/10);

sigma2 = 1 / (mu*noise_linear); % noise

%bf param
% BFgain = 60;
f=3e9;
c=3e8; %speed of light
lambdaBF = c/f;
k = 2*pi/lambdaBF;

d=lambdaBF/2;
n=16;

theta0 = 0;
theta = 0;

delta = sind(theta) - sind(theta0);

if abs(delta) < 1e-12
    BFgain = n^2;
else
    BFgain = (sin(n*(k*d*delta)/2)/sin((k*d*delta)/2))^2;
end

% T range (SINR threshold)
threshold_dB = linspace(-10,20,50);
Tvec = 10.^(threshold_dB./10);


pc = zeros(size(Tvec));

for i = 1:length(Tvec)
    T = Tvec(i);
    pc(i) = pc_fun(T, lambda, alpha, mu, sigma2, BFgain);
end

p_c=pc;
plot(threshold_dB, p_c, '*');
xlabel('Threshold T (dB)');
ylabel('Coverage Probability');
grid on;

save("ExpressionResults5__3_16n","p_c","threshold_dB");


function pc = pc_fun(T, lambda, alpha, mu, sigma2, BFgain)

    beta = beta_fun(T, alpha, mu, BFgain);

    integrand_v = @(v) exp( ...
        -pi*lambda*v.*beta ...
        -mu*T/BFgain*sigma2.*v.^(alpha/2) );

    pc = pi*lambda * integral(integrand_v, 0, Inf);

end

function beta = beta_fun(T, alpha, mu, BFgain)

    % Rayleigh fading: g ~ Exp(1)
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



