clear all;
close all;

cd('C:\Users\Yumi\OneDrive - The University of Manchester\Documents\3rd year Individual Project\Model');

alpha = 2.5;
sigma_db = 1000000000;
SNR_linear = 10^(sigma_db/10);

monteCarloRun = 1000;

squarelenth = 2e3;
area = squarelenth^2;

u = 1;

f = 1000;
c = 3e8;
lambdaBF = c/f;
k = 2*pi/lambdaBF;

d = lambdaBF/2;
n = 8;

theta0 = 0;
theta = 0;

delta = sind(theta) - sind(theta0);

if delta == 0
    BFgain = n^2;
else
    BFgain = (sin(n*(k*d*delta)/2) / sin((k*d*delta)/2))^2;
end

nbrPoints = 25;
lambda_vec =  linspace(1e-4,1e-2,nbrPoints);

threshold_dB = -10;
threshold = 10^(threshold_dB/10);

coverageProbability=[];

for i = 1:nbrPoints

    lambda = lambda_vec(i);
    count = 0;

    for j = 1:monteCarloRun

        nbrBS = poissrnd(lambda * area);

        bsx = (rand(1, nbrBS) - 0.5) * squarelenth;
        bsy = (rand(1, nbrBS) - 0.5) * squarelenth;

        distance_all = sqrt(bsx.^2 + bsy.^2);
        distance = distance_all;

        [closestBSloc, indexClosestBS] = min(distance_all);
        distance(indexClosestBS) = [];

        % fading
        fading_exp = exprnd(1/u, 1, nbrBS);

        fading = fading_exp(indexClosestBS);
        fadingInterferingBS = fading_exp;
        fadingInterferingBS(indexClosestBS) = [];

        % path loss
        pathloss = closestBSloc^(-alpha);
        interferingBSPathLoss = distance.^(-alpha);

        % noise
        sigma = 1 / (u * SNR_linear);

        % interference
        interference = sum(fadingInterferingBS .* interferingBSPathLoss);

        SINR = (fading * pathloss * BFgain) / (sigma + interference);

        if SINR > threshold
            count = count + 1;
        end
    end

    % coverageProbability(i) = count / monteCarloRun;
    avg = count/monteCarloRun;
    coverageProbability= [coverageProbability, avg];
    fprintf('Simulating lambda = %.2e\n', lambda);
end

% ===== plot =====
plot(lambda_vec, coverageProbability, '-^');
xlabel('Base station density');
ylabel('Coverage probability');
grid on;

%set(gca, 'XScale', 'log');
ylim([0 1.1]);

save("MonteCarloResultsBF_BS_4.mat",  "coverageProbability", "lambda_vec", "threshold_dB");