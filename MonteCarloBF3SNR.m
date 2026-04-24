clear all;
close all;

hold on;

cd('C:\Users\Yumi\OneDrive - The University of Manchester\Documents\3rd year Individual Project\Model');

alpha = 4;
monteCarloRun = 100;

lambda = 1e-0;
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

nbrPoints = 50;
SNR_dB_vec = linspace(-10,100,nbrPoints);
SNR_linear_vec = 10.^(SNR_dB_vec./10);

threshold_dB = -10;
threshold = 10^(threshold_dB/10);

coverageProbability=[];

for i = 1:nbrPoints

    SNR_linear = SNR_linear_vec(i);
    sigma = 1 / (u * SNR_linear);

    count = 0;

    for j = 1:monteCarloRun


        nbrBS = poissrnd(lambda * area);

        bsx = (rand(1, nbrBS) - 0.5) * squarelenth;
        bsy = (rand(1, nbrBS) - 0.5) * squarelenth;

        distance_all = sqrt(bsx.^2 + bsy.^2);

        [closestBSloc, indexClosestBS] = min(distance_all);
        distance = distance_all;
        distance(indexClosestBS) = [];


        fading_exp = exprnd(1/u, 1, nbrBS);

        fading = fading_exp(indexClosestBS);
        fadingInterferingBS = fading_exp;
        fadingInterferingBS(indexClosestBS) = [];


        pathloss = closestBSloc^(-alpha);
        interferingBSPathLoss = distance.^(-alpha);


        interference = sum(fadingInterferingBS .* interferingBSPathLoss);


        SINR = (fading * pathloss * BFgain) / (sigma + interference);

        if SINR > threshold
            count = count + 1;
        end
    end
    avg = count/monteCarloRun;
    coverageProbability= [coverageProbability, avg];

    fprintf('Simulating SNR = %.1f dB\n', SNR_dB_vec(i));
end

% ===== plot =====
plot(SNR_dB_vec, coverageProbability, '-^');
xlabel('SNR (dB)');
ylabel('Coverage probability');
grid on;

ylim([0 1.1]);

save("MonteCarloResultsBF_SNR_4.mat", "coverageProbability", "SNR_dB_vec", "threshold_dB");