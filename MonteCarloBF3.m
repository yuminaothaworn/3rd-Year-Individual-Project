%match wit beamforming5 -> only serving BS has BF
%{
- correccted expression relating to noise
%}
% hold on;
clear all;
close all;
cd('C:\Users\Yumi\OneDrive - The University of Manchester\Documents\3rd year Individual Project\Model');

alpha= 4; %path loss exponent
sigma_db = 1000000000;
%sigma_db = inf;
SNR_linear = 10^(sigma_db/10);
%SNR_linear = 10;
%sigma=10; %noise
monteCarloRun=10000;
lambda = 1e-4;
squarelenth = 2e3;
area = squarelenth^2;
% this is never defined power_db = 0;
u=1;


%bf param

f=1000;
c=3e8; %speed of light
lambdaBF = c/f;
k = 2*pi/lambdaBF;

d=lambdaBF/2;
n=8;

theta0 = 20;
theta = 0;

nbrPoints = 50;
% used as power threshold for coverage proability
threshold_dB=linspace(-10,20,nbrPoints); %number of points on x axis, i.e. resolution of graph
threshold = 10.^(threshold_dB./10); %convert threshold to linear function

delta = sind(theta) - sind(theta0);
%BFgain = (sin(n*(k*d*delta)/2)/sin((k*d*delta)/2))^2;

if delta ==0;
    BFgain = n^2;
else
    BFgain = (sin(n*(k*d*delta)/2)/sin((k*d*delta)/2))^2;
end

coverageProbability=[];

%run loop for each value of threshold
for i=1:nbrPoints
    count = 0;
    for j=1:monteCarloRun

        %generate number of BS according to PPP
        nbrBS=poissrnd(lambda*area);
        
        %generate BS coords centered around user at origin
        bsx=(rand(1,nbrBS)-0.5)*squarelenth;
        bsy=(rand(1,nbrBS)-0.5)*squarelenth;

        %distance from user at origin
        distance_all = (bsx.^2 + bsy.^2).^0.5;
        distance = distance_all;

        %find location of the closest BS
        [closestBSloc,indexClosestBS]=min(distance_all);

        %remove element of tagged BS
        distance(indexClosestBS) = [];
        

        

        %fading
        %fading = exprnd(1);
        %fadingInterferingBS=exprnd(1,1,nbrBS);
        fading_exp = exprnd(1/u,1,nbrBS);

        %Remove element of tagged BS
        fading = fading_exp(indexClosestBS);
        fadingInterferingBS = fading_exp;
        fadingInterferingBS(indexClosestBS) = [];

        %path loss --> r^(-aplpha)
        pathloss = closestBSloc^(-alpha);
        interferingBSPathLoss = distance.^(-alpha);

        %snr=1/(u*sigma^2), u=1--> rayleigh fading with mean 1
        sigma = 1 / (u*SNR_linear);

        
        
        %sigma = closestBSloc^(-alpha) / SNR_linear;
        %{
        noise is usually fixed value
        calculated based on SNR
        relative signal power --> transmite power/noise (ignore
        interference)
        noise is due to circuit --> amplifyer TX & Rx
        independent from signal propagation
        %}

        %Interference=sum of fading of all other BS multiplied by distance
        %of each BS
        %sum of multiplication of path loss and distance of all BS - tagged
        %BS                                                     
        interference = sum(fadingInterferingBS.*interferingBSPathLoss);
        %multiply the distance of the tagged BS with the gi because of
        %factoring --> want to remove it from the total interference value

        %SINR = (fading*pathloss)/(noise + Ir)
        SINR= (fading*pathloss*BFgain)/((sigma)+interference);

        if SINR > threshold(i)
                count = count +1;
        end
    end
    avg = count/monteCarloRun;
    coverageProbability= [coverageProbability, avg];
    fprintf('Simulating threshold = %.1f dB\n', threshold_dB(i));
end

%plot results
x=threshold_dB; %row of matrix
y=coverageProbability; %column of matrix
plot(x,y,'-^');
xlabel('threshold power (dB)');
ylabel('coverage probability');
ylim([0 1.1]);


save("MonteCarloResultsBF3_1_20deg.mat","coverageProbability","threshold_dB")