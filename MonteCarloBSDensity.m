clear all;
close all;
cd('C:\Users\Yumi\OneDrive - The University of Manchester\Documents\3rd year Individual Project\Model');

alpha= 2.5; %path loss exponent

%SNR_linear = 10;
%sigma=10; %noise
monteCarloRun=10000;
% lambda = 1e-4;
squarelenth = 2e3;
area = squarelenth^2;
% this is never defined power_db = 0;
u=1;


% used as power threshold for coverage proability
threshold_dB=-10; %number of points on x axis, i.e. resolution of graph
threshold = 10^(threshold_dB/10); %convert threshold to linear function

sigma_db = 20;
%sigma_db = inf;
SNR_linear = 10^(sigma_db/10);

nbrPoints = 25;
lambda = linspace(1e-4,1e-2,nbrPoints);


coverageProbability=[];

%run loop for each value of threshold
for i=1:nbrPoints
    count = 0;
    for j=1:monteCarloRun

        %generate number of BS according to PPP
        nbrBS=poissrnd(lambda(i)*area);
        
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
        sigma = 1 / (u*SNR_linear); %signal power = 1
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
        SINR= (fading*pathloss)./((sigma)+interference);

        if SINR > threshold
                count = count +1;
        end
    end
    avg = count/monteCarloRun;
    coverageProbability= [coverageProbability, avg];
    fprintf('Simulating lambda = %.3e \n', lambda(i));
end

%plot results
x=lambda; %row of matrix
y=coverageProbability; %column of matrix
plot(x,y,'-^');
xlabel('BS density');
ylabel('coverage probability');


save("MonteCarloResultsBS_5_1.mat","coverageProbability","lambda")