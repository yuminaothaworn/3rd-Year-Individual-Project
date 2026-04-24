%5th version with all variables correctly and appropriately named

clear all;
close all;

%path loss exponent
alpha=2.5;

SNR_db = 100000;
SNR_linear = 10^(SNR_db/10);

monteCarloRun = 10;

%base station density
lambda = 1e-4;

squareLength = 2e5;
area = squareLength^2;

%fading exponent
u=1;

nbrPoints = 50;
threshold_dB = linspace(-10,20,nbrPoints);
threshold=10.^(threshold_dB./10);

coverageProbability=[];

%run loop for each value of threshold
for i=1:nbrPoints
    count = 0;
    for j=1:monteCarloRun

        %generate number of BS according to PPP
        nbrBS=poissrnd(lambda*area);
        
        %generate BS coords centered around user at origin
        bsx=(rand(1,nbrBS)-0.5)*squareLength;
        bsy=(rand(1,nbrBS)-0.5)*squareLength;

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
        SINR= (fading*pathloss)/((sigma)+interference);

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


save("MonteCarloResultsx_2_3.mat","coverageProbability","threshold_dB")