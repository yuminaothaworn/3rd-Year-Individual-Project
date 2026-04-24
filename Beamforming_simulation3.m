%final version for beamforming
close all;
clear all;
cd('C:\Users\Yumi\OneDrive - The University of Manchester\Documents\3rd year Individual Project\Model');


f=1000;
c=3e8; %sped of light
lambda = c/f;
k = 2*pi/lambda;

d=lambda/2;
n=2;

theta0 = 0;
theta = linspace(-90,90,4000);

%psi = (k.*d.*(sind(theta)-sind(theta0)));

gain = (sin(n.*(k.*d.*(sind(theta)-sind(theta0)))./2)./sin((k.*d.*(sind(theta)-sind(theta0)))./2)).^2;

plot(theta, gain);
xlabel('\theta (deg)');
ylabel('Gain');

save("BeamformingGain_n2.mat","gain","theta")