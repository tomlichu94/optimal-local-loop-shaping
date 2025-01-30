clear all
close all
clc
addpath('Functions','Simulink','Data')
load mainPlantData;
%%%%%%%%%%%%%%%%%%%%% Loading Simulink parameters %%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 10; % multiplier to increase sampling time
PlantData.Ts = PlantData.Ts*Nx;
PlantData.Tu = PlantData.Tu*Nx;
Tu = PlantData.Tu; % baseline plant sampling time
Ps = PlantData.Pn; % continuous-time plant
Ps = tf(Ps); % define CT plnat as a transfer function
Pz = c2d(Ps,Tu,'zoh'); % discretize the plant
z = tf('z',Tu);
Q = (1-z^-1);
Qs1 = d2c(Q,'tustin'); % through bi-linear transformation Q(s) = 2s/(s+2/T)
% zero at the origin, pole at -2/T.
% Qs2= d2c(Q,'zoh');
% Qs3 = d2c(Q,'foh');
figure
bode(Q,Qs1)
legend('Qz','Tustin')

S1 = feedback(Pz,1);
S2 = feedback(Pz*Q,1);
figure
bode(S1,S2)
legend('P','PQ')