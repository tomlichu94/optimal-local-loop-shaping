close all
clear all
clc

addpath('Functions','Simulink','Data')
load mainPlantData;

%%%%%%%%%%%%%%%%%%%%% Loading Simulink parameters %%%%%%%%%%%%%%%%%%%%%%%%%
PlantData.Ts = PlantData.Ts*Nx;
PlantData.Tu = PlantData.Tu*Nx;
Tu = PlantData.Tu; % baseline plant sampling time
Ps = PlantData.Pn; % continuous-time plant
Ps = tf(Ps); % define CT plnat as a transfer function
Pz_delay = c2d(Ps,Tu,'zoh'); % discretize the plant
z = tf('z',Tu); % discrete time based on fast sampling
s = tf('s'); % continous time
w_d = 0.8; % fast measurement of disturbance in radians
w_n = w_d/(Tu);

% Generalized Q(z) form is: Q(z) = Q0*Q_FIR or Q0*Q_IIR
Q0_num = [1 -1];
Q1_num = [1 1];
Q0_den = [1 0]; 
Q0 = tf(Q0_num,Q0_den,Tu); % Q0 = 1-z^-1 = (z-1)/z
Q1 = tf(Q1_num,Q0_den,Tu);
Q2_num = conv(Q0_num,Q1_num);
Q2_den = conv(Q0_den,Q0_den);
Q2 = tf(Q2_num,Q2_den,Tu);


%% from Automatica
rho = 0.9;
Q0_num = q0(w_d,rho);
Q0_den = zeros(size(Q0_num));
Q0_den(1) = 1;
Q3 = tf(Q0_num,Q0_den,Tu);


figure
bode(Q0)
hold on
bode(Q1)
bode(Q2)
bode(Q3)
legend('1-z','1+z','Combined','B0')
xline(w_n)