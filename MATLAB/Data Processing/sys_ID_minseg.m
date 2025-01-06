% clear all
close all
clc

%% run SID to identify the model, MinsegMotor
Ts = 1/500; % sampled time

% data collected using quad encoder
t_sim = in_voltage.time;
u_sim = in_voltage.signals.values;
y_sim = squeeze(double(out_encoder.signals.values));

%% pre-process data using filter from Spiral Servo Writing in Hard Disk 
% Drives Using Iterative Learning Based Tracking Control
% Q_filter = (1+z^-1)(1+z)/4 = (z+2-z^-1)/4, can only do this offline
n_enc = height(y_sim);
y_filt = double.empty(n_enc,0);
y_filt(1) = y_sim(1);
for i = 2:(n_enc-1)
    y_filt(i) = 0.25*(y_sim(i+1)+2*y_sim(i)-y_sim(i-1));
end
y_filt = y_filt';
t_filt = t_sim(1:end-1);
u_filt = u_sim(1:end-1);
figure
plot(t_sim,y_sim,t_filt,y_filt)
legend('Unprocessed','Post-processed')

figure
plot(t_sim,u_sim,t_filt,u_filt)
legend('Unprocessed','Post-processed')
xlim([0.5 0.52])
title('Input')

%% n4sid (sysID)
clear A1 B1 C1 D1 A2 B2 C2 D2
nx = 4;
sys_pre = n4sid(u_sim,y_sim,nx,'Ts',Ts);
sys_post = n4sid(u_filt,y_filt,nx,'Ts',Ts);

A1 = sys_pre.A;
B1 = sys_pre.B;
C1 = sys_pre.C;
D1 = sys_pre.D;
A2 = sys_post.A;
B2 = sys_post.B;
C2 = sys_post.C;
D2 = sys_post.D;

sys_pre = ss(A1,B1,C1,D1);
sys_post = ss(A2,B2,C2,D2);

figure
bode(sys_pre)
hold on
bode(sys_post)
legend('Pre-processing','Post-processing')