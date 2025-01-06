clear all
close all
clc

addpath('Data for Signal Reconstruction/PRBS')
load sys_ID_data.mat;

%% run SID to identify the model, MinsegMotor
Ts = 1/500; % sampled time
Fs = 1/Ts;
f_in = 8; % input in Hz
% data collected using quad encoder
t_sim = in_voltage.time;
u_sim = squeeze(in_voltage.signals.values);
y_sim = squeeze(double(out_encoder.signals.values));
T_fin = 10;
%% pre-process data using filter from Spiral Servo Writing in Hard Disk 
y_filt = zero_phase_low_pass(y_sim); % apply zero-phase error low pass filter
y_filt = y_filt';

t_filt = t_sim(1:end-1);
u_filt = zero_phase_low_pass(u_sim);
u_filt = u_filt';
% u_filt = u_sim(1:end-1);

figure
plot(t_sim,y_sim,t_filt,y_filt)
legend('Unprocessed','Post-processed')

figure
plot(t_sim,u_sim,t_filt,u_filt)
legend('Unprocessed','Post-processed')
% xlim([0.5 0.52])
title('Input')

%% n4sid (sysID)
close all
clear A1 B1 C1 D1 A2 B2 C2 D2
nx = 1;
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

[num_pre den_pre] = ss2tf(A1,B1,C1,D1);
[num_post den_post] = ss2tf(A2,B2,C2,D2);

z = tf('z',Ts);
G_tf_post = tf(num_post,den_post,Ts);
G_inv = z^(-1)*tf(den_post,num_post,Ts);

sys_pre = ss(A1,B1,C1,D1,Ts);
sys_post = ss(A2,B2,C2,D2,Ts);
zero(sys_pre)
zero(sys_post)

opts = bodeoptions;
opts.FreqUnits = 'Hz';

figure
bodeplot(sys_pre,opts)
hold on
bodeplot(sys_post,opts)
legend('Pre-processing','Post-processing')

%% validating this simulation
[y_lsim_post t_lsim_post] = lsim(sys_post,u_filt,t_filt);
[y_lsim_pre t_lsim_pre] = lsim(sys_pre,u_sim,t_sim);

figure
plot(t_sim,y_sim,t_filt,y_filt,t_lsim_pre,y_lsim_pre)
legend('Unprocessed','Processed','SysID-Pre')

figure
plot(t_sim,y_sim,t_filt,y_filt,t_lsim_post,y_lsim_post)
legend('Unprocessed','Processed','SysID-Post')

%% simulink feedfoward
t_lsim_post = 0:Ts:T_fin;
% Amp = 2;
% u_in = sin(2*pi*f_in*t_lsim);
y_des = 200*sin(2*pi*f_in*t_lsim_post);
figure
plot(y_des)
[u_des t_lsim] = lsim(G_inv,y_des,t_lsim_post);
figure
plot(t_lsim,u_des)
title('Desired Input')

data_out = specCal(u_sim,Fs);
y = data_out.amp;
f_lin = data_out.f;

figure
semilogy(f_lin,y)