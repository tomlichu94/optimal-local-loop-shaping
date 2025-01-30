clear all
close all
clc

% run this master file to get SID and signal reconstruction using multirate
% FIR signal estimation.
%% SID
addpath('C:\Users\tpjch\OneDrive\Documents\MATLAB\RASPlib\Multirate Sampling/Data for Signal Reconstruction/PRBS')
load sys_ID_data.mat;
T_fs = 1/500; % sampled time
F_fs = 1/T_fs;
f_in = 8; % input in Hz
% data collected using quad encoder
t_sim = in_voltage.time;
u_sim = squeeze(in_voltage.signals.values);
y_sim = squeeze(double(out_encoder.signals.values));
T_fin = 10;
y_filt = zero_phase_low_pass(y_sim); % apply zero-phase error low pass filter
y_filt = y_filt';
t_filt = t_sim(1:end-1);
u_filt = zero_phase_low_pass(u_sim);
u_filt = u_filt';
nx = 1;
sys_post = n4sid(u_sim,y_sim,nx,'Ts',T_fs);
A = sys_post.A;
B = sys_post.B;
C = sys_post.C;
D = sys_post.D;
[num_post den_post] = ss2tf(A,B,C,D);

z = tf('z',T_fs);
G_tf_post = tf(num_post,den_post,T_fs);
G_inv = z^(-1)*tf(den_post,num_post,T_fs);

% run this file when before running Simulink
%% model parameters
% DC motor from MinSegMotor
% sampled at 250 Hz
% input is 3.9V with a sin wave at 8 Hz.
% noise is introduced into the "clean" measured data
% Fast sampling here defined at 250 Hz, slow sampling at 50 Hz

a_g = 0.9;
L_t = 50;
T_ss = T_fs*L_t;
[w_k_IIR B_para] = W_coeff_IIR(L_t,f_in,a_g,T_fs);
[w_k] = W_coeff_FIR(L_t,f_in,T_fs);

%% run this after running the experimentexporting data

y_encoder = squeeze(double(out_encoder.signals.values));
t_encoder = out_encoder.time;
y_w_FIR = out_W_FIR.signals.values;
t_w = out_W_FIR.time;
u_w_FIR = squeeze(in_W_FIR.signals.values);
t_u = squeeze(in_W_FIR.time);

figure
stairs(t_encoder,y_encoder)
hold on
stairs(t_w,y_w_FIR)
stairs(t_u,u_w_FIR)
xlim([3 5])
legend('Encoder','Reconstructed - FIR','Slow-Sampled')
