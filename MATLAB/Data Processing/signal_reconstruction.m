clc
close all
clear all
%% model parameters
% DC motor from MinSegMotor
% sampled at 250 Hz
% input is 3.9V with a sin wave at 8 Hz.
% noise is introduced into the "clean" measured data
% Fast sampling here defined at 250 Hz, slow sampling at 50 Hz

Fs = 100;
T_fs = 1/Fs;
T_fin = 10;
a_g = 0.9;
L_t = 10;
T_ss = T_fs*L_t;
f_d = 8;
f_in = 8;
[w_k_iir B_para] = W_coeff_IIR(L_t,f_d,a_g,T_fs);
[w_k_fir] = W_coeff_FIR(L_t,f_d,T_fs);

%% run hardware, then export the data exporting data
addpath('Experimental Runs')
load run_7.mat
y_encoder = squeeze(out_encoder.signals.values)';
t_encoder = out_encoder.time';
d_ss = squeeze(in_W.signals.values)';
t_ss = in_W.time';
y_w = out_W.signals.values';
t_w = out_W.time';
u_w_FIR = squeeze(in_W.signals.values)';
t_u = squeeze(in_W.time)';
figure
stairs(t_encoder,y_encoder)
hold on
stairs(t_w,y_w(1,:))
stairs(t_w,y_w(2,:))
xlim([3 5])
legend('Encoder','FIR','IIR')
hold off

%% test with built-in function
close all
addpath('Functions')
[dest_fir dest_iir] = signal_recovery(w_k_fir,w_k_iir,B_para,L_t,d_ss);

figure
stairs(t_w,y_w(1,:))
hold on
stairs(t_w,dest_fir)
xlim([3 5])
legend('Simulink FIR','MATLAB FIR')
hold off

figure
stairs(t_w,y_w(2,:))
hold on
stairs(t_w,dest_iir)
xlim([3 5])
legend('Simulink IIR','MATLAB IIR')
hold off

figure
stairs(t_w,y_encoder)
hold on
stairs(t_w,dest_iir)
stairs(t_w,y_w(2,:))
xlim([3 5])
legend('Simulink Enconder','MATLAB IIR','Simulink IIR')
hold off
