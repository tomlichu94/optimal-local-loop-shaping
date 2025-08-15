clc
close all
clear
%% model parameters
% DC motor from MinSegMotor
% sampled at 250 Hz
% input is 3.9V with a sin wave at 8 Hz.
% Fast sampling here defined at 250 Hz, slow sampling at 50 Hz
addpath('Functions')

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
close all
addpath('Experimental Runs\Fractional Recovery')
y_encoder = squeeze(out_encoder.signals.values)';
t_encoder = out_encoder.time';
t_fs = out_W.time';
t_ss = squeeze(in_W.time)';
d_ss = squeeze(in_W.signals.values)'; % slow sampled signal
peak_mean_fs = mean(findpeaks(y_encoder(2:end))); % find average peak value
peak_mean_ss = mean(findpeaks(d_ss(2:end)));
y_norm_fs = y_encoder/peak_mean_fs;
y_norm_ss = d_ss/peak_mean_ss;
y_w = out_W.signals.values'/peak_mean_fs; % normalized