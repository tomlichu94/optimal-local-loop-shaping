close all
clear
%% model parameters
% DC motor from MinSegMotor
% sampled at 250 Hz
% input is 3.9V with a sin wave at 8 Hz.
% Fast sampling here defined at 250 Hz, slow sampling at 50 Hz
clc
addpath('Functions')
Fs = 100;
T_fs = 1/Fs;
T_fin = 10;
a_g = 0.9;
L = 25/2;
R = 2;
RL = R*L;
T_ss = T_fs*L;
T_cs = T_fs/R;
f_ny = 1/(2*T_ss);
f_d = 7;
f_in = 7;
[w_k_iir B_para] = w_kiir_frac(f_d, T_fs, a_g, R, L);
[w_k_fir] = w_kfir_frac(f_d, T_fs, R, L);

%% run hardware, then export the data exporting data
close all
addpath('Experimental Runs\Fractional Recovery')
y_encoder = squeeze(out_encoder.signals.values)';
t_encoder = out_encoder.time';
t_lin_cs = out_W.time';
t_ss = squeeze(in_W.time)';
d_ss = squeeze(in_W.signals.values)'; % slow sampled signal
peak_mean_fs = mean(findpeaks(y_encoder(2:end))); % find average peak value
peak_mean_ss = mean(findpeaks(d_ss(2:end)));
y_norm_cs = y_encoder/peak_mean_fs;
y_norm_ss = d_ss/peak_mean_ss;
y_w1 = out_W1.signals.values'/peak_mean_fs; % normalized
y_w2 = out_W2.signals.values'/peak_mean_fs; % normalized

t1 = t_lin_cs(1:R:end);
t2 = t_lin_cs(2:R:end);
y1 = y_w1(1:R:end);
y2 = y_w2(2:R:end);

y_cs = zeros(size(t_lin_cs));
y_cs(1:R:end) = y1;
y_cs(2:R:end) = y2;

figure
stairs(t_lin_cs,y_norm_cs)
hold on
stairs(t1,y1,'x')
stairs(t2,y2,'o')
legend('Ground Truth','Recovery 1','Recovery 2')

figure
stairs(t_lin_cs,y_norm_cs)
hold on
stairs(t_lin_cs,y_cs,':x')
legend('Ground Truth','Combined Sampling')
