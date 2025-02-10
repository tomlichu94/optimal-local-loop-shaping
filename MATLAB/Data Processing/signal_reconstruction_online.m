clc
close all
clear all
%% model parameters
% DC motor from MinSegMotor
% sampled at 250 Hz
% input is 3.9V with a sin wave at 8 Hz.
% noise is introduced into the "clean" measured data
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
addpath('Experimental Runs')
load run_8.mat
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


figure
s = stairs(t_encoder,y_norm_fs);
s.Color = [0.4 0.4 0.4];
s.LineWidth = 1.3;
hold on
s = stairs(t_ss,y_norm_ss);
s.LineWidth = 1.3;
s.Color = [0 0 0.65];
s.Marker = '*';

s = stairs(t_fs,y_w(1,:));
s.LineWidth = 1.4;
s.LineStyle = '-.';
s.Marker = 'x';
s.MarkerSize = 8;
s.Color = [0.9290 0.6940 0.1250];

s = stairs(t_fs,y_w(2,:));
s.LineWidth = 1.4;
s.Marker = 'o';
s.MarkerSize = 7;
s.Color = [1 0 0];
s.LineStyle = ':';
legend('Fast Sampled Signal','Slow Sampled Signal','FIR MMP','IIR MMP')
hold off
ylabel('Normalized Enconder Count')
xlabel('Time (sec)')
xlim([3.45 3.75])
ylim([-2 1.3])

% find the rms error after 1 second
idx_err = 101;
y_err = abs(y_norm_fs(idx_err:end)-y_w(:,idx_err:end));
y_rms = rms(y_err,2);

%% test with built-in function, built-in function is two steps ahead
[dest_fir dest_iir] = signal_recovery(w_k_fir,w_k_iir,B_para,L_t,y_norm_ss);
figure
stairs(t_fs,y_w(1,:))
hold on
stairs(t_fs,dest_fir)
xlim([3 5])
legend('Simulink FIR','MATLAB FIR')
hold off

figure
stairs(t_fs,y_w(2,:))
hold on
stairs(t_fs,dest_iir)
xlim([3 5])
legend('Simulink IIR','MATLAB IIR')
hold off

figure
stairs(t_fs,y_norm_fs)
hold on
stairs(t_fs,dest_iir)
stairs(t_fs,y_w(2,:))
xlim([3 5])
legend('Simulink Enconder','MATLAB IIR','Simulink IIR')
hold off