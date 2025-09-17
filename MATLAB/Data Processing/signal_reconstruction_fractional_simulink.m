close all
clear
clc
%% model parameters
% DC motor from MinSegMotor
% sampled at 250 Hz
% input is 3.9V with a sin wave at 8 Hz.
% Fast sampling here defined at 250 Hz, slow sampling at 50 Hz
addpath('Functions')
Fs = 100;
T_fs = 1/Fs;
T_fin = 20;
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

%% online recovery, no noise. Run hardware, then export the data exporting data
addpath('Experimental Runs\Fractional Recovery')
load run_1.mat
y_encoder = squeeze(out_encoder.signals.values)';
y_encoder = double(y_encoder);
t_encoder = out_encoder.time';
t_lin_cs = out_W1.time';
t_ss = squeeze(in_W.time)';
d_ss = squeeze(in_W.signals.values)'; % slow sampled signal
peak_mean_fs = mean(findpeaks(y_encoder(2:end))); % find average peak value
peak_mean_ss = mean(findpeaks(d_ss(2:end)));
y_norm_fs = y_encoder/peak_mean_fs;
y_norm_ss = d_ss/peak_mean_ss;
y_w1 = out_W1.signals.values'/peak_mean_fs; % normalized
y_w2 = out_W2.signals.values'/peak_mean_fs; % normalized

%% plotting
close all

% recovery time
t1 = t_lin_cs(1:R:end); % recovery 1
t2 = t_lin_cs(2:R:end); % recovery 2

% recovery signal
y_fir_1 = y_w1(1,1:R:end); % FIR recovery 1
y_fir_2 = y_w2(1,2:R:end); % FIR recovyer 2
y_iir_1 = y_w1(2,1:R:end); % IIR recovery 1
y_iir_2 = y_w2(2,2:R:end); % IIR recovery 2

% combined sampling of two recovery
y_fir_cs = zeros(1, length(y_w1));
y_iir_cs = zeros(1, length(y_w1));
y_fir_cs(1:R:end) = y_fir_1; % fir combined sampling
y_fir_cs(2:R:end) = y_fir_2;
y_iir_cs(1:R:end) = y_iir_1; % iir combined sampling
y_iir_cs(2:R:end) = y_iir_2; 

% plot of fast/slow sample with recovery 1 and 2
x_lim = [16.12, 16.37];
y_lim = [-2, 1.5];
figure
s = stairs(t_lin_cs,y_norm_fs);
s.Color = [0.4 0.4 0.4];
s.LineWidth = 1;
hold on
s = stairs(t_ss,y_norm_ss);
s.LineWidth = 1;
s.Color = [0 0 0.65];
s.Marker = '*';
s.MarkerSize = 7;
s = stairs(t1,y_fir_1,'x');
s.LineWidth = 0.8;
s.MarkerSize = 8;
s.Color = [0.9290 0.6940 0.1250];
s = stairs(t2,y_fir_2,'o');
s.LineWidth = 0.8;
s.MarkerSize = 6;
s.Color = [1 0 0];
legend('Fast-Sampled Signal','Slow-Sampled Signal','Recovery #1','Recovery #2','location','Best')
title('Fractional Signal Recovery - FIR')
xlim(x_lim)
ylim(y_lim)

figure
s = stairs(t_lin_cs,y_norm_fs);
s.Color = [0.4 0.4 0.4];
s.LineWidth = 1;
hold on
s = stairs(t_ss,y_norm_ss);
s.LineWidth = 1;
s.Color = [0 0 0.65];
s.Marker = '*';
s.MarkerSize = 7;
s = stairs(t1,y_iir_1,'x');
s.LineWidth = 0.8;
s.MarkerSize = 8;
s.Color = [0.9290 0.6940 0.1250];
s = stairs(t2,y_iir_2,'o');
s.LineWidth = 0.8;
s.MarkerSize = 6;
s.Color = [1 0 0];
legend('Fast-Sampled Signal','Slow-Sampled Signal','Recovery #1','Recovery #2','location','Best')
title('Fractional Signal Recovery - IIR')
xlim(x_lim)
ylim(y_lim)

% plot of fast/slow sample with combined sampling
figure
s = stairs(t_lin_cs,y_norm_fs);
s.Color = [0.4 0.4 0.4];
s.LineWidth = 1;
hold on
s = stairs(t_ss,y_norm_ss);
s.LineWidth = 1;
s.Color = [0 0 0.65];
s.Marker = '*';
s.MarkerSize = 7;
s = stairs(t_lin_cs,y_fir_cs);
s.LineWidth = 0.8;
s.LineStyle = '-.';
s.Marker = 'x';
s.MarkerSize = 8;
s.Color = [0.9290 0.6940 0.1250];
s = stairs(t_lin_cs,y_iir_cs);
s.LineWidth = 0.8;
s.Marker = 'o';
s.MarkerSize = 6;
s.Color = [1 0 0];
s.LineStyle = ':';
legend('Fast Sampled Signal','Slow Sampled Signal','FIR MMP','IIR MMP','location','best')
hold off
ylabel('Normalized Enconder Count')
xlabel('Time (sec)')
xlim(x_lim)
ylim(y_lim)

%% offline recovery, with noise, requires data from previous runs
% post processing data
snr = 5; % signal to noise ratio
y_fs_noisy = awgn(y_norm_fs,snr,'measured'); % create noise for data
y_ss_noisy = y_fs_noisy(1:RL:end); % slow sampled data
y_fir_noisy = multi_phase_recovery_fir(y_ss_noisy, f_in, T_fs, T_fin, R, L);
y_iir_noisy = multi_phase_recovery_iir(y_ss_noisy, f_in, T_fs, T_fin, a_g, R, L);

% plot comparison of true encoder, noisy encoder, FIR, and IIR recovery
figure
s = plot(t_encoder, y_norm_fs); % plot of true encoder
s.Color = [0.4 0.4 0.4];
s.LineWidth = 1;
hold on
s = stairs(t_encoder, y_fs_noisy); % plot of noisy encoder
s.LineWidth = 1;
s.Color = [0 0 0.65];
s = stairs(y_fir_noisy(1,:), y_fir_noisy(2,:)); % plot of FIR-MMP
s.LineWidth = 1;
s.LineStyle = '-.';
s.Marker = 'x';
s.MarkerSize = 8;
s.Color = [0.9290 0.6940 0.1250];
s = stairs(y_iir_noisy(1,:), y_iir_noisy(2,:)); % plot of IIR-MMP
s.LineWidth = 1;
s.Marker = 'o';
s.MarkerSize = 7;
s.Color = [1 0 0];
s.LineStyle = ':';
legend('True Signal','Noisy Signal','FIR MMP','IIR MMP','Location','southeast')
ax = gca;
ax.FontSize= 12;
xlim(x_lim)
ylim([-2.7, 1.7])
xlabel('Time (sec)')
ylabel('Normalized Encoder Count')
hold off

% calculating RMS of noisy data
% plotting RMS
% remove matched measurements (i.e. d_fs[nL] = d_ss[n])
% finding absolute error
err_FIR = y_fir_noisy(2,:) - y_norm_fs;
err_IIR = y_iir_noisy(2,:) - y_norm_fs;

% plotting error
figure
s = stairs(t_lin_cs,err_FIR);
s.LineWidth = 0.7;
s.LineStyle = '-';
hold on
s = stairs(t_lin_cs, err_IIR);
s.LineWidth = 1;
s.Color = [1 0 0];
s.LineStyle = ':';
hold off
legend('FIR MMP','IIR MMP','Location','northeast')
ax = gca;
ax.FontSize= 12;
ylabel('Error')
xlabel('Time (sec)')
ylim([-2 2])