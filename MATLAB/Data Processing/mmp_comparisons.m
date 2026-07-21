clear all
close all
clc

addpath('Functions')
%% Sampling
L = 5;
f_fs = 120;
f_ss = f_fs/L;
t_fs = 1/f_fs;
t_ss = 1/f_ss;

ct_idx = 100;
t_cs = t_fs/ct_idx;

t_end = 20;

f_sample = [23, 29, 33];
a_g = 0.9;
mem_length = length(f_sample)*4; % need this many for mmp
% lambda = 1-1/mem_length;
lambda = 0.9;
mem_length = 1/(1-lambda);
snr = 5;
delta = 1;

% Number of slow-sampling intervals
N_ss = floor(t_end/t_ss) + 1;

% Corresponding number of fast and continuous samples
N_fs = (N_ss-1)*L + 1;
N_cs = (N_fs-1)*ct_idx + 1;

% Time vectors
t_lin_ss = (0:N_ss-1)*t_ss;
t_lin_fs = (0:N_fs-1)*t_fs;
t_lin_cs = (0:N_cs-1)*t_cs;

A_func = {
    @(t) 1.0 + 0.3*sin(2*pi*0.2*t)
    @(t) 0.8 + 0.4*cos(2*pi*0.4*t)
    @(t) 0.5 + 0.1*sin(2*pi*0.1*t)
};

% A_func = [1,1,1];

phase = [0, pi/4, pi/2];

[y a] = variable_amplitude_sine(t_lin_cs, f_sample, A_func, phase);
y_fs = y(1:ct_idx:end);
y_fs_noisy = awgn(y_fs,snr,'measured'); % create noise for data
% y_ss = y(1:L*ct_idx:end); % clean data
y_ss = y_fs_noisy(1:L:end); % noisy data


% using MMP
y_mmp = signal_recovery_iir(y_ss, f_sample, t_fs, a_g, L);
N = min([length(y_fs), length(y_mmp), length(t_lin_fs)]);
y_fs     = y_fs(1:N);
y_mmp    = y_mmp(1:N);
t_lin_fs = t_lin_fs(1:N);
err = y_fs - y_mmp;

figure
stairs(t_lin_fs,y_fs)
hold on
stairs(t_lin_fs,y_mmp)
legend('Ground Truth', 'MMP')

figure
stairs(t_lin_fs, err)
title('Residuals')

% using RLS
[y_rls, theta_rls] = rls_sine_reconstruction(y_ss, t_lin_ss, t_lin_fs, f_sample, lambda, delta);

N = min([length(t_lin_fs), length(y_fs), length(y_mmp), length(y_rls)]);

t_compare = t_lin_fs(1:N);
y_true    = y_fs(1:N);
y_mmp     = y_mmp(1:N);
y_rls     = y_rls(1:N);

err_mmp = y_true - y_mmp;
err_rls = y_true - y_rls;

figure
plot(t_compare, y_true, 'LineWidth', 1.2)
hold on
plot(t_compare, y_mmp, '--', 'LineWidth', 1.2)
plot(t_compare, y_rls, ':', 'LineWidth', 1.4)
grid on

xlabel('Time (s)')
ylabel('Signal')
legend('Ground truth', 'IIR-MMP', 'RLS', ...
       'Location', 'best')
title('Signal Reconstruction Comparison')

figure
plot(t_compare, err_mmp, 'LineWidth', 1.2)
hold on
plot(t_compare, err_rls, 'LineWidth', 1.2)
grid on

xlabel('Time (s)')
ylabel('Reconstruction error')
legend('IIR-MMP error', 'RLS error', ...
       'Location', 'best')
title('Reconstruction Residuals')

t_transient = 1;
idx_eval = t_compare >= t_transient;

rmse_mmp = sqrt(mean(err_mmp(idx_eval).^2));
rmse_rls = sqrt(mean(err_rls(idx_eval).^2));

nrmse_mmp = rmse_mmp / rms(y_true(idx_eval));
nrmse_rls = rmse_rls / rms(y_true(idx_eval));

fprintf('MMP RMSE:       %.6f\n', rmse_mmp);
fprintf('RLS RMSE:       %.6f\n', rmse_rls);
fprintf('MMP NRMSE:      %.2f %%\n', 100*nrmse_mmp);
fprintf('RLS NRMSE:      %.2f %%\n', 100*nrmse_rls);