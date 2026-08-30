%% Fractional signal reconstruction from Simulink data
close all;
clearvars;
clc;

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'Functions'));

%% Model parameters
% DC motor data sampled at 100 Hz with a 7 Hz sinusoidal input.
Fs = 100;
T_fs = 1 / Fs;
T_fin = 20;
L = 25 / 2;
[N_L, D_L] = rat(L);
T_ss = T_fs * L;
f_ny = 1 / (2 * T_ss);
f_d = 7;
f_in = 7;

% The recorded model contains two interleaved recovery channels.
assert(D_L == 2, ...
    'This script expects a reconstruction factor with denominator 2.');

%% Online recovery without noise
run_file = fullfile(script_dir, 'Experimental Runs', ...
    'Fractional Recovery', 'run_1.mat');
assert(isfile(run_file), 'Data file not found: %s', run_file);
run_data = load(run_file);

required_fields = {'out_encoder', 'out_W1', 'out_W2', 'in_W'};
assert(all(isfield(run_data, required_fields)), ...
    'The run file does not contain all required Simulink outputs.');

out_encoder = run_data.out_encoder;
out_W1 = run_data.out_W1;
out_W2 = run_data.out_W2;
in_W = run_data.in_W;

y_encoder = double(squeeze(out_encoder.signals.values)).';
t_encoder = squeeze(out_encoder.time).';
t_lin_cs = squeeze(out_W1.time).';
t_ss = squeeze(in_W.time).';
d_ss = squeeze(in_W.signals.values).';

encoder_peaks = findpeaks(y_encoder(2:end));
slow_peaks = findpeaks(d_ss(2:end));
assert(~isempty(encoder_peaks) && ~isempty(slow_peaks), ...
    'Peak-based normalization requires peaks in both input signals.');
peak_mean_fs = mean(encoder_peaks);
peak_mean_ss = mean(slow_peaks);
assert(all(isfinite([peak_mean_fs, peak_mean_ss])) && ...
    all([peak_mean_fs, peak_mean_ss] ~= 0), ...
    'Peak normalization values must be finite and nonzero.');

y_norm_fs = y_encoder/peak_mean_fs;
y_norm_ss = d_ss/peak_mean_ss;
y_w1 = squeeze(out_W1.signals.values).'/peak_mean_fs;
y_w2 = squeeze(out_W2.signals.values).'/peak_mean_fs;

assert(numel(t_lin_cs) == numel(y_norm_fs), ...
    'The combined time vector and encoder signal must have equal lengths.');
assert(size(y_w1, 1) >= 2 && size(y_w2, 1) >= 2, ...
    'Each recovery output must contain FIR and IIR channels.');

%% plotting
close all

% recovery time
t1 = t_lin_cs(1:D_L:end); % recovery 1
t2 = t_lin_cs(2:D_L:end); % recovery 2

% recovery signal
y_fir_1 = y_w1(1,1:D_L:end); % FIR recovery 1
y_fir_2 = y_w2(1,2:D_L:end); % FIR recovyer 2
y_iir_1 = y_w1(2,1:D_L:end); % IIR recovery 1
y_iir_2 = y_w2(2,2:D_L:end); % IIR recovery 2

% combined sampling of two recovery
y_fir_cs = zeros(1, length(y_w1));
y_iir_cs = zeros(1, length(y_w1));
y_fir_cs(1:D_L:end) = y_fir_1; % fir combined sampling
y_fir_cs(2:D_L:end) = y_fir_2;
y_iir_cs(1:D_L:end) = y_iir_1; % iir combined sampling
y_iir_cs(2:D_L:end) = y_iir_2; 

% plot of fast/slow sample with recovery 1 and 2, FIR
x_lim = [16.1, 16.2];
y_lim = [-1.5, 1.5];
figure
s = stairs(t_lin_cs,y_norm_fs);
s.Color = [0.4 0.4 0.4];
s.LineWidth = 1;
hold on
s = stairs(t_ss,y_norm_ss);
s.LineWidth = 1;
s.Color = [0 0 0.65];
s.Marker = '*';
s.MarkerSize = 9;
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

% plot of fast/slow sample with recovery 1 and 2, IIR
figure
s = stairs(t_lin_cs,y_norm_fs);
s.Color = [0.4 0.4 0.4];
s.LineWidth = 1;
hold on
s = stairs(t_ss,y_norm_ss);
s.LineWidth = 1;
s.Color = [0 0 0.65];
s.Marker = '*';
s.MarkerSize = 9;
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

%% Offline recovery for fixed SNR and varying alpha
fixed_snr_db = 5;
y_fs_noisy_fixed = awgn(y_norm_fs, fixed_snr_db, 'measured');
y_ss_noisy_fixed = y_fs_noisy_fixed(1:N_L:end);
y_fir_noisy = multi_phase_recovery_fir(y_ss_noisy_fixed, f_in, T_fs, T_fin, L);
alpha_values = [0.3, 0.5, 0.9];
y_iir_noisy_a_g = zeros(numel(alpha_values), size(y_fir_noisy, 2));
for alpha_idx = 1:numel(alpha_values)
    recovery = multi_phase_recovery_iir(y_ss_noisy_fixed, f_in, ...
        T_fs, T_fin, alpha_values(alpha_idx), L);
    y_iir_noisy_a_g(alpha_idx, :) = recovery(2, :);
end

assert(size(y_fir_noisy, 2) == numel(y_norm_fs), ...
    'Offline recovery and reference signal lengths do not match.');

%% Monte Carlo reconstruction error versus SNR and alpha
snr_db_values = [0, 3, 5, 10, 15, 20, 25, 30];
Ntrial = 50;

rng('default');
rng(1,'twister');

rms_FIR_trials = zeros(Ntrial, numel(snr_db_values));
rms_IIR_trials = zeros(Ntrial, numel(snr_db_values), ...
    numel(alpha_values));

% Evaluate only reconstructed intersamples after the startup transient.
evaluation_mask = true(size(y_norm_fs));
evaluation_mask(1:N_L:end) = false;
evaluation_mask = evaluation_mask & (t_lin_cs >= 4);

for trial = 1:Ntrial
    for snr_idx = 1:numel(snr_db_values)

        y_fs_noisy = awgn(y_norm_fs, snr_db_values(snr_idx), 'measured');
        y_ss_noisy = y_fs_noisy(1:N_L:end);

        y_fir_tmp = multi_phase_recovery_fir(y_ss_noisy, f_in, T_fs, T_fin, L);
        y_fir_tmp = squeeze(y_fir_tmp(2,:));

        err_FIR = y_fir_tmp - y_norm_fs;
        rms_FIR_trials(trial, snr_idx) = rms(err_FIR(evaluation_mask));

        for alpha_idx = 1:numel(alpha_values)

            y_iir_tmp = multi_phase_recovery_iir( ...
                y_ss_noisy, f_in, T_fs, T_fin, ...
                alpha_values(alpha_idx), L);

            y_iir_tmp = squeeze(y_iir_tmp(2,:));

            err_IIR = y_iir_tmp - y_norm_fs;
            rms_IIR_trials(trial, snr_idx, alpha_idx) = ...
                rms(err_IIR(evaluation_mask));
        end
    end
end

rms_FIR_mean = mean(rms_FIR_trials,1);
rms_FIR_std  = std(rms_FIR_trials,0,1);

rms_IIR_mean = squeeze(mean(rms_IIR_trials,1));  % SNR x alpha
rms_IIR_std  = squeeze(std(rms_IIR_trials,0,1));

figure; hold on; grid on;

p = errorbar(snr_db_values, rms_FIR_mean, rms_FIR_std);
p.LineStyle = '-';
p.LineWidth = 1.3;
p.Marker = 'x';
p.MarkerSize = 4;
p.CapSize = 8;

markers = {'o','s','d'};

for alpha_idx = 1:numel(alpha_values)
    p = errorbar(snr_db_values, rms_IIR_mean(:, alpha_idx), ...
        rms_IIR_std(:, alpha_idx));
    p.LineStyle = '-';
    p.LineWidth = 1.3;
    p.Marker = markers{alpha_idx};
    p.MarkerSize = 4;
    p.MarkerFaceColor = 'w';
    p.CapSize = 8;
end

xlabel('SNR (dB)')
ylabel('RMS reconstruction error')
xlim([snr_db_values(1)-1, snr_db_values(end)+1])
legend_entries = [{'W-FIR'}, ...
    arrayfun(@(alpha) sprintf('W-IIR, \\alpha = %.1f', alpha), ...
    alpha_values, 'UniformOutput', false)];

legend(legend_entries, 'Location', 'best')

%% plotting for signal outputs
% plot of fast/slow sample with combined sampling, FIR and IIR
% first five seconds
x_lim = [0, 5];
y_lim = [-2, 2];
figure
s = stairs(t_lin_cs,y_norm_fs);
s.Color = [0.4 0.4 0.4];
s.LineWidth = 1;
hold on
s = stairs(t_ss(2:end),y_norm_ss(1:end-1));
s.LineWidth = 1;
s.Color = [0 0 0.65];
s = stairs(t_lin_cs,y_fir_cs);
s.LineWidth = 1;
s.LineStyle = '-.';
s.Color = [0.9290 0.6940 0.1250];
s = stairs(t_lin_cs,y_iir_cs);
s.LineWidth = 1;
s.LineStyle = ':';
ax = gca;
ax.FontSize= 12;
s.Color = [1 0 0];
legend('Fast sampled signal','Slow sampled signal','W-FIR','W-IIR','location','best')
hold off
ylabel('Normalized Encoder Count')
xlabel('Time (sec)')
xlim(x_lim)
ylim(y_lim)

% offline plot comparison of true encoder, noisy encoder, FIR, and IIR recovery
figure
s = stairs(t_encoder, y_norm_fs); % plot of true encoder
s.Color = [0.4 0.4 0.4];
s.LineWidth = 1.2;
hold on
s = stairs(t_ss, y_ss_noisy_fixed); % plot of slow-sampled input
s.LineWidth = 1;
s.Color = [0 0 0.65];
s = stairs(y_fir_noisy(1,:), y_fir_noisy(2,:)); % plot of FIR-MMP
s.LineWidth = 1;
s.LineStyle = '-.';
s.Color = [0.9290 0.6940 0.1250];
s = stairs(t_lin_cs, y_iir_noisy_a_g(3,:)); % plot of IIR-MMP
s.LineWidth = 1;
s.Color = [1 0 0];
s.LineStyle = ':';
legend('True fast sampled signal','Slow sampled noisy signal','W-FIR','W-IIR','Location','southeast')
ax = gca;
ax.FontSize= 12;
xlim(x_lim)
ylim(y_lim)
xlabel('Time (sec)')
ylabel('Normalized Encoder Count')
hold off

% small window during transient
x_lim = [15.7, 16.2];
y_lim = [-2, 2];
figure
s = stairs(t_lin_cs,y_norm_fs);
s.Color = [0.4 0.4 0.4];
s.LineWidth = 1;
hold on
s = stairs(t_ss(2:end),y_norm_ss(1:end-1));
s.LineWidth = 1;
s.Color = [0 0 0.65];
s.Marker = '*';
s.MarkerSize = 10;
s = stairs(t_lin_cs,y_fir_cs);
s.LineWidth = 1;
s.LineStyle = '-.';
s.Marker = 'x';
s.MarkerSize = 8;
s.Color = [0.9290 0.6940 0.1250];
s = stairs(t_lin_cs,y_iir_cs);
s.LineWidth = 1;
s.Marker = 'o';
s.MarkerSize = 7;
s.LineStyle = ':';
ax = gca;
ax.FontSize= 12;
s.Color = [1 0 0];
legend('Fast sampled signal','Slow sampled signal','W-FIR','W-IIR','location','best')
hold off
ylabel('Normalized Encoder Count')
xlabel('Time (sec)')
xlim(x_lim)
ylim(y_lim)

% offline plot comparison of true encoder, noisy encoder, FIR, and IIR recovery
figure
s = stairs(t_encoder, y_norm_fs); % plot of true encoder
s.Color = [0.4 0.4 0.4];
s.LineWidth = 1.2;
hold on
s = stairs(t_ss, y_ss_noisy_fixed); % plot of slow-sampled input
s.LineWidth = 1;
s.Color = [0 0 0.65];
s.Marker = '*';
s.MarkerSize = 10;
s = stairs(y_fir_noisy(1,:), y_fir_noisy(2,:)); % plot of FIR-MMP
s.LineWidth = 1;
s.LineStyle = '-.';
s.Marker = 'x';
s.MarkerSize = 8;
s.Color = [0.9290 0.6940 0.1250];
s = stairs(t_lin_cs, y_iir_noisy_a_g(3,:)); % plot of IIR-MMP
s.LineWidth = 1;
s.Marker = 'o';
s.MarkerSize = 7;
s.Color = [1 0 0];
s.LineStyle = ':';
legend('True fast sampled signal','Slow sampled noisy signal','W-FIR','W-IIR','Location','southeast')
ax = gca;
ax.FontSize= 12;
xlim(x_lim)
ylim(y_lim)
xlabel('Time (sec)')
ylabel('Normalized Encoder Count')
hold off

% calculating RMS of noisy data
% plotting RMS
% remove matched measurements (i.e. d_fs[nL] = d_ss[n])
% finding absolute error
err_FIR = y_fir_noisy(2,:) - y_norm_fs;
err_IIR = y_iir_noisy_a_g - y_norm_fs;
err_FIR_clean = y_fir_cs - y_norm_fs;
err_IIR_clean = y_iir_cs - y_norm_fs;
% plotting error
figure
s = stairs(t_lin_cs,err_FIR);
s.LineWidth = 0.7;
s.LineStyle = '-';
hold on
s = stairs(t_lin_cs, err_IIR(3,:));
s.LineWidth = 1;
s.Color = [1 0 0];
s.LineStyle = ':';
hold off
legend('W-FIR','W-IIR','Location','northeast')
ax = gca;
ax.FontSize= 12;
ylabel('Error')
xlabel('Time (sec)')
ylim([-2 2])

% box chart error and RMS error
idx_ss = find(t_lin_cs >= 5);

rms_FIR = rms(err_FIR(idx_ss));
rms_IIR = rms(err_IIR(:,idx_ss),2);
err_all = [err_FIR(idx_ss);err_IIR(:,idx_ss)]';
rms_all = [rms_FIR;rms_IIR]';
rms_FIR_noise_free = rms(err_FIR_clean(idx_ss));
rms_IIR_noise_free = rms(err_IIR_clean(idx_ss));
fprintf('Noise-free online RMS error (FIR): %.6f\n', rms_FIR_noise_free);
fprintf('Noise-free online RMS error (IIR): %.6f\n', rms_IIR_noise_free);
figure
boxchart(err_all)
xaxisproperties = get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
xticklabels({'FIR','$\alpha = 0.3$','$\alpha = 0.5$','$\alpha = 0.9$'})
hold on
p = plot(rms_all,'-o');
p.LineWidth = 1.5;
txt = cell(1,4);
for i = 1:4
    str = sprintf('%.3f',rms_all(i));
    txt{i} = cellstr(str);
end
text((1:4)+0.25,rms_all+0.1,txt);
legend('','RMS Error')
ylabel('Error')

% bode plot options and setup
wk_fir = w_kfir_frac(f_in, T_fs, L);
k_tot = 1:(N_L-1); % index of intersample points
color_all = {[0 0 0], [0.9290 0.6940 0.1250], [0 1 0], [1 0 0]};
line_style = {'-','--','--','-.',':'};
w_in_Hz_end = 1/(2*T_fs); % Nyq freq of fast sampling
w_in_Hz = logspace(-2,log10(w_in_Hz_end),10001); % range of freq (Hz)
w_in_rad = w_in_Hz*2*pi; % convert to rad/s
font_size = 12; % font size
x_lim_PQ = [1 w_in_Hz(end)]; % PQ x_lim
y_lim = [-50 15]; % y_lim
k_idx = 1; % pick a transfer function from the intersample

% calc bode plot values
[W_FIR_num, W_FIR_den] = w_tf_fir(wk_fir); % coeff for TF
W_k_FIR = tf(W_FIR_num(k_idx,:),W_FIR_den(k_idx,:),T_fs); % FIR-MMP TF
[mag_W_FIR, phi_W_FIR, w_out] = bode(W_k_FIR,w_in_rad); % mag and deg values
mag_W_FIR = 20*log10(mag_W_FIR(:));
phi_W_FIR = wrapTo180(phi_W_FIR(:));

% plot FIR Bode plot
figure
h = gobjects(1, numel(alpha_values) + 1);
h(1) = plot(w_in_Hz,mag_W_FIR);
hold on
h(1).Color = color_all{1};
h(1).LineStyle = line_style{1};
h(1).LineWidth = 1.5;
xlabel('Hz')
ylabel('Magnitude (dB)')
xlim(x_lim_PQ)
wk_iir = zeros(numel(alpha_values), size(wk_fir, 1), size(wk_fir, 2));
B_para_all = cell(1, numel(alpha_values));

% IIR bode plot options
for i = 1:numel(alpha_values)
    [wk_iir(i,:,:), B_para_all{i}] = w_kiir_frac( ...
        f_in, T_fs, alpha_values(i), L);
    [W_IIR_num, W_IIR_den] = w_tf_iir( ...
        squeeze(wk_iir(i,:,:)), B_para_all{i});
    W_k_IIR = tf(W_IIR_num(k_idx,:),W_IIR_den(k_idx,:),T_fs); % IIR TF for k=1
    mag_W_IIR = bode(W_k_IIR, w_in_rad);
    mag_W_IIR = 20*log10(mag_W_IIR(:));
    h(i+1) = plot(w_in_Hz,mag_W_IIR);
    h(i+1).Color = color_all{i+1};
    h(i+1).LineStyle = '-';
    h(i+1).LineWidth = 1.5;
end
Legend_bode = cell(1,5);
Legend_bode{1} = strcat('FIR');
for iter = 1:numel(alpha_values)
  Legend_bode{iter+1} = strcat('\alpha =', num2str(alpha_values(iter)));
end
Legend_bode{5}='k = 2,...24';
% plotting IIR-MMP Bode for k = 2...L-1
k_length = k_tot(end);
% plotting intersamples of IIR-MMP for a = 0.9
[W_IIR_num, W_IIR_den] = w_tf_iir( ...
    squeeze(wk_iir(3,:,:)), B_para_all{3});
for k = 2:k_length
        W_k_IIR = tf(W_IIR_num(k,:), W_IIR_den(k,:), T_fs);
        mag_W_IIR = bode(W_k_IIR, w_in_rad);
        mag_W_IIR = 20*log10(mag_W_IIR(:));
        gray_line = plot(w_in_Hz, mag_W_IIR);
        gray_line.Color = [0.3 0.3 0.3];
        gray_line.LineStyle = ':';
        gray_line.LineWidth = 1.15;
end
% plot disturbance frequency line
for i = 1:numel(f_d)
    x_line = xline_default(f_d(i));
    x_line.Color = [0 0 0];
    x_line.LineWidth = 1;
    x_line.LineStyle = '--';
end
x_line.Label = strcat(num2str(f_in),{' '},'Hz');
x_line.LabelHorizontalAlignment = 'center';
x_line.LabelVerticalAlignment = 'bottom';
x_line = xline_default(f_ny);
x_line.Label = strcat(num2str(f_ny), {' '}, 'Hz');
x_line.LabelHorizontalAlignment = 'center';
x_line.LabelVerticalAlignment = 'bottom';
x_line.LineWidth = 1;
hold off
legend(Legend_bode)
xlim([0,20])

% plotting pzmap
figure
W_k_IIR = tf(1,W_IIR_den(k_idx,:),T_fs); % IIR-MMP TF poles only
pzmap(W_k_IIR,'k') % map of IIR-MMP poles
hold on

% range of zeros to plot
iter_zeros = 3:9:k_length;
[W_IIR_num, ~] = w_tf_iir(squeeze(wk_iir(3,:,:)), B_para_all{3});
for k = iter_zeros
    zero_system = tf(W_IIR_num(k,:), 1, T_fs);
    pzmap(zero_system)
end
pole_zero_lines = findobj(gca, 'Type', 'line');
set(pole_zero_lines, 'MarkerSize', 8, 'LineWidth', 1.5);
Legend_pz = cell(1,length(iter_zeros)+1);
Legend_pz{1} = strcat('Poles (\alpha = ', ...
    num2str(alpha_values(3)), ')');
n_leg = 2;
for iter=iter_zeros
  Legend_pz{n_leg}=strcat('k = ', num2str(k_tot(iter)));
  n_leg = n_leg+1;
end
legend(Legend_pz,'Location','Best')

% pzmap plot of FIR and IIR poles only
figure
hold on
pole_colors = {'k', 'b', 'r'};
for alpha_idx = 1:numel(alpha_values)
    W_k_IIR = tf(1, W_IIR_den_for_alpha( ...
        wk_iir, B_para_all, alpha_idx, k_idx), T_fs);
    pzmap(W_k_IIR, pole_colors{alpha_idx})
end
hold off
pole_zero_lines = findobj(gca, 'Type', 'line');
set(pole_zero_lines, 'MarkerSize', 8, 'LineWidth', 1.5);
Legend_pz = arrayfun(@(alpha) sprintf('Poles (\\alpha = %.1f)', alpha), ...
    alpha_values, 'UniformOutput', false);
legend(Legend_pz,'Location','Best')

function denominator = W_IIR_den_for_alpha( ...
        wk_iir, B_para_all, alpha_idx, k_idx)
    [~, all_denominators] = w_tf_iir( ...
        squeeze(wk_iir(alpha_idx, :, :)), B_para_all{alpha_idx});
    denominator = all_denominators(k_idx, :);
end
