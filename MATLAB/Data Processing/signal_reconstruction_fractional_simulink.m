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
a_g = 0.5;
L = 25/2;
[N_L, D_L] = rat(L);
T_ss = T_fs*L;
T_cs = T_fs/D_L;
f_ny = 1/(2*T_ss);
f_d = 7;
f_in = 7;
[wk_iir B_para] = w_kiir_frac(f_d, T_fs, a_g, L);
[wk_fir] = w_kfir_frac(f_d, T_fs, L);

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

% plot of fast/slow sample with combined sampling, FIR and IIR
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
ax.FontSize= 12;s.Color = [1 0 0];
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
y_ss_noisy = y_fs_noisy(1:N_L:end); % slow sampled data
y_fir_noisy = multi_phase_recovery_fir(y_ss_noisy, f_in, T_fs, T_fin, L);
a_g = [0.3, 0.5, 0.9]; % alpha to adjust bandwidth of IIR-MMP
length_cs = length(y_fir_noisy);
y_iir_noisy = zeros(length(a_g),2,length_cs);
for i = 1:length(a_g)
    y_iir_noisy(i,:,:) = multi_phase_recovery_iir(y_ss_noisy, f_in, T_fs, T_fin, a_g(i), L);
end
y_iir_noisy = squeeze(y_iir_noisy(:,2,:));

%% plotting online signal recovery
% plot comparison of true encoder, noisy encoder, FIR, and IIR recovery
figure
s = stairs(t_encoder, y_norm_fs); % plot of true encoder
s.Color = [0.4 0.4 0.4];
s.LineWidth = 1;
hold on
s = stairs(t_ss, y_ss_noisy); % plot of noisy encoder
s.Marker = '*';
s.LineWidth = 1;
s.Color = [0 0 0.65];
s = stairs(y_fir_noisy(1,:), y_fir_noisy(2,:)); % plot of FIR-MMP
s.LineWidth = 1;
s.LineStyle = '-.';
s.Marker = 'x';
s.MarkerSize = 8;
s.Color = [0.9290 0.6940 0.1250];
s = stairs(t_lin_cs, y_iir_noisy(3,:)); % plot of IIR-MMP
s.LineWidth = 1;
s.Marker = 'o';
s.MarkerSize = 7;
s.Color = [1 0 0];
s.LineStyle = ':';
legend('True Signal','Noisy Signal','FIR MMP','IIR MMP','Location','southeast')
ax = gca;
ax.FontSize= 12;
xlim(x_lim)
ylim([-2.3, 1.5])
xlabel('Time (sec)')
ylabel('Normalized Encoder Count')
hold off

% calculating RMS of noisy data
% plotting RMS
% remove matched measurements (i.e. d_fs[nL] = d_ss[n])
% finding absolute error
err_FIR = y_fir_noisy(2,:) - y_norm_fs;
err_IIR = y_iir_noisy - y_norm_fs;
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
legend('FIR MMP','IIR MMP','Location','northeast')
ax = gca;
ax.FontSize= 12;
ylabel('Error')
xlabel('Time (sec)')
ylim([-2 2])

% box chart error and RMS error
rms_FIR = rms(err_FIR);
rms_IIR = rms(err_IIR,2);
err_all = [err_FIR;err_IIR]';
rms_all = [rms_FIR;rms_IIR]';
rms_FIR_clean = rms(err_FIR_clean);
rms_IIR_clean = rms(err_IIR_clean);
figure
boxchart(err_all)
xaxisproperties= get(gca, 'XAxis');
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
h(1) = plot(w_in_Hz,mag_W_FIR);
hold on
h(1).Color = color_all{1};
h(1).LineStyle = line_style{1};
h(1).LineWidth = 1.5;
xlabel('Hz')
ylabel('Magnitude (dB)')
xlim(x_lim_PQ)
wk_iir = zeros(length(a_g),size(wk_fir,1),size(wk_fir,2));

% IIR bode plot options
for i = 1:length(a_g)
    [wk_iir(i,:,:), B_para] = w_kiir_frac(f_in, T_fs, a_g(i), L); % iir coeff
    [W_IIR_num, W_IIR_den] = w_tf_iir(squeeze(wk_iir(i,:,:)),B_para);
    W_k_IIR = tf(W_IIR_num(k_idx,:),W_IIR_den(k_idx,:),T_fs); % IIR TF for k=1
    [mag_W_IIR, phi_W_IIR, w_out] = bode(W_k_IIR,w_in_rad); % IIR mag/deg vals
    mag_W_IIR = 20*log10(mag_W_IIR(:));
    h(i+1) = plot(w_in_Hz,mag_W_IIR);
    h(i+1).Color = color_all{i+1};
    h(i+1).LineStyle = '-';
    h(i+1).LineWidth = 1.5;
end
Legend_bode = cell(1,5);
Legend_bode{1} = strcat('FIR');
for iter=1:length(a_g)
  Legend_bode{iter+1}=strcat('\alpha =', num2str(a_g(iter)));
end
Legend_bode{5}='k = 2,...24';
% plotting IIR-MMP Bode for k = 2...L-1
k_length = k_tot(end);
% plotting intersamples of IIR-MMP for a = 0.9
[W_IIR_num, W_IIR_den] = w_tf_iir(squeeze(wk_iir(3,:,:)),B_para); % IIR TF coeff
for k = 1:k_length
        W_k_IIR(k) = tf(W_IIR_num(k,:),W_IIR_den(k,:),T_fs);
        [mag_W_IIR, phi_W_IIR, w_out] = bode(W_k_IIR(k),w_in_rad);
        mag_W_IIR = 20*log10(mag_W_IIR(:));
        phi_W_IIR = wrapTo180(phi_W_IIR(:));
        h(k) = plot(w_in_Hz,mag_W_IIR);
        h(k).Color = [0.3 0.3 0.3];
        h(k).LineStyle = ':';
        h(k).LineWidth = 1.15;
        xlabel('Hz')
        ylabel('Magnitude (dB)')
        xlim(x_lim_PQ)
end
% plot disturbance frequency line
for i = 1:width(f_d)
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
a = findobj(gca,'type','line');
for i = 1:length(a)
    set(a(i),'markersize',8)
    set(a(i),'linewidth',1.5)
end
hold on
for i = 1:length(a_g)
    [W_IIR_num, W_IIR_den] = w_tf_iir(squeeze(wk_iir(i,:,:)),B_para);
    for k = 1:k_length
    W_k_IIR(i,k) = tf(W_IIR_num(k,:),1,T_fs); % IIR-MMP TF zeros only
    end
end

% range of zeros to plot
iter_zeros = 3:9:k_length;
for k = iter_zeros
    pzmap(W_k_IIR(3,k)) % plot zeros
    for i = 1:length(a)
    a = findobj(gca,'type','line');
    set(a(i),'markersize',8)
    set(a(i),'linewidth',1.5)
    end
end
Legend_pz = cell(1,length(iter_zeros)+1);
Legend_pz{1} = strcat('Poles (\alpha = ',num2str(a_g(3)),')');
n_leg = 2;
for iter=iter_zeros
  Legend_pz{n_leg}=strcat('k = ', num2str(k_tot(iter)));
  n_leg = n_leg+1;
end
legend(Legend_pz,'Location','Best')

% pzmap plot of FIR and IIR poles only
figure
[wk_iir, B_para] = w_kiir_frac(f_in, T_fs, a_g(1), L); % coefficients for IIR-MMP
[W_IIR_num, W_IIR_den] = w_tf_iir(wk_iir,B_para);
[~, W_IIR_den] = w_tf_iir(wk_iir,B_para);
W_k_IIR = tf(1,W_IIR_den(1,:),T_fs); % IIR-MMP TF poles only
pzmap(W_k_IIR,'k') % map of IIR-MMP poles
a = findobj(gca,'type','line');
    for i = 1:length(a)
        set(a(i),'markersize',8)
        set(a(i),'linewidth',1.5)
    end
hold on
[wk_iir, B_para] = w_kiir_frac(f_in, T_fs, a_g(2), L); % coefficients for IIR-MMP
[W_IIR_num, W_IIR_den] = w_tf_iir(wk_iir,B_para);
[~, W_IIR_den] = w_tf_iir(wk_iir,B_para);
W_k_IIR = tf(1,W_IIR_den(1,:),T_fs); % IIR-MMP TF poles only
pzmap(W_k_IIR,'b') % map of IIR-MMP poles
a = findobj(gca,'type','line');
    for i = 1:length(a)
        set(a(i),'markersize',8)
        set(a(i),'linewidth',1.5)
    end
[wk_iir, B_para] = w_kiir_frac(f_in, T_fs, a_g(3), L); % coefficients for IIR-MMP
[W_IIR_num, W_IIR_den] = w_tf_iir(wk_iir,B_para);
W_k_IIR = tf(1,W_IIR_den(1,:),T_fs); % IIR-MMP TF poles only
pzmap(W_k_IIR,'r') % map of IIR-MMP poles
a = findobj(gca,'type','line');
    for i = 1:length(a)
        set(a(i),'markersize',8)
        set(a(i),'linewidth',1.5)
    end
Legend_pz = cell(1,3);
for iter=1:3
Legend_pz{iter} = strcat('Poles (\alpha = ',num2str(a_g(iter)),')');
end
legend(Legend_pz,'Location','Best')