clc
close all
clear
%% model parameters
% DC motor from MinSegMotor
% sampled at 100 Hz
% input is 2.5V with a sin wave at 8 Hz.
% noise is introduced into the "clean" measured data
% Fast sampling here defined at 100 Hz, slow sampling at 10 Hz
addpath('Functions','Experimental Runs')
load run_8.mat

% load data from encoder
in_V = double(squeeze(in_voltage.signals.values)); % voltage
out_p = double(squeeze(out_encoder.signals.values)); % encoder position

% define parameters for signal, slow/fast sampling
F_fs = 100; % fast sampling rate, Hz
T_fs = 1/F_fs; % fast sampling time
L = 10; % factor of slow sampling
T_ss = L/F_fs; % slow sampling time
f_d = 8; % disturbance frequency

% post processing data
snr = 5; % signal to noise ratio
peak_mean = mean(findpeaks(out_p(2:end))); % find average peak value
d_fs = out_p'/peak_mean; % normalize encoder to 1
d_fs_noisy = awgn(d_fs,snr,'measured'); % create noise for data
d_ss = d_fs_noisy(1:L:end); % slow sampled data

% signal reconstruction
t_slow = in_W.time'; % simulation time, slow-sampled
t_fast = out_encoder.time'; % simulation time, fast sampled
a_g = [0.3, 0.5 0.9]; % bandwidth sizes
d_est_fir = zeros(1,length(d_fs)); % preallocate size of FIR and IIR signal reconstruction
d_est_iir = zeros(length(a_g),length(d_fs));
[wk_fir] = w_k_fir(L,f_d,T_fs); % coefficients for FIR-MMP
for k = 1:length(a_g) % find coeff for IIR-MMP for various alpha values
[wk_iir, Bpara] = w_k_iir(L,f_d,a_g(k),T_fs); % coefficients for IIR-MMP
[d_est_fir, d_est_iir(k,:)] = signal_recovery(wk_fir,wk_iir,Bpara,L,d_ss);
end

%% plotting signal reconstruction
% plot comparison of true encoder, noisy encoder, FIR, and IIR recovery
x_lim = [3.45 3.75];
y_lim = [-3 2];
figure
s = stairs(t_fast,d_fs); % plot of true encoder
s.Color = [0.4 0.4 0.4];
s.LineWidth = 1;
hold on
s = stairs(t_fast,d_fs_noisy); % plot of noisy encoder
s.LineWidth = 1;
s.Color = [0 0 0.65];
s = stairs(t_fast,d_est_fir); % plot of FIR-MMP
s.LineWidth = 1.3;
s.LineStyle = '-.';
s.Marker = 'x';
s.MarkerSize = 8;
s.Color = [0.9290 0.6940 0.1250];
s = stairs(t_fast,d_est_iir(3,:)); % plot of IIR-MMP
s.LineWidth = 1.3;
s.Marker = 'o';
s.MarkerSize = 7;
s.Color = [1 0 0];
s.LineStyle = ':';
legend('True Signal','Noisy Signal','FIR MMP','IIR MMP','Location','southeast')
ax = gca;
ax.FontSize= 12;
xlim(x_lim)
ylim(y_lim)
xlabel('Time (sec)')
ylabel('Normalized Encoder Count')
hold off

%% plotting RMS
% remove matched measurements (i.e. d_fs[nL] = d_ss[n])
d_fs_inter = d_fs;
d_fir_inter = d_est_fir; 
d_iir_inter = d_est_iir;
d_fs_inter(1:L:end) = []; 
d_fir_inter(1:L:end) = [];
d_iir_inter(:,1:L:end) = [];

% finding absolute error
err_FIR = abs(d_fs_inter - d_fir_inter);
err_IIR = zeros(length(a_g),length(d_iir_inter));
for k = 1:length(a_g)
err_IIR(k,:) = abs(d_fs_inter - d_iir_inter(k,:));
end

% plotting error
figure
s = stairs(err_FIR);
s.LineWidth = 1.5;
s.LineStyle = '-';
s.Marker = 'x';
s.Color = [0.9290 0.6940 0.1250];
hold on
s = stairs(err_IIR(3,:));
s.LineWidth = 1.3;
s.Color = [1 0 0];
s.LineStyle = ':';
s.Marker = 'o';
hold off
legend('FIR MMP','IIR MMP','Location','northeast')
ax = gca;
ax.FontSize= 12;
ylabel('Absolute Error')
xlabel('Time (sec)')
ylim([0 2])

% spectral of error
d_IIR_spec = specCal(d_est_iir(3,:),F_fs);
d_FIR_spec = specCal(d_est_fir,F_fs);

% spectral of error
figure()
h(1) = semilogy(d_FIR_spec.f,d_FIR_spec.amp);
hold on
h(2) = semilogy(d_IIR_spec.f,d_IIR_spec.amp);
hold on
legend('FIR','IIR')
xlim([d_FIR_spec.f(1) d_FIR_spec.f(end)])

% finding RMS error
rms_FIR = rms(err_FIR);
rms_IIR = rms(err_IIR,2);
err_all = [err_FIR;err_IIR]';
rms_all = [rms_FIR;rms_IIR]';

% box chart of RMS error
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

%% bode plots
k_tot = 1:(L-1); % index of intersample points

% plotting options
color_all = {[0 0 0], [0.9290 0.6940 0.1250], [0 1 0], [1 0 0]};
line_style = {'-','--','--','-.',':'};
w_in_Hz_end = 1/(2*T_fs); % Nyq freq of fast sampling
w_in_Hz = logspace(-2,log10(w_in_Hz_end),10001); % range of freq (Hz)
w_in_rad = w_in_Hz*2*pi; % convert to rad/s
font_size = 12; % font size
x_lim_PQ = [1 w_in_Hz(end)]; % PQ x_lim
y_lim = [-50 15]; % y_lim

% calc bode plot values
[W_FIR_num, W_FIR_den] = w_tf_fir(wk_fir); % coeff for TF
W_k_FIR = tf(W_FIR_num(9,:),W_FIR_den(9,:),T_fs); % FIR-MMP TF
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
% plot IIR Bode plot 
for i = 1:length(a_g)
    [wk_iir(i,:,:), B_para] = w_k_iir(L,f_d,a_g(i),T_fs); % IIR-MMP coeff
    [W_IIR_num, W_IIR_den] = w_tf_iir(squeeze(wk_iir(i,:,:)),B_para);
    W_k_IIR = tf(W_IIR_num(1,:),W_IIR_den(1,:),T_fs); % IIR TF for k=1
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
Legend_bode{5}='k = 2,...9';
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
x_line = xline(f_d(i));
x_line.Color = [0 0 0];
x_line.LineWidth = 1;
x_line.LineStyle = '--';
end
x_line.Label = '8 Hz';
x_line.LabelHorizontalAlignment = 'center';
x_line.LabelVerticalAlignment = 'bottom';
hold off
legend(Legend_bode)
xlim([0,20])
%% plotting pzmap
figure
W_k_IIR = tf(1,W_IIR_den(1,:),T_fs); % IIR-MMP TF poles only
pzmap(W_k_IIR,'k') % map of IIR-MMP poles
hold on
for i = 1:length(a_g)
    [W_IIR_num, W_IIR_den] = w_tf_iir(squeeze(wk_iir(i,:,:)),B_para);
    for k = 1:k_length
    % W_k_IIR(i,k) = tf(W_IIR_num(k,:),W_IIR_den(k,:),T_fs);
    W_k_IIR(i,k) = tf(W_IIR_num(k,:),1,T_fs); % IIR-MMP TF zeros only
    end
end
for k = 3:3:k_length
    pzmap(W_k_IIR(3,k)) % plot zeros
end
Legend_pz = cell(1,4);
Legend_pz{1} = strcat('Poles (\alpha = ',num2str(a_g(3)),')');
n_leg = 2;
for iter=3:3:k_length
  Legend_pz{n_leg}=strcat('k = ', num2str(k_tot(iter)));
  n_leg = n_leg+1;
end
legend(Legend_pz,'Location','Best')