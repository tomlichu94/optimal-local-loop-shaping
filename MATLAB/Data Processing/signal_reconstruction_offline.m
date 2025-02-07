clc
close all
clear all
%% model parameters
% DC motor from MinSegMotor
% sampled at 100 Hz
% input is 2.5V with a sin wave at 8 Hz.
% noise is introduced into the "clean" measured data
% Fast sampling here defined at 100 Hz, slow sampling at 10 Hz
F_fs = 100;
snr = 3;
addpath('Functions','Experimental Runs')
load run_8.mat
in_V = double(squeeze(in_voltage.signals.values));
out_p = double(squeeze(out_encoder.signals.values));
peak_p = findpeaks(out_p(2:end));
peak_mean = mean(peak_p);
out_p = out_p/peak_mean;
t_sim = in_voltage.time;
% t_sim = t_sim(2:end); % starts 2 steps after
out_p_noisy = awgn(out_p,5,'measured'); %awgn(signal,snr,input_mag,'Linear' or 'db')
% out_p_noisy = out_p; % for ACC paper, no noise
% t_sim = out_encoder.time;
L = 10;
d_fs = out_p_noisy';
d_ss = d_fs(1:L:end);
f_d = 8;
T_fs = 1/F_fs;
T_ss = L/F_fs;
t_slow = in_W.time;
t_fast = t_sim;
length_dest = length(d_fs);
a_g = [0.3, 0.5 0.9];
n_t = length(a_g);
d_est_fir = zeros(1,length_dest);
d_est_iir = zeros(n_t,length_dest);
[w_k_fir] = W_coeff_FIR(L,f_d,T_fs);
for k = 1:n_t
[w_k_fir Bpara] = W_coeff_IIR(L,f_d,a_g(k),T_fs);
[d_est_fir d_est_iir(k,:)] = signal_recovery(w_k_fir,w_k_fir,Bpara,L,d_ss);
end

%% plotting
x_lim = [3.45 3.75];
y_lim = [-3 2];
% figure 1
figure
s = stairs(t_sim,out_p);
s.Color = [0.4 0.4 0.4];
s.LineWidth = 1;
hold on
s = stairs(t_sim,out_p_noisy);
s.LineWidth = 1;
s.Color = [0 0 0.65];
s = stairs(t_fast,d_est_fir);
s.LineWidth = 1.3;
s.LineStyle = '-.';
s.Marker = 'x';
s.MarkerSize = 8;
s.Color = [0.9290 0.6940 0.1250];
s = stairs(t_fast,d_est_iir(2,:));
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

%% calculating RMS
% remove the matched sampled points
out_p_inter = out_p';
out_p_inter(1:L:end) = [];
out_fir = d_est_fir;
out_fir(1:L:end) = [];
out_iir = d_est_iir;
out_iir(:,1:L:end) = [];

err_IIR = zeros(n_t,length(out_iir));
for k = 1:n_t
err_IIR(k,:) = abs(out_p_inter - out_iir(k,:));
end
err_FIR = abs(out_p_inter - out_fir);
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

d_IIR_spec = specCal(d_est_iir(2,:),F_fs);
d_FIR_spec = specCal(d_est_fir,F_fs);

figure()
h(1) = semilogy(d_FIR_spec.f,d_FIR_spec.amp);
hold on
h(2) = semilogy(d_IIR_spec.f,d_IIR_spec.amp);
hold on
legend('FIR','IIR')
xlim([d_FIR_spec.f(1) d_FIR_spec.f(end)])

err_FIR(1:L:end) = [];
err_IIR(:,1:L:end) = [];
rms_FIR = rms(err_FIR);
rms_IIR = rms(err_IIR,2);
err_all = [err_FIR;err_IIR]';
rms_all = [rms_FIR;rms_IIR]';

figure
boxchart(err_all)
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
xticklabels({'FIR','$\alpha = 0.3$','$\alpha = 0.5$','$\alpha = 0.9$'})
hold on
p = plot(rms_all,'-o');
p.LineWidth = 1.5;
for i = 1:4
str = sprintf('%.3f',rms_all(i));
cell{i} = cellstr(str);
end
text((1:4)+0.25,rms_all+0.1,cell);
% text(x, y, sprintf('Text %f more text', variable))
legend('','RMS Error')
ylabel('Error')

% bode plots
L = 10;
k_tot = 1:(L-1);

% plotting options
color_all = {[0 0 0], [0.9290 0.6940 0.1250], [0 1 0], [1 0 0]};
w_in_Hz_end = 1/(2*T_fs); % end of plotting for Nyq freq of the fast
w_in_Hz = logspace(0,log10(w_in_Hz_end),6001);
w_in_rad = w_in_Hz*2*pi;
w_start = 1; % starting frequency for x_lim
l_width = 1.5; % linewidth
font_size = 12; % font size
x_lim_PQ = [w_start w_in_Hz(end)]; % PQ x_lim
y_lim = [-50 15]; % y_lim


[W_FIR_num, W_FIR_den] = W_tf_fir(w_k_fir);
W_k_FIR = tf(W_FIR_num(1,:),W_FIR_den(1,:),T_fs);
[mag_W_FIR phi_W_FIR w_out] = bode(W_k_FIR,w_in_rad);
mag_W_FIR = 20*log10(mag_W_FIR(:));
phi_W_FIR = wrapTo180(phi_W_FIR(:));

figure
h(1) = plot(w_in_Hz,mag_W_FIR);
hold on
h(1).Color = color_all{1};
h(1).LineStyle = line_style{1};
h(1).LineWidth = 1.5;
xlabel('Hz')
ylabel('Magnitude (dB)')
xlim(x_lim_PQ)
title('W Predictor Bode')
%%
for i = 1:length(a_g)
    [w_k_IIR B_para] = W_coeff_IIR(L,f_d,a_g(i),T_fs);
    [W_IIR_num W_IIR_den] = W_TF_IIR(w_k_IIR,B_para);
    W_k_IIR = tf(W_IIR_num(1,:),W_IIR_den(1,:),T_fs);
    [mag_W_IIR phi_W_IIR w_out] = bode(W_k_IIR,w_in_rad);
    mag_W_IIR = 20*log10(mag_W_IIR(:));
    h(i+1) = plot(w_in_Hz,mag_W_IIR);
    h(i+1).Color = color_all{i+1};
    h(i+1).LineStyle = '-';
    h(i+1).LineWidth = 1.5;
end
Legend=cell(length(a_g)+1,1);
Legend{1} = strcat('FIR');
for iter=1:length(a_g)
  Legend{iter+1}=strcat('\alpha =', num2str(a_g(iter)));
end

%%
% plotting interamples
k_length = size(k_tot,2);
[w_k_IIR B_para] = W_coeff_IIR(L,f_d,a_g(1),T_fs);
[W_IIR_num W_IIR_den] = W_TF_IIR(w_k_IIR,B_para);
for k = 1:k_length
        W_k_IIR = tf(W_IIR_num(k_tot(k),:),W_IIR_den(k_tot(k),:),T_fs);
        [mag_W_IIR phi_W_IIR w_out] = bode(W_k_IIR,w_in_rad);
        mag_W_IIR = 20*log10(mag_W_IIR(:));
        phi_W_IIR = wrapTo180(phi_W_IIR(:));
        h(k) = plot(w_in_Hz,mag_W_IIR);
        h(k).Color = [0.3 0.3 0.3];
        h(k).LineStyle = ':';
        h(k).LineWidth = 1;
        % h.Color = color_all{i+1};
        % h.LineStyle = line_style{i+1};
        % h.LineWidth = 2-k/10;
        xlabel('Hz')
        ylabel('Magnitude (dB)')
        xlim(x_lim_PQ)
end
Nyq_Hz = 1/(2*T_ss);
x_line = xline(Nyq_Hz);
x_line.Color = [0 0 0];
x_line.LineWidth = 1.5;
x_line.FontSize = 10;
x_line.FontWeight = 'bold';
for i = 1:m_d
x_line = xline(f_d(i));
x_line.Color = [0 0 0];
x_line.LineWidth = 1;
x_line.LineStyle = '--';
end
hold off

figure
for i = 1:length(a_g)
    [w_k_IIR B_para] = W_coeff_IIR(L,f_d,a_g(i),T_fs);
    [W_IIR_num W_IIR_den] = W_TF_IIR(w_k_IIR,B_para);
    for k = 1:k_length
    % for k = 1:1
        W_k_IIR = tf(W_IIR_num(k_tot(k),:),W_IIR_den(k_tot(k),:),T_fs);
        pzmap(W_k_IIR)
        hold on
    end
end
n_leg = 1;
for iter=1:k_length
  Legend{n_leg}=strcat('k = ', num2str(k_tot(iter)));
  n_leg = n_leg+1;
end
legend(Legend,'Location','Best')

W_k_IIR
abs(zero(W_k_IIR))

a = cell2mat(W_k_IIR.den);
b = cell2mat(W_k_IIR.num);

[A B C D] = tf2ss(b,a);