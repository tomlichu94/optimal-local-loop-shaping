clear all
close all
clc
%%%%%%%%%% signal reconstruction for multiple disturbances in simulink
m_d = 1; % define number of disturbances
L = 5; % scale of the slow sampling, T_ss = L*T_fs
batches = 100; % number of slow sampled measurements

T_fs = 1/100; % sampled time
T_ss = T_fs*L; % slow sampling time
f_nyq = 1/(2*T_ss); % Nyquist frequency
f_d = f_nyq*(1+rand(1,m_d)); % disturbances beyond Nyquist
A_d = 0.5+2*rand(1,m_d); % amplitudes
phi_off = rand(1,m_d)*pi; % random phase offsets up to pi

%% disturbance generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end_time = batches*L*T_fs;
d_fs_time = 0:T_fs:end_time;
d_ss_time = 0:T_ss:end_time;
dc_time = 0:T_fs/100:end_time;
[dc,d_fs,d_ss] = disturb_gen(f_d, A_d, dc_time, d_fs_time, d_ss_time,phi_off);
% d_ss = awgn(d_ss,5); % add noise

%% signal reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_g = 0.9;
[w_k_IIR Bpara] = W_coeff_IIR(L,f_d,a_g,T_fs);
[w_k] = W_coeff_FIR(L,f_d,T_fs);
[d_est_FIR, d_est_IIR] = signal_recovery(w_k,w_k_IIR,Bpara,L,d_ss);
%% plotting
f = figure();
% f.Position = [60 60 800 400]; 
plot(dc_time,dc,'Linewidth',1)
xlim([1 1.3])

f = figure();
% f.Position = [60 60 800 400]; 
plot(dc_time,dc,'Linewidth',1)
hold on
s = stairs(d_ss_time,d_ss);
s.LineWidth = 1;
s.Marker = 'o';
s.MarkerSize = 8;
xlim([1 1.3])
hold off
% legend('Continuous Time Signal','Slow Sampled')

f = figure();
% f.Position = [60 60 800 400]; 
plot(dc_time,dc,'Linewidth',1)
hold on
s = stairs(d_fs_time,d_est_FIR);
s.LineWidth = 1;
s.LineStyle = '-';
s.Marker = 'x';
s.MarkerSize = 8;
s.Color = [1 0 0];
% ylim([-2 2])
xlim([1 1.3])
hold off