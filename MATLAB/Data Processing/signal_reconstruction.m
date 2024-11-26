% clc
% close all
% clear all
%% model parameters
% DC motor from MinSegMotor
% sampled at 250 Hz
% input is 3.9V with a sin wave at 8 Hz.
% noise is introduced into the "clean" measured data
% Fast sampling here defined at 250 Hz, slow sampling at 50 Hz
% add this to the path C:\Users\tpjch\OneDrive\Documents\MATLAB\RASPlib
Fs = 500;
T_fs = 1/Fs;
a_g = 0.9;
L_t = 50;
T_ss = T_fs*L_t;
f_d = 8;
% [w_k_IIR B_para] = W_coeff_IIR(L_t,f_d,a_g,T_fs);
[w_k] = W_coeff_FIR(L_t,f_d,T_fs);
return
%% exporting data
close all
y_encoder = squeeze(double(out_encoder.signals.values));
t_encoder = out_encoder.time;
y_w_FIR = out_W_FIR.signals.values;
t_w = out_W_FIR.time;
u_w_FIR = squeeze(in_W_FIR.signals.values);
t_u = squeeze(in_W_FIR.time);

figure
stairs(t_encoder,y_encoder)
hold on
stairs(t_w,y_w_FIR)
stairs(t_u,u_w_FIR)
xlim([3 5])
legend('Encoder','Reconstructed - FIR','Slow-Sampled')

%% from spiral servo writing in HDD paper
close all
z = tf('z',T_fs);
num = [1 2 -1];
den = [0 4 0];
G_zpet = tf(num,den,T_fs);

figure
bode(G_zpet)
legend('G_zpet')
% 
% y_zpet = lsim()
n_enc = height(y_encoder);
y_filt(1) = y_encoder(1);
for i = 2:(n_enc-1)
    y_filt(i) = 0.25*(y_encoder(i+1)+2*y_encoder(i)-y_encoder(i-1));
end
y_filt = y_filt';
t_filt = t_encoder(1:end-1);
figure
stairs(t_encoder,y_encoder)
hold on
stairs(t_filt,y_filt)
legend('Encoder','Low-Pass')

y_encoder_norm = y_filt-mean(y_filt);

figure
stairs(t_filt,y_encoder_norm)
legend('Zerod Encoder')