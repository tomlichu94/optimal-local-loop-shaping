clear all
close all
clc
addpath('Functions')

%% signal properties
f_hz = 8; % signal freq
f_amp = 1; % signal amp
f_off = 0.4; % signal phase offset
t_ss = 1/10; % slow sampling rate, Nyq freq = 5 hz
L = 5/2; % fast sampling rate scale
t_fs = t_ss/L; % fast sampling rate
[R, den] = rat(1/L);
RL = L*R; % scaling to make R*L an integer

%% simulink function
% simulink environment
t_end = 10; 
ts_sim = t_fs/R; % match the upsampling rate of multi-phase signal recovery
t_sim = 0:ts_sim:t_end;

% generate signal to be reconstructed
[sig_ct, sig_ss, sig_fs] = sin_generator(f_hz,f_amp,f_off,t_fs,t_ss,t_end);

% plotting slow and fast sampled signal
x_lim = [5.1, 5.5];
close all
figure
plot(sig_ct(1,:),sig_ct(2,:))
hold on
stairs(sig_ss(1,:),sig_ss(2,:),'-o')
stairs(sig_fs(1,:),sig_fs(2,:),':x')
legend('CT','SS','FS')
xlim(x_lim)

%% signal recovery properties
w_kfir = w_kfir_frac(f_hz, t_fs, R, L);

%%%%%%%%%%%%%%%%%%%%%% run simulink test function %%%%%%%%%%%%%%%%%%%%%%%%%
% fixed variables  
    i_sim_ds = 1;
    n_w = size(w_kfir,1); % n_w = 2m disturbances (index starts at 0,1,...,2m-1)
    phi = zeros(1, n_w);
    n_fs = 0; % Fast sampling count
    n_ss = 0; % Slow sampling count
d_temp = zeros(RL-1,1);
for i_sim = 1:length(t_sim)

% downsample slow signal slow signal
    if mod(i_sim,RL) == 0
        i_sim_ds = i_sim_ds + 1;
    end
    % ss_frac_idx = i_sim_ds-1;
    % d_ss = sig_ss(2,ss_frac_idx);
    d_ss = sig_ss(2,i_sim_ds);

    % Increment fast sampling count
    n_fs = n_fs + 1;
    k = mod((n_fs-1),RL); % check if on intersample
    % Update slow sampling count when k == 0
    if k == 0
        n_ss = n_ss + 1; % Number of slow samples
    end
    % d_ss = sig_ss(2,n_ss*2-1);
    % Case 1: n_ss < n_w (not enough data points for FIR)
    if n_ss < n_w
        if k == 0 % Fast and slow align
            phi = [d_ss, phi(1:end-1)]; % Shift and update phi
            d_fir = d_ss;
        else % Not enough data for intersample estimation
            d_fir = 0;
        end
    % Case 2: n_ss >= n_w (enough data points for FIR)
    elseif n_ss
        if k == 0
            phi = [d_ss, phi(1:end-1)];
            d_fir = d_ss;
        else % FIR reconstruction
            d_fir = phi * w_kfir(:, k); 
        end
    end
    d_fs(i_sim) = d_fir;
end   
figure
stairs(sig_fs(1,:),sig_fs(2,:));
hold on
stairs(t_sim,d_fs,'o')
%%%%%%%%%%%%%%%%%%%%%% end simulink test function %%%%%%%%%%%%%%%%%%%%%%%%%

%% functions
function [sig_ct, sig_ss, sig_fs] = sin_generator(f_hz,A_dist,p_off,t_fs,t_ss,t_end)
    t_ct = t_fs/100;
    t_lin_ct = 0:t_ct:t_end;
    t_lin_ss = 0:t_ss:t_end;
    t_lin_fs = 0:t_fs:t_end;
    sig_ct = [t_lin_ct; zeros(1,length(t_lin_ct))];
    sig_ss = [t_lin_ss; zeros(1,length(t_lin_ss))];
    sig_fs = [t_lin_fs; zeros(1,length(t_lin_fs))];
    
    for i = 1:length(f_hz)
        sig_ct(2,:) = sig_ct(2,:) + A_dist(i)*sin(2*pi*f_hz(i)*t_lin_ct+p_off(i));
        sig_ss(2,:) = sig_ss(2,:) + A_dist(i)*sin(2*pi*f_hz(i)*t_lin_ss+p_off(i));
        sig_fs(2,:) = sig_fs(2,:) + A_dist(i)*sin(2*pi*f_hz(i)*t_lin_fs+p_off(i));
    end  
end