clear all
close all
clc

%% fractional signal reconstruction generalized

% signal recovery parameters
L = 5/3;
[R, den] = rat(1/L);
RL = L*R; % scaling to make R*L an integer
t_ss = 0.1;
t_fs = t_ss/L; % fast rate
f_ny_ss = 1/(2*t_ss); % slow sampling Nyquist frequency
f_ny_fs = 1/(2*t_fs); % fast sampling Nyquist freqquency
f_ny_comb = 1/(2*t_ss/RL); % Nyquist frequency of the combined rate

% simulation time and sin generator
t_end = 10;
% generate signals for reconstruction'%
% m_d = 7;
% f_hz = sort(f_ny_ss+(f_ny_comb-f_ny_ss)*rand(1,m_d)); % generate random 
% A_hz = rand(1,m_d);
% p_off = rand(1,m_d);
f_hz = 8;
A_hz = 1;
p_off = 0;

[sig_ct, sig_ss, sig_fs] = sin_generator(f_hz,A_hz,p_off,t_fs,t_ss,t_end);
% plotting the signal
figure
plot(sig_ct(1,:),sig_ct(2,:))
hold on
stairs(sig_ss(1,:),sig_ss(2,:))
stairs(sig_fs(1,:),sig_fs(2,:))
legend('CT','SS','FS')

input_signal = awgn(sig_ss(2,:),10,'measured');

y_est = multi_phase_recovery(input_signal, f_hz, t_fs, t_end, R, L);

figure
stairs(sig_ct(1,:),sig_ct(2,:))
hold on
stairs(y_est(1,:), y_est(2,:),'o')
% xlim([R^2 R^2+2*RL]*t_ss)
xlabel('Time (sec)')
ylabel('Output')
legend('Ground Truth','Fractional Signal Recovery')
title('Signal Recovery')

%% test for IIR
a_g = 0.9;
y_est_iir = multi_phase_recovery_iir(input_signal, f_hz, t_fs, t_end, a_g, R, L);
figure
stairs(sig_ct(1,:),sig_ct(2,:))
hold on
stairs(y_est_iir(1,:), y_est_iir(2,:),'o')
% xlim([R^2 R^2+2*RL]*t_ss)
xlabel('Time (sec)')
ylabel('Output')
legend('Ground Truth','IIR Fractional Signal Recovery')
title('Signal Recovery')


%% MATLAB functions
% signal recovery using iir
function out = signal_recovery_iir(input_signal, f_hz, t_fs, a_g, R, L)
% SIGNAL_RECOVERY_IIR - Performs IIR-based signal recovery on slow-sampled input
%
% Syntax:
%   out = signal_recovery_iir(input_signal, f_hz, t_fs, a_g, R, L)
%
% Inputs:
%   input_signal : Row vector of slow-sampled signal values
%   f_hz         : Vector of harmonic frequencies in Hz
%   t_fs         : Fast sampling time
%   a_g          : Group delay parameter for IIR filter design
%   R            : Downsampling factor
%   L            : Upsampling factor
%
% Output:
%   out          : Reconstructed signal at fast sampling rate

    RL = R * L;
    [w_kiir, Bpara] = w_kiir_frac(f_hz, t_fs, a_g, R, L); % extract W_k coefficients
    length_ss = length(input_signal);
    set_ss = 1:R:length_ss;
    length_set = length(set_ss) - 1;
    length_fs = length_set*RL + 1;
    out = zeros(1, length_fs); % preallocate recovered signal array
    n_w = height(w_kiir);
    n_ss = 1; % slow sample index counter
    % Create circular buffer for IIR reconstruction
    buffer_indices = 2:RL:((n_w - 1) * RL + 2);
    max_buffer_size = n_w * RL + 1;
    d_buff = zeros(1, max_buffer_size); % fast signal buffer
    phi = zeros(1, n_w); % store past slow samples

    for n_fs = 1:length_fs
        d_temp = input_signal((n_ss - 1) * R + 1);
        k = mod(n_fs - 1, RL); % intersample point index

        if n_ss < n_w + 1
            % Case 1: not enough slow samples
            if k == 0
                phi = [d_temp, phi(1:end - 1)];
                out(n_fs) = 0;
                n_ss = n_ss + 1;
            else
                out(n_fs) = 0;
            end
        elseif n_ss >= n_w && n_ss < n_w + R
            % Case 2: enough for FIR
            if k == 0
                phi = [d_temp, phi(1:end - 1)];
                out(n_fs) = d_temp;
                n_ss = n_ss + 1;
            else
                out(n_fs) = phi * w_kiir(:, k);
            end
        else
            % Case 3: IIR reconstruction
            if k == 0
                phi = [d_temp, phi(1:end - 1)];
                out(n_fs) = d_temp;
                n_ss = n_ss + 1;
            else
                d_temp_buff = flip(d_buff);
                out(n_fs) = phi * w_kiir(:, k) - ...
                            d_temp_buff(buffer_indices) * flip(Bpara(2:end))';
            end
        end

        % Update circular buffer
        d_buff = [out(n_fs), d_buff(1:end - 1)];
    end
end

% coefficients for iir
function [w_kiir, Bpara] = w_kiir_frac(f_d, T_s, a_g, R, L)
    % T_s is the fast sampling rate
    RL = R*L;
    if mod(RL,1) ~= 0
        error('RL must be an integer')
    end

    k_max = RL-1;
    m_d = numel(f_d); % num of frequencies
    
    Apara = 1; % defined as the coefficients for A(z) = 1 + a_1 z^-1 + ... + z^-2m
    Bpara = 1;
    for i = 1:m_d
        omega = 2 * pi * f_d(i) * T_s;
        Apara = conv(Apara,[1 -2*cos(omega) 1]);
        Bpara = conv(Bpara,[1 -2*a_g*cos(RL*omega) a_g^2]);
    end
    m_a = numel(Apara); % num of coeff for Apara
    n_w = 2* m_d - 1; % number of W coefficients starting from 0

    % preallocate dimensions
    m_k1 = 2* m_d * RL; % rows of M_kt
    m_k2 = 2* m_d * (RL-1); % columns of M_kt

    % construct M_kt matrix (Toeplitz-like)
    M_kt = zeros(m_k1, m_k2);
    for i = 1:m_k2
        M_kt(i:i + m_a - 1, i) = Apara(:);
    end

    % preallocate E_k and M_k using 3D arrays, (3rd dim is for kth sample)
    M_k = zeros(m_k1, m_k1, k_max);
    E_k = zeros(m_k1, 2*m_d, k_max);
    
    % determine 1 location for each col of E_k
    for ki = 1:k_max
        for j = 1:2*m_d
            row_idx = ki + RL * (j - 1);
            E_k(row_idx, j, ki) = 1;
        end
        M_k(:,1:m_k2, ki) = M_kt;
        M_k(:, m_k2 + 1:end, ki) = E_k(:, :, ki);
    end

    % the signal reconstruction is M_k * [f_vec w_vec]' = [a_vec 0]'
    % let us express this as Ax = b
    a_sol = -[Apara(2:end)'; zeros(m_k1 - m_a + 1,1)];
    a_size = size(a_sol);
    b_sol = zeros(a_size);
    for i = 1:2*m_d
        b_sol(i*RL) = Bpara(i+1);
    end

    % finding w_k coefficients
    x_sol = zeros(m_k1, k_max);
    for i = 1:k_max
        x_sol(:,i) = M_k(:,:,i) \ (a_sol+b_sol);
    end
    w_kiir = x_sol(end - n_w:end, :);
end

% multi phase signal recovery for IIR
function [out] = multi_phase_recovery_iir(input_signal, f_hz, t_fs, t_end, a_g, R, L)
    length_ss = length(input_signal);         % e.g., 51
    out = zeros(2, (length_ss-1)*R*L+1);                % preallocate output
    for i = 1:R
        % shifted_input = [input_signal(i:end), zeros(1, i-1)];
        y_fcn = signal_recovery_iir(input_signal(i:end), f_hz, t_fs, a_g, R, L);
        base_idx = (i-1)*R*L+1; % starting index
        idx = base_idx:R:base_idx+R*(length(y_fcn)-1); % indexed range
        out(2,idx) = y_fcn;
    end
    out(1,:) = 0:t_fs/R:t_end;
end

% multi phase signal recovery for FIR
function [out] = multi_phase_recovery(input_signal, f_hz, t_fs, t_end, R, L)
    length_ss = length(input_signal);  % Length of the input signal
    out = zeros(2,(length_ss-1)*R*L+1);
    for i = 1:R
        y_fcn = signal_recovery(input_signal(i:end), f_hz, t_fs, R, L);
        base_idx = (i-1)*R*L+1; % starting index
        idx = base_idx:R:base_idx+R*(length(y_fcn)-1); % indexed range
        out(2,idx) = y_fcn;
    end
    out(1,:) = 0:t_fs/R:t_end;
end

% returns w_k coefficient
function [w_k] = w_k_frac(f_d, T_s, R, L)
    % Inputs:
    %   f_d : vector of digital frequencies
    %   T_s : fast sampling period
    %   RL  : rate conversion factor (assumed integer)
    RL = R*L;
    if mod(RL,1) ~= 0
        error('RL must be an integer');
    end

    k_max = RL - 1;
    m_d = numel(f_d);

    % Compute A(z) coefficients
    Apara = 1;
    for i = 1:m_d
        omega = 2 * pi * f_d(i) * T_s;
        Apara = conv(Apara, [1 -2*cos(omega) 1]);
    end
    m_a = numel(Apara);
    n_w = 2*m_d - 1;

    % Preallocate dimensions
    m_k1 = 2*m_d * RL;
    m_k2 = 2*m_d * (RL - 1);

    % Construct M_kt matrix (Toeplitz-like)
    M_kt = zeros(m_k1, m_k2);
    for i = 1:m_k2
        M_kt(i:i + m_a - 1, i) = Apara(:);
    end

    % Preallocate E_k and M_k using 3D arrays
    E_k = zeros(m_k1, 2*m_d, k_max);
    M_k = zeros(m_k1, m_k1, k_max);

    % Determine 1 location for each col of E_k
    for ki = 1:k_max
        for j = 1:2*m_d
            row_idx = ki + RL*(j - 1);
            E_k(row_idx, j, ki) = 1;
        end
        M_k(:, 1:m_k2, ki) = M_kt;
        M_k(:, m_k2 + 1:end, ki) = E_k(:, :, ki);
    end

    % Construct b vector
    a_sol = [-Apara(2:end)'; zeros(m_k1 - m_a + 1, 1)];

    % Solve least squares for each k
    x_sol = zeros(m_k1, k_max);
    for ki = 1:k_max
        x_sol(:, ki) = M_k(:, :, ki) \ a_sol;  % faster than pinv
    end

    % Extract w_k (last n_w + 1 rows of x_sol)
    w_k = x_sol(end - n_w:end, :);
end

% signal recovery
function [out] = signal_recovery(input_signal, f_hz, t_fs, R, L)
    persistent phi
    w_k = w_k_frac(f_hz, t_fs, R, L);
    length_fs = floor((length(input_signal)-1)*L)+1;
    out = zeros(1,length_fs);
    n_w = height(w_k); 
    n_ss = 1;
    if isempty(phi)
        phi = zeros(1, n_w); % Preallocate phi using the same type as d_ss
    end
    for n_fs = 1:length_fs
        k = mod((n_fs-1),R*L); % check if on intersample
        if n_ss < n_w+1
            if k == 0 % Fast and slow align
                phi = [input_signal((n_ss-1)*R+1), phi(1:end-1)]; % Shift and update phi
                out(n_fs) = 0;
                n_ss = n_ss+1;
            else % Not enough data for intersample estimation
                out(n_fs) = 0;
            end
        % Case 2: n_ss >= n_w (enough data points for FIR)
        else
            if k == 0
                phi = [input_signal((n_ss-1)*R+1), phi(1:end-1)];
                out(n_fs) = input_signal((n_ss-1)*R+1);
                n_ss = n_ss+1;
            else % FIR reconstruction
                out(n_fs) = phi * w_k(:, k); 
            end
        end
    end
end
    
% sin generator
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

function [w_k, Bpara] = w_kiir_check(f_d,a_g,T_s,L)
% T_s is the fast sampling rate
k = 1:(L-1);
m_d = max(size(f_d));
Apara = 1; % defined as the coefficients for A(z) = 1 + a_1 z^-1 + ... + z^-2m
Bpara = 1;
for i = 1:m_d
    Apara = conv(Apara,[1 -2*cos(f_d(i)*2*pi*T_s) 1]);
    Bpara = conv(Bpara,[1 -2*a_g*cos(L*f_d(i)*2*pi*T_s) a_g^2]);
end
% conv returns as [1 z^-1 z^-2 ... z^-2m_d], returning 2m+1 terms
%%%%%%%%% defining M_k and M_kt (t for tilde) %%%%%%%%%
n_w = 2*m_d-1;
m_kr = 2*m_d*L; % rows of M_kt
m_kc = 2*m_d*(L-1); % columns of M_kt
m_a = 2*m_d+1; % number of coefficients for A(z)
M_kt = zeros(m_kr,m_kc);
for i = 1:m_kc
    M_kt(i:(i+m_a-1),i) = Apara';
end

% create row vectors for e_k with the kth row = 1
E_k = zeros(m_kr,2*m_d);
for i = 1:max(k)
    for j = 1:2*m_d
    E_k(i+L*(j-1),j,i) = 1;
    end
end

% build M_k matrix, where there are k matrices nested within M_k
M_k = zeros(m_kr,m_kr);
for i = 1:max(k)
M_k(:,1:m_kc,i) = M_kt;
M_k(:,(m_kc+1):end,i) = E_k(:,:,i);
end

% the signal reconstruction is M_k * [f_vec w_vec]' = [a_vec 0]'
% let us express this as Ax = b
a_sol = -[Apara(2:end)'; zeros(m_kr-m_a+1,1)];
a_size = size(a_sol);
b_sol = zeros(a_size);
for i = 1:2*m_d
    b_sol(i*L) = Bpara(i+1);
end

for i = 1:max(k)
    x_sol(:,i) = pinv(M_k(:,:,i))*(a_sol+b_sol);
    f_max = 2*m_d*(L-1); % highest coefficient for F
    w_max = m_kr - f_max-1; % highest coefficient for W based on f_max
    f_k = x_sol(1:(m_kr-(n_w+1)),:);
    w_k = x_sol((m_kr-n_w):end,:);
end
end