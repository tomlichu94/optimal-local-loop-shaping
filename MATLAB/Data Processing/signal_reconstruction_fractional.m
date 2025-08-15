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

% input_signal = awgn(sig_ss(2,:),10,'measured');
% input_signal = sig_ss(2,:);

y_est = multi_phase_recovery_fir(input_signal, f_hz, t_fs, t_end, R, L);

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
function out = signal_recovery_iir(input_signal, f_hz, t_s, a_g, R, L)
% IIR-signal recovery using fractional-speeds, every Rth slow sample is
% used since L is a fraction and L = num/R.
% (e.g. fast sampling is 5/2 times faster than slow sampling)

% Inputs:
%   input_signal : row vec of slow sampled signal only
%   f_hz         : signal frequency to be recovered
%   t_s          : fast sampling time
%   a_g          : bandwidth of the IIR signal recovery
%   R            : if t_ss = t_s * L, where t_ss is the slow sampling time
%                  and L = num/den, den = R (num and den are integers)
%   L            : upsampling factor
%
% Output:
%   out:         : output recoverd signal upsampled by L

    RL = R * L;
    % find W_k coefficiets
    [w_kiir, Bpara] = w_kiir_frac(f_hz, t_s, a_g, R, L); % extract W_k coefficients
    
    % determine signal lengths for fast sampling
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

    % signal recovery algorithm
    for n_fs = 1:length_fs
        d_temp = input_signal((n_ss - 1) * R + 1);
        k = mod(n_fs - 1, RL); % intersample point index
        % Case 1: not enough slow samples
        if n_ss < n_w + 1
            if k == 0
                phi = [d_temp, phi(1:end - 1)];
                out(n_fs) = 0;
                n_ss = n_ss + 1;
            else
                out(n_fs) = 0;
            end
        % Case 2: enough for FIR
        elseif n_ss >= n_w && n_ss < n_w + R
            if k == 0
                phi = [d_temp, phi(1:end - 1)];
                out(n_fs) = d_temp;
                n_ss = n_ss + 1;
            else
                out(n_fs) = phi * w_kiir(:, k);
            end
        % Case 3: IIR reconstruction
        else
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
function [w_kiir, Bpara] = w_kiir_frac(f_hz, t_s, a_g, R, L)
% IIR-MMP coefficents ...
% ... (e.g. fast sampling is 5/2 times faster than slow sampling)
% Inputs:
%   f_hz         : signal frequency to be recovered
%   t_s         : fast sampling time
%   a_g          : bandwidth of the IIR signal recovery
%   R            : if t_ss = T_fs * L, ...
%                  ... and L = num/den, den = R (num and den are integers)
%   L            : upsampling factor
%
% Output:
%   w_kiir       : outputs coefficients for signal recovery ...
%                  ... 2m_d x (RL-1), where m is the number of frequencies
%   Bpara        : denominator coefficients of the IIR-MMP

    RL = R*L;
    if mod(RL,1) ~= 0
        error('RL must be an integer')
    end

    k_max = RL-1; % max num of intersamples
    m_d = numel(f_hz); % num of frequencies
    
    % finding coefficeints for A(z) and B(z), used in diophantine Eqn, ...
    % ... where F(z)A(z) + z^{-k} W_{k}(z) - B(z) = 0
    Apara = 1;
    Bpara = 1;
    for i = 1:m_d
        omega = 2 * pi * f_hz(i) * t_s;
        Apara = conv(Apara,[1 -2*cos(omega) 1]);
        Bpara = conv(Bpara,[1 -2*a_g*cos(RL*omega) a_g^2]); % output
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
    
    % determine location of 1 for each col of E_k
    for ki = 1:k_max
        for j = 1:2*m_d
            row_idx = ki + RL * (j - 1);
            E_k(row_idx, j, ki) = 1;
        end
        M_k(:,1:m_k2, ki) = M_kt;
        M_k(:, m_k2 + 1:end, ki) = E_k(:, :, ki);
    end

    % the signal reconstruction is M_k * [f_vec w_vec]' = [a_vec 0]'
    % let us express this as Ax = b, where the solution is x = inv(A)b
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
    w_kiir = x_sol(end - n_w:end, :); % output
end

% multi phase signal recovery for IIR
function [out] = multi_phase_recovery_iir(input_signal, f_hz, t_s, t_end, a_g, R, L)
% IIR signal recovery, with R recoveries in tandem, increasing output ...
% ... sampling rate by t_ss/(RL) 
% R signal recoveries performend, each delayed by t_ss

% Inputs:
%   input_signal : row vec of slow sampled signal only
%   f_hz         : signal frequency to be recovered
%   t_s          : fast sampling time
%   t_end        : end simulation time
%   a_g          : bandwidth of the IIR signal recovery
%   R            : if t_ss = t_s * L, where t_ss is the slow sampling time
%                  and L = num/den, den = R (num and den are integers)
%   L            : upsampling factor
%
% Output:
%   out:         : output recoverd signal upsampled by RL ...
%                : ... out(1,:) is the time, out(2,:) is the signal

    % preallocate output
    length_ss = length(input_signal);     % length of slow sampled output
    out = zeros(2, (length_ss-1)*R*L+1);

    % perform multiple signal recovery, offset by one slow sample
    for i = 1:R
        y_fcn = signal_recovery_iir(input_signal(i:end), f_hz, t_s, a_g, R, L);
        base_idx = (i-1)*R*L+1; % starting index
        idx = base_idx:R:base_idx+R*(length(y_fcn)-1); % indexed range
        out(2,idx) = y_fcn;
    end
    out(1,:) = 0:t_s/R:t_end;
end

% multi phase signal recovery for FIR
function [out] = multi_phase_recovery_fir(input_signal, f_hz, t_s, t_end, R, L)
% FIR signal recovery, with R recoveries in tandem, increasing output ...
% ... sampling rate by t_ss/(RL) 
% R signal recoveries performend, each delayed by t_ss

% Inputs:
%   input_signal : row vec of slow sampled signal only
%   f_hz         : signal frequency to be recovered
%   t_s          : fast sampling time
%   t_end        : end simulation time
%   R            : if t_ss = t_s * L, where t_ss is the slow sampling time
%                  and L = num/den, den = R (num and den are integers)
%   L            : upsampling factor
%
% Output:
%   out:         : output recoverd signal upsampled by RL ...
%                : ... out(1,:) is the time, out(2,:) is the signal
    % preallocate output
    length_ss = length(input_signal);  % Length of the input signal
    out = zeros(2,(length_ss-1)*R*L+1);

    % perform multiple signal recovery, offset by one slow sample
    for i = 1:R
        y_fcn = signal_recovery_fir(input_signal(i:end), f_hz, t_s, R, L);
        base_idx = (i-1)*R*L+1; % starting index
        idx = base_idx:R:base_idx+R*(length(y_fcn)-1); % indexed range
        out(2,idx) = y_fcn;
    end
    out(1,:) = 0:t_s/R:t_end;
end

% returns w_k coefficient
function [w_kfir] = w_kfir_frac(f_hz, t_s, R, L)
% FIR-MMP coefficents ...
% ... (e.g. fast sampling is 5/2 times faster than slow sampling)
% Inputs:
%   f_hz         : signal frequency to be recovered
%   t_s         : fast sampling time
%   R            : if t_ss = T_fs * L, ...
%                  ... and L = num/den, den = R (num and den are integers)
%   L            : upsampling factor
%
% Output:
%   w_kfir       : outputs coefficients for signal recovery ...
%                  ... 2m_d x (RL-1), where m is the number of frequencies

    RL = R*L;
    if mod(RL,1) ~= 0
        error('RL must be an integer');
    end

    k_max = RL - 1;
    m_d = numel(f_hz);

    % Compute A(z) coefficients
    Apara = 1;
    for i = 1:m_d
        omega = 2 * pi * f_hz(i) * t_s;
        Apara = conv(Apara, [1 -2*cos(omega) 1]);
    end
    m_a = numel(Apara);
    n_w = 2*m_d - 1;

    % Preallocate dimensions
    m_k1 = 2*m_d * RL;
    m_k2 = 2*m_d * (RL - 1);

    % Construct M_kt matrix
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

    % output coefficients
    w_kfir = x_sol(end - n_w:end, :);
end

% signal recovery
function [out] = signal_recovery_fir(input_signal, f_hz, t_s, R, L)
% FIR-signal recovery using fractional-speeds, every Rth slow sample is
% used since L is a fraction and L = num/R.
% (e.g. fast sampling is 5/2 times faster than slow sampling)

% Inputs:
%   input_signal : row vec of slow sampled signal only
%   f_hz         : signal frequency to be recovered
%   t_s          : fast sampling time
%   a_g          : bandwidth of the IIR signal recovery
%   R            : if t_ss = t_s * L, where t_ss is the slow sampling time
%                  and L = num/den, den = R (num and den are integers)
%   L            : upsampling factor
%
% Output:
%   out:         : output recoverd signal upsampled by L

    % preallocate output
    length_fs = floor((length(input_signal)-1)*L)+1;
    out = zeros(1,length_fs);

    % finding coefficients for FIR-MMP
    w_k = w_kfir_frac(f_hz, t_s, R, L);
    
    % signal recovery algorithm
    n_w = height(w_k); 
    phi = zeros(1, n_w); % store past slow samples
    n_ss = 1; % start slow sample count
    for n_fs = 1:length_fs
        k = mod((n_fs-1),R*L); % find intersample count
        % case 1: insufficient number of historical measurements
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