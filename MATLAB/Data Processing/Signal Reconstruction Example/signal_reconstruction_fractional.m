clear all
close all
clc

%% fractional signal reconstruction generalized

% signal recovery parameters
L = 9/4;
[R, den] = rat(1/L);
RL = L*R; % scaling to make R*L an integer
t_ss = 0.1;
t_fs = t_ss/L; % fast rate
f_ny_ss = 1/(2*t_ss); % slow sampling Nyquist frequency
f_ny_fs = 1/(2*t_fs); % fast sampling Nyquist freqquency
f_ny_comb = 1/(2*t_ss/RL); % Nyquist frequency of the combined rate

% simulation time and sin generator
t_end = 100;
% generate signals for reconstruction'%
m_d = 7;
f_hz = sort(f_ny_ss+(f_ny_comb-f_ny_ss)*rand(1,m_d)); % generate random 
A_hz = rand(1,m_d);
p_off = rand(1,m_d);

[sig_ct, sig_ss, sig_fs] = sin_generator(f_hz,A_hz,p_off,t_fs,t_ss,t_end);
% plotting the signal
figure
plot(sig_ct(1,:),sig_ct(2,:))
hold on
stairs(sig_ss(1,:),sig_ss(2,:))
stairs(sig_fs(1,:),sig_fs(2,:))
legend('CT','SS','FS')

y_est = multi_phase_recovery(sig_ss(2,:), f_hz, t_fs, t_end, R, L);

figure
stairs(sig_ct(1,:),sig_ct(2,:))
hold on
stairs(y_est(1,:), y_est(2,:),'o')
xlim([R^2 R^2+2*RL]*t_ss)
xlabel('Time (sec)')
ylabel('Output')
legend('Ground Truth','Fractional Signal Recovery')
title('Signal Recovery')

%% MATLAB functions
% multi phase signal recovery
function [out] = multi_phase_recovery(input_signal, f_hz, t_fs, t_end, R, L)
    length_ss = length(input_signal);  % Length of the input signal
    out = zeros(2,(length_ss-1)*R*L);
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
    % Efficient version of w_k_frac
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

    % Precompute identity-type matrices for each k
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