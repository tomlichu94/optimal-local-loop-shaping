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