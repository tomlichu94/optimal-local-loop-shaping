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