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