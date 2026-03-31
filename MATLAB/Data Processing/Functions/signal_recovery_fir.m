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