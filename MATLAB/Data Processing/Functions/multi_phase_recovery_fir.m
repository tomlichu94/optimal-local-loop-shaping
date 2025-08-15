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