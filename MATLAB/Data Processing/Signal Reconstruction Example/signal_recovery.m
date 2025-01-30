function [d_fir,d_iir] = signal_recovery(w_k_fir,w_k_iir, B_para, L, d_ss)
    % Assigning initial values   
    n_w = size(w_k_iir,1); % n_w = 2m disturbances (index starts at 0,1,...,2m-1)
    buffer_indices = 2:L:((n_w-1)*L+2);
    max_buffer_size = (n_w*L+1); % Limit to 32,768 or less
    phi = zeros(1, n_w); % Preallocate phi using the same type as d_ss
    n_ss = 0; % Slow sampling count
    d_buff = zeros(1, max_buffer_size); % Circular buffer
    d_dim = size(d_ss);
    fs_total = (d_dim(2)-1)*L+1;
    d_fir = zeros(1,fs_total);
    d_iir = zeros(1,fs_total);
    for n_fs = 1:fs_total
        k = mod((n_fs-1),L); % check if on intersample
        % Update slow sampling count when k == 0
        if k == 0
            n_ss = n_ss + 1; % Number of slow samples
        end
        d_temp = d_ss(n_ss);
        % Case 1: n_ss < n_w (not enough data points for FIR)
        if n_ss < n_w
            if k == 0 % Fast and slow align
                phi = [d_temp, phi(1:end-1)]; % Shift and update phi
                d_iir(n_fs) = d_temp; 
                d_fir(n_fs) = d_temp;
            else % Not enough data for intersample estimation
                d_iir(n_fs) = 0; 
                d_fir(n_fs) = 0;
            end
    % Case 2: n_ss >= n_w (enough data points for FIR)
        elseif n_ss >= n_w && n_ss < (n_w + 1)
            if k == 0
                phi = [d_temp, phi(1:end-1)];
                d_fir(n_fs) = d_temp;
                d_iir(n_fs) = d_temp;
            else % FIR reconstruction
                d_fir(n_fs) = phi * w_k_fir(:, k); 
                d_iir(n_fs) = 0;
            end
    % Case 3: n_ss >= n_w + 1 (enough data points for IIR)
        else
            if k == 0
                phi = [d_temp, phi(1:end-1)];
                d_fir(n_fs) = d_temp;
                d_iir(n_fs) = d_temp;        
            else % IIR reconstruction
                d_fir(n_fs) = phi * w_k_fir(:, k); 
                d_temp = flip(d_buff);
                d_iir(n_fs) = phi*w_k_iir(:,k) - ...
                              d_temp(buffer_indices)*flip(B_para(2:end))';
            end
        end
    % % update buffer
        d_buff = [d_iir(n_fs), d_buff(1:end-1)];
    end
end