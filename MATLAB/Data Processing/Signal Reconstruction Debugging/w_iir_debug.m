function d_est = w_iir_debug(w_k, B_para, d_ss, L)
    % debugging function
    % Keeps data of variables
    persistent phi n_fs n_ss d_IIR d_IIR_idx;
    n_w = size(w_k, 1); % n_w = 2m disturbances (index starts at 0,1,...,2m-1)

    % Assigning initial values
    if isempty(phi)
        % phi = zeros(1, n_w, 'like', d_ss); % Preallocate phi using the same type as d_ss
        phi = zeros(1,n_w);
    end
    if isempty(n_fs)
        n_fs = 0; % Fast sampling count
    end
    if isempty(n_ss)
        n_ss = 0; % Slow sampling count
    end
    if isempty(d_IIR)
        % Use a smaller buffer size, such as storing only required history
        % max_buffer_size = min(2 * n_w * L, 32768); % Limit to 32,768 or less
        max_buffer_size = 100;
        % d_IIR = zeros(1, max_buffer_size, 'like', d_ss); % Circular buffer
        d_IIR = zeros(1,100);
        d_IIR_idx = 1; % Index for circular buffer
    end

    % Increment fast sampling count
    n_fs = n_fs + 1;

    % Update slow sampling count when required
    if mod((n_fs - 1), L) == 0
        n_ss = n_ss + 1; % Number of slow samples
    end

    % Calculate the kth intersample point (k = 1,...,L-1)
    k = n_fs - (n_ss - 1) * L - 1;

    % Case 1: n_ss < n_w (not enough data points for FIR)
    if n_ss < n_w
        if k == 0
            phi = [d_ss, phi(1:end-1)]; % Shift and update phi
            d_est = d_ss; % Fast and slow align
        else
            d_est = 0; % Not enough data for intersample estimation
        end
    % Case 2: n_ss >= n_w (enough data points for FIR)
    elseif n_ss < (n_w + 1)
        if k == 0
            phi = [d_ss, phi(1:end-1)];
            d_est = d_ss;
        else
            d_est = phi * w_k(:, k); % FIR reconstruction
        end
    % Case 3: n_ss >= n_w + 1 (enough data points for IIR)
    else
        if k == 0
            phi = [d_ss, phi(1:end-1)];
            d_est = d_ss;
        else
            % Circular buffer handling
            % buffer_start = mod(d_IIR_idx - n_w * L, length(d_IIR)) + 1;
            % buffer_indices = buffer_start:L:d_IIR_idx - 1;
            % d_est = phi * w_k(:, k) - ...
            %     flip(d_IIR(buffer_indices)) * B_para(2:end);
                 % buffer_start = mod(d_IIR_idx - n_w * L, length(d_IIR)) + 1;
            buffer_indices = (L+1+k):L:(k+1);
            d_est = phi * w_k(:, k) - ...
                flip(d_IIR(buffer_indices)) * B_para(2:end)';
        end
    end
    % Update the circular buffer
    d_IIR(d_IIR_idx) = d_est;
    d_IIR_idx = mod(d_IIR_idx, length(d_IIR)) + 1;
end
