function [y_hat_fs, theta_hist, y_components] = ...
    rls_sine_reconstruction(y_ss, t_ss, t_fs, frequencies, lambda, delta)
%RLS_MULTISINE_RECONSTRUCTION
% Reconstruct multiple known-frequency sinusoids from slow samples using RLS.
%
% Model:
%   y(t) = sum_i [
%       a_i*sin(2*pi*f_i*t) + b_i*cos(2*pi*f_i*t)
%   ]
%
% Inputs:
%   y_ss        - Slow-sampled measurements
%   t_ss        - Slow-rate time vector
%   t_fs        - Fast-rate reconstruction time vector
%   frequencies - Known frequency vector [f1, f2, ..., fM] in Hz
%   lambda      - Forgetting factor, 0 < lambda <= 1
%   delta       - Initial covariance scaling
%
% Outputs:
%   y_hat_fs     - Reconstructed combined fast-rate signal
%   theta_hist   - RLS coefficient history, 2M x length(y_ss)
%   y_components - Reconstructed component signals, M x length(t_fs)
%
% Parameter ordering:
%   theta = [a1; b1; a2; b2; ...; aM; bM]

    % Preserve requested output orientation
    output_is_row = isrow(t_fs);

    % Convert inputs to columns
    y_ss = y_ss(:);
    t_ss = t_ss(:);
    t_fs_column = t_fs(:);
    frequencies = frequencies(:);

    Ns = length(y_ss);
    Nf = length(t_fs_column);
    M = length(frequencies);
    number_parameters = 2*M;

    if length(t_ss) ~= Ns
        error('y_ss and t_ss must have the same length.');
    end

    if isempty(frequencies)
        error('frequencies cannot be empty.');
    end

    if lambda <= 0 || lambda > 1
        error('lambda must satisfy 0 < lambda <= 1.');
    end

    if delta <= 0
        error('delta must be positive.');
    end

    % Initial parameter estimates and covariance
    theta = zeros(number_parameters,1);
    P = delta*eye(number_parameters);

    theta_hist = zeros(number_parameters,Ns);

    %% Recursive least-squares estimation
    for k = 1:Ns

        phi = zeros(number_parameters,1);

        for i = 1:M
            omega_i = 2*pi*frequencies(i);

            phi(2*i - 1) = sin(omega_i*t_ss(k));
            phi(2*i)     = cos(omega_i*t_ss(k));
        end

        % RLS gain
        denominator = lambda + phi.'*P*phi;
        K = P*phi/denominator;

        % Prediction error
        prediction = phi.'*theta;
        error_k = y_ss(k) - prediction;

        % Parameter and covariance updates
        theta = theta + K*error_k;
        P = (P - K*phi.'*P)/lambda;

        % Improve numerical symmetry
        P = 0.5*(P + P.');

        theta_hist(:,k) = theta;
    end

    %% Expand coefficient estimates to the fast time vector
    theta_fs = zeros(number_parameters,Nf);

    for j = 1:number_parameters
        theta_fs(j,:) = interp1( ...
            t_ss, ...
            theta_hist(j,:), ...
            t_fs_column, ...
            'previous', ...
            'extrap');
    end

    %% Reconstruct individual frequency components
    y_components = zeros(M,Nf);

    for i = 1:M
        omega_i = 2*pi*frequencies(i);

        a_i = theta_fs(2*i - 1,:);
        b_i = theta_fs(2*i,:);

        y_components(i,:) = ...
            a_i.*sin(omega_i*t_fs_column.') + ...
            b_i.*cos(omega_i*t_fs_column.');
    end

    % Combined reconstruction
    y_hat_fs = sum(y_components,1);

    if ~output_is_row
        y_hat_fs = y_hat_fs.';
        y_components = y_components.';
    end
end