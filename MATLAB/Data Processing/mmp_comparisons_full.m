%% RLS versus MMP multisine reconstruction benchmark
clear
close all
clc

addpath('Functions')

%% Configuration
cfg.f_fs = 100;                 % Fast sampling frequency [Hz]
cfg.L = 5;                      % Slow-to-fast sampling ratio
cfg.t_end = 50;                 % Simulation duration [s]

% True signal model
cfg.f_true = [23, 27, 34];      % Physical frequencies [Hz]
cfg.A0 = [1.0, 0.8, 0.6];       % Mean amplitudes
cfg.mod_depth = [0.4, 0.3, 0.2];
cfg.f_amp = [0.10, 0.20, 0.05]; % Amplitude modulation frequencies [Hz]
cfg.phase = [0, pi/4, -pi/3];   % Signal phases [rad]

% Noise
cfg.SNR_dB = 20;

% Reconstruction parameters
cfg.lambda = 0.98;              % RLS forgetting factor
cfg.delta = 1e4;                % RLS initial covariance
cfg.alpha = 0.90;               % IIR-MMP tuning parameter

% Analysis
cfg.t_eval_start = 20;           % Exclude startup for steady-state metrics
cfg.N_trials = 30;              % Monte Carlo trials

f_ss = cfg.f_fs/cfg.L;

fprintf('Fast sampling frequency: %.2f Hz\n', cfg.f_fs)
fprintf('Slow sampling frequency: %.2f Hz\n', f_ss)
fprintf('Slow Nyquist frequency:  %.2f Hz\n\n', f_ss/2)

check_frequency_aliases(cfg.f_true, f_ss);


%% ================================================================
%  1. Baseline comparison
%  ================================================================

data = generate_test_data(cfg, 1);

rls_result = evaluate_rls(data, cfg, cfg.f_true);
mmp_result = evaluate_mmp(data, cfg, cfg.f_true);

N = min([
    length(data.t_fast), ...
    length(data.y_fast), ...
    length(rls_result.y_hat), ...
    length(mmp_result.y_hat)
]);

t_plot = data.t_fast(1:N);
y_true = data.y_fast(1:N);
y_rls = rls_result.y_hat(1:N);
y_mmp = mmp_result.y_hat(1:N);

err_rls = y_true - y_rls;
err_mmp = y_true - y_mmp;

fprintf('Baseline results\n')
fprintf('----------------\n')
fprintf('RLS full RMSE:         %.6f\n', ...
    rls_result.metrics.rmse_full)
fprintf('MMP full RMSE:         %.6f\n', ...
    mmp_result.metrics.rmse_full)
fprintf('RLS steady RMSE:       %.6f\n', ...
    rls_result.metrics.rmse_steady)
fprintf('MMP steady RMSE:       %.6f\n', ...
    mmp_result.metrics.rmse_steady)
fprintf('RLS steady NRMSE:      %.2f %%\n', ...
    100*rls_result.metrics.nrmse_steady)
fprintf('MMP steady NRMSE:      %.2f %%\n\n', ...
    100*mmp_result.metrics.nrmse_steady)

figure
tiledlayout(2,1)

nexttile
plot(t_plot, y_true, 'LineWidth', 1.2)
hold on
plot(t_plot, y_rls, '--', 'LineWidth', 1.2)
plot(t_plot, y_mmp, ':', 'LineWidth', 1.4)
grid on
xline(cfg.t_eval_start, '--', 'Steady-state evaluation');
xlabel('Time (s)')
ylabel('Signal')
legend('Ground truth', 'RLS', 'IIR-MMP', ...
    'Location', 'best')
title('Baseline Reconstruction')

nexttile
plot(t_plot, err_rls, 'LineWidth', 1.1)
hold on
plot(t_plot, err_mmp, 'LineWidth', 1.1)
grid on
xline(cfg.t_eval_start, '--', 'Steady-state evaluation');
xlabel('Time (s)')
ylabel('Reconstruction error')
legend('RLS error', 'IIR-MMP error', ...
    'Location', 'best')
title('Baseline Residuals')


%% ================================================================
%  2. RLS delta sweep
%
%  Delta should primarily affect the startup transient.
%  ================================================================

delta_values = [1, 1e2, 1e4, 1e6];

rmse_delta_full = nan(length(delta_values), cfg.N_trials);
rmse_delta_steady = nan(length(delta_values), cfg.N_trials);

for trial = 1:cfg.N_trials

    data = generate_test_data(cfg, trial);

    for i_delta = 1:length(delta_values)

        cfg_temp = cfg;
        cfg_temp.delta = delta_values(i_delta);

        result = evaluate_rls(data, cfg_temp, cfg.f_true);

        rmse_delta_full(i_delta,trial) = ...
            result.metrics.rmse_full;

        rmse_delta_steady(i_delta,trial) = ...
            result.metrics.rmse_steady;
    end
end

mean_delta_full = mean(rmse_delta_full, 2, 'omitnan');
std_delta_full = std(rmse_delta_full, 0, 2, 'omitnan');

mean_delta_steady = mean(rmse_delta_steady, 2, 'omitnan');
std_delta_steady = std(rmse_delta_steady, 0, 2, 'omitnan');

figure
errorbar(delta_values, mean_delta_full, std_delta_full, ...
    'o-', 'LineWidth', 1.2)
hold on
errorbar(delta_values, mean_delta_steady, std_delta_steady, ...
    's-', 'LineWidth', 1.2)
set(gca, 'XScale', 'log')
grid on
xlabel('\delta')
ylabel('RMSE')
legend('Full-record RMSE', 'Steady-state RMSE', ...
    'Location', 'best')
title('Effect of RLS Initial Covariance')

[~, best_delta_index] = min(mean_delta_full);
cfg.delta = delta_values(best_delta_index);

fprintf('Selected delta: %.3g\n\n', cfg.delta)


%% ================================================================
%  3. Tune RLS lambda and MMP alpha
%  ================================================================

lambda_values = [0.90, 0.94, 0.96, 0.98, 0.99, 0.995, 1.0];
alpha_values = [0, 0.3, 0.5, 0.7, 0.9, 0.95];

rmse_lambda = nan(length(lambda_values), cfg.N_trials);
rmse_alpha = nan(length(alpha_values), cfg.N_trials);

for trial = 1:cfg.N_trials

    % The same noise realization is used for all parameter values.
    data = generate_test_data(cfg, trial);

    for i_lambda = 1:length(lambda_values)

        cfg_temp = cfg;
        cfg_temp.lambda = lambda_values(i_lambda);

        result = evaluate_rls(data, cfg_temp, cfg.f_true);

        rmse_lambda(i_lambda,trial) = ...
            result.metrics.rmse_steady;
    end

    for i_alpha = 1:length(alpha_values)

        cfg_temp = cfg;
        cfg_temp.alpha = alpha_values(i_alpha);

        result = evaluate_mmp(data, cfg_temp, cfg.f_true);

        rmse_alpha(i_alpha,trial) = ...
            result.metrics.rmse_steady;
    end
end

mean_lambda = mean(rmse_lambda, 2, 'omitnan');
std_lambda = std(rmse_lambda, 0, 2, 'omitnan');

mean_alpha = mean(rmse_alpha, 2, 'omitnan');
std_alpha = std(rmse_alpha, 0, 2, 'omitnan');

figure
tiledlayout(1,2)

nexttile
errorbar(lambda_values, mean_lambda, std_lambda, ...
    'o-', 'LineWidth', 1.2)
grid on
xlabel('\lambda')
ylabel('Steady-state RMSE')
title('RLS Forgetting Factor')

nexttile
errorbar(alpha_values, mean_alpha, std_alpha, ...
    'o-', 'LineWidth', 1.2)
grid on
xlabel('\alpha')
ylabel('Steady-state RMSE')
title('IIR-MMP Parameter')

[~, best_lambda_index] = min(mean_lambda);
[~, best_alpha_index] = min(mean_alpha);

cfg.lambda = lambda_values(best_lambda_index);
cfg.alpha = alpha_values(best_alpha_index);

fprintf('Selected lambda: %.4f\n', cfg.lambda)
fprintf('Selected alpha:  %.4f\n\n', cfg.alpha)


%% ================================================================
%  4. SNR sweep
%  ================================================================

SNR_values = [Inf, 30, 20, 10, 5, 0];

rmse_rls_snr = nan(length(SNR_values), cfg.N_trials);
rmse_mmp_snr = nan(length(SNR_values), cfg.N_trials);

for i_snr = 1:length(SNR_values)

    cfg_temp = cfg;
    cfg_temp.SNR_dB = SNR_values(i_snr);

    for trial = 1:cfg.N_trials

        data = generate_test_data(cfg_temp, trial);

        rls_result = evaluate_rls( ...
            data, cfg_temp, cfg_temp.f_true);

        mmp_result = evaluate_mmp( ...
            data, cfg_temp, cfg_temp.f_true);

        rmse_rls_snr(i_snr,trial) = ...
            rls_result.metrics.rmse_steady;

        rmse_mmp_snr(i_snr,trial) = ...
            mmp_result.metrics.rmse_steady;
    end
end

mean_rls_snr = mean(rmse_rls_snr, 2, 'omitnan');
std_rls_snr = std(rmse_rls_snr, 0, 2, 'omitnan');

mean_mmp_snr = mean(rmse_mmp_snr, 2, 'omitnan');
std_mmp_snr = std(rmse_mmp_snr, 0, 2, 'omitnan');

figure

x_snr = 1:length(SNR_values);

errorbar(x_snr, mean_rls_snr, std_rls_snr, ...
    'o-', 'LineWidth', 1.2)
hold on
errorbar(x_snr, mean_mmp_snr, std_mmp_snr, ...
    's-', 'LineWidth', 1.2)

xticks(x_snr)
xticklabels(compose('%g', SNR_values))
grid on

xlabel('SNR (dB)')
ylabel('Steady-state RMSE')
legend('RLS', 'IIR-MMP', 'Location', 'best')
title('Noise Robustness')


%% ================================================================
%  5. Amplitude modulation frequency sweep
%  ================================================================

f_amp_values = [0, 0.05, 0.1, 0.2, 0.5, 1.0];

rmse_rls_famp = nan(length(f_amp_values), cfg.N_trials);
rmse_mmp_famp = nan(length(f_amp_values), cfg.N_trials);

% Different components change at slightly different rates.
amp_rate_ratios = [1.0, 0.8, 1.2];

for i_amp = 1:length(f_amp_values)

    cfg_temp = cfg;
    cfg_temp.f_amp = ...
        f_amp_values(i_amp)*amp_rate_ratios;

    for trial = 1:cfg.N_trials

        data = generate_test_data(cfg_temp, trial);

        rls_result = evaluate_rls( ...
            data, cfg_temp, cfg_temp.f_true);

        mmp_result = evaluate_mmp( ...
            data, cfg_temp, cfg_temp.f_true);

        rmse_rls_famp(i_amp,trial) = ...
            rls_result.metrics.rmse_steady;

        rmse_mmp_famp(i_amp,trial) = ...
            mmp_result.metrics.rmse_steady;
    end
end

mean_rls_famp = mean(rmse_rls_famp, 2, 'omitnan');
std_rls_famp = std(rmse_rls_famp, 0, 2, 'omitnan');

mean_mmp_famp = mean(rmse_mmp_famp, 2, 'omitnan');
std_mmp_famp = std(rmse_mmp_famp, 0, 2, 'omitnan');

figure
errorbar(f_amp_values, mean_rls_famp, std_rls_famp, ...
    'o-', 'LineWidth', 1.2)
hold on
errorbar(f_amp_values, mean_mmp_famp, std_mmp_famp, ...
    's-', 'LineWidth', 1.2)
grid on

xlabel('Nominal amplitude modulation frequency (Hz)')
ylabel('Steady-state RMSE')
legend('RLS', 'IIR-MMP', 'Location', 'best')
title('Amplitude Tracking Performance')


%% ================================================================
%  6. Frequency-model mismatch sweep
%
%  The true signal uses cfg.f_true, while the reconstruction methods
%  receive f_model.
%  ================================================================

frequency_error_values = [0, 0.01, 0.05, 0.1, 0.25, 0.5];

rmse_rls_mismatch = ...
    nan(length(frequency_error_values), cfg.N_trials);

rmse_mmp_mismatch = ...
    nan(length(frequency_error_values), cfg.N_trials);

% Alternate the direction of the frequency errors.
frequency_error_sign = ones(size(cfg.f_true));
frequency_error_sign(2:2:end) = -1;

for i_error = 1:length(frequency_error_values)

    delta_f = frequency_error_values(i_error);

    f_model = cfg.f_true + ...
        delta_f*frequency_error_sign;

    for trial = 1:cfg.N_trials

        data = generate_test_data(cfg, trial);

        rls_result = evaluate_rls(data, cfg, f_model);
        mmp_result = evaluate_mmp(data, cfg, f_model);

        rmse_rls_mismatch(i_error,trial) = ...
            rls_result.metrics.rmse_steady;

        rmse_mmp_mismatch(i_error,trial) = ...
            mmp_result.metrics.rmse_steady;
    end
end

mean_rls_mismatch = ...
    mean(rmse_rls_mismatch, 2, 'omitnan');

std_rls_mismatch = ...
    std(rmse_rls_mismatch, 0, 2, 'omitnan');

mean_mmp_mismatch = ...
    mean(rmse_mmp_mismatch, 2, 'omitnan');

std_mmp_mismatch = ...
    std(rmse_mmp_mismatch, 0, 2, 'omitnan');

figure
errorbar(frequency_error_values, ...
    mean_rls_mismatch, std_rls_mismatch, ...
    'o-', 'LineWidth', 1.2)
hold on
errorbar(frequency_error_values, ...
    mean_mmp_mismatch, std_mmp_mismatch, ...
    's-', 'LineWidth', 1.2)
grid on

xlabel('Frequency-model error (Hz)')
ylabel('Steady-state RMSE')
legend('RLS', 'IIR-MMP', 'Location', 'best')
title('Sensitivity to Frequency-Model Error')


%% ================================================================
%  7. Number-of-frequencies sweep
%  ================================================================

number_frequency_values = [1, 2, 3, 5, 7];

% These frequencies have distinct aliases at f_ss = 20 Hz.
frequency_bank = 21:29;

rmse_rls_num_freq = ...
    nan(length(number_frequency_values), cfg.N_trials);

rmse_mmp_num_freq = ...
    nan(length(number_frequency_values), cfg.N_trials);

regressor_condition = ...
    nan(length(number_frequency_values), 1);

for i_num = 1:length(number_frequency_values)

    M = number_frequency_values(i_num);

    cfg_temp = cfg;

    cfg_temp.f_true = frequency_bank(1:M);
    cfg_temp.A0 = linspace(1.0, 0.5, M);
    cfg_temp.mod_depth = linspace(0.4, 0.15, M);
    cfg_temp.f_amp = linspace(0.05, 0.20, M);
    cfg_temp.phase = linspace(0, pi/2, M);

    check_frequency_aliases( ...
        cfg_temp.f_true, cfg_temp.f_fs/cfg_temp.L);

    for trial = 1:cfg.N_trials

        data = generate_test_data(cfg_temp, trial);

        if trial == 1
            regressor_condition(i_num) = ...
                calculate_regressor_condition( ...
                    data.t_slow, cfg_temp.f_true);
        end

        rls_result = evaluate_rls( ...
            data, cfg_temp, cfg_temp.f_true);

        mmp_result = evaluate_mmp( ...
            data, cfg_temp, cfg_temp.f_true);

        rmse_rls_num_freq(i_num,trial) = ...
            rls_result.metrics.rmse_steady;

        rmse_mmp_num_freq(i_num,trial) = ...
            mmp_result.metrics.rmse_steady;
    end
end

mean_rls_num_freq = ...
    mean(rmse_rls_num_freq, 2, 'omitnan');

std_rls_num_freq = ...
    std(rmse_rls_num_freq, 0, 2, 'omitnan');

mean_mmp_num_freq = ...
    mean(rmse_mmp_num_freq, 2, 'omitnan');

std_mmp_num_freq = ...
    std(rmse_mmp_num_freq, 0, 2, 'omitnan');

figure
errorbar(number_frequency_values, ...
    mean_rls_num_freq, std_rls_num_freq, ...
    'o-', 'LineWidth', 1.2)
hold on
errorbar(number_frequency_values, ...
    mean_mmp_num_freq, std_mmp_num_freq, ...
    's-', 'LineWidth', 1.2)
grid on

xlabel('Number of frequency components')
ylabel('Steady-state RMSE')
legend('RLS', 'IIR-MMP', 'Location', 'best')
title('Effect of Signal Complexity')

figure
semilogy(number_frequency_values, ...
    regressor_condition, 'o-', 'LineWidth', 1.2)
grid on

xlabel('Number of frequency components')
ylabel('Condition number of \Phi')
title('RLS Regressor Conditioning')


%% ================================================================
%  8. Create summary tables
%  ================================================================

snr_summary = table( ...
    SNR_values(:), ...
    mean_rls_snr, ...
    std_rls_snr, ...
    mean_mmp_snr, ...
    std_mmp_snr, ...
    'VariableNames', {
        'SNR_dB', ...
        'RLS_RMSE_Mean', ...
        'RLS_RMSE_Std', ...
        'MMP_RMSE_Mean', ...
        'MMP_RMSE_Std'
    });

amplitude_summary = table( ...
    f_amp_values(:), ...
    mean_rls_famp, ...
    std_rls_famp, ...
    mean_mmp_famp, ...
    std_mmp_famp, ...
    'VariableNames', {
        'AmplitudeFrequency_Hz', ...
        'RLS_RMSE_Mean', ...
        'RLS_RMSE_Std', ...
        'MMP_RMSE_Mean', ...
        'MMP_RMSE_Std'
    });

mismatch_summary = table( ...
    frequency_error_values(:), ...
    mean_rls_mismatch, ...
    std_rls_mismatch, ...
    mean_mmp_mismatch, ...
    std_mmp_mismatch, ...
    'VariableNames', {
        'FrequencyError_Hz', ...
        'RLS_RMSE_Mean', ...
        'RLS_RMSE_Std', ...
        'MMP_RMSE_Mean', ...
        'MMP_RMSE_Std'
    });

frequency_count_summary = table( ...
    number_frequency_values(:), ...
    mean_rls_num_freq, ...
    std_rls_num_freq, ...
    mean_mmp_num_freq, ...
    std_mmp_num_freq, ...
    regressor_condition, ...
    'VariableNames', {
        'NumberFrequencies', ...
        'RLS_RMSE_Mean', ...
        'RLS_RMSE_Std', ...
        'MMP_RMSE_Mean', ...
        'MMP_RMSE_Std', ...
        'RegressorCondition'
    });

disp('SNR summary')
disp(snr_summary)

disp('Amplitude variation summary')
disp(amplitude_summary)

disp('Frequency mismatch summary')
disp(mismatch_summary)

disp('Number of frequencies summary')
disp(frequency_count_summary)


%% ================================================================
%  Local functions
%  ================================================================

function data = generate_test_data(cfg, random_seed)
%GENERATE_TEST_DATA Generate aligned fast- and slow-sampled signals.

    validate_signal_configuration(cfg)

    T_fast = 1/cfg.f_fs;
    T_slow = cfg.L*T_fast;

    % Define the slow-sample count first so the fast and slow records
    % terminate at exactly corresponding samples.
    N_slow = floor(cfg.t_end/T_slow) + 1;
    N_fast = (N_slow - 1)*cfg.L + 1;

    t_fast = (0:N_fast-1)*T_fast;
    t_slow = t_fast(1:cfg.L:end);

    [y_fast, A_fast, y_components] = ...
        variable_amplitude_multisine( ...
            t_fast, ...
            cfg.f_true, ...
            cfg.A0, ...
            cfg.mod_depth, ...
            cfg.f_amp, ...
            cfg.phase);

    y_slow_clean = y_fast(1:cfg.L:end);

    rng(random_seed, 'twister')

    if isinf(cfg.SNR_dB)
        y_slow_noisy = y_slow_clean;
    else
        signal_power = mean(y_slow_clean.^2);

        noise_power = signal_power / ...
            10^(cfg.SNR_dB/10);

        measurement_noise = ...
            sqrt(noise_power)*randn(size(y_slow_clean));

        y_slow_noisy = ...
            y_slow_clean + measurement_noise;
    end

    data.t_fast = t_fast;
    data.t_slow = t_slow;

    data.y_fast = y_fast;
    data.y_slow_clean = y_slow_clean;
    data.y_slow_noisy = y_slow_noisy;

    data.A_fast = A_fast;
    data.y_components = y_components;
end


function result = evaluate_rls(data, cfg, frequencies)
%EVALUATE_RLS Run RLS and calculate reconstruction metrics.

    [y_hat, theta_history, y_components] = ...
        rls_multisine_reconstruction( ...
            data.y_slow_noisy, ...
            data.t_slow, ...
            data.t_fast, ...
            frequencies, ...
            cfg.lambda, ...
            cfg.delta);

    result.y_hat = y_hat;
    result.theta_history = theta_history;
    result.y_components = y_components;

    result.metrics = calculate_metrics( ...
        data.y_fast, ...
        y_hat, ...
        data.t_fast, ...
        cfg.t_eval_start);
end


function result = evaluate_mmp(data, cfg, frequencies)
%EVALUATE_MMP Run IIR-MMP and calculate reconstruction metrics.

    T_fast = 1/cfg.f_fs;

    % This assumes signal_recovery_iir accepts a frequency vector.
    y_hat = signal_recovery_iir( ...
        data.y_slow_noisy, ...
        frequencies, ...
        T_fast, ...
        cfg.alpha, ...
        cfg.L);

    y_hat = y_hat(:).';

    result.y_hat = y_hat;

    result.metrics = calculate_metrics( ...
        data.y_fast, ...
        y_hat, ...
        data.t_fast, ...
        cfg.t_eval_start);
end


function metrics = calculate_metrics( ...
    y_true, y_hat, t, evaluation_start)
%CALCULATE_METRICS Calculate full and steady-state errors.

    y_true = y_true(:).';
    y_hat = y_hat(:).';
    t = t(:).';

    N = min([
        length(y_true), ...
        length(y_hat), ...
        length(t)
    ]);

    y_true = y_true(1:N);
    y_hat = y_hat(1:N);
    t = t(1:N);

    valid = isfinite(y_true) & isfinite(y_hat);

    full_index = valid;
    steady_index = valid & t >= evaluation_start;

    if ~any(full_index)
        error('No valid samples are available for full-record metrics.')
    end

    if ~any(steady_index)
        error('No valid samples are available after t_eval_start.')
    end

    error_signal = y_true - y_hat;

    metrics.rmse_full = sqrt(mean( ...
        error_signal(full_index).^2));

    metrics.rmse_steady = sqrt(mean( ...
        error_signal(steady_index).^2));

    reference_rms = sqrt(mean( ...
        y_true(steady_index).^2));

    metrics.nrmse_steady = ...
        metrics.rmse_steady/max(reference_rms, eps);

    metrics.max_error_steady = max( ...
        abs(error_signal(steady_index)));
end


function [y, A, y_components] = ...
    variable_amplitude_multisine( ...
        t, frequencies, A0, modulation_depth, ...
        amplitude_frequencies, phase)
%VARIABLE_AMPLITUDE_MULTISINE Generate a variable-amplitude multisine.
%
% A_i(t) = A0_i * [1 + m_i*sin(2*pi*f_amp_i*t)]
%
% y(t) = sum_i A_i(t)*sin(2*pi*f_i*t + phase_i)

    output_is_column = iscolumn(t);

    t = t(:).';
    frequencies = frequencies(:);

    A0 = A0(:);
    modulation_depth = modulation_depth(:);
    amplitude_frequencies = amplitude_frequencies(:);
    phase = phase(:);

    M = length(frequencies);

    if length(A0) ~= M || ...
       length(modulation_depth) ~= M || ...
       length(amplitude_frequencies) ~= M || ...
       length(phase) ~= M

        error(['frequencies, A0, modulation_depth, ', ...
               'amplitude_frequencies, and phase must ', ...
               'have the same number of elements.'])
    end

    A = A0 .* ( ...
        1 + modulation_depth .* ...
        sin(2*pi*amplitude_frequencies*t));

    signal_angle = ...
        2*pi*frequencies*t + phase;

    y_components = A .* sin(signal_angle);
    y = sum(y_components, 1);

    if output_is_column
        y = y.';
        A = A.';
        y_components = y_components.';
    end
end


function [y_hat_fast, theta_history, y_components] = ...
    rls_multisine_reconstruction( ...
        y_slow, t_slow, t_fast, frequencies, ...
        lambda, delta)
%RLS_MULTISINE_RECONSTRUCTION
% Estimate sine and cosine coefficients for multiple known frequencies.
%
% Parameter ordering:
% theta = [a1; b1; a2; b2; ...; aM; bM]
%
% Signal model:
% y(t) = sum_i [
%   a_i*sin(2*pi*f_i*t) + b_i*cos(2*pi*f_i*t)
% ]

    output_is_row = isrow(t_fast);

    y_slow = y_slow(:);
    t_slow = t_slow(:);
    t_fast_column = t_fast(:);
    frequencies = frequencies(:);

    N_slow = length(y_slow);
    M = length(frequencies);
    N_parameters = 2*M;

    if length(t_slow) ~= N_slow
        error('y_slow and t_slow must have the same length.')
    end

    if lambda <= 0 || lambda > 1
        error('lambda must satisfy 0 < lambda <= 1.')
    end

    if delta <= 0
        error('delta must be positive.')
    end

    theta = zeros(N_parameters, 1);
    P = delta*eye(N_parameters);

    theta_history = zeros(N_parameters, N_slow);

    for k = 1:N_slow

        frequency_angle = ...
            2*pi*frequencies*t_slow(k);

        phi = zeros(N_parameters, 1);

        phi(1:2:end) = sin(frequency_angle);
        phi(2:2:end) = cos(frequency_angle);

        denominator = ...
            lambda + phi.'*P*phi;

        K = P*phi/denominator;

        prediction = phi.'*theta;
        prediction_error = y_slow(k) - prediction;

        theta = theta + K*prediction_error;

        P = (P - K*phi.'*P)/lambda;

        % Reduce numerical asymmetry.
        P = 0.5*(P + P.');

        theta_history(:,k) = theta;
    end

    % Hold each parameter estimate between slow samples.
    theta_fast = interp1( ...
        t_slow, ...
        theta_history.', ...
        t_fast_column, ...
        'previous', ...
        'extrap').';

    fast_angle = ...
        2*pi*frequencies*t_fast_column.';

    sine_coefficients = theta_fast(1:2:end,:);
    cosine_coefficients = theta_fast(2:2:end,:);

    y_components = ...
        sine_coefficients.*sin(fast_angle) + ...
        cosine_coefficients.*cos(fast_angle);

    y_hat_fast = sum(y_components, 1);

    if ~output_is_row
        y_hat_fast = y_hat_fast.';
        y_components = y_components.';
    end
end


function condition_number = ...
    calculate_regressor_condition(t_slow, frequencies)
%CALCULATE_REGRESSOR_CONDITION Condition number of the RLS regressor.

    t_slow = t_slow(:);
    frequencies = frequencies(:);

    M = length(frequencies);
    Phi = zeros(length(t_slow), 2*M);

    for i = 1:M
        Phi(:,2*i-1) = ...
            sin(2*pi*frequencies(i)*t_slow);

        Phi(:,2*i) = ...
            cos(2*pi*frequencies(i)*t_slow);
    end

    condition_number = cond(Phi);
end


function check_frequency_aliases(frequencies, f_slow)
%CHECK_FREQUENCY_ALIASES Warn about duplicate or degenerate aliases.

    frequencies = frequencies(:);

    signed_alias = mod( ...
        frequencies + f_slow/2, f_slow) - f_slow/2;

    alias_magnitude = abs(signed_alias);

    fprintf('Frequency and alias mapping\n')
    fprintf('---------------------------\n')

    for i = 1:length(frequencies)
        fprintf('%8.3f Hz -> %8.3f Hz\n', ...
            frequencies(i), signed_alias(i))
    end

    tolerance = 1e-8;

    for i = 1:length(alias_magnitude)
        for j = i+1:length(alias_magnitude)

            if abs(alias_magnitude(i) - ...
                   alias_magnitude(j)) < tolerance

                warning(['Frequencies %.4f Hz and %.4f Hz have ', ...
                         'the same slow-rate alias magnitude. ', ...
                         'RLS cannot uniquely separate them.'], ...
                         frequencies(i), frequencies(j))
            end
        end
    end

    if any(alias_magnitude < tolerance)
        warning(['At least one frequency aliases to DC. ', ...
                 'Its sine and cosine coefficients are not ', ...
                 'both independently identifiable.'])
    end

    if any(abs(alias_magnitude - f_slow/2) < tolerance)
        warning(['At least one frequency aliases to the slow ', ...
                 'Nyquist frequency. The RLS basis may be ', ...
                 'rank deficient.'])
    end

    fprintf('\n')
end


function validate_signal_configuration(cfg)
%VALIDATE_SIGNAL_CONFIGURATION Check signal-vector dimensions.

    M = length(cfg.f_true);

    parameter_lengths = [
        length(cfg.A0), ...
        length(cfg.mod_depth), ...
        length(cfg.f_amp), ...
        length(cfg.phase)
    ];

    if any(parameter_lengths ~= M)
        error(['f_true, A0, mod_depth, f_amp, and phase ', ...
               'must have the same number of elements.'])
    end

    if any(cfg.mod_depth < 0) || any(cfg.mod_depth >= 1)
        error('mod_depth should satisfy 0 <= mod_depth < 1.')
    end

    if cfg.L < 1 || cfg.L ~= round(cfg.L)
        error('This script assumes L is a positive integer.')
    end

    if cfg.t_eval_start >= cfg.t_end
        error('t_eval_start must be smaller than t_end.')
    end
end