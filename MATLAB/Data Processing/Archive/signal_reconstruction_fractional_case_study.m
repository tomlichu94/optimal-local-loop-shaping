clear all
close all
clc
addpath('Functions')

%% note on this function

% this signal reconstruction is for the specific case where L = 2.5, so
% R = 2 to make the sampling feasible by integers. We must scale the delay
% by R*L = 5, meaning there are 4 intersample points between two slow
% measurements. Because L is not an integer, the smallest integer we can
% scale up to from 2.5 is 5. Therefore, we can only use odd indexes such as
% d_ss[rn+1] and d_ss[r(n+1)+1] for signal recovery (e.g. d_ss[1] and 
% d_ss[3], which as 4 intersample points). We can also use d_ss[rn] and
% d_ss[r(n+1)] for the even indexes and do a separate signal recovery,
% (e.g. d_ss[2] and d_ss[4]). 
% 
% Combining these two recovered signals, we can double the sampling rate 
% to be 5x faster than the slow sampled rate. We also only need to solve
% for W_k for k = 1...4 once because we can just shift the slow sampled
% signal by one time step to make it an even or odd index.

%% sampling parameters
T_ss = 1/10; % slow sampling rate
L = 2.5; % ratio of fast to slow sampling
T_fs = T_ss/L; % fast rate
R = 2; % want R*L to be the smallest intger (2 in this case)
L_max = R*L; % # intersamples-1 between two slow-sampled points
T_comb = T_ss/(R*L); % combined sample rate

% signal recovery
f_target = 7; % siganl frequency
T_end = 3; % signal end time

% sin generator
[y_ct, y_ss, y_frac, t_ct, t_ss, t_fs] = sin_generator(T_fs,T_ss,f_target,T_end);
[y_ct, y_ss, y_comb, t_ct, t_ss, t_comb] = sin_generator(T_comb,T_ss,f_target,T_end);

%% signal reconstruction parameters
[w_k, M_k] = w_k_frac_odd(f_target,T_fs,R*L);
% [w_k_even, M_k_even] = w_k_frac_even(f_target,T_fs,R*L); % even index
% [w_k_test] = w_k_fir(L*R,f_target,T_fs); % validate w_k
length_odd = length(y_frac);
length_slow = length(y_ss);

%% odd indexing reconstruction
y_est_odd = zeros(1,length_odd);
n_ss = 1;
n_w = 2; 
phi = zeros(1, n_w); % Preallocate phi using the same type as d_ss
for n_fs = 1:length_odd
    k = mod((n_fs-1),L_max); % check if on intersample
    if n_ss < n_w+1
        if k == 0 % Fast and slow align
            phi = [y_ss((n_ss-1)*R+1), phi(1:end-1)]; % Shift and update phi
            y_est_odd(n_fs) = 0;
            n_ss = n_ss+1;
        else % Not enough data for intersample estimation
            y_est_odd(n_fs) = 0;
        end
    % Case 2: n_ss >= n_w (enough data points for FIR)
    else
        if k == 0
            phi = [y_ss((n_ss-1)*R+1), phi(1:end-1)];
            y_est_odd(n_fs) = y_ss((n_ss-1)*R+1);
            n_ss = n_ss+1;
        else % FIR reconstruction
            y_est_odd(n_fs) = phi * w_k(:, k); 
        end
    end
end

%% even indexing reconstruction using shifted slow sample
length_even = length_odd-R-1;
y_ss_even = y_ss(2:end);
y_est_even = zeros(1,length_even);
n_ss = 1;
n_w = 2; 
phi = zeros(1, n_w); % Preallocate phi using the same type as d_ss
for n_fs = 1:length_even
    k = mod((n_fs-1),L_max); % check if on intersample
    if n_ss < n_w+1
        if k == 0 % Fast and slow align
            phi = [y_ss_even((n_ss-1)*R+1), phi(1:end-1)]; % Shift and update phi
            y_est_even(n_fs) = 0;
            n_ss = n_ss+1;
        else % Not enough data for intersample estimation
            y_est_even(n_fs) = 0;
        end
    % Case 2: n_ss >= n_w (enough data points for FIR)
    else
        if k == 0
            phi = [y_ss_even((n_ss-1)*R+1), phi(1:end-1)];
            y_est_even(n_fs) = y_ss_even((n_ss-1)*R+1);
            n_ss = n_ss+1;
        else % FIR reconstruction
            y_est_even(n_fs) = phi * w_k(:, k); 
        end
    end
end

%% even indices shifted slow sampling
% clear n_fs n_ss phi
% y_ss_even = y_ss(2:end);
% t_ss_even = t_lin_ss(2:end);
% slow_length_even = length(y_ss_even);
% % fast_length_even = L*R*(slow_length_even-2)/2+1;
% fast_length_even = 23;
% y_est_even = zeros(1,fast_length_even);
% n_ss = 1;
% n_w = 2; 
% phi = zeros(1, n_w); % Preallocate phi using the same type as d_ss
% for n_fs = 1:fast_length_even
%     k = mod((n_fs-1),L_max); % check if on intersample
%     if n_ss < n_w+1
%         if k == 0 % Fast and slow align
%             phi = [y_ss_even((n_ss-1)*R+1), phi(1:end-1)]; % Shift and update phi
%             y_est_even(n_fs) = 0;
%             n_ss = n_ss+1;
%         else % Not enough data for intersample estimation
%             y_est_even(n_fs) = 0;
%         end
%     % Case 2: n_ss >= n_w (enough data points for FIR)
%     else
%         if k == 0
%             phi = [y_ss_even((n_ss-1)*R+1), phi(1:end-1)];
%             y_est_even(n_fs) = y_ss_even((n_ss-1)*R+1);
%             n_ss = n_ss+1;
%         else % FIR reconstruction
%             y_est_even(n_fs) = phi*w_k_even(:, k); 
%         end
%     end
% end 

% this is for the even indices
% clear n_fs n_ss n_count phi
% slow_length_even = slow_length-1; % account for starting 1 slow step later
% fast_length_even = slow_length_even*L-R-1;
% y_est_even = zeros(1,fast_length_even);
% n_ss = 1;
% n_count = 1;
% n_w = 2; 
% phi = zeros(1, n_w); % Preallocate phi using the same type as d_ss
% for n_fs = 1:fast_length_even
%     k = mod((n_fs-1),L_max); % check if on intersample
%     if n_count < n_w+1
%         if k == 0 % Fast and slow align
%             phi = [y_ss(n_ss*R), phi(1:end-1)]; % Shift and update phi
%             y_est_even(n_fs) = 0;
%         else % Not enough data for intersample estimation
%             y_est_even(n_fs) = 0;
%         end
%     % Case 2: n_ss >= n_w (enough data points for FIR)
%     else
%         if k == 0
%             phi = [y_ss(n_ss*R), phi(1:end-1)];
%             y_est_even(n_fs) = y_ss(n_ss*R);
%         else % FIR reconstruction
%             y_est_even(n_fs) = phi * w_k_even(:, k); 
%         end
%     end
%     if k == 0
%         n_ss = n_ss + 1; % Number of slow samples
%         n_count = n_count + 1;
%     end
% end

%% Plotting
t_lin_fs_even = T_ss:T_fs:T_end; % skipped the first slow sampling point
length_combined = R*L*(length_slow-1)+1;
y_est_fs = zeros(1,length_combined);
y_est_fs(1:R:end) = y_est_odd;
y_est_fs((R*L+1):R:end) = y_est_even;

figure
plot(t_ct,y_ct)
hold on
s(1) = stairs(t_ss,y_ss,'diamond-');
s(2) = stairs(t_fs,y_est_odd,'x');
s(3) = stairs(t_lin_fs_even,y_est_even,'o');
for i = 1:3
    s(i).MarkerSize = 8;
    s(i).LineWidth = 1;
end
% stairs(t_lin_ss,y_ss,'o')
legend('Ground truth','Slow-Sampling','Odd Indexing','Even Indexing')

figure
plot(t_ct,y_ct)
hold on
s(1) = stairs(t_comb,y_comb,'o');
s(2) = stairs(t_comb,y_est_fs,'square-');
for i = 1:2
    s(i).MarkerSize = 8;
    s(i).LineWidth = 1;
end
legend('Ground truth','Fast Sampling','Combined Signal Recovery')

%% this is for MATLAB's even indexing - not needed
function [w_k, M_k] = w_k_frac_even(f_d,T_s,rL)
% test function create M_k matrix and solve
    m_d = length(f_d);
    Apara = 1; % defined as the coefficients for A(z) = 1 + a_1 z^-1 + ... + z^-2m
    for i = 1:m_d
        Apara = conv(Apara,[1 -2*cos(f_d(i)*2*pi*T_s) 1]);
    end
% Apara = flip(Apara); % flips from [z^n ... z 1] to [1 z^-1 ... z^-n)]
    n_w = 2*m_d-1;
    m_k1 = 2*m_d*rL; % rows of M_kt
    m_k2 = 2*m_d*(rL-1); % columns of M_kt
    m_a = 2*m_d+1; % number of coefficients for A(z)
    M_kt = zeros(m_k1,m_k2);
    b = zeros(m_k1,1);
    b(1:2) = -Apara(2:3)';
    w_k = zeros(2*n_w,rL-1);
    for i = 1:m_k2
        M_kt(i:(i+m_a-1),i) = Apara';
    end
% create row vectors for e_k with the kth row = 1
    E_k = zeros(m_k1,2*m_d,rL-1);
    for i = 1:(rL-1)
        for j = 1:2*m_d
            E_k(i+1+rL*(j-1),j,i) = 1;
        end
    end
    M_k = zeros(m_k1,m_k1);
    for i = 1:(rL-1)
        M_k(:,1:m_k2,i) = M_kt;
        M_k(:,(m_k2+1):end,i) = E_k(:,:,i);
        w_k1 = pinv(M_k(:,:,i))*b;
        w_k1 = w_k1(end-1:end);
        w_k(:,i) = w_k1;
    end
end

%% this is for MATLAB's odd indexing 
function [w_k, M_k] = w_k_frac_odd(f_d,T_s,rL)
% test function create M_k matrix and solve
    m_d = length(f_d);
    Apara = 1; % defined as the coefficients for A(z) = 1 + a_1 z^-1 + ... + z^-2m
    for i = 1:m_d
        Apara = conv(Apara,[1 -2*cos(f_d(i)*2*pi*T_s) 1]);
    end
% Apara = flip(Apara); % flips from [z^n ... z 1] to [1 z^-1 ... z^-n)]
    n_w = 2*m_d-1;
    m_k1 = 2*m_d*rL; % rows of M_kt
    m_k2 = 2*m_d*(rL-1); % columns of M_kt
    m_a = 2*m_d+1; % number of coefficients for A(z)
    M_kt = zeros(m_k1,m_k2);
    b = zeros(m_k1,1);
    b(1:2) = -Apara(2:3)';
    w_k = zeros(2*n_w,rL-1);
    for i = 1:m_k2
        M_kt(i:(i+m_a-1),i) = Apara';
    end
% create row vectors for e_k with the kth row = 1
    E_k = zeros(m_k1,2*m_d,rL-1);
    for i = 1:(rL-1)
        for j = 1:2*m_d
            E_k(i+rL*(j-1),j,i) = 1;
        end
    end
    M_k = zeros(m_k1,m_k1);
    for i = 1:(rL-1)
        M_k(:,1:m_k2,i) = M_kt;
        M_k(:,(m_k2+1):end,i) = E_k(:,:,i);
        w_k1 = pinv(M_k(:,:,i))*b;
        w_k1 = w_k1(end-1:end);
        w_k(:,i) = w_k1;
    end
end

%% other functions
function [y_ct, y_ss, y_fs, t_lin_ct, t_lin_ss, t_lin_fs] = sin_generator(t_fs, t_ss, f_hz, t_end)
    t_ct = t_fs/100;
    t_lin_ct = 0:t_ct:t_end;
    t_lin_ss = 0:t_ss:t_end;
    t_lin_fs = 0:t_fs:t_end;
    y_ct = sin(2*pi*f_hz*t_lin_ct);
    y_fs = sin(2*pi*f_hz*t_lin_fs);
    y_ss = sin(2*pi*f_hz*t_lin_ss);
end
