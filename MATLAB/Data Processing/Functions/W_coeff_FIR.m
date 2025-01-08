function [w_k] = W_coeff_FIR(L,f_d,T_s)
% T_s is the fast sampling rate
k = 1:(L-1);
m_d = max(size(f_d));
Apara = 1; % defined as the coefficients for A(z) = 1 + a_1 z^-1 + ... + z^-2m
for i = 1:m_d
    Apara = conv(Apara,[1 -2*cos(f_d(i)*2*pi*T_s) 1]);
end

%%%%%%%%% defining M_k and M_kt (t for tilde) %%%%%%%%%
n_w = 2*m_d-1;
m_k1 = 2*m_d*L; % rows of M_kt
m_k2 = 2*m_d*(L-1); % columns of M_kt
m_a = 2*m_d+1; % number of coefficients for A(z)
M_kt = zeros(m_k1,m_k2);
for i = 1:m_k2
    M_kt(i:(i+m_a-1),i) = Apara';
end

% create row vectors for e_k with the kth row = 1
E_k = zeros(m_k1,2*m_d);
for i = 1:max(k)
    for j = 1:2*m_d
    E_k(i+L*(j-1),j,i) = 1;
    end
end

% build M_k matrix, where there are k matrices nested within M_k
M_k = zeros(m_k1,m_k1);
for i = 1:max(k)
M_k(:,1:m_k2,i) = M_kt;
M_k(:,(m_k2+1):end,i) = E_k(:,:,i);
end

% the signal reconstruction is M_k * [f_vec w_vec]' = [a_vec 0]'
% let us express this as Ax = b
a_sol = [-Apara(2:end)'; zeros(m_k1-m_a+1,1)];

for i = 1:max(k)
    x_sol(:,i) = pinv(M_k(:,:,i))*a_sol;
    f_max = 2*m_d*(L-1); % highest coefficient for F
    w_max = m_k1 - f_max-1; % highest coefficient for W based on f_max
    f_k = x_sol(1:(m_k1-(n_w+1)),:);
    w_k = x_sol((m_k1-n_w):end,:);
end
% if w_max >= n_w
%     fprintf('Reconstruction is feasible\n')
% else
%     fprintf('Reconstruction unfeasible\n')
%     return
% end