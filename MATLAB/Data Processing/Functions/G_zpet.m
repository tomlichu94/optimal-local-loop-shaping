% function [num_zpet den_zpet] = G_zpet(num,den,Ts)
% function to find the inverse G with zero phase error tracking

% the plant
% G(z) = (b_m z^m + ... + b_1 z + b_0)/(z^n + a_{n-1} z^{n-1} + ... + a_0)
% order G(z) = n, relative degree, d = n-m

% test code
den = poly([0.5 1 0]);
num = poly([-0.5 0.8 1.2]);
% end test code

n_den = max(size(den));
m_num = max(size(num));

% state casuality
if n_den > m_num
    fprintf('Strictly casual')
elseif n_den == m_num
    fprintf('Casual')
else
    fprintf('Noncasual, function halted')
    return
end

d_deg = n_den - m_num; % relative degree
G_zeros = roots(num); % find zeros of G
NMP_zeros = G_zeros(abs(G_zeros)>1)';
MP_zeros = G_zeros(abs(G_zeros)<=1)';

A_c = den;
B_neg = poly(NMP_zeros);
B_pos = poly(MP_zeros);