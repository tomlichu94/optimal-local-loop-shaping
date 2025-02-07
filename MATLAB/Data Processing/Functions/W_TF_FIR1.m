function [W_num W_den] = W_TF_FIR2(w_k)
% transfer function for W_TF = d[nL+k]/d[nL], based on the fast sampling
m_d = height(w_k)/2;
k_all = width(w_k);
L = k_all+1;
W_num = zeros(k_all,(2*m_d-1)*L+1);
W_den = [ones(k_all,1), zeros(k_all,(2*m_d-1)*L)];
for k = 1:k_all
    W_num(k,1:L:end) = w_k(:,k)';
end
