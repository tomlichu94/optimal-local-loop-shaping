function [W_num W_den] = W_TF_IIR(w_k,B_para)
% transfer function for W_TF = d[nL+k]/d[nL], based on the fast sampling
m_d = height(w_k)/2;
k_all = width(w_k);
L = k_all+1;
W_num = zeros(k_all,2*m_d*L+1);
W_den = zeros(k_all,2*m_d*L+1);
for k = 1:k_all
    W_num(k,1:L:2*m_d*L-1) = w_k(:,k)';
    W_den(k,1:L:end) = B_para;
end
% W_TF = tf(W_num,W_den,T_s);