function [d_est_IIR] = reconst_signal_iir(w_k,Bpara,d_ss,batches,L)
k = 1:(L-1);
n_w = height(w_k)-1;
d_mat1 = zeros(n_w,n_w+1);
d_mat2 = [];
n_dss = length(d_ss); % find length of d_ss
m_d = (n_w+1)/2;
% puts this matrix similar to Hankel matrix form
for i = (n_w+1):n_dss
    d_mat2(i-n_w,:) = flip(d_ss((i-n_w):(i)));
end
d_mat = [d_mat1; d_mat2];
d_est_IIR = [];
d_est_IIR(1) = d_mat(1,:)*w_k(:,1); % first row
for i = 1:(batches-1)
     d_est_IIR(i*L+1) = d_ss(i+1); % slow and fast match
     for p = 1:k(end) % intersample loop
         if i <= 2*m_d
         d_est_IIR(i*L+1+k(p)) = d_mat(i+1,:)*w_k(:,k(p));
         else
         d_est_IIR(i*L+1+k(p)) = d_mat(i+1,:)*w_k(:,k(p))-...
                                 flip(d_est_IIR(((i-2*m_d)*L+1+k(p)):L:((i-1)*L+1+k(p))))*Bpara(2:end)';             
         end
     end
 end
d_est_IIR(batches*L+1) = d_ss(end);