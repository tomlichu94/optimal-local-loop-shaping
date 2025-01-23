function [d_est] = reconst_signal_fir(w_k,d_ss,batches,L)
k = 1:(L-1);
n_w = size(w_k,1)-1;
% batch one is incomplete due to not enough data available
% d_mat1 = zeros(n_w+1,n_w+1);
% for i = 1:(n_w)
%     d_mat1 = d_mat1+diag(ones(n_w+2-i,1)*d_ss(i),(1-i));
% end
% d_mat1(end,:) = []; % remove last row, which is not needed;
d_mat1 = zeros(n_w,n_w+1);
d_mat2 = [];
n_dss = length(d_ss); % find length of d_ss
for i = (n_w+1):n_dss % information of d_ss similar to hankel matrix
    d_mat2(i-n_w,:) = flip(d_ss((i-n_w):(i)));
end
d_mat = [d_mat1; d_mat2];
d_est = [];
d_est(1) = d_mat(1,:)*w_k(:,1); % first row
for i = 1:(batches-1)
     d_est(i*L+1) = d_ss(i+1); % slow and fast match
     for p = 1:k(end) % intersample loop
         d_est(i*L+1+k(p)) = d_mat(i+1,:)*w_k(:,k(p));
     end
 end
d_est(batches*L+1) = d_ss(end);