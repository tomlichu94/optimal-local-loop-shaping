function [Fz_num, Fz_den] = F_notch(w_d,alpha)
% outputs F_nf = \prod_i^n_nf (1-2*cos(w_d(i))*z^(-1)+z^(-2))/...
%                 ((1-2*alpha*cos(w_d(i))*z^(-1)+alpha^2*z^(-2)))
fprintf('Input frequency must be in rad\n')
m_d = size(w_d,2);
Fz_num = 1;
Fz_den = 1;
for i = 1:m_d
    Fz_num = conv(Fz_num,[1 -2*cos(w_d(i)) 1]);
    Fz_den = conv(Fz_den,[1 -2*alpha*cos(w_d(i)) alpha^2]);
end
end