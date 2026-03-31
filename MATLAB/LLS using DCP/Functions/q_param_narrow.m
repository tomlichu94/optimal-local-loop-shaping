function [Qnarrow q0] = q_param_narrow(Pz,w_d,Ts)
% based on (7) from the LLS for Rejecting Band-Limited Disturbances in NMP
% Sys with Appliction to Laser Beam Steering for AM

[w_rad Pz_mag Pz_re Pz_im] = mag_phase_DT(Pz,w_d,Ts); % find output mag and phase at w_i
n = max(size(w_d));
m = 2*n-1;
z = tf('z',Ts);

q_vec = [];
Mk_o = cos(kron(w_rad',(1:m)));
Mk_e = sin(kron(w_rad',(1:m)));
Mk = zeros(2*n,2*n);
for k = 1:n
    Mk(2*k-1,1) = 1;
    Mk(2*k,2:end) = Mk_e(k,:);
    Mk(2*k-1,2:end) = Mk_o(k,:);
    q_vec(2*k-1) = Pz_re(k)/(Pz_mag(k)^2);
    q_vec(2*k) = Pz_im(k)/(Pz_mag(k)^2);
end

q0 = Mk\q_vec';
m = max(size(q0));
Qnarrow = 0;
for k = 1:m
    Qnarrow = q0(k)*z^(-k+1)+Qnarrow;
    Qnarrow = minreal(Qnarrow);
end

% legacy code, replaced with more efficient operations
% Mk_o{1} = 1;
% Mk_e{1} = 0;
% syms W
% for k = 1:m
%     Mk_o{k+1} = cos(k*W);
%     Mk_e{k+1} = sin(k*W);
% end
% Mk = [];
% q_vec = [];
% for k = 1:n
%     Mk(2*k-1,:) = double(subs(cell2sym(Mk_o),W,w_rad(k)));
%     Mk(2*k,:) = double(subs(cell2sym(Mk_e),W,w_rad(k)));
%     q_vec(2*k-1) = Pz_re(k)/(Pz_mag(k)^2);
%     q_vec(2*k) = Pz_im(k)/(Pz_mag(k)^2);
% endend
