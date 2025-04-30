function [w_rad Pz_mag Pz_re Pz_im] = mag_phase_DT(Pz,w,Ts)
% input frequency is in rad/s. Use Ts == 0 if in rad/s. If not, specify Ts
if Ts == 0;
    w_rad = w;
else
    w_rad = w*(2*pi*Ts); % convert to rad/s
end

syms Z
Pnum = poly2sym(cell2mat(Pz.num),Z);
Pden = poly2sym(cell2mat(Pz.den),Z);
Pz_sym = simplifyFraction(Pnum/Pden); % symbolic transfer function

n = max(size(w_rad));
Pz_out = [];
Pz_mag = [];
for k = 1:n
    Pz_out(k) = double(subs(Pz_sym,Z,exp(w_rad(k)*j))); % sub in P(Z = e^jw),
    Pz_mag(k) = norm(Pz_out(k),2); % magnitude
end
Pz_re = real(Pz_out);
Pz_im = imag(Pz_out);

end