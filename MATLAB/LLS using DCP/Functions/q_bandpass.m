function [Qbp Q_bp1] = q_bandpass(Pz,w,B_bw,Ts)
% based on (16) from the LLS for Rejecting Band-Limited Disturbances in NMP
% Sys with Appliction to Laser Beam Steering for AM
[w_rad Pz_mag Pz_re Pz_im] = mag_phase_DT(Pz,w,Ts); % find output mag and phase at w_i
[Qnarrow] = q_param_narrow(Pz,w,Ts);

B_bw_rad = B_bw*2*pi*Ts;
z = tf('z',Ts);
n = max(size(w_rad));
Q_bp1 = 1;
for i = 1:n    
    k1 = -cos(w_rad(i));
    k2 = (1-tan(B_bw_rad(i)/2))/(1+tan(B_bw_rad(i)/2));
    Q_bp1 = 1/2*(1+k2)*(1+2*k1*z^-1+z^-2)/(1+k1*(1+k2)*z^-1+k2*z^-2)*Q_bp1;
    Q_bp1 = minreal(Q_bp1);
end
Q_bp = 1-Q_bp1;
Qbp = minreal(Q_bp*Qnarrow);
