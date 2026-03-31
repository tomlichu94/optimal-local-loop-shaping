clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% notes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for generalized number of disturbances and selected L_t
% and L_t sampling multiplier (T_{ss} = L_t*T_{fs})
% test showing the signal reconstruction built in SIMULINK with CVX using IIR
%%%%%%%%%%%%%%%%%%%%%%%%%% load plant data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('Functions')
cvx_precision low
cvx_save_prefs

%%%%%%%%%%%%%%%%% Plant parameters of the galvo scanner %%%%%%%%%%%%%%%%%%%
% continuous time plant with delay
Ps = tf(3.937e9, [1, 565.5, 3.198e5], 'InputDelay', 1e-5); 

% discrete time plant
f_hz = 5280; % sampling frequency
Tu = 1/f_hz; % sampling time
Pz = c2d(Ps,Tu,'zoh'); % discretized the plant
Pz = absorbDelay(Pz); % nest delay term into transfer function

% store num and den coefficients
Pz_num = cell2mat(Pz.num); 
Pz_den = cell2mat(Pz.den);

%%%%%%%%%%%%%%%%%%%% Q(z), filter Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_order = 30; % this is the value of "r" from Eq. (14), (33)
a_q = 0.90; % QIIR \alpha from Eq. (32). Adjusts notch bandwidth

%%%%%%%%%%%%%%%%%%%% narrow-band disturbance frequency %%%%%%%%%%%%%%%%%%%%
% generates random disturbance freq
% will generate and pick multiplier based on the number of disturbances
% spaced roughly equidistance. e.g. m_d = 3, L_t = 2. Will pick multiplier
% between 1 and 2. 1st between 1-1.33, 2nd between 1.333-1.67, 3rd between
% 1.67-2
L_t = 2; % fast sampling output of the plant
Ts = Tu/L_t; % slow sampling speed of the feedback sensor
m_d = 1; % number of disturbances
tempW = zeros(1,m_d); % temp variable
for i = 1:(m_d) % for loop generates random disturbance frequencies
    tempW(i) = 1+rand(1)*((L_t-1)/m_d)+(i-1)*(L_t-1)/m_d; 
end
w_d = tempW*pi/L_t; % fast measurement of disturbance in radians
f_d = w_d/(2*pi*Tu); % disturbance in Hz

%% %%%%%%%%%%%%%%%%%%%%% discrete set of frequencies %%%%%%%%%%%%%%%%%%%%%%
% based on S. Boyds' paper on approximating DT systems as CT systems
range_mult = 30;

% generate bandpass frequencies
w_lin = w_lin_spacing(max_order,range_mult,tempW,L_t);  % see Fig. 3
w_lin_Hz = w_lin/(2*pi*Tu); % in Hz

% %% %%%%%%%%%%%%%%%%%%%%%%%%% hdd model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % defining CT system with delay
% Ps_pade = pade(Ps,4); % pade approximation for delay term
% Ps_sys = ss(Ps_pade); % state space of the CT system

%% %%%%%%%%%%%%%%%%%%%%%%%%% PQ Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = tf('z',Tu); % discrete time based on fast sampling

% parameters for bandpass constraint 1 - PQ = 1
phi = tf(); 
phi(1) = 1;
for i = 1:max_order % Eq. (15), phi = [1, z^-1, ..., z^-r]
    phi(i+1) = z^(-i); 
end
P_phi = Pz*phi; % = [P*phi(1), ... ,P*phi(r)], part of Eq. (17)

% parameters for bandstop constraint 
PQ_d = freqresp(P_phi,exp(1j*w_d)); % subs in z = e^{jw_d}, w_d is the disturbance frequency

%%  %%%%%%%%%% QFIR SDP implementation, see Eq (21) - (29) %%%%%%%%%%%%%%%%

% state space realization of P(z), Eq. (23)
q_dim = size(Pz_den,2); % dimension for state space of P(z)
% A matrix
A_p = zeros(q_dim-1,q_dim-1); 
A_p(1:end-1,2:end) = eye(q_dim-2);
A_p(end,:) = -flip(Pz_den(2:end));
% B matrix
B_p = zeros(q_dim-1,1);
B_p(end) = 1;
% C matrix
C_p = flip(Pz_num(2:end));
D_p = Pz_num(1);

% state space representation of Q(z), Eq. (26)
% A matrix
A_q = zeros(max_order,max_order);
A_q(1:(end-1),2:end) = eye(max_order-1,max_order-1);
% B matrix
B_q = zeros(max_order,1);
B_q(end) = 1;

% real bounded lemma set up, Eq. (29), \tilde{A}  = [A_p, zero_mat; B_qC_p, A_q] 
zero_mat = zeros(size(A_p,1),size(A_q,2)); 
M_size = size(A_p)+size(A_q); % size of M in Eq. (20) 

% DCP using SDP
clear q
cvx_begin sdp
        % q -> Eq. (14), gamma, M -> Eq. (20)
        variables q((max_order+1),1) gamma_sdp
        variable M(M_size) symmetric
        
        % state space of Q(z) for C and D, stores Q(z) coeff int C and D
        C_q = [flip(q(2:end))]'; 
        D_q = [q(1)];

        % \tilde{A, B, C, D} formulation, see Eq. (29)
        A_t = [A_p zero_mat; B_q*C_p A_q];
        B_t = [B_p; B_q*D_p];
        C_t = -[D_q*C_p, C_q];
        D_t = 1-D_q*D_p;

        % this is the matrix in Eq. (20)
        L_cvx = [M, M*A_t, M*B_t, zeros(M_size(1),1);...
             A_t'*M, M, zeros(M_size(1),1), C_t';...
             B_t'*M, zeros(1,M_size(1)), -gamma_sdp, D_t';...
             zeros(1,M_size(1)), C_t, D_t, -gamma_sdp]; 

        % SDP formulation
        minimize gamma_sdp % objective function

        % constraint functions
        subject to
            % bandstop constraint
            for i = 1:m_d
               1-PQ_d(1,:,i)*q == 0;
            end

            % real-bounded lemma constraints, Eq. (20)
            L_cvx <= 0;
            gamma_sdp >= 0;
cvx_end

% Q filter construction
Qcvx_num = q'; % coefficients of Eq. (14)
Qcvx_den = [1 zeros(1,max_order)];
PQcvx_FIR_num = conv(Pz_num,Qcvx_num); %
PQcvx_FIR_den = conv(Pz_den,Qcvx_den);

% coefficients of P(z)Q(z)
PQcvx_num = conv(Pz_num,Qcvx_num);
PQcvx_den = conv(Pz_den,Qcvx_den);

% transfer function for P(z)Q(z)
PQ_sdp = tf(PQcvx_num,PQcvx_den,Tu);

% transfer function of 1-P(z)Q(z), Eq. (10)
T_sdp_fir = minreal(1-PQ_sdp);

%%  %%%%%%%%%% QIIR SDP implementation, see Eq (31) - (37) %%%%%%%%%%%%%%%%

[Vz_num, Vz_den] = F_notch(w_d,a_q); % coefficients of Eq. (32)
Vz = tf(Vz_num,Vz_den,Tu); % Eq. (32)

% H(z) formulation, Eq. (35), coefficients
H_num = conv(Pz_num,(Vz_num-Vz_den)); % order from z^-6 z^-5 ... z^-1 1
H_den = conv(Pz_den,Vz_den); % order from z^6 z^5 ... z^1 1

% H(z), Eq. (35)
Hz = tf(H_num,H_den,Tu);

% Getting the dimension of Eq. (33)
k_dim = size(H_num,2);

% state space realization of H(z), similar to FIR fomrulation
% A matrix
A_h = zeros(k_dim-1,k_dim-1);
A_h(1:end-1,2:end) = eye(k_dim-2);
A_h(end,:) = -flip(H_den(2:end));
% B matrix
B_h = zeros(k_dim-1,1);
B_h(end) = 1;
% C matrix
C_h = flip(H_num(2:end));
% D matrix
D_h = [H_num(1)];

% bandstop constraint values
Hz_phi = Hz*phi;
HK_d = freqresp(Hz_phi,exp(1j*w_d));

% state space for K(z), Eq. (33)
% A matrix
A_k = zeros(max_order,max_order);
A_k(1:(end-1),2:end) = eye(max_order-1,max_order-1);
% B matrix
B_k = zeros(max_order,1);
B_k(end) = 1;

% real bounded lemma set up, Eq. (37), \bar{A}  = [A_H, zero_mat; B_KC_H, A_K] 
zero_mat = zeros(size(A_h,1),size(A_k,2));
M_size = size(A_h)+size(A_k);
clear gamma_sdp M L_cvx quad_k
cvx_begin sdp
        % q -> Eq. (14), gamma, M -> Eq. (20)
        variables k((max_order+1),1) gamma_sdp
        variable M(M_size) symmetric

        % state space for K(z), Eq. (33)
        C_k = [flip(k(2:end))]';
        D_k = [k(1)];
        
        % state space for 1+H(z)K(z), \bar{A, B, C, D}, see Eq. (37)
        A_b = [A_h zero_mat; B_k*C_h A_k]; 
        B_b = [B_h; B_k*D_h];
        C_b = [D_k*C_h, C_k];
        D_b = D_k*D_h+1;

        % matrix from Eq. (2)
        L_cvx = [M, M*A_b, M*B_b, zeros(M_size(1),1);... 
            A_b'*M, M, zeros(M_size(1),1), C_b';...
            B_b'*M, zeros(1,M_size(1)), -gamma_sdp, D_b';...
            zeros(1,M_size(1)), C_b, D_b, -gamma_sdp];  

        % objective function
        minimize gamma_sdp 
        
        % constraint function
        subject to
            % bandstop constraint
            for i = 1:m_d
               1+HK_d(1,:,i)*k == 0;
            end
            
            % real-bounded lemma constraints, Eq. (20)
            L_cvx <= 0;
            gamma_sdp >= 0;
cvx_end

% QIIR filter, QIIR = Q0*(1-Fz)*Kcvx_IIR
Fz1 = (1-Vz);
Fz1_num = cell2mat(Fz1.numerator);
Fz1_den = cell2mat(Fz1.denominator);

% transfer function of (1-F(z))K(z)
Qcvx_num = conv(Fz1_num,k'); % conv((1-Fz),K)
Qcvx_den = conv(Fz1_den,[1 zeros(1,max_order)]);
PQcvx_IIR_num = conv(Pz_num,Qcvx_num);
PQcvx_IIR_den = conv(Pz_den,Qcvx_den);

% transfer function of P(z)Q_IIR(z)
PQcvx_IIR = tf(PQcvx_IIR_num,PQcvx_IIR_den,Tu);

% transfer function of 1-P(z)Q(z), Eq. (10) 
T_cvx_IIR = 1 - PQcvx_IIR;

%% bandpass filter baseline from previous work
w_hz_d = w_d/(2*pi*Tu);
B_bw = w_hz_d*0.2; % bandwidth
[Q_bp, Q0_bp] = q_bandpass(Pz,w_hz_d,B_bw,Tu);
PQ_bp = minreal(Pz*Q_bp); % transfer function of P(z)Q(z)
T_bp = 1-PQ_bp; % 1-P(z)Q(z), Eq. (10) 


%% =============== Plotting bode plots ==================================
% defining the frequency range to be plotted and axes limits
w_in_Hz_end = 1/(2*Tu); % end of plotting for Nyq freq of the fast
w_in_Hz = 1:1:w_in_Hz_end;
w_in_rad = w_in_Hz*2*pi;
l_width = 1.2; % linewidth
font_size = 11; % font size
x_lim_loop = [400 w_in_Hz_end]; % 1-PQ x_lim
y_lim = [-50 20]; % y_lim

% ==================== plotting style ==================================
color_cvx = {[0.9290 0.6940 0.1250], [0 1 0], [1 0 0], [0 0 1]};
line_style = {'-','--','--',':'};
n_all = 3;
Nyq_Hz = 1/(2*Ts);

% =================== Q-filter Bode ====================================
[mag_T_FIR, phi_T_FIR, ~] = bode(T_sdp_fir,w_in_rad);
[mag_PQ_FIR, phi_PQ_FIR, ~] = bode(PQ_sdp,w_in_rad);
mag_T_FIR = 20*log10(mag_T_FIR(:));
mag_PQ_FIR = 20*log10(mag_PQ_FIR(:));
phi_T_FIR = wrapTo180(phi_T_FIR(:));
phi_PQ_FIR = wrapTo180(phi_PQ_FIR(:));

[mag_T_IIR, phi_T_IIR, ~] = bode(T_cvx_IIR,w_in_rad);
[mag_PQ_IIR, phi_PQ_IIR, ~] = bode(PQcvx_IIR,w_in_rad);
mag_T_IIR = 20*log10(mag_T_IIR(:));
mag_PQ_IIR = 20*log10(mag_PQ_IIR(:));
phi_T_IIR = wrapTo180(phi_T_IIR(:));
phi_PQ_IIR = wrapTo180(phi_PQ_IIR(:));

[mag_T_bp, phi_T_bp, ~] = bode(T_bp,w_in_rad);
[mag_PQ_bp, phi_PQ_bp, ~] = bode(PQ_bp,w_in_rad);
mag_T_bp = 20*log10(mag_T_bp(:));
mag_PQ_bp = 20*log10(mag_PQ_bp(:));
phi_T_bp = wrapTo180(phi_T_bp(:));
phi_PQ_bp = wrapTo180(phi_PQ_bp(:));

% ======================== PQ plots =====================================
figure()
subplot(2,1,1)
h = semilogx(w_in_Hz,mag_PQ_bp, w_in_Hz,mag_PQ_FIR,w_in_Hz,mag_PQ_IIR);
set(h,'linewidth',l_width);
hold on
x_line = xline(Nyq_Hz);
x_line.Color = [0 0 0];
x_line.LineWidth = 1.5;
x_line.FontSize = 10;
x_line.FontWeight = 'bold';
for i = 1:m_d
    x_line = xline(f_d(i));
    x_line.Color = [0 0 0];
    x_line.LineWidth = 1;
    x_line.LineStyle = '--';
end
xlim(x_lim_loop);
ylim(y_lim);
xlabel('Hz')
ylabel('dB')
title('Magnitude of PQ')
ax = gca;
ax.FontSize= font_size;
for i = 1:n_all
    h(i).Color = color_cvx{i};
    h(i).LineStyle = line_style{i};
end
hold off

% ============================ 1-PQ plots ===========================
subplot(2,1,2)
h = semilogx(w_in_Hz,mag_T_bp, w_in_Hz,mag_T_FIR,w_in_Hz,mag_T_IIR);
set(h,'linewidth',l_width);
hold on
x_line = xline(Nyq_Hz);
x_line.Label = sprintf('%.f Hz',Nyq_Hz);
x_line.Color = [0 0 0];
x_line.LineWidth = 1.5;
x_line.FontWeight = 'bold';
x_line.LabelVerticalAlignment = 'bottom';
x_line.LabelHorizontalAlignment = 'left';
for i = 1:size(w_d,2)
    x_line = xline(f_d(i));
    x_line.Color = [0 0 0];
    x_line.LineWidth = 1;
    x_line.Label = sprintf('%.f Hz',f_d(i));
    x_line.LabelHorizontalAlignment = 'left';
    x_line.LabelOrientation = 'aligned';
    x_line.LabelVerticalAlignment = 'Bottom';
    x_line.FontSize = 10;
    x_line.LineStyle = '--';
    x_line.FontWeight = 'bold';
end
xlim(x_lim_loop);
ylim(y_lim);
xlabel('Hz')
ylabel('dB')
title('Magnitude 1-PQ')
ax = gca;
ax.FontSize= font_size;
for i = 1:n_all
    h(i).Color = color_cvx{i};
    h(i).LineStyle = line_style{i};
end
legend('Bandpass','SDP:FIR','SDP:IIR','location','southwest')