%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% notes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is a modal analysis on the T(z) on looking at the Bode's integral
% theorem. We look at the block diagram from Hui's multirate paper

% if time allows, comapre to the 2024 ACC paper
%%%%%%%%%%%%%%%%%%%%%%%%%% load plant data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear
% close all
% clc

addpath('Functions','Simulink','Data')
load mainPlantData;

%%%%%%%%%%%%%%%%%%%%% Loading Simulink parameters %%%%%%%%%%%%%%%%%%%%%%%%%
% the plant is originally a high-precision fast-sampling system. To minic
% the case of slow sensor speed, the plant sampling frequency is manually reduced
% (i.e. increase the plant/sensor sampling time).
Nx = 10; % multiplier to increase sampling time
PlantData.Ts = PlantData.Ts*Nx;
PlantData.Tu = PlantData.Tu*Nx;
Tu = PlantData.Tu; % baseline plant sampling time
Ps = PlantData.Pn; % continuous-time plant
Ps = tf(Ps); % define CT plnat as a transfer function
Pz_delay = c2d(Ps,Tu,'zoh'); % discretize the plant
z = tf('z',Tu); % discrete time based on fast sampling
s = tf('s'); % continous time

%%%%%%%%%%%%%%%%%%%% input from the user %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_t = 4; % sampling rate multilpier
k_t = L_t - 1;
Ts = Tu*L_t; % slow sampling time
batches = 220; % the number of cycles 
max_order = 30; % max filter order
max_amp = 1; % max disturb amplitude

PQ_max = 3; % max value of PQ for the quadratic constraint
beta = PQ_max^2; % FIR SDP, set max value for quad, play around with quadratic
f_stop = 600; % SOCP stop constraint past this frequency
a_g_IIR = 0.9; % predictor alpha
alpha = 0.9; % QIIR alpha

%%%%%%%%%%%%%%%%%%%% narrow-band disturbance frequency %%%%%%%%%%%%%%%%%%%%
% generates random multiplier for beyond Nyq freq of slow samlper
% will generate and pick multiplier based on the number of disturbances
% spaced roughly equidistance. e.g. m_d = 3, L_t = 2. Will pick multiplier
% between 1 and 2. 1st between 1-1.33, 2nd between 1.333-1.67, 3rd between
% 1.67-2
% m_d = 3;
% tempW = zeros(1,m_d);
% for i = 1:(m_d) 
%     tempW(i) = 1+rand(1)*((L_t-1)/m_d)+(i-1)*(L_t-1)/m_d; 
% end

% general test results
% tempW = [1.337 1.739 2.312 2.618]; % okay, bad robustness
% tempW = [1.32 1.67 1.93 2.18]; % good run
% tempW = [1.4482 1.7411 2.0070 2.8114]; % okay, mid robustness
% tempW = [0.23 0.56 0.8];
% tempW = [1.7178 2.0072 2.9787 3.5073];
tempW = [1.32 1.67 2.39 3.41];

m_d = size(tempW,2); % number of disturbances
w_d = tempW*pi/L_t; % fast measurement of disturbance in radians
f_d = w_d/(2*pi*Tu); % disturbance in Hz
p_off = rand(1,m_d)*pi; % random phase shift from 0 - 1
A_amp = rand(1,m_d)*max_amp; % random amplitude from 0 - max_ampl

% ============== Predictor Coefficients and TF, W_k ======================
[w_k_fir] = w_kfir_frac(f_d, Tu, 1, L_t); % predictor coefficients
[w_k_iir, B_para] = w_kiir_frac(f_d, Tu, a_g_IIR, 1, L_t);
[W_FIR_num, W_FIR_den] = W_TF_FIR(w_k_fir);
[W_IIR_num, W_IIR_den] = W_TF_IIR(w_k_iir,B_para);
W_k_FIR(1,k_t) = tf(0);
W_k_iir(1,k_t) = tf(0);

for k = 1:k_t
    W_k_FIR(k) = tf(W_FIR_num(k,:),W_FIR_den(k,:),Tu);
    W_k_iir(k) = tf(W_IIR_num(k,:),W_IIR_den(k,:),Tu);
end
%%
% ============== discrete set of frequencies ==========================
range_mult = 30;
w_lin = w_lin_spacing(max_order,range_mult,tempW,L_t);
w_lin_Hz = w_lin/(2*pi*Tu);
stop_indx = find(w_lin_Hz>f_stop);
stop_indx = stop_indx(1);
fprintf('Stop Constraint %u Hz \n',w_lin_Hz(stop_indx))

% %%%%%%%%%%%%%%%%%%%%%%%%% hdd model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ps_pade = pade(Ps,4); % pade approximation for delay term
Ps_sys = ss(Ps_pade);
Pz = absorbDelay(Pz_delay);
k_c = 0.1;
kp= 1/13320;
ki= 1/33300;
kd= 1/2775;
Cz = k_c*(kp + ki/(z-1) + kd*(z-1)/z); % PID controller
% check stability
Pz = feedback(Pz*Cz,1);

Pz_num = cell2mat(Pz.num);
Pz_den = cell2mat(Pz.den);
Pz_sys = ss(Pz); % used in Simulink

% %%%%%%%%%%%%%%%%%%%%%%%%% Q Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generalized Q(z) form is: Q(z) = Q0*Q_FIR or Q0*Q_IIR
Q0_num = 1;
Q0_den = 1;

% %%%%%%%%%%%%%%%%%%%%%%%%% PQ Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = tf(); % create phi = [1 z^-1 ... z^-r]
phi(1) = 1;
for i = 1:max_order
    phi(i+1) = z^(-i);
end
PQ0_num = conv(Pz_num,Q0_num);
PQ0_den = conv(Pz_den,Q0_den);
PQ0 = tf(PQ0_num,PQ0_den,Tu); % transfer function of Pz*Q0
P_phi = PQ0*phi;
PQ_d = freqresp(P_phi,exp(1j*w_d)); % sub in disturbance vals into P_phi

% SOCP implementation
%%%%%%%%%%%%%%%%%%%% Quad implementation for PQ %%%%%%%%%%%%%%%%%%%%%%%%%%%
phi_val = freqresp(phi,exp(1j*w_lin(1:length(w_lin))));
phi_val = squeeze(phi_val);
phi_r = real(phi_val);
phi_i = -imag(phi_val);
Pz_val = freqresp(PQ0,exp(1j*w_lin(1:length(w_lin))));
Pz_val = squeeze(Pz_val)';
Pz_2 = real(Pz_val).^2 + imag(Pz_val).^2;
quad = zeros(max_order+1,max_order+1,length(w_lin));
for i = 1:length(w_lin)
    quad(:,:,i) = Pz_2(i)*phi_r(:,i)*phi_r(:,i)'+...
                  Pz_2(i)*phi_i(:,i)*phi_i(:,i)';
end
cvx_begin
        variables q_vec((max_order+1),1) b(length(w_lin),1)
        beta_sum = sum(b);
        for i = 1:length(w_lin)
           quad_q(i) = quad_form(q_vec,quad(:,:,i));
        end
        minimize beta_sum
        subject to
            for i = 1:m_d
               1-PQ_d(1,:,i)*q_vec == 0;
            end            
            for i = 1:length(w_lin)
                quad_q(i) <= b(i);
            end
cvx_end
Qcvx_num = conv(Q0_num,q_vec');
Qcvx_den = conv(Q0_den,[1 zeros(1,max_order)]);
Qcvx_socp = tf(Qcvx_num,Qcvx_den,Tu);
Qcvx_socp = minreal(Qcvx_socp);
PQ_socp = minreal(Pz*Qcvx_socp*W_k_FIR(1));
T_socp = minreal(1-PQ_socp*W_k_FIR(1));

% QFIR SDP implementation
%%%%%%%%%%%%%%%%%%%%%%%%% FIR implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
% state space realization of P(z)
k_fir = size(PQ0_den,2);
A_p = zeros(k_fir-1,k_fir-1);
A_p(1:end-1,2:end) = eye(k_fir-2);
A_p(end,:) = -flip(PQ0_den(2:end));
B_p = zeros(k_fir-1,1);
B_p(end) = 1;
C_p = flip(PQ0_num(2:end));
D_p = PQ0_num(1);

% bounded real lemma setup
A_q = zeros(max_order,max_order);
A_q(1:(end-1),2:end) = eye(max_order-1,max_order-1);
B_q = zeros(max_order,1);
B_q(end) = 1;
zero_mat = zeros(size(A_p,1),size(A_q,2));
M_size = size(A_p)+size(A_q);
clear quad_q q
cvx_begin sdp
        variables q((max_order+1),1) rho
        variable M(M_size) symmetric
        C_q = [flip(q(2:end))]';
        D_q = [q(1)];
        A_t = [A_p zero_mat; B_q*C_p A_q];
        B_t = [B_p; B_q*D_p];
        %%%% 1-PQ formulation
        C_t = -[D_q*C_p, C_q];
        D_t = 1-D_q*D_p;
        L_cvx = [M, M*A_t, M*B_t, zeros(M_size(1),1);...
             A_t'*M, M, zeros(M_size(1),1), C_t';...
             B_t'*M, zeros(1,M_size(1)), -rho, D_t';...
             zeros(1,M_size(1)), C_t, D_t, -rho]; 
        for i = 1:stop_indx
             quad_q(i) = quad_form(q,quad(:,:,i));
        end
        minimize rho
        subject to
            for i = 1:m_d
               1-PQ_d(1,:,i)*q == 0;
            end
            for i = 1:stop_indx
                quad_q(i) <= beta
            end
            L_cvx <= 0;
            rho >= 0;
cvx_end

% Q filter construction
Qcvx_num = conv(Q0_num,q');
Qcvx_den = conv(Q0_den,[1 zeros(1,max_order)]);
Qcvx_FIR = tf(Qcvx_num,Qcvx_den,Tu);
PQcvx_FIR_num = conv(Pz_num,Qcvx_num);
PQcvx_FIR_den = conv(Pz_den,Qcvx_den);
PQcvx_FIR = tf(PQcvx_FIR_num,PQcvx_FIR_den,Tu)*W_k_FIR(1);
T_cvx_FIR = 1 - PQcvx_FIR*W_k_FIR(1);

% QIIR SDP implementation
%%%%%%%%%%%%%%%%%%%%%%%%%% IIR implementation %%%%%%%%%%%%%%%%%%%%%%%%%%
[Fz_num, Fz_den] = F_notch(w_d,alpha);
Fz = tf(Fz_num,Fz_den,Tu);

%%%%%%%%%%%%%%%%%%%%%% Notes on H(z) formulation %%%%%%%%%%%%%%%%%%%%%%%%%%
H_num = conv(PQ0_num,(Fz_num-Fz_den)); % order from z^-6 z^-5 ... z^-1 1
H_den = conv(PQ0_den,Fz_den); % order from z^6 z^5 ... z^1 1
Hz = tf(H_num,H_den,Tu);
k_iir = size(H_num,2);

% quadratic constraint
Hz_val = freqresp(Hz,exp(1j*w_lin(1:stop_indx)));
Hz_2 = real(Hz_val).^2 + imag(Hz_val).^2;
quad_IIR = zeros(max_order+1,max_order+1,stop_indx);
for i = 1:stop_indx
    quad_IIR(:,:,i) = Hz_2(i)*phi_r(:,i)*phi_r(:,i)'+...
                      Hz_2(i)*phi_i(:,i)*phi_i(:,i)';
end

% state space realization of H(z)
A_h = zeros(k_iir-1,k_iir-1);
A_h(1:end-1,2:end) = eye(k_iir-2);
A_h(end,:) = -flip(H_den(2:end));
C_h = flip(H_num(2:end));
B_h = zeros(k_iir-1,1);
B_h(end) = 1;
D_h = [H_num(1)];

Hz_phi = Hz*phi;
HK_d = freqresp(Hz_phi,exp(1j*w_d));
% bounded real lemma setup
A_k = zeros(max_order,max_order);
A_k(1:(end-1),2:end) = eye(max_order-1,max_order-1);
B_k = zeros(max_order,1);
B_k(end) = 1;
zero_mat = zeros(size(A_h,1),size(A_k,2));
M_size = size(A_h)+size(A_k);
clear rho M L_cvx quad_k
cvx_begin sdp
        variables k((max_order+1),1) rho
        variable M(M_size) symmetric
        C_k = [flip(k(2:end))]';
        D_k = [k(1)];
        A_t = [A_h zero_mat; B_k*C_h A_k];
        B_t = [B_h; B_k*D_h];
        %%%%% 1+HK %%%%%%%%%%%
        C_t = [D_k*C_h, C_k];
        D_t = D_k*D_h+1;
        %%%%% HK formulation, where HK = -PQ
        L_cvx = [M, M*A_t, M*B_t, zeros(M_size(1),1);... 
            A_t'*M, M, zeros(M_size(1),1), C_t';...
            B_t'*M, zeros(1,M_size(1)), -rho, D_t';...
            zeros(1,M_size(1)), C_t, D_t, -rho];  
        for i = 1:stop_indx
           quad_k(i) = quad_form(k,quad_IIR(:,:,i));
        end
        minimize rho 
        subject to
            for i = 1:m_d
               1+HK_d(1,:,i)*k == 0
            end
            for i = 1:stop_indx
                quad_k(i) <= beta;
            end
            L_cvx<=0;
            rho >= 0;
cvx_end

% QIIR filter, QIIR = Q0*(1-Fz)*Kcvx_IIR
Fz1 = (1-Fz);
Fz1_num = cell2mat(Fz1.numerator);
Fz1_den = cell2mat(Fz1.denominator);
Q0Fz1_num = conv(Q0_num,Fz1_num); % conv(Q0,1-Fz)
Q0Fz1_den = conv(Q0_den,Fz1_den);
Qcvx_num = conv(Q0Fz1_num,k'); % conv(Q0(1-Fz),K)
Qcvx_den = conv(Q0Fz1_den,[1 zeros(1,max_order)]);
Qcvx_IIR = tf(Qcvx_num,Qcvx_den,Tu); % for Simulink
PQcvx_IIR_num = conv(Pz_num,Qcvx_num);
PQcvx_IIR_den = conv(Pz_den,Qcvx_den);
PQcvx_IIR = tf(PQcvx_IIR_num,PQcvx_IIR_den,Tu)*W_k_FIR(1);
T_cvx_IIR = 1 - PQcvx_IIR*W_k_FIR(1);

% bandpass baseline
w_hz_d = w_d/(2*pi*Tu);
B_bw = w_hz_d*0.2; % bandwidth
[Q_bp, Q0_bp] = q_bandpass(Pz,w_hz_d,B_bw,Tu);
PQ_bp = minreal(Pz*Q_bp);
T_bp = 1-PQ_bp;

% Storing Qcvx terms
Qcvx(1) = Qcvx_socp;
Qcvx(2) = Qcvx_FIR;
Qcvx(3) = Qcvx_IIR;

%% =============== Plotting ==================================
% defining the frequency range to be plotted and axes limits

w_in_Hz_end = 1/(2*Tu); % end of plotting for Nyq freq of the fast
w_in_Hz = 1:1:w_in_Hz_end;
w_in_rad = w_in_Hz*2*pi;
Nyq_Hz = 1/(2*Ts);
w_start = 1; % starting frequency for x_lim
l_width = 1.2; % linewidth
font_size = 11; % font size
x_lim_loop = [400 w_in_Hz_end]; % 1-PQ x_lim
y_lim = [-50 20]; % y_lim

% ==================== plotting style ==================================
color_cvx = {[0.9290 0.6940 0.1250], [0 1 0], [1 0 0], [0 0 1]};
color_base = {[0 0 0],[0.4940 0.1840 0.5560],[0.3010 0.7450 0.9330]};
marker_style = {'none','diamond','o','x','square'};
line_style = {'-','--','--',':'};
n_all = 4;

% =================== W Predictor Bode =================================
[mag_W_IIR1, phi_W_IIR1, ~] = bode(W_k_IIR(1),w_in_rad);
[mag_W_IIR2, phi_W_IIR2, ~] = bode(W_k_IIR(2),w_in_rad);
mag_W_IIR1 = 20*log10(mag_W_IIR1(:));
phi_W_IIR1 = wrapTo180(phi_W_IIR1(:));
mag_W_IIR2 = 20*log10(mag_W_IIR2(:));
phi_W_IIR2 = wrapTo180(phi_W_IIR2(:));

figure()
subplot(2,1,1)
semilogx(w_in_Hz,mag_W_IIR1,w_in_Hz,mag_W_IIR2);
xlabel('Hz')
ylabel('Magnitude (dB)')
xlim([1 x_lim_loop(end)])
title('W Predictor Bode')
subplot(2,1,2)
semilogx(w_in_Hz,phi_W_IIR1,w_in_Hz,phi_W_IIR2);
hold off
xlabel('Hz')
ylabel('Phase (deg.)')
xlim([1 x_lim_loop(end)])
legend('IIR 1','IIR 2')

% =================== Q-filter Bode ====================================
[mag_T_FIR, phi_T_FIR, ~] = bode(T_cvx_FIR,w_in_rad);
[mag_PQ_FIR, phi_PQ_FIR, ~] = bode(PQcvx_FIR,w_in_rad);
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

[mag_T_socp, phi_T_socp, ~] = bode(T_socp,w_in_rad);
[mag_PQ_socp, phi_PQ_socp, ~] = bode(PQ_socp,w_in_rad);
mag_T_socp = 20*log10(mag_T_socp(:));
mag_PQ_socp = 20*log10(mag_PQ_socp(:));
phi_T_socp = wrapTo180(phi_T_socp(:));
phi_PQ_socp = wrapTo180(phi_PQ_socp(:));
% ======================== PQ plots =====================================
figure()
subplot(2,1,1)
h = semilogx(w_in_Hz,mag_PQ_bp,w_in_Hz,mag_PQ_socp, ...
              w_in_Hz,mag_PQ_FIR,w_in_Hz,mag_PQ_IIR);
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
h = semilogx(w_in_Hz,mag_T_bp,w_in_Hz,mag_T_socp, ...
              w_in_Hz,mag_T_FIR,w_in_Hz,mag_T_IIR);
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
legend('Bandpass','SOCP:FIR', 'SDP:FIR','SDP:IIR','location','southwest')

%% robustness analysis
% derived T = PQW
% T_fir(i) = f(Qfir,Wfir_{mmp,i})
% T_iir(i) = f(Qiir,Wfir_{mmp,i})
close all

PQ_fir_num = conv(Pz_num, cell2mat(Qcvx(2).num));
PQ_fir_den = conv(Pz_den, cell2mat(Qcvx(2).den));
PQ_iir_num = conv(Pz_num, cell2mat(Qcvx(3).num));
PQ_iir_den = conv(Pz_den, cell2mat(Qcvx(3).den));

T_fir(1,3) = tf(0);
T_iir(1,3) = tf(0);
% % using PQW
% for i = 1:3
%     PQW_fir_num = conv(PQ_fir_num, cell2mat(W_k_IIR(i).num));
%     PQW_fir_den = conv(PQ_fir_den, cell2mat(W_k_IIR(i).den));
%     PQW_iir_num = conv(PQ_iir_num, cell2mat(W_k_IIR(i).num));
%     PQW_iir_den = conv(PQ_iir_den, cell2mat(W_k_IIR(i).den));
%     T_fir(i) = tf(PQW_fir_num, PQW_fir_den, Tu);
%     T_iir(i) = tf(PQW_iir_num, PQW_iir_den, Tu);
% end
% using PQ
for i = 1:3
    T_fir(i) = tf(PQ_fir_num, PQ_fir_den, Tu);
    T_iir(i) = tf(PQ_iir_num, PQ_iir_den, Tu);
end

for i = 1:3
    figure
    pzmap(T_fir(i))
    legend(strcat('Q FIR, k=', num2str(i)));
    title('Q-FIR')
    figure
    pzmap(T_iir(i))
    legend(strcat('Q IIR, k=', num2str(i)));
    title('Q-IIR')
end
T_int_fir = zeros(2,3);
T_int_iir = zeros(2,3);
% finding NMP zeros
for i = 1:3
    T_zeros = abs(zero(T_fir(i)));
    T_zeros = log(T_zeros(T_zeros>1));
    T_int_fir(1,i) = 2*pi*sum(T_zeros);
    T_zeros = abs(zero(T_iir(i)));
    T_zeros = log(T_zeros(T_zeros>1));
    T_int_iir(1,i) = 2*pi*sum(T_zeros);
end
Tinv(1:6) = tf(0);
% bode of inv(T)
for i = 1:3
    Tinv(i) = inv(T_fir(i));
    Tinv(i+3) = inv(T_iir(i));
end

% bode of T
[mag1, phi, ~] = bode(T_fir,w_in_rad);
mag1 = squeeze(mag1);
for i = 1:3
    T_int_fir(2,i) = trapz(w_in_rad,mag1(i,:));
end
T_int_fir

[mag2, phi, ~] = bode(T_iir,w_in_rad);
mag2 = squeeze(mag2);
for i = 1:3
    T_int_iir(2,i) = trapz(w_in_rad,mag2(i,:));
end
T_int_iir

% plotting 1/T
Tinv = minreal(Tinv);
[mag, phi, ~] = bode(Tinv,w_in_rad);
mag = 20*log10(squeeze(mag));
phi = wrapTo180(squeeze(phi));

legend_bode = cell(1,3);
for i = 1:3
    legend_bode{i}=strcat('k =', num2str(i));
end
figure
hold on
for i = 1:3
    h(i) = semilogx(w_in_Hz,mag(i,:));
    h(i).Color = color_cvx{i};
    h(i).LineStyle = line_style{i};
    h(i).LineWidth = 1.5;
end
hold off
ax = gca;
ax.FontSize= font_size;
xlim([1 2640])
legend(legend_bode)
xlabel('Hz')
ylabel('Magnitude (dB)')
title('W-FIR')

figure
hold on
for i = 1:3
    h(i) = semilogx(w_in_Hz,mag(i+3,:));
    h(i).Color = color_cvx{i};
    h(i).LineStyle = line_style{i};
    h(i).LineWidth = 1.5;
end
hold off
ax = gca;
ax.FontSize= font_size;
legend(legend_bode)
xlim([1 2640])
xlabel('Hz')
ylabel('Magnitude (dB)')
title('W-IIR')

% troubleshooting
S_fir = minreal(1-T_fir);
S_iir = minreal(1-T_iir);

figure
    [mag, phi, ~] = bode(S_fir(1),w_in_rad);
    mag = 20*log10(squeeze(mag));
    phi = wrapTo180(squeeze(phi));
    h = semilogx(w_in_Hz,mag);
    h.Color = color_cvx{1};
    h.LineStyle = line_style{1};
    h.LineWidth = 1.5;
    hold on
    [mag, phi, ~] = bode(T_fir(1),w_in_rad);
    mag = 20*log10(squeeze(mag));
    phi = wrapTo180(squeeze(phi));
    h = semilogx(w_in_Hz,mag);
    h.Color = color_cvx{2};
    h.LineStyle = line_style{2};
    h.LineWidth = 1.5;
    hold off
    ax = gca;
    ax.FontSize= font_size;
    xlim([1 2640])
    legend('S','T')
    xlabel('Hz')
    ylabel('Magnitude (dB)')
    title('S vs T for Q-FIR')
    x_line = xline(Nyq_Hz);
    x_line.Color = [0 0 0];
    x_line.LineWidth = 1;
    x_line.FontSize = 10;
    x_line.FontWeight = 'bold';
    for i = 1:m_d
        x_line = xline(f_d(i));
        x_line.Color = [0 0 0];
        x_line.LineWidth = 1;
        x_line.LineStyle = '--';
    end

figure
    [mag, phi, ~] = bode(S_iir(1),w_in_rad);
    mag = 20*log10(squeeze(mag));
    phi = wrapTo180(squeeze(phi));
    h = semilogx(w_in_Hz,mag);
    h.Color = color_cvx{1};
    h.LineStyle = line_style{1};
    h.LineWidth = 1.5;
    hold on
    [mag, phi, ~] = bode(T_iir(1),w_in_rad);
    mag = 20*log10(squeeze(mag));
    phi = wrapTo180(squeeze(phi));
    h = semilogx(w_in_Hz,mag);
    h.Color = color_cvx{2};
    h.LineStyle = line_style{2};
    h.LineWidth = 1.5;
    hold off
    ax = gca;
    ax.FontSize= font_size;
    xlim([1 2640])
    legend('S','T')
    xlabel('Hz')
    ylabel('Magnitude (dB)')
    title('S vs T for Q-IIR')
    x_line = xline(Nyq_Hz);
    x_line.Color = [0 0 0];
    x_line.LineWidth = 1;
    x_line.FontSize = 10;
    x_line.FontWeight = 'bold';
    for i = 1:m_d
        x_line = xline(f_d(i));
        x_line.Color = [0 0 0];
        x_line.LineWidth = 1;
        x_line.LineStyle = '--';
    end