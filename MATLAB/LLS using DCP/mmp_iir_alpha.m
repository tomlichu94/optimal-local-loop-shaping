clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% notes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for generalized number of disturbances and selected L_t
% and L_t sampling multiplier (T_{ss} = L_t*T_{fs})
% test showing the signal reconstruction built in SIMULINK with CVX using IIR
%%%%%%%%%%%%%%%%%%%%%%%%%% load plant data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
max_order = 30; % max filter order
max_amp = 1; % max disturb amplitude

PQ_max = 3; % max value of PQ for the quadratic constraint
beta = PQ_max^2; % FIR SDP, set max value for quad, play around with quadratic
f_stop = 600; % SOCP stop constraint past this frequency
a_g_IIR = 0.95; % predictor alpha
alpha = [0.3 0.5 0.9 0.99];

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
tempW = [1.32 1.67 2.39 3.41];

m_d = size(tempW,2); % number of disturbances
w_d = tempW*pi/L_t; % fast measurement of disturbance in radians
f_d = w_d/(2*pi*Tu); % disturbance in Hz

%% ============== Predictor Coefficients and TF, W_k ======================
[w_k] = W_coeff_FIR(L_t,f_d,Tu); % predictor coefficients
[w_k_IIR, B_para] = W_coeff_IIR(L_t,f_d,a_g_IIR,Tu);
[W_FIR_num, W_FIR_den] = W_TF_FIR(w_k);
[W_IIR_num, W_IIR_den] = W_TF_IIR(w_k_IIR,B_para);
for k = 1:k_t
    W_k_FIR(k) = tf(W_FIR_num(k,:),W_FIR_den(k,:),Tu);
    W_k_IIR(k) = tf(W_IIR_num(k,:),W_IIR_den(k,:),Tu);
end

%% ============== discrete set of frequencies ==========================
range_mult = 30;
w_lin = w_lin_spacing(max_order,range_mult,tempW,L_t);
w_lin_Hz = w_lin/(2*pi*Tu);
stop_indx = find(w_lin_Hz>f_stop);
stop_indx = stop_indx(1);
fprintf('Stop Constraint %u Hz \n',w_lin_Hz(stop_indx))

%% %%%%%%%%%%%%%%%%%%%%%%%%% hdd model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ps_pade = pade(Ps,4); % pade approximation for delay term
Ps_sys = ss(Ps_pade);
Pz = absorbDelay(Pz_delay);
Pz_num = cell2mat(Pz.num);
Pz_den = cell2mat(Pz.den);
%% %%%%%%%%%%%%%%%%%%%%%%%%% Q Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generalized Q(z) form is: Q(z) = Q0*Q_FIR or Q0*Q_IIR
%%%%% using Q0 = 1-z^-1
Q0_num = [1 -1];
Q0_den = [1 0]; 
Q0 = tf(Q0_num,Q0_den,Tu); % Q0 = 1-z^-1 = (z-1)/z
% Q0_num = 1;
% Q0_den = 1;
%%%%% using Q0 = 1+z^-1
% Q0_num = [1 -1];
% Q1_num = [1 1];
% Q0_den = [1 0]; 
% Q2_num = conv(Q0_num,Q1_num);
% Q2_den = conv(Q0_den,Q0_den);
% Q0 = tf(Q2_num,Q2_den,Tu);
%%%%% using Q0 from Automatica
% rho = 0.8;
% Q0_num = q0(w_d,rho);
% Q0_den = zeros(size(Q0_num));
% Q0_den(1) = 1;
% Q0 = tf(Q0_num,Q0_den,Tu);

%% %%%%%%%%%%%%%%%%%%%%%%%%% PQ Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%% Stabilizing controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_c = 0.1;
kp= 1/13320;
ki= 1/33300;
kd= 1/2775;
Cz = k_c*(kp + ki/(z-1) + kd*(z-1)/z); % PID controller
% check stability
Tz = feedback(Pz*Cz,1);
if isstable(Tz) == 1
    fprintf('CL system is stable\n')
else
    fprintf('Unstable CL, program stopped\n')
    return
end

%% Setting up variables for SDP
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
%%%%%%%%%%%%%%%%%%%%%%%%%% IIR implementation %%%%%%%%%%%%%%%%%%%%%%%%%%
%% for loop starts here
%% QFIR SDP implementation
M_size = size(A_p)+size(A_q);
zero_mat = zeros(size(A_p,1),size(A_q,2));
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

%% QIIR SDP implementation
for u = 1:length(alpha)
[Fz_num, Fz_den] = F_notch(w_d,alpha(u));
Fz = tf(Fz_num,Fz_den,Tu);
H_num = conv(PQ0_num,(Fz_num-Fz_den)); % order from z^-6 z^-5 ... z^-1 1
H_den = conv(PQ0_den,Fz_den); % order from z^6 z^5 ... z^1 1
Hz = tf(H_num,H_den,Tu);
k_iir = size(H_num,2);
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
Fz1 = (1-Fz);
Fz1_num = cell2mat(Fz1.numerator);
Fz1_den = cell2mat(Fz1.denominator);
Q0Fz1_num = conv(Q0_num,Fz1_num); % conv(Q0,1-Fz)
Q0Fz1_den = conv(Q0_den,Fz1_den);
Qcvx_num = conv(Q0Fz1_num,k'); % conv(Q0(1-Fz),K)
Qcvx_den = conv(Q0Fz1_den,[1 zeros(1,max_order)]);
PQcvx_IIR_num = conv(Pz_num,Qcvx_num);
PQcvx_IIR_den = conv(Pz_den,Qcvx_den);
PQcvx_IIR(u) = tf(PQcvx_IIR_num,PQcvx_IIR_den,Tu);
T_cvx_IIR(u) = minreal(1 - PQcvx_IIR(u));
opt_value(u) = cvx_optval;
end

%% plotting
close all
color_cvx = {[0.9290 0.6940 0.1250], [0 1 0], [1 0 0], [0 0 1]};
color_base = {[0 0 0],[0.4940 0.1840 0.5560],[0.3010 0.7450 0.9330]};
marker_style = {'none','diamond','o','x','square'};
line_style = {'-','--','--',':'};
w_in_Hz_end = 1/(2*Tu); % end of plotting for Nyq freq of the fast
w_in_Hz = 1:1:w_in_Hz_end;
w_in_rad = w_in_Hz*2*pi;
Nyq_Hz = 1/(2*Ts);
x_lim_loop = [1 2400]; % 1-PQ x_lim
font_size = 11; % font size
y_lim = [-50 20]; % y_lim
l_width = 1.3;

figure()
subplot(2,1,1)
for i = 1:length(alpha)
    [mag_PQ_IIR, phi_PQ_IIR, ~] = bode(PQcvx_IIR(i),w_in_rad);
    mag_PQ_IIR = 20*log10(mag_PQ_IIR(:));
    h(i) = semilogx(w_in_Hz,mag_PQ_IIR);
    hold on
end
x_line = xline(Nyq_Hz);
x_line.Color = [0 0 0];
x_line.LineWidth = 1.5;
x_line.FontSize = 11;
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
for i = 1:length(alpha)
    h(i).Color = color_cvx{i};
    h(i).LineStyle = line_style{i};
    h(i).LineWidth = 1.3;
end
hold off
subplot(2,1,2)
for i = 1:length(alpha)
    [mag_T_IIR, phi_T_IIR, ~] = bode(T_cvx_IIR(i),w_in_rad);
    mag_T_IIR = 20*log10(mag_T_IIR(:));
    h(i) = semilogx(w_in_Hz,mag_T_IIR);
    hold on
end
x_line = xline(Nyq_Hz);
x_line.Color = [0 0 0];
x_line.LineWidth = 1.5;
x_line.FontSize = 11;
x_line.FontWeight = 'bold';
for i = 1:size(w_d,2)
    x_line = xline(f_d(i));
    x_line.Color = [0 0 0];
    x_line.LineWidth = 1;
    x_line.Label = sprintf('%.f Hz',f_d(i));
    x_line.LabelHorizontalAlignment = 'left';
    x_line.LabelOrientation = 'aligned';
    x_line.LabelVerticalAlignment = 'Bottom';
    x_line.FontSize = 11;
    x_line.LineStyle = '--';
    x_line.FontWeight = 'bold';
end
xlim(x_lim_loop);
ylim(y_lim);
xlabel('Hz')
ylabel('dB')
title('Magnitude of 1-PQ')
ax = gca;
ax.FontSize= font_size;
for i = 1:length(alpha)
    h(i).Color = color_cvx{i};
    h(i).LineStyle = line_style{i};
    h(i).LineWidth = l_width;
end
Legend_alpha = cell(1,4);
for iter=1:4
  Legend_alpha{iter}=strcat('\alpha = ', num2str(alpha(iter)));
end
legend(Legend_alpha,'Location','Best')