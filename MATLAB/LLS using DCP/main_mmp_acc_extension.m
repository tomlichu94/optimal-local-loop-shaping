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
batches = 220; % the number of cycles 
max_order = 30; % max filter order
max_amp = 1; % max disturb amplitude

PQ_max = 3; % max value of PQ for the quadratic constraint
beta = PQ_max^2; % FIR SDP, set max value for quad, play around with quadratic
f_stop = 600; % SOCP stop constraint past this frequency
a_g_IIR = 0.95; % predictor alpha
alpha = 0.95; % QIIR alpha

%%%%%%%%%%%%%%%%%%%%% simulation run time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tsim = batches*150*Ts; %   simulation time
simStepSize = Ts/80; % simulation step size
Ts_CT_approx = Ts/20; % approximating continuous time sys
[SensorNoise, ForceDist] = SetDisturbance(Tsim,Ts,Tu); % generate noise

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
Pz_sys = ss(Pz); % used in Simulink
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

%% SOCP implementation
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
Q_socp = tf(Qcvx_num,Qcvx_den,Tu);
Q_socp = minreal(Q_socp);
PQ_socp = minreal(Pz*Q_socp);
T_socp = minreal(1-PQ_socp);

%% QFIR SDP implementation
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
PQcvx_FIR = tf(PQcvx_FIR_num,PQcvx_FIR_den,Tu);
T_cvx_FIR = 1 - PQcvx_FIR;

%% QIIR SDP implementation
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
PQcvx_IIR = tf(PQcvx_IIR_num,PQcvx_IIR_den,Tu);
T_cvx_IIR = 1 - PQcvx_IIR;

%% bandpass baseline
w_hz_d = w_d/(2*pi*Tu);
B_bw = w_hz_d*0.2; % bandwidth
[Q_bp, Q0_bp] = q_bandpass(Pz,w_hz_d,B_bw,Tu);
PQ_bp = minreal(Pz*Q_bp);
T_bp = 1-PQ_bp;

SimModelBase = 'main_mmp_baseline';
[SimOutBase] = sim(SimModelBase,'StopTime','Tsim','FixedStep','simStepSize');
d_out = SimOutBase.dist_fast.signals.values;
y_base_fast = SimOutBase.y_fast.signals.values;
y_base_slow = SimOutBase.y_slow.signals.values;
t_base_fast = SimOutBase.y_fast.time;
t_base_slow = SimOutBase.y_slow.time;

%% Storing Qcvx terms
Qcvx(1) = Q_socp;
Qcvx(2) = Qcvx_FIR;
Qcvx(3) = Qcvx_IIR;

%% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run Simulink %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SimModel = 'main_mmp';
for i = 1:3
    Qcvx_sim = Qcvx(i);
    [SimOut] = sim(SimModel,'StopTime','Tsim','FixedStep','simStepSize');
    d_out(:,i) = SimOut.dist_fast.signals.values;
    y_fast(:,:,i) = SimOut.y_fast.signals.values;
    y_slow(:,:,i) = SimOut.y_slow.signals.values;
    y_mmp(:,:,i) = SimOut.mmp.signals.values;
    t_fast = SimOut.y_fast.time;
    t_slow = SimOut.y_slow.time;
end
% y0: Q_bp
% y1: Q_socp
% y2: Q_fir
% y3: Q_iir
% y#(1):w_fir, y#(2):w_iir
d_out = d_out';
y0_tu = y_base_fast(:,1)'; % fast output for q_bp using w_iir
y1_tu = y_fast(:,:,1)'; % output for q_socp
y2_tu = y_fast(:,:,2)'; % output for q_fir
y3_tu = y_fast(:,:,3)'; % output for q_iir

y0_ts = y_base_slow(:,1)'; % slow output for q_bp
y1_ts = y_slow(:,:,1)'; % slow output for q_socp
y2_ts = y_slow(:,:,2)'; % output for q_fir
y3_ts = y_slow(:,:,3)'; % output for q_iir

mmp0 = y_base_fast(:,2)'; % output of mmp for Q_bp using w_iir
mmp1 = y_mmp(:,:,1)'; % output of mmp for q_socp
mmp2 = y_mmp(:,:,2)'; % output of mmp for q_fir
mmp3 = y_mmp(:,:,3)'; % output of mmp for q_iir

%% =============== Plotting bode plots ==================================
% defining the frequency range to be plotted and axes limits
close all
clc
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

% ================== Error ==========================================
err_mmp_base = abs(d_out(1,:)-mmp0);
err_mmp_fir(1,:) = abs(d_out(1,:)-mmp1(1,:)); % q-socp, w-fir
err_mmp_fir(2,:) = abs(d_out(2,:)-mmp2(1,:)); % q-fir, w-fir
err_mmp_fir(3,:) = abs(d_out(3,:)-mmp3(1,:)); % q-iir, w-fir
err_mmp_iir(1,:) = abs(d_out(1,:)-mmp1(2,:)); % q-socp, w-iir
err_mmp_iir(2,:) = abs(d_out(2,:)-mmp2(2,:)); % q-fir, w-iir
err_mmp_iir(3,:) = abs(d_out(3,:)-mmp3(2,:)); % q-iir, w-iir

rms_base = rms(err_mmp_base);
rms_fir(1) = rms(err_mmp_fir(1,:));
rms_fir(2) = rms(err_mmp_fir(2,:));
rms_fir(3) = rms(err_mmp_fir(3,:));
rms_iir(1) = rms(err_mmp_iir(1,:));
rms_iir(2) = rms(err_mmp_iir(2,:));
rms_iir(3) = rms(err_mmp_iir(3,:));

fprintf('RMS of Q-BP, W-FIR %d\n', rms(y_base_fast(:,4)))
fprintf('RMS of Q-BP, W-IIR: %d\n',rms(y0_tu));
fprintf('RMS of Q-SOCP, W-FIR: %d\n',rms(y1_tu(1,:)));
fprintf('RMS of Q-SOCP, W-IIR: %d\n',rms(y1_tu(2,:)));
fprintf('RMS of Q-FIR, W-FIR: %d\n',rms(y2_tu(1,:)));
fprintf('RMS of Q-FIR, W-IIR: %d\n',rms(y2_tu(2,:)));
fprintf('RMS of Q-IIR, W-FIR: %d\n',rms(y3_tu(1,:)));
fprintf('RMS of Q-IIR, W-IIR: %d\n',rms(y3_tu(2,:)));

% ==================== Plotting outputs ==============================
% ==================== Spectral of noise ==============================
% spec_Noise = specCal(SensorNoise.Data,1/SensorNoise.Ts);
% figure()
% semilogy(spec_Noise.f,spec_Noise.amp);
% title('Sensor Noise Spectral')
% xlabel('Hz')
% ylabel('Magnitude')
% 
% spec_F_dist = specCal(ForceDist.Data,1/ForceDist.Ts);
% figure()
% semilogy(spec_F_dist.f,spec_F_dist.amp);
% title('Force Disturbance Spectral')
% xlabel('Hz')
% ylabel('Magnitude')

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting - Simulink %%%%%%%%%%%%%%%%%%%%%%%%%
% y0: Q_bp
% y1: Q_socp
% y2: Q_fir
% y3: Q_iir
% y#(1,:):w_fir, y#(2,:):w_iir
y_out_lim = [-6 6];
x_lim = [1.58 1.6];
size_mark = 6;
l_width2 = 1;

figure()
h(1) = stairs(t_fast,y0_tu);
hold on
h(2) = stairs(t_fast,y1_tu(1,:));
h(3) = stairs(t_fast,y1_tu(2,:));
hold off
title('Output Position, Fast')
legend('Bandpass - W:IIR','SOCP - W:FIR','SOCP - W:IIR','location','southeast')
ax = gca;
ax.FontSize= font_size;
xlim(x_lim)
ylim([-10 10]);
xlabel('Time (sec)')
ylabel('Magnitude')
set(h,'linewidth',l_width);

for i = 1:3
    h(i).Color = color_base{i};
    h(i).MarkerSize = size_mark;
    h(i).LineWidth = l_width2;
end

figure()
h(1) = stairs(t_slow,y2_ts(1,:));
hold on
h(2) = stairs(t_slow,y2_ts(2,:));
h(3) = stairs(t_slow,y3_ts(1,:));
h(4) = stairs(t_slow,y3_ts(2,:));
hold off
legend('SDP:FIR - W:FIR','SDP:FIR - W:IIR','SDP:IIR - W:FIR','SDP:IIR - W:IIR','location','southeast')
title('Output Position, Slow')
ax = gca;
ax.FontSize= font_size;
xlim(x_lim)
ylim(y_out_lim);
xlabel('Time (sec)')
ylabel('Magnitude')
set(h,'linewidth',l_width);
for i = 1:n_all
    h(i).Color = color_cvx{i};
    h(i).LineStyle = line_style{i};
    h(i).Marker = marker_style{i};
    h(i).MarkerSize = size_mark;
    h(i).LineWidth = l_width2;
end

figure()
h(1) = stairs(t_fast,y2_tu(1,:));
hold on
h(2) = stairs(t_fast,y2_tu(2,:));
h(3) = stairs(t_fast,y3_tu(1,:));
h(4) = stairs(t_fast,y3_tu(2,:));
hold off
title('Output Position, Fast')
legend('SDP:FIR - W:FIR','SDP:FIR - W:IIR','SDP:IIR - W:FIR','SDP:IIR - W:IIR','location','southeast')
ax = gca;
ax.FontSize= font_size;
xlim(x_lim)
ylim(y_out_lim);
xlabel('Time (sec)')
ylabel('Magnitude')
set(h,'linewidth',l_width);
for i = 1:n_all
    h(i).Color = color_cvx{i};
    h(i).LineStyle = line_style{i};
    h(i).Marker = marker_style{i};
    h(i).MarkerSize = size_mark;
    h(i).LineWidth = l_width2;
end

spec_Qsocp_fir = specCal(y1_tu(1,:),1/Tu);
spec_Qsocp_iir = specCal(y1_tu(2,:),1/Tu);
spec_QFIR_FIR = specCal(y2_tu(1,:),1/Tu);
spec_QFIR_IIR = specCal(y2_tu(2,:),1/Tu);
spec_QIIR_FIR = specCal(y3_tu(1,:),1/Tu);
spec_QIIR_IIR = specCal(y3_tu(2,:),1/Tu);

figure()
h(1) = semilogy(spec_Qsocp_fir.f,spec_Qsocp_fir.amp);
hold on
h(2) = semilogy(spec_Qsocp_iir.f,spec_Qsocp_iir.amp);
ax = gca;
ax.FontSize= font_size;
x_line = xline(Nyq_Hz);
x_line.Color = [0 0 0];
x_line.LineWidth = 1.5;
x_line.Label = sprintf('%.f Hz',Nyq_Hz);
x_line.LabelOrientation = 'aligned';
x_line.LabelVerticalAlignment = 'bottom';
x_line.LabelHorizontalAlignment = 'center';
x_line.FontSize = 10;
x_line.FontWeight = 'bold';
for i = 1:m_d
x_line = xline(f_d(i));
x_line.Color = [0 0 0];
x_line.LineWidth = 1;
x_line.Label = sprintf('%.f Hz',f_d(i));
x_line.LabelOrientation = 'aligned';
x_line.LabelVerticalAlignment = 'bottom';
x_line.LabelHorizontalAlignment = 'center';
x_line.FontSize = 10;
x_line.FontWeight = 'bold';
x_line.LineStyle = '--';
end
xlim([spec_Qsocp_fir.f(1) spec_Qsocp_fir.f(end)])
ylabel('Magnitude')
xlabel('Hz')
legend('W:FIR','W:IIR','location','southeast')
title('Y output for SOCP filters')
for i = 1:2
    h(i).Color = color_base{i+1};
    h(i).LineStyle = line_style{i};
end
set(h,'linewidth',1);
hold off

figure()
h(1) = semilogy(spec_QFIR_FIR.f,spec_QFIR_FIR.amp);
hold on
h(2) = semilogy(spec_QFIR_IIR.f,spec_QFIR_IIR.amp);
ax = gca;
ax.FontSize= font_size;
x_line = xline(Nyq_Hz);
x_line.Color = [0 0 0];
x_line.LineWidth = 1.5;
x_line.Label = sprintf('%.f Hz',Nyq_Hz);
x_line.LabelOrientation = 'aligned';
x_line.LabelVerticalAlignment = 'bottom';
x_line.LabelHorizontalAlignment = 'center';
x_line.FontSize = 10;
x_line.FontWeight = 'bold';
for i = 1:m_d
x_line = xline(f_d(i));
x_line.Color = [0 0 0];
x_line.LineWidth = 1;
x_line.Label = sprintf('%.f Hz',f_d(i));
x_line.LabelOrientation = 'aligned';
x_line.LabelVerticalAlignment = 'bottom';
x_line.LabelHorizontalAlignment = 'center';
x_line.FontSize = 10;
x_line.FontWeight = 'bold';
x_line.LineStyle = '--';
end
xlim([spec_QIIR_FIR.f(1) spec_QIIR_FIR.f(end)])
ylabel('Magnitude')
xlabel('Hz')
legend('W:FIR','W:IIR','location','southeast')
title('Y output for SDP:FIR filters')
for i = 1:2
    h(i).Color = color_cvx{i};
    h(i).LineStyle = line_style{i};
end
set(h,'linewidth',1);
hold off

figure()
h(3) = semilogy(spec_QIIR_FIR.f,spec_QIIR_FIR.amp);
hold on
h(4) = semilogy(spec_QIIR_IIR.f,spec_QIIR_IIR.amp);
ax = gca;
ax.FontSize= font_size;
x_line = xline(Nyq_Hz);
x_line.Color = [0 0 0];
x_line.LineWidth = 1.5;
x_line.Label = sprintf('%.f Hz',Nyq_Hz);
x_line.LabelOrientation = 'aligned';
x_line.LabelVerticalAlignment = 'bottom';
x_line.LabelHorizontalAlignment = 'center';
x_line.FontSize = 10;
x_line.FontWeight = 'bold';
for i = 1:m_d
x_line = xline(f_d(i));
x_line.Color = [0 0 0];
x_line.LineWidth = 1;
x_line.Label = sprintf('%.f Hz',f_d(i));
x_line.LabelOrientation = 'aligned';
x_line.LabelVerticalAlignment = 'bottom';
x_line.LabelHorizontalAlignment = 'center';
x_line.FontSize = 10;
x_line.FontWeight = 'bold';
x_line.LineStyle = '--';
end
xlim([spec_QIIR_FIR.f(1) spec_QIIR_FIR.f(end)])
ylabel('Magnitude')
xlabel('Hz')
legend('W:FIR','W:IIR','location','southeast')
title('Y output for SDP:IIR filters')
for i = 3:4
    h(i).Color = color_cvx{i};
    h(i).LineStyle = line_style{i};
end
set(h,'linewidth',1);
hold off

% mmp0 = y_base_fast(:,2)'; % output of mmp for Q_bp using w_iir
% mmp1 = y_mmp(:,:,1)'; % output of mmp for q_socp
% mmp2 = y_mmp(:,:,2)'; % output of mmp for q_fir
% mmp3 = y_mmp(:,:,3)'; % output of mmp for q_iir

% MMP output
spec_y_MMP_QF_FIR = specCal(mmp2(1,:),1/Tu);
spec_y_MMP_QF_IIR = specCal(mmp2(2,:),1/Tu);
spec_y_MMP_QI_FIR = specCal(mmp3(1,:),1/Tu);
spec_y_MMP_QI_IIR = specCal(mmp3(2,:),1/Tu);

figure
h(1) = semilogy(spec_y_MMP_QF_FIR.f,spec_y_MMP_QF_FIR.amp);
hold on
h(2) = semilogy(spec_y_MMP_QF_IIR.f,spec_y_MMP_QF_IIR.amp);
title('MMP Output')
for i = 1:2
    % h(i).Color = color_all{i+2};
    h(i).LineStyle = line_style{i+2};
end
set(h,'linewidth',1);
ax = gca;
ax.FontSize= font_size;
x_line = xline(Nyq_Hz);
x_line.Color = [0 0 0];
x_line.LineWidth = 1.5;
x_line.Label = sprintf('%.f Hz',Nyq_Hz);
x_line.LabelOrientation = 'aligned';
x_line.LabelVerticalAlignment = 'bottom';
x_line.LabelHorizontalAlignment = 'center';
x_line.FontSize = 10;
x_line.FontWeight = 'bold';
for i = 1:m_d
x_line = xline(f_d(i));
x_line.Color = [0 0 0];
x_line.LineWidth = 1;
x_line.Label = sprintf('%.f Hz',f_d(i));
x_line.LabelOrientation = 'aligned';
x_line.LabelVerticalAlignment = 'bottom';
x_line.LabelHorizontalAlignment = 'center';
x_line.FontSize = 10;
x_line.FontWeight = 'bold';
x_line.LineStyle = '--';
end
legend('W:FIR','W:IIR','location','southeast')
xlim([spec_QIIR_FIR.f(1) spec_QIIR_FIR.f(end)])
ylabel('Magnitude')
xlabel('Hz')

% robust stability analysis, T and inv(T)
% states that |T(z) w| < 1 for all z=e^j\omega to be stable with
% uncertainty. w is characterized as the magnitude of the uncertainty
Ps = tf(Ps); % define CT plnat as a transfer function
Pz_delay = c2d(Ps,Tu,'zoh'); % discretize the plant
Pz_num = cell2mat(Pz_delay.num);
Pz_den = cell2mat(Pz_delay.den);
Pz_nd = tf(Pz_num,Pz_den,Tu);
T_1 = feedback(Pz*Cz,1);
T_2 = feedback(Pz,Cz,-1);

T_all(1) = minreal(T_1 + T_2*Qcvx_FIR*W_k_FIR(1)); % Q:FIR, W:FIR
T_all(2) = minreal(T_1 + T_2*Qcvx_FIR*W_k_IIR(1)); % Q:FIR, W:IIR
T_all(3) = minreal(T_1 + T_2*Qcvx_IIR*W_k_FIR(1)); % Q:IIR, W:FIR
T_all(4) = minreal(T_1 + T_2*Qcvx_IIR*W_k_IIR(1)); % Q:IIR, W:IIR
T_all(5) = minreal(T_1 + T_2*Qcvx_FIR*W_k_FIR(2)); % Q:FIR, W:FIR
T_all(6) = minreal(T_1 + T_2*Qcvx_FIR*W_k_IIR(2)); % Q:FIR, W:IIR
T_all(7) = minreal(T_1 + T_2*Qcvx_IIR*W_k_FIR(2)); % Q:IIR, W:FIR
T_all(8) = minreal(T_1 + T_2*Qcvx_IIR*W_k_IIR(2)); % Q:IIR, W:IIR
T_all(9) = minreal(T_1 + T_2*Qcvx_FIR); % Q:FIR
T_all(10) = minreal(T_1 + T_2*Qcvx_IIR); % Q:IIR


% close all
% bode of inv(T)
Tinv = tf();
for i = 1:10
Tinv(i) = inv(T_all(i));
end
[mag, phi, ~] = bode(Tinv,w_in_rad);
mag = 20*log10(squeeze(mag));
phi = wrapTo180(squeeze(phi));
figure
hold on
for i = 1:4
    h(i) = semilogx(w_in_Hz,mag(i,:));
    h(i).Color = color_cvx{i};
    h(i).LineStyle = line_style{i};
    h(i).LineWidth = 1.5;
end
hold off
ax = gca;
ax.FontSize= font_size;
xlim([1 2640])
xlabel('Hz')
ylabel('Magnitude (dB)')
legend('SDP:FIR - W:FIR','SDP:FIR - W:IIR','SDP:IIR - W:FIR','SDP:IIR - W:IIR','location','southeast')
title('W(1), inv(T)')

figure
hold on
for i = 1:4
    h(i+4) = semilogx(w_in_Hz,mag(i+4,:));
    h(i+4).Color = color_cvx{i};
    h(i+4).LineStyle = line_style{i};
    h(i+4).LineWidth = 1.5;
end
hold off
ax = gca;
ax.FontSize= font_size;
xlim([1 2640])
xlabel('Hz')
ylabel('Magnitude (dB)')
legend('SDP:FIR - W:FIR','SDP:FIR - W:IIR','SDP:IIR - W:FIR','SDP:IIR - W:IIR','location','southeast')
title('W(2), inv(T)')

figure
hold on
for i = 1:2
h(i+8) = semilogx(w_in_Hz,mag(i+4,:));
h(i+8).LineStyle = line_style{i};
h(i+8).LineWidth = 1.5;
end
hold off
ax = gca;
ax.FontSize= font_size;
xlim([1 2640])
xlabel('Hz')
ylabel('Magnitude (dB)')
legend('SDP:FIR','SDP:IIR','location','southeast')
title('1/|T(z)|')
% % bode of T
% [mag, phi, ~] = bode(T_all,w_in_rad);
% mag = 20*log10(squeeze(mag));
% phi = wrapTo180(squeeze(phi));
% figure
% hold on
% for i = 1:4
% h(i) = semilogx(w_in_Hz,mag(i,:));
% h(i).Color = color_all{i};
% h(i).LineStyle = line_style{i};
% h(i).LineWidth = 2.5-0.1*i;
% end
% ax = gca;
% ax.FontSize= font_size;
% xlim([1 2640])
% xlabel('Hz')
% ylabel('Magnitude (dB)')
% legend('Q:FIR - W:FIR','Q:FIR - W:IIR','Q:IIR - W:FIR','Q:IIR - W:IIR','location','southeast')
% title('T')

%% other plots
opts = bodeoptions;
opts.FreqUnits = 'Hz';
opts.PhaseWrapping = 'on';
figure
bodeplot(Qcvx_FIR,Qcvx_IIR,opts)
legend('Q:FIR','Q:IIR','location','southeast')
figure
bodeplot(W_k_IIR(1),W_k_IIR(2),opts)
legend('W_1','W_2','location','southeast')

figure
bodeplot(PQcvx_FIR,PQcvx_IIR,opts)
legend('PQ_FIR','PQ_IIR','location','southeast')

figure
bodeplot(PQcvx_FIR*W_k_FIR(1),PQcvx_IIR*W_k_IIR(1),opts)
legend('PWQ_FIR','PWQ_IIR','location','southeast')

S1 = feedback(Pz*Cz,1);
figure
bodeplot(S1,opts);
title('PC/(1+PC)')

%% calculating max uncertainty allowed
T_all(1) = minreal(T_1 + T_2*Qcvx_FIR*W_k_FIR(1)); % Q:FIR, W:FIR
T_all(2) = minreal(T_1 + T_2*Qcvx_FIR*W_k_IIR(1)); % Q:FIR, W:IIR
T_all(3) = minreal(T_1 + T_2*Qcvx_IIR*W_k_FIR(1)); % Q:IIR, W:FIR
T_all(4) = minreal(T_1 + T_2*Qcvx_IIR*W_k_IIR(1)); % Q:IIR, W:IIR
T_all(5) = minreal(T_1 + T_2*Qcvx_FIR*W_k_FIR(2)); % Q:FIR, W:FIR
T_all(6) = minreal(T_1 + T_2*Qcvx_FIR*W_k_IIR(2)); % Q:FIR, W:IIR
T_all(7) = minreal(T_1 + T_2*Qcvx_IIR*W_k_FIR(2)); % Q:IIR, W:FIR
T_all(8) = minreal(T_1 + T_2*Qcvx_IIR*W_k_IIR(2)); % Q:IIR, W:IIR
h_inf = [];
for i_it = 1:8
    h_inf(i_it) = hinfnorm(T_all(i_it));
end
h_inf = 1./h_inf;
h_inf_db = 20*log10(h_inf);