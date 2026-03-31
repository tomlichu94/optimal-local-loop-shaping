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

PQ_max = 8; % max value of PQ for the quadratic constraint
beta = PQ_max^2; % FIR SDP, set max value for quad, play around with quadratic
f_stop = 600; % SOCP stop constraint past this frequency
a_g_IIR = 0.90; % predictor alpha

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
% L = 3;
% tempW = [1.337 1.739 2.312 2.618]; % okay, bad robustness
tempW = [1.32 1.67 1.93 2.18]; % good run
% tempW = [1.4482 1.7411 2.0070 2.8114]; % okay, mid robustness
% tempW = [0.23 0.56 0.8];
% L_t = 3;
% tempW = [1.7178 2.0072 2.9787 3.5073];

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
phi_val = freqresp(phi,exp(1j*w_lin(1:stop_indx)));
phi_val = squeeze(phi_val);
phi_r = real(phi_val);
phi_i = -imag(phi_val);
Pz_val = freqresp(PQ0,exp(1j*w_lin(1:stop_indx)));
Pz_val = squeeze(Pz_val)';
Pz_2 = real(Pz_val).^2 + imag(Pz_val).^2;
quad = zeros(max_order+1,max_order+1,stop_indx);
for i = 1:stop_indx
    quad(:,:,i) = Pz_2(i)*phi_r(:,i)*phi_r(:,i)'+...
                      Pz_2(i)*phi_i(:,i)*phi_i(:,i)';
end
cvx_begin
        variables q_vec((max_order+1),1) b(stop_indx,1)
        beta_sum = sum(b);
        for i = 1:stop_indx
           quad_q(i) = quad_form(q_vec,quad(:,:,i));
        end
        minimize beta_sum
        subject to
            for i = 1:m_d
               1-PQ_d(1,:,i)*q_vec == 0;
            end            
            for i = 1:stop_indx
                quad_q(i) <= b(i)
            end
cvx_end
Qcvx_num = conv(Q0_num,q_vec');
Qcvx_den = conv(Q0_den,[1 zeros(1,max_order)]);
Qcvx_socp = tf(Qcvx_num,Qcvx_den,Tu);
Qcvx_socp = minreal(Qcvx_socp);
PQ_cvx_quad = minreal(Pz*Qcvx_socp);
T_cvx_quad = minreal(1-PQ_cvx_quad);

%% =============== Plotting bode plots ==================================
% defining the frequency range to be plotted and axes limits
w_in_Hz_end = 1/(2*Tu); % end of plotting for Nyq freq of the fast
w_in_Hz = 1:1:w_in_Hz_end;
w_in_rad = w_in_Hz*2*pi;
Nyq_Hz = 1/(2*Ts);
w_start = 1; % starting frequency for x_lim
l_width = 1.2; % linewidth
font_size = 11; % font size
x_lim_loop = [100 w_in_Hz_end]; % 1-PQ x_lim
y_lim = [-50 20]; % y_lim

% ==================== plotting style ==================================
color_all = {[0 0 0], [0.9290 0.6940 0.1250], [0 1 0], [1 0 0], [0 0 1]};
marker_style = {'none','diamond','o','x','square'};
line_style = {'-','--','--',':'};
n_all = 4;

% =================== Q-filter Bode ====================================
[mag_T_socp, phi_T_socp, ~] = bode(T_cvx_quad,w_in_rad);
[mag_PQ_socp, phi_PQ_socp, ~] = bode(PQ_cvx_quad,w_in_rad);
mag_T_socp = 20*log10(mag_T_socp(:));
mag_PQ_socp = 20*log10(mag_PQ_socp(:));
phi_T_socp = wrapTo180(phi_T_socp(:));
phi_PQ_socp = wrapTo180(phi_PQ_socp(:));
% ======================== PQ plots =====================================
figure()
subplot(2,1,1)
h = semilogx(w_in_Hz,mag_PQ_socp);
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
h(1).Color = [0 1 0];
hold off
% ============================ 1-PQ plots ===========================
subplot(2,1,2)
h = semilogx(w_in_Hz,mag_T_socp);
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
h(1).Color = [0 1 0];
legend('Q:SOCP','location','southwest')