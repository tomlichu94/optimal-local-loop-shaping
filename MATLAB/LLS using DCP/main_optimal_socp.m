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

%% SOCP implementation
%%%%%%%%%%%%%%%%%%%% Quadratic Constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%
% sub in non-disutrance frequencies, where z = e^{jw}
% such that phi_val = phi(e^{jw}) = [1, e^{-jw}, e^{-2*jw} ... e^{-r*jw}]
phi_val = freqresp(phi,exp(1j*w_lin(1:length(w_lin)))); 
phi_val = squeeze(phi_val); % remove extra dimension

% separate real and imag terms
phi_r = real(phi_val);
phi_i = -imag(phi_val);

% finding values of the P(z) using z = e^{jw}
Pz_val = freqresp(P_phi,exp(1j*w_lin(1:length(w_lin))));
Pz_val = squeeze(Pz_val)';
Pz_2 = real(Pz_val).^2 + imag(Pz_val).^2; % squared values of P(z)

% create constraint array to store the values in it
% quad is defined in Eq. (18) as q^T*quad*q, where q is our unknown
quad = zeros(max_order+1,max_order+1,length(w_lin));
for i = 1:length(w_lin)
    quad(:,:,i) = Pz_2(i)*phi_r(:,i)*phi_r(:,i)'+...
                  Pz_2(i)*phi_i(:,i)*phi_i(:,i)';
end

% DCP for SOCP
cvx_begin
        % q_vec is coefficients of Eq. (14)
        % rho is from Eq. (18)
        variables q_vec((max_order+1),1) rho(length(w_lin),1) 
        rho_sum = sum(rho);
        
        % formulating quadratic with unknown variable, see Eq. (18)
        for i = 1:length(w_lin)
           quad_q(i) = quad_form(q_vec,quad(:,:,i));
        end

        % objective function, see Eq. (19)
        minimize rho_sum

        % constraint functions
        subject to
            % bandstop constraint
            for i = 1:m_d
               1-PQ_d(1,:,i)*q_vec == 0;
            end
            % bandpass constraint
            for i = 1:length(w_lin)
                quad_q(i) <= rho(i);
            end
cvx_end

Qcvx_num = q_vec'; % coefficients of Eq. (14)
Qcvx_den = [1 zeros(1,max_order)];

% coefficients of PQ
PQcvx_num = conv(Pz_num,Qcvx_num);
PQcvx_den = conv(Pz_den,Qcvx_den);

% transfer function for P(z)Q(z)
PQ_socp = tf(PQcvx_num,PQcvx_den,Tu);

% transfer function of 1-P(z)Q(z), Eq. (10)
T_socp = minreal(1-PQ_socp);

%% bandpass filter baseline
w_hz_d = w_d/(2*pi*Tu);
B_bw = w_hz_d*0.2; % bandwidth
[Q_bp, Q0_bp] = q_bandpass(Pz,w_hz_d,B_bw,Tu);
PQ_bp = minreal(Pz*Q_bp);
T_bp = 1-PQ_bp;


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
n_all = 2;
Nyq_Hz = 1/(2*Ts);

% =================== Q-filter Bode ====================================
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
h = semilogx(w_in_Hz,mag_PQ_bp,w_in_Hz,mag_PQ_socp);
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
h = semilogx(w_in_Hz,mag_T_bp,w_in_Hz,mag_T_socp);
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
legend('Bandpass','SOCP:FIR','southwest')