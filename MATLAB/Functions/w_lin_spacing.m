function w_dt = w_lin_spacing(max_order,range_mult,tempW,L_t)
% range_mult = 30;
m_d = length(tempW); % number of disturbances
w_d_range = [0 tempW L_t]*(pi/L_t); % range of disturb in rad for T_{fs}
n_it_d = floor(w_d_range*range_mult*max_order/pi)+1; % the nth point of disturbance
n_it_d(end) = max_order*range_mult;
% generate set of non-disturbance frequencies
if m_d == 1
    w_dt(1:n_it_d(2)) = linspace(0,w_d_range(2)-pi/max_order, ...
                                  n_it_d(2)-n_it_d(1)+1);
    w_dt(n_it_d(end-1):n_it_d(end)) = linspace(w_d_range(end-1)+ ...
                              pi/max_order,pi,n_it_d(end)-n_it_d(end-1)+1);
else
    w_dt(1:n_it_d(2)) = linspace(0,w_d_range(2)-pi/max_order, ...
                                  n_it_d(2)-n_it_d(1)+1);
    for i = 2:m_d
        w_dt(n_it_d(i):n_it_d(i+1)) = linspace(w_d_range(i)+pi/max_order, ...
                      w_d_range(i+1)-pi/max_order,n_it_d(i+1)-n_it_d(i)+1);
    end
    w_dt(n_it_d(end-1):n_it_d(end)) = linspace(w_d_range(end-1)+ ...
                            pi/max_order,pi,n_it_d(end)-n_it_d(end-1)+1);
end
end