function w_dt = w_lin_spacing(max_order,range_mult,tempW,L_t)
% w_lin_spacing spaces the discrete set of frequencies such that the
% constraints are spaced to be away from the disturbance frequencies
% (bandpass frequencies are generated here)

% number of disturbances
m_d = length(tempW); 

% array with [starting frequency, disturbance frequencies, end frequency]
w_d_range = [0 tempW L_t]*(pi/L_t); 

% identify the nth value where the disturbance is at
n_it = floor(w_d_range*range_mult*max_order/pi)+1;
n_it(end) = max_order*range_mult;

% generate set of frequencies without including frequencies near the
% disturbance
if m_d == 1
    w_dt(1:n_it(2)) = linspace(0,w_d_range(2)-pi/max_order, ...
                                  n_it(2)-n_it(1)+1);
    w_dt(n_it(end-1):n_it(end)) = linspace(w_d_range(end-1)+ ...
                              pi/max_order,pi,n_it(end)-n_it(end-1)+1);
else
    w_dt(1:n_it(2)) = linspace(0,w_d_range(2)-pi/max_order, ...
                                  n_it(2)-n_it(1)+1);
    for i = 2:m_d
        w_dt(n_it(i):n_it(i+1)) = linspace(w_d_range(i)+pi/max_order, ...
                      w_d_range(i+1)-pi/max_order,n_it(i+1)-n_it(i)+1);
    end
    w_dt(n_it(end-1):n_it(end)) = linspace(w_d_range(end-1)+ ...
                            pi/max_order,pi,n_it(end)-n_it(end-1)+1);
end
end