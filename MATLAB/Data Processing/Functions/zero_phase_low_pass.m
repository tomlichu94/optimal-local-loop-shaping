function filtered_out = zero_phase_low_pass(data_in)
% function for applying zero phase error with low pass filtering. This
% filter can only be applied a posteriori.
% Q_filter = (1+z^-1)(1+z)/4 = (z+2-z^-1)/4, can only do this offline,
% where z^{-1} is a one-step delay operator.

n_size = size(data_in);
n_it = max(n_size);
filtered_out(1) = data_in(1);
for i = 2:(n_it-1)
    filtered_out(i) = 0.25*(data_in(i+1)+2*data_in(i)-data_in(i-1));
end