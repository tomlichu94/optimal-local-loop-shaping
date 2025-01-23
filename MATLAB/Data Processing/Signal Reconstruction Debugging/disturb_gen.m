function [dc,d_fs,d_ss] = disturb_gen(f_dist, A_dist, dc_t, d_fs_t, d_ss_t,p_off)
%MULTISIN Generate sensor measurements of harmonic signal.
%   Detailed explanation goes here
dc = multisin(f_dist,A_dist,dc_t,p_off);
d_fs = multisin(f_dist,A_dist,d_fs_t,p_off);
d_ss = multisin(f_dist,A_dist,d_ss_t,p_off);
    function y = multisin(f_dist,A,time_series,p_off)
        y = zeros(1,length(time_series));
        for i = 1:length(f_dist)
            y = y + A(i)*sin(2*pi*f_dist(i)*time_series+p_off(i));
        end
    end
end