function [] = mmp_bode_plot(G, w_in)
% Function used to autoplot Bode plots quickly
% Variables used
%   G - transfer function for bode plot
%   w_in - frequency range, keep in rad/s

[mag phi ~] = bode(G, w_in);
mag = 20*log10(squeeze(mag));
phi = wrapTo180(squeeze(phi));
db_lim = max(abs(mag));

% plotting
figure()
    subplot(2,1,1)
    h = semilogx(w_in,mag);
    hold on
    h.Color = 'k';
    h.LineWidth = 1;
    xlabel('rad/s')
    ylabel('Magnitude (dB)')
    xlim([w_in(1), w_in(end)])
    ylim([min(mag)-0.1*db_lim, max(mag)+0.1*db_lim]);
    
    subplot(2,1,2)
    h = semilogx(w_in,phi);
    hold on
    h.Color = 'k';
    h.LineWidth = 1;
    xlabel('Hz')
    ylabel('Phase (Degrees)')
    xlim([w_in(1), w_in(end)])
    ylim([-200, 200]);
end