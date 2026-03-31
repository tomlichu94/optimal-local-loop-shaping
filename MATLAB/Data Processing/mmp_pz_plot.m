function [] = mmp_pz_plot(G)
% plotting pzmap
figure
pzmap(G) % map of IIR-MMP poles
a = findobj(gca,'type','line');
for i = 1:length(a)
    set(a(i),'markersize',8)
    set(a(i),'linewidth',1.5)
end
end