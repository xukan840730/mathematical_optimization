%% plot evolution of PDDs at a point over time.
%  requires presence of arrays generated in cum_time_mapper.m
close all
year = 1850:2839;
lim(1:1000)=200.;
clear p
figure
hold on
p(1)=plot(year,posdd(:,97,79), 'linewidth', 2, 'color', 'k','LineStyle','--');
p(2)=plot(year,posdd(:,14,84), 'linewidth', 2, 'color', 'r','LineStyle','--');
p(3)=plot(year,posdd(:,12,80), 'linewidth', 2, 'color', 'b','LineStyle','--');
p(4)=plot(year(152:990),posddcom((1:839),97,79), 'linewidth', 2, 'color', 'k');
p(5)=plot(year(152:990),posddcom((1:839),14,84), 'linewidth', 2, 'color', 'r');
p(6)=plot(year(152:990),posddcom((1:839),12,80), 'linewidth', 2, 'color', 'b');
xlabel('Time (years)', 'FontSize', 14)
ylabel('PDD/year', 'FontSize', 14)
legend(p,'Ward Hunt (high)','Larsen B (high)','Wilkins (high)','Ward Hunt (com)','Larsen B (com)','Wilkins (com)');
legend('Location', 'NorthWest')
plot(year,lim(1:990), 'linewidth', 2, 'color', 'g')
axis image
axis normal
set(gca, 'FontSize', 14)
print -depsc2 figures/temporalvalid
hold off
figure
hold on
clear p
p(1)=plot(year,posdd(:,5,53), 'linewidth', 2, 'color', 'k','LineStyle','--');
p(2)=plot(year,posdd(:,5,84), 'linewidth', 2, 'color', 'r','LineStyle','--');
p(3)=plot(year,posdd(:,11,20), 'linewidth', 2, 'color', 'g','LineStyle','--');
p(4)=plot(year,posdd(:,94,91), 'linewidth', 2, 'color', 'b','LineStyle','--');
p(5)=plot(year(152:990),posddcom((1:839),5,53), 'linewidth', 2, 'color', 'k');
xlabel('Time (years)', 'FontSize', 14)
ylabel('PDD/year', 'FontSize', 14)
legend(p, 'Ross (high)', 'Ronne-Filchner (high)', 'Amery (high)', 'GIS (high)', 'All (com)');
legend('Location', 'NorthWest')
axis image
axis normal
set(gca, 'FontSize', 14)
print -depsc2 figures/temporalinterior
hold off