clear
%plot barotropic stream fucntion for Atlantic
load('100free.mat');

ncload('~/srvccrc/data01/z3263455/escm27/pdetop5/100free/tavg.03001.01.01.nc','zw');
ncload('~/srvccrc/data01/z3263455/escm27/pdetop5/100free/tavg.03001.01.01.nc','v');
ncload('~/srvccrc/data01/z3263455/escm27/pdetop5/100free/tavg.03001.01.01.nc','mskhr');

tmp = find(gmoc == 0);
gmoc(tmp)=NaN;


%line contour plot with labels put on mauaully
v1 = [-80:5:0];
v2 = [0:5:80];

figure
title('1.8^\circ \times 3.6');
%this plots negative and positive values seperately
[C1,h1] = contour(yu(1:end),-zw(1:end),gmoc(1:end,1:end)./1e6,v1,':k');
hold on;
[C2,h2] = contour(yu(1:end),-zw(1:end),gmoc(1:end,1:end)./1e6,v2,'-k');
set(h1,'LineWidth',2);
set(h2,'LineWidth',2);
shading flat;
colormap(gray);
%title('100free: am=2e9 athkdf=4e6');
clabel(C1,h1,v1);
clabel(C2,h2,v2);
%clabel(C1,'manual');
%clabel(C2,'manual');
hold on;
caxis manual;
shading flat;
set(gca,'TickDir','out');
set(gca,'TickDir','out');
%set(gca,'XLim',[ 80.0]);
%set(gca,'XTick',[0 20 40 60 80]);
%set(gca,'XTickLabel',{'0', '20N', '40N','60N','80N'})
