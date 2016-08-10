%%Load lgm and modern ocean temperatures
clear all
ncload('lgm.nc','O_temp');
lgm_temp=O_temp;
i=find(lgm_temp > 1.e20);
lgm_temp(i) = nan;

ncload('modern.nc','O_temp');
modern_temp=O_temp;
i=find(modern_temp > 1.e20);
modern_temp(i) = nan;

clear O_temp

ncload('modern.nc','zw');

%%Get difference in ocean temperatures (moden-lgm)
diff=modern_temp-lgm_temp;

%%Plot temperature difference by layer

for n=1:18
  temp(:,:)=diff(n,:,:);
  subplot(9,2,n);
  [c,h]=contour(temp);
  clabel(c,h);
  %pcolor(temp), shading flat;
  set(gca,'xtick',[],'ytick',[]);
  text(1,75,num2str(zw(n)));
  caxis([-1 10]);  
  colorbar;
end





