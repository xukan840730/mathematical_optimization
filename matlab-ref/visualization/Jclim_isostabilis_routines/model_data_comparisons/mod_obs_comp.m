%% script to interpolate NCEP SAT to UVic grid for comparison with UVic output
clear all
close all

%% load data

ncload('air.mon.mean.nc','air');
ncload('NSHWINDtavg.01951.01.01.nc','sat');

%% extract 1950-2007 from UVic ESCM data

UVic_sat=sat(1:58,:,:)-273.15;
clear sat

%% extract 1950-2007 from NCEP data

NCEP_tmp=air((25:720),:,:);
clear air

%% set annual mean array size
NCEP_arr_dim = size(NCEP_tmp);
NCEP_sat= zeros((NCEP_arr_dim(1)/12),NCEP_arr_dim(2),NCEP_arr_dim(3));

clear NCEP_arr_dim

%% accumulate annual mean from monthly means

for year=1:58;
  for month=1:12;
    nrec=(year-1)*12 + month;
    NCEP_sat(year,:,:) = NCEP_sat(year,:,:) + NCEP_tmp(nrec,:,:);
  end
  if month == 12;
    NCEP_sat(year,:,:) = NCEP_sat(year,:,:)/12;
  end
end

clear nrec NCEP_tmp month

%% interpolate to UVic grid

UVic_arr_dim = size(UVic_sat);
NCEP_arr_dim = size(NCEP_sat);

xi=linspace(1,NCEP_arr_dim(3),100);
yi=linspace(1,NCEP_arr_dim(2),100);
xo=1:1:NCEP_arr_dim(3);
yo=1:1:NCEP_arr_dim(2);

[X,Y] = meshgrid(xo,yo);
[XI,YI] = meshgrid(xi,yi);

NCEP_sat_interp = zeros(58,UVic_arr_dim(3),UVic_arr_dim(2));
for year = 1:58;
    interp_in_array(:,:)=NCEP_sat(year,:,:);  
    NCEP_sat_interp(year,:,:)=interp2(X,Y,interp_in_array,XI,YI,'*spline');
end

clear xo yo xi yi X Y XI YI NCEP_sat UVic_arr_dim NCEP_arr_dim

%% create differences between UVic and NCEP timeslices

NCEP_minus_UVic = zeros(58,100,100);

for year = 1:58;
  NCEP_slice(:,:)=NCEP_sat_interp(year,:,:);
  NCEP_slice=flipud(NCEP_slice);
  
  UVic_slice(:,:)=UVic_sat(year,:,:);
 
  for i=1:100
    for j=1:100
      UVic_minus_NCEP(year,i,j) = UVic_slice(i,j)-NCEP_slice(i,j);     
    end
  end
  
end

clear year i j

%% plot
close all
lat_array=linspace(-89.1,89.1);lon_array=linspace(1.8,358.2);
[Plg,Plt]=meshgrid(lon_array,lat_array);

hold on
m_proj('stereographic','longitude', 300,'latitude',80,'radius', 30);
hold on
m_pcolor(Plg,Plt,NCEP_slice);shading flat;
m_coast('color',[0 0 0]);
m_grid('xticklabels',[],'yticklabels',[]);

h=colorbar('h');
title('NCEP Greenland/Ellesmere Annual Average SAT, 2005 (C)','FontSize', 20);

print -djpeg90 ../../NCEP_NH_sat.jpg
hold off

close


hold on
m_proj('stereographic','longitude', 300,'latitude',80,'radius', 30);
hold on
m_pcolor(Plg,Plt,UVic_slice);shading flat;
m_coast('color',[0 0 0]);
m_grid('xticklabels',[],'yticklabels',[]);

h=colorbar('h');
title('UVic Greenland/Ellesmere Annual Average SAT, 2005 (C)','FontSize', 20);

print -djpeg90 ../../UVic_NH_sat.jpg
hold off

%% movie

arr_sz = size(UVic_minus_NCEP);
lat_array=linspace(-89.1,89.1);lon_array=linspace(1.8,358.2);
[Plg,Plt]=meshgrid(lon_array,lat_array);
%m_proj('miller','longitude', 180);
m_proj('stereographic','longitude', 300,'latitude',80,'radius', 30);

for n=1:arr_sz(1);
  plot_var(:,:)=UVic_minus_NCEP(n,:,:); 
  
  m_pcolor(Plg,Plt,plot_var);shading flat;
  m_coast('color',[0 0 0]);
  %m_grid('xaxis','bottom');
  m_grid('xticklabels',[],'yticklabels',[]);
  h=colorbar('h');
  
  date = 1947+n;
  year = 'SAT diff: UVic ESCM - NCEP, year= ';
  date_text = num2str(1947+n);
  tag = strcat(year,date_text);
  title(tag, 'Fontsize', 20)
  
  M(n) = getframe(gcf);

end
movie(M)
movie2avi(M,'NCEP-UVic_NH.avi')

close all

