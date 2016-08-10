%% script to interpolate RECON SAT to UVic grid for comparison with UVic output
clear all
close all

%% load data

RECON_tmp=load('monaghan_recon_nsfc_temp_1960_2005_ann.txt');
ncload('NSHWINDtavg.01951.01.01.nc','sat');

%% extract 1960-2005 from UVic ESCM data (temp in Kelvin)

UVic_sat=sat(10:55,:,:);
clear sat

%% extract 1960-2005 from RECON data (temp in Kelvin)
 %(time,lat,lon)
    RECON_tmp2 = zeros(46,180,360);
    RECON_sat = zeros(46,180,360);    
for n=1:64800;
    n_int=int32(n);
    i=int32(RECON_tmp(n_int,1));
    j=int32(RECON_tmp(n_int,2));
    RECON_tmp2(:,j,i) = RECON_tmp(n_int,3:48);
end

%%

for n=1:46;
    tmp(:,:)=RECON_tmp2(n,:,:);
    tmp=flipud(tmp);
    RECON_sat(n,:,:)=tmp;    
end
%%

RECON_arr_dim = size(RECON_sat);

ocnmsk=(RECON_sat>0.);  %ocnmsk = 1 over land/ice shelf
for year=1:46    
  for j=2:(RECON_arr_dim(2))
    for i=1:RECON_arr_dim(3)
      if RECON_sat(year,j,i) < 0.
        RECON_sat(year,j,i) = RECON_sat(year,j-1,i);
      end
    end
  end
end


%% interpolate to UVic grid

UVic_arr_dim = size(UVic_sat);
RECON_arr_dim = size(RECON_sat);

xi=linspace(1,RECON_arr_dim(3),100);
yi=linspace(1,RECON_arr_dim(2),100);
xo=1:1:RECON_arr_dim(3);
yo=1:1:RECON_arr_dim(2);

[X,Y] = meshgrid(xo,yo);
[XI,YI] = meshgrid(xi,yi);

RECON_sat_interp = zeros(46,UVic_arr_dim(3),UVic_arr_dim(2));
ocnmsk_interp    = zeros(46,UVic_arr_dim(3),UVic_arr_dim(2));
for year = 1:46;    

    interp_in_array(:,:)=RECON_sat(year,:,:);  
    RECON_sat_interp(year,:,:)=interp2(X,Y,interp_in_array,XI,YI,'linear');

    interp_in_array(:,:)=ocnmsk(year,:,:);      
    ocnmsk_interp(year,:,:)=interp2(X,Y,interp_in_array,XI,YI,'linear');
end

%% create differences between UVic and RECON timeslices

UVic_minus_RECON = zeros(46,100,100);

for year = 1:46;
  RECON_slice(:,:)=RECON_sat_interp(year,:,:);
  ocnmsk_slice(:,:)=ocnmsk_interp(year,:,:);
  UVic_slice(:,:)=UVic_sat(year,:,:);
 
  for i=1:100;
    for j=1:100;
        if ocnmsk_slice(i,j) == 1;
        UVic_minus_RECON(year,i,j) = UVic_slice(i,j)-RECON_slice(i,j);
        end
    end
  end
  
end

%%
close all

lat_array=linspace(-89.1,89.1);lon_array=linspace(1.8,358.2);
[Plg,Plt]=meshgrid(lon_array,lat_array);

for i=1:100;
  for j=1:100;
      if ocnmsk_slice(i,j) == 0;
      RECON_slice(i,j) = 283.15;
      end
  end
end

hold on
m_proj('stereographic','longitude', 180,'latitude',-90,'radius', 30);
hold on
m_pcolor(Plg,Plt,(RECON_slice-273.15));shading flat;
m_coast('color',[0 0 0]);
m_grid('xticklabels',[],'yticklabels',[]);

h=colorbar('h');
title('RECON Annual Average SAT, 2005 (C)','FontSize', 20);

print -djpeg90 ../../RECON_SH_sat.jpg
hold off



%% movie

arr_sz = size(UVic_minus_RECON);
lat_array=linspace(-89.1,89.1);lon_array=linspace(1.8,358.2);
[Plg,Plt]=meshgrid(lon_array,lat_array);
%m_proj('miller','longitude', 180);
m_proj('stereographic','longitude', 180,'latitude',-90,'radius', 30);

for n=1:arr_sz(1);
  plot_var(:,:)=UVic_minus_RECON(n,:,:); 
  
  m_pcolor(Plg,Plt,plot_var);shading flat;
  m_coast('color',[0 0 0]);
  %m_grid('xaxis','bottom');
  m_grid('xticklabels',[],'yticklabels',[]);
  h=colorbar('h');
  
  year = 'SAT diff: UVic ESCM - RECON, year= ';
  date_text = num2str(1959+n);
  tag = strcat(year,date_text);
  title(tag, 'Fontsize', 20)
  
  M(n) = getframe(gcf);

end
movie(M)
movie2avi(M,'RECON-UVic_SH.avi')

close all

