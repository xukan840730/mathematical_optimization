%% map 1,2,3 = global, antarctica, greenland
map = 3
year = 2500 %%Note: make sure to change plot_var below to correct slice of data

if map==1
mname = 'Globe';
lat_array=linspace(-90,90);lon_array=linspace(0,360);  
[Plg,Plt]=meshgrid(lon_array,lat_array);                   %define array grid
m_proj('equidistant cylindrical','lon',180,'lat',[-89.5 89.5]); 
end
if map==2
mname = 'Antarctica';   
lat_array=linspace(-89.1,89.1);lon_array=linspace(0,360);
[Plg,Plt]=meshgrid(lon_array,lat_array);
m_proj('stereographic','longitude', 0,'latitude',-90,'radius', 30);
end
if map==3
mname = 'Greenland';   
lat_array=linspace(-89.1,89.1);lon_array=linspace(0,360);
[Plg,Plt]=meshgrid(lon_array,lat_array);
m_proj('stereographic','longitude', 300,'latitude',75,'radius', 20);
end
close all
figure

plot_var(:,:)=mask(11,:,:); 
m_pcolor(Plg,Plt,plot_var), shading flat


%% plot recent collapses as red squares
%Ayles
[X,Y]=m_ll2xy((-77.558+360),83.025);
a(1)=line(X,Y);
%Ward Hunt
[X,Y]=m_ll2xy((-75+360),83.1);
a(2)=line(X,Y);
%larsen A/B/Inlet
[X,Y]=m_ll2xy((-60+360),-65);
a(3)=line(X,Y);    
%Wilkins
[X,Y]=m_ll2xy((-73+360),-70.25);
a(4)=line(X,Y);
%set gridlines, color axis  
set(a,'marker','square','markersize',15,'markerfacecolor','r','markeredgecolor','r')  
m_grid('xticklabels',[],'yticklabels',[]);
caxis([1 4]);

fname = strcat(mname,num2str(year));
print('-depsc2', '-r300', '-tiff',fname) 
 
  
  