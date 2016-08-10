%%
clear all
cd /Network/Servers/cl00.seos.uvic.ca/Volumes/Casa/Users2/jer/Desktop/School/Analysis/isostabilis/'Isostabilis output'/Inline_NCEP_output/SHwind/

%%
% %load commitment run data... must remove every 99th year, and method depends on file input
 ncload('iscatcom.nc');

 posddcom = zeros(999,100,100);
 posddcom(1:9) = pdd(1:9);
 for c = 0:9 
   posddcom((c*100+1-c):(c*100+99-c),:,:) = pdd((c*100+12):((c+1)*100+10),:,:);  
 end

%%
%Load entire SHwind record: 1850-2850.  Trim 1 year/century (PDDs can't be
%summed correctly for first year of run after restarted).
 ncload('iscat2.nc');
% iso15 = zeros(990,100,100);
% iso30 = zeros(990,100,100);
% iso45 = zeros(990,100,100);
 posdd = zeros(990,100,100);
 for c = 0:9
%   iso15((c*100+1-c):(c*100+99-c),:,:) = isostab15((c*100+2):((c+1)*100),:,:);
%   iso30((c*100+1-c):(c*100+99-c),:,:) = isostab30((c*100+2):((c+1)*100),:,:);  
%   iso45((c*100+1-c):(c*100+99-c),:,:) = isostab45((c*100+2):((c+1)*100),:,:);  
   posdd((c*100+1-c):(c*100+99-c),:,:) = pdd((c*100+2):((c+1)*100),:,:);  
 end

clear isostab15 isostab30 isostab45 pdd xt yt

%%
ncload('tmsk.nc','tmsk');
temp = tmsk;
kmt(1:100,1:100) = temp(2:101,2:101);
clear temp

%%

array_test=posddcom;  %array to apply test to
test_param=200;  %test value
array_dim = size(array_test);
mask=zeros(array_dim(1), array_dim(2), array_dim(3));

for t=1:array_dim(1);
    array_slice(:,:)=array_test(t,:,:);
    %more than melt period over land    
    i = find(kmt < 0.5 & array_slice > test_param);
    mask(t,i) = 3;   
    %more than melt period over ocean
    i = find(kmt > 0.5 & array_slice > test_param);
    mask(t,i) = 3;  
    %less than melt period over land
    i = find(kmt < 0.5 & array_slice < test_param);
    mask(t,i) = 1;
    %less than melt period over ocean
    i = find(kmt > 0.5 & array_slice < test_param);
    mask(t,i) = 1;
end


%%
slice = [1 150 250 990];
for time=1:4;
  if time==1
    date = '1850'
  elseif time==2  
    date = '2000' 
  elseif time==3 
    date = '2100' 
  elseif time==4 
    date = '2850' 
  end
  
  for sheet=1:3;  
    if (sheet==1)
      lat_array=linspace(-90,90);lon_array=linspace(0,360);  
      [Plg,Plt]=meshgrid(lon_array,lat_array);                   %define array grid
      m_proj('equidistant cylindrical','lon',180,'lat',[-89.5 89.5]);
    elseif (sheet==2)
      lat_array=linspace(-89.1,89.1);lon_array=linspace(0,360);
      [Plg,Plt]=meshgrid(lon_array,lat_array);
      m_proj('stereographic','longitude', 0,'latitude',-90,'radius', 30);
    elseif (sheet==3)
      lat_array=linspace(-89.1,89.1);lon_array=linspace(0,360);
      [Plg,Plt]=meshgrid(lon_array,lat_array);
      m_proj('stereographic','longitude', 300,'latitude',80,'radius', 30);
    end

    figure                                         
    plot_var(:,:)=mask(slice(time),:,:);   
    m_pcolor(Plg,Plt,plot_var);shading flat;
    m_coast('color','k','linewidth',1)

    caxis([1 4]);
    %Ward Hunt/Ayles
    [X,Y]=m_ll2xy((-75+360),83.1);
    a(1)=line(X,Y);

    %larsen A/B/Inlet
    [X,Y]=m_ll2xy((-60+360),-65);
    a(2)=line(X,Y);
    
    %Wilkins
    [X,Y]=m_ll2xy((-73+360),-70.25);
    a(3)=line(X,Y);
  
    %Ross
    [X,Y]=m_ll2xy((-170+360),-82);
    a(4)=line(X,Y);
  
    %Ronne-Filchner
    [X,Y]=m_ll2xy((-58+360),-82);
    a(5)=line(X,Y);
    
    %GIS
    [X,Y]=m_ll2xy((-34+360),78);
    a(6)=line(X,Y); 

    %Amery
    %[X,Y]=m_ll2xy(70,-70);
    %a(6)=line(X,Y);    
  


    set(a(1:3),'marker','pentagram','markersize',12,'markerfacecolor','r','markeredgecolor','r')  
    set(a(4:6),'marker','diamond','markersize',10,'markerfacecolor','r','markeredgecolor','r')  
    
    if (sheet==1)
      m_grid('xtick',6,'ytick',6);
      print('-depsc2',strcat('figures/Global',date))
    elseif (sheet==2)
      m_grid('XaxisLocation', 'top', 'xtick',12,'tickdir','out','ytick',[-70 -80],'linest','-');
      print('-depsc2',strcat('figures/Antarctica',date))
    elseif (sheet==3)
      m_grid('XaxisLocation', 'bottom', 'xtick',12,'tickdir','out','YaxisLocation', 'middle','ytick',[60 70 80],'linest','-');
      print('-depsc2',strcat('figures/Greenland',date))
    end    
  end
end

close all



