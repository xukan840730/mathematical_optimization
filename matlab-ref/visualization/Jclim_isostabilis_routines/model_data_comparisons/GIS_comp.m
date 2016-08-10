%% Open Box et. al. (2004) GIS surface mass balance datasets for 1991-2000 and get average values.
cd /Network/Servers/cl00.seos.uvic.ca/Volumes/Casa/Users2/jer/Desktop/School/Data/BoxGIS
clear all
 
 smb=zeros(55,101);
 tmp=zeros(55,101);
 count=0.
 for year=1991:2000;
   count=count+1.
   txtyr=num2str(year);
   string=[txtyr,'_smb_101x55_24km.asc']
   input=textread(string,'%f');
   for row=0:100;
     n=row*55+1;
     tmp(:,row+1)=input(n:n+54);
   end
   clear input
   smb=smb+tmp;

 end

 smb=smb/count;
   
 smb=fliplr(smb);
 smb=rot90(smb);
 
 %% Open UVic-produced GIS mass balance
 ncload('budgsnow_Greenland.nc');
 avg=mean(output)
 i=find(smb==0.);
 smb(i)=nan;
 i=find(output==0.);
 output(i)=nan;
 
 %plot Box and Uvic side by side with same color scale
 close all
 hold on
 figure
 pcolor(smb), shading flat;
 caxis([-3. 1.9]);
 colorbar; 
 print -djpeg90 
 cax=caxis;
 figure
 pcolor(output), shading flat;
 caxis(cax);
 colorbar;
 print -djpeg90  
 hold off
 
 
