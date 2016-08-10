%% workspace to test ocean melt parameterization
clear all
cd /Network/Servers/cl00.seos.uvic.ca/Volumes/Casa/Users2/jer/Desktop/School/Analysis/oceanmelt

nx=100;
ny=100;

%array dimensions are: depth, latitude, longitude

ncload('lgm.nc','O_temp');
lgmtemp=O_temp;
clear O_temp;

ncload('modern.nc','O_temp');
modtemp=O_temp;
clear O_temp;

i=find(lgmtemp > 1.e20);
lgmtemp(i) = nan;
modtemp(i) = nan;

ncload('ifrac.nc','output_field');
ifrac=output_field;
clear output_field;

ncload('kmt.nc','kmt');
ncload('tavg.nc','yt','zt');

difftemp = modtemp-lgmtemp;

mn(:,:) =isnan(difftemp(1,:,:));
i=find((ifrac > 0.5) & (mn < 0.5));

%set up volume matrix (cubic meters)
dx=400238.893464630;
dy=200119.446732315;
for j=1:ny;
  for k=1:19;
    volume(k,j,1:nx) = dx*cos(yt(j)*pi/180.)*dy*zt(k);
  end 
end

mask=zeros(nx,ny);

mind=4
maxd=9

%% RIS region

tvdiff=0.;
tvLGM=0.;
tvmod=0.;
vol=0.;

for i=47:61;
j=1;
  %search northwards for ice edge
  while (ifrac(j,i) > 0.5);
    j=j+1;
  end;
  %while between mind and maxd ocean levels
  while (kmt(j,i) <= maxd);
    if (kmt(j,i) >= mind);
      mask(j,i)=1;
      for k=mind:kmt(j,i);
        %accumulate T*V value (Cm^3)
        tvdiff=tvdiff+difftemp(k,j,i)*volume(k,j,i);
        tvLGM=tvLGM+lgmtemp(k,j,i)*volume(k,j,i);
        tvmod=tvmod+modtemp(k,j,i)*volume(k,j,i);        
        vol=vol+volume(k,j,i);
      end
    end
    j=j+1;
  end
end
RISlgm=tvLGM/vol
RISmod=tvmod/vol
RISdt=tvdiff/vol

%% RFIS region

tv=0.;
vol=0.;
for i=84:92;
j=1;
  %search northwards for ice edge
  while (ifrac(j,i) > 0.5);
    j=j+1;
  end
  %while between mind and maxd ocean levels
  while (kmt(j,i) <= maxd);
    if (kmt(j,i) >= mind);
      mask(j,i)=2;
      for k=mind:kmt(j,i);
        %accumulate T*V value (C*m^3)
        tv=tv+difftemp(k,j,i)*volume(k,j,i);
        vol=vol+volume(k,j,i);
      end
    end
    j=j+1;
  end
end

RFISdt=tv/vol    

%%
pcolor(mask);

%%
interior=[0.0 0.1 2.0]
expshelf=[0.0 5.0 10.0]
opnocean=[2.0 5.0 10.0]
dT=[0 1.5 3]

pinterior=polyfit(dT,interior,2)
pexpshelf=polyfit(dT,expshelf,2)
popnocean=polyfit(dT,opnocean,2)

test=(0:0.1:2)
f1=polyval(pinterior,test);
f2=polyval(pexpshelf,test);
f3=polyval(popnocean,test);
close all
hold on
plot(test,f1, 'r')
plot(test,f2, 'k')
plot(test,f3, 'b')
hold off

%%

dT=[0 1.5 3]

% i1=[0.0 0.1 2.0]
% i2=[0.0 0.1 4.0]
% i3=[0.0 0.1 6.0]

% i1=[0.0 5.0 10.0]
% i2=[0.0 5.0 15.0]
% i3=[0.0 5.0 20.0]

i1=[2.0 5.0 10.0]
i2=[2.0 5.0 15.0]
i3=[2.0 5.0 20.0]

p1=polyfit(dT,i1,2)
p2=polyfit(dT,i2,2)
p3=polyfit(dT,i3,2)

test=(0:0.1:2)
f1=polyval(p1,test);
f2=polyval(p2,test);
f3=polyval(p3,test);
close all
hold on
plot(test,f1, 'r')
plot(test,f2, 'k')
plot(test,f3, 'b')
hold off




