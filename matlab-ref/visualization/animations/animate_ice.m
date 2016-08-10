function animate_ice(hin,hsin,fname)

%% make animation of evolving ice geometries

%% interpolate in time
%interpolation factor
ifac=1;
%original array length in time
arrdimin=size(hin);
arrdim=arrdimin;
%expand array length in time
arrdim(3)=arrdimin(3)*ifac;

%set time to time dimension length
time=arrdim(3)

h=zeros(arrdim);
hs=zeros(arrdim);
size(h)
wt=linspace(1.,1./ifac,ifac)
for no=1:arrdimin(3)-1
  %set start point minus 1 for filling expanded array
   n_expstrt=(no-1)*ifac;
   for ni=1:ifac;
      %get time indice of expanded matrix to fill
      ind=n_expstrt+ni;
      %fill expanded array with interpolated values
      h(:,:,ind)=hin(:,:,no)*wt(ni)+hin(:,:,no+1)*(1.-wt(ni));
      hs(:,:,ind)=hsin(:,:,no)*wt(ni)+hsin(:,:,no+1)*(1.-wt(ni));      
    end
end
%% generate animation

%set all initial locations where elev is lower than sealev to 0.
hs=max(0.,hs);

%basic pan using camera and light pan/zoom arrays
azimuth=linspace(85,95,time);
elevation=linspace(55,65,time);
zoom=linspace(.1,.1,time);
xlight=linspace(40,80,time);
dx=linspace(-1,1,time);
%get maximum and minimum hs values
minh=min(min(min(hs)));
maxh=max(max(max(hs)));

%set number of color contours to total range in elevation
ccontours=floor(maxh-minh);
cmap(:,1)=linspace(200,255,ccontours);
cmap(:,2)=linspace(255,255,ccontours);
cmap(:,3)=linspace(255,255,ccontours);
cmap(1,1:3)=[6 113 148]; %ocean blue
cmap(2,1:3)=[158 128 110]; %land brown
cmap=cmap/255.;
%% Set figure size
close all
set(0,'Units','pixels');
scnsize=get(0,'Screensize')
fig1=figure;
position=get(fig1,'Position')
set(fig1,'Position',scnsize*.75)
%% Generate movie frames
for n=1:time;
  temph(:,:)=(h(:,:,n));
  temphs(:,:)=(hs(:,:,n));  
  temph=rot90(temph');
  temphs=rot90(temphs'); 
  
  %set array that will be used to determine colors
  temphscolor=temphs;
  
  %find bare land and drop to near bottom of scale
  clear i
  i=find(temph<=0.1);
  temphscolor(i)=minh+1.5; 
  %then drop ocean off bottom of scale
  clear i
  i=find(temphs==0.);
%   temphscolor(i)=minh-1; 
%   zexag=5./1000.;
%   %figh=surf(temphs*zexag,temphscolor);
%   
%   caxis([minh maxh]);
%   colormap(cmap);
%   %camlookat(figh)
%   %camdolly(dx(n),dx(n),dx(n),'targetmode')
%   %campan(dx(n)*90,dx(n)*90)
%   camproj('perspective')
% 
%   view(azimuth(n),elevation(n))
%   %camzoom(zoom(n));
%   
%   light('Position',[xlight(n),-1,0.3],'Style','Infinite')
%   lighting phong
%   shading interp;
%   grid off
%   axis off
%   axis equal tight;
%   clear temph temphs
  
temphs(i)=nan;
pcolor(temphs),shading flat
axis square
  command=strcat('print -dpng movies/',num2str(n));
  eval(command);
end

!animate *.png