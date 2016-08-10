
%% Antarctic shelves that have collapsed
clear all
close all


figure
m_proj('stereographic','lon',310,'lat',-65,'rad',15);
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_grid;

%note: coordinate pair: [lon lat]

%larsen A/B/Inlet
[X,Y]=m_ll2xy(-61,-65);
a(1)=line(X,Y);
t(1)=text(X,Y,'Larsen A/B/Inlet, 1995-2002','vertical','top');

%Wordie
[X,Y]=m_ll2xy(-67.75,-69.25);
a(2)=line(X,Y);
t(2)=text(X,Y,'Wordie, 1980s','vertical','top');

%Muller/Jones
%cant find coordinates?
%1970-2003

%Wilkins
[X,Y]=m_ll2xy(-73,-70.25);
a(3)=line(X,Y);
t(3)=text(X,Y,'Wilkins, 2008','vertical','top');

set(a,'marker','square','markersize',10,'markerfacecolor','r','markeredgecolor','r')
set(t,'fontsize',15)
print ('-djpeg100', ('Antarctic_collapse'))
 
%% Arctic shelves that have collapsed
clear all
close all

figure
m_proj('stereographic','lon',280,'lat',85,'rad',5);
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_grid;

%Ayles
[X,Y]=m_ll2xy(-77.558,83.025);
a(1)=line(X,Y,'marker','square','markersize',10);
t(1)=text(X,Y,'Ayles, 2005','vertical','top');

%Ward Hunt
[X,Y]=m_ll2xy(-75,83.1);
a(2)=line(X,Y,'marker','square','markersize',10);
t(2)=text(X,Y,'Ward Hunt, 2002','vertical','bottom');

set(a,'marker','square','markersize',10,'markerfacecolor','r','markeredgecolor','r')
set(t,'fontsize',15)
print ('-djpeg100', ('Arctic_collapse'))

