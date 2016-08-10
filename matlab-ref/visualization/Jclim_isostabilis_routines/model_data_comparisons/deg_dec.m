%% Clear workspace
clear all;

%% Uvic trend
ncload('SH2000-2100sat.nc','sat');

dtdecUVic = temp_trend(sat);

clear sat;

a = 'UVic trend done'

%% CCCMA trend

ncload('pcmdi.ipcc4.cccma_cgcm3_1.sresa2.run1.monthly.tas_a1_sresa2_1_cgcm3.1_t47_2001_2100.nc','tas');
dim = size(tas);
sat = zeros(dim(1)/12,dim(2),dim(3));

sat = mon2yearavg(tas);
dtdec = temp_trend(sat);
dtdeccccma = interp2uvic(dtdec);

clear P dtdec sat dim data temp tas e s i j y;

a = 'CCCMA trend done'

%% GDFL trend

ncload('pcmdi.ipcc4.gfdl_cm2_0.sresa2.run1.monthly.tas_A1.200101-210012.nc','tas');

dim = size(tas);
sat = zeros(dim(1)/12,dim(2),dim(3));

sat = mon2yearavg(tas);
dtdec = temp_trend(sat);
dtdecgdfl = interp2uvic(dtdec);

clear dtdec sat dim tas;

a = 'GDFL trend done'

%% ECHAM5 trend

ncload('pcmdi.ipcc4.mpi_echam5.sresa2.run1.monthly.tas_A1.nc','tas');

dim = size(tas);
sat = zeros(dim(1)/12,dim(2),dim(3));

sat = mon2yearavg(tas);
dtdec = temp_trend(sat);
dtdececham5 = interp2uvic(dtdec);

clear dtdec sat dim tas;

a = 'ECHAM5 trend done'

%% HADCM3 trend

ncload('pcmdi.ipcc4.ukmo_hadcm3.sresa2.run1.monthly.tas_A1.nc','tas');

dim = size(tas);
sat = zeros(dim(1)/12,dim(2),dim(3));

sat = mon2yearavg(tas);
dtdec = temp_trend(sat);
dtdechadcm3 = interp2uvic(dtdec);

clear dtdec sat dim tas;

a = 'HADCM3 trend done'

%% CCSM3 trend

ncload('pcmdi.ipcc4.ncar_ccsm3_0.sresa2.run1.monthly.tas_A1.nc','tas');

dim = size(tas);
sat = zeros(dim(1)/12,dim(2),dim(3));

sat = mon2yearavg(tas);
dtdec = temp_trend(sat);
dtdecccsm3 = interp2uvic(dtdec);

clear dtdec sat dim tas;

a = 'CCSM3 trend done'

%% ensemble trend

dtdecens = (dtdecgdfl + dtdeccccma + dtdececham5 + dtdechadcm3 + dtdecccsm3)/5;
uvic_min_ens = dtdecUVic - dtdecens;

%% Plot

% figure
% subplot(2,1,1);
% [cs,h]=contourf(uvic_min_ens);
% clabel(cs,h);
% colorbar
% subplot(2,1,2);
% [cs,h]=contourf(dtdecUVic);
% clabel(cs,h);
% colorbar

figure
%GL 
lat_array=linspace(-89.1,89.1);lon_array=linspace(0,360);
[Plg,Plt]=meshgrid(lon_array,lat_array);
m_proj('stereographic','longitude', 300,'latitude',80,'radius', 30);

%Ant
% lat_array=linspace(-89.1,89.1);lon_array=linspace(0,360);
% [Plg,Plt]=meshgrid(lon_array,lat_array);
% m_proj('stereographic','longitude', 0,'latitude',-90,'radius', 30);


m_contourf(Plg,Plt,uvic_min_ens);
m_coast;
colorbar
%m_grid;

