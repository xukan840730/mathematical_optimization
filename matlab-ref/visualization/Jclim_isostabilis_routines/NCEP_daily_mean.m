%% routine calculates NCEP-UVic bias for each day of daily average SAT, over
%% 1971-2001.
%1. Load NCEP data, interpolate to UVic grid, obtain daily average
%2. Load UVic data
%3. Subtract NCEP-UVic data, for each day at each grid point to obtain bias
%field.
%4. Write bias field to netcdf

clear all

%% 
%load consecutive .nc files, accumulate

SAT_accum = zeros(365,73,144);

count = 0;
for file = 1971:2001;
  count = count + 1
  fstring=int2str(file);
  fname = strcat('air.sig995.',fstring,'.nc');
  ncload(fname,'air');
  %unpack integer values using NCEP formula
  air = air*0.01;
  air = air+512.81;
  SAT_accum = SAT_accum + (air(1:365,:,:));
end

SAT_accum(:,:,:) = (SAT_accum(:,:,:)/count) - 273.15;

clear count file fname fstring

%%
%interpolate SAT_accum to the UVic grid.

NCEP_arr_dim = size(SAT_accum);          %(time,latitude,longitude)

xi=linspace(1,size(SAT_accum,3),100);
yi=linspace(1,size(SAT_accum,2),100);
xo=1:1:size(SAT_accum,3);
yo=1:1:size(SAT_accum,2);

[X,Y] = meshgrid(xo,yo);
[XI,YI] = meshgrid(xi,yi);

SAT_accum_intrp = zeros(size(SAT_accum,1),100,100);
for n=1:size(SAT_accum,1);
  temp(:,:) = SAT_accum(n,:,:);
  temp = flipud(temp);
  SAT_accum_intrp(n,:,:)=interp2(X,Y,temp,XI,YI);
end
clear temp X Y XI YI xi yi xo yo NCEP_arr_dim SAT_accum;

%%
%load UVic SAT data

ncload ('UVic_SAT.nc','output_field');
UVic_SAT = output_field;
clear output_field

%add lapse rate (5 C/km, UVic standard)

ncload('elev.nc','elev')
lapse_adjust = elev(2:101,2:101)*5.e-3;
clear elev

%apply lapse rate adjustment to SAT data.  Use UVic lapse rate (5K/km).
clear temp
for n = 1:size(UVic_SAT,1);
   temp(:,:) = UVic_SAT(n,:,:);
   UVic_SAT(n,:,:) = temp - lapse_adjust;
end
clear temp lapse_adjust n

%%
%subtract: NCEP - UVic

SAT_bias = SAT_accum_intrp - UVic_SAT;

%%
%write to NetCDF

array = SAT_bias;
fname = 'NCEP_UVic_bias.nc';
netcdf_writer





