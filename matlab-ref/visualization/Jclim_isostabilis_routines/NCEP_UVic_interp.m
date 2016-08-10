function[NCEP_mon_mean_sat] = NCEP_UVic_interp(start_year,num_years)

%% Function to interpolate NCEP monthly mean fields to UVic grid, get monthly
%% averages over a number of years.

%% interpolate NCEP to UVic grid

ncload('../air.mon.mean.nc','air');
NCEP_arr_dim = size(air);          %(time,latitude,longitude)

xi=linspace(1,NCEP_arr_dim(3),100);
yi=linspace(1,NCEP_arr_dim(2),100);
xo=1:1:NCEP_arr_dim(3);
yo=1:1:NCEP_arr_dim(2);

[X,Y] = meshgrid(xo,yo);
[XI,YI] = meshgrid(xi,yi);

NCEP_sat_interp = zeros(NCEP_arr_dim(1),100,100);
for n=1:NCEP_arr_dim(1);
  temp(:,:) = air(n,:,:);
  NCEP_sat_interp(n,:,:)=interp2(X,Y,temp,XI,YI);
end
clear temp;

%%  accumulate monthly averages

%preallocate multi-year avereaged monthly mean array
NCEP_mon_mean_sat = zeros(12,100,100);

%accumulate monthly means.  Specify start year and number of years to
%accumulate monthly averages.

start = ((start_year-1948)*12)+1;
finish= start + (num_years*12-1);

count = 1;
for n = start:finish 
  NCEP_mon_mean_sat(count,:,:) = NCEP_mon_mean_sat(count,:,:) + NCEP_sat_interp(n,:,:);
  count = count + 1;
  if (count == 13)
    count = 1;
  end
end

%flip NCEP array (to correspond to UVIC array) and get average monthly
%means

for n=1:12
  temp(:,:) = NCEP_mon_mean_sat(n,:,:);
  temp = flipud(temp);
  NCEP_mon_mean_sat(n,:,:) = temp(:,:)/num_years;
end

clear temp