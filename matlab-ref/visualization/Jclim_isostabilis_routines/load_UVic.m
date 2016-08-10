function[UVic_mon_mean_sat,UVic_mon_mean_sat_base,kmt] = load_UVic(start_year,num_years)
%% load uvic data

%Define size of monthly mean array. 
%Array size is hard-coded as: 100x100x(500 yearsx12 months/year).

UVic_mon_mean_sat = zeros(6000,100,100);

ncload('mon_sat.01950.01.01.nc','mon_mean_t')
UVic_mon_mean_sat(1:1200,:,:) = mon_mean_t(1:1200,:,:);
clear mon_mean_t

ncload('mon_sat.02050.01.01.nc','mon_mean_t')
UVic_mon_mean_sat(1201:2400,:,:) = mon_mean_t(1:1200,:,:);
clear mon_mean_t

ncload('mon_sat.02150.01.01.nc','mon_mean_t')
UVic_mon_mean_sat(2401:3600,:,:) = mon_mean_t(1:1200,:,:);
clear mon_mean_t

ncload('mon_sat.02250.01.01.nc','mon_mean_t')
UVic_mon_mean_sat(3601:4800,:,:) = mon_mean_t(1:1200,:,:);
clear mon_mean_t

ncload('mon_sat.02350.01.01.nc','mon_mean_t')
UVic_mon_mean_sat(4801:6000,:,:) = mon_mean_t(1:1200,:,:);
clear mon_mean_t

ncload('../kmt.nc','kmt')
temp = kmt;
clear kmt
kmt = temp(2:101,2:101);

ncload('../elev.nc','elev')
lapse_adjust = elev(2:101,2:101)*5.e-3;
clear elev

%%

%apply lapse rate adjustment to SAT data.  Use UVic lapse rate (5K/km).
clear temp
for n = 1:size(UVic_mon_mean_sat,1);
   temp(:,:) = UVic_mon_mean_sat(n,:,:);
   UVic_mon_mean_sat(n,:,:) = temp - lapse_adjust;
end
clear temp

%% get base climate

%accumulate monthly means.  Specify start year and number of years to
%accumulate monthly averages.

UVic_mon_mean_sat_base = zeros(12,100,100);

start = ((start_year-1950)*12)+1;
finish= start + (num_years*12-1);

count = 1;
for n = start:finish 
  UVic_mon_mean_sat_base(count,:,:) = UVic_mon_mean_sat_base(count,:,:) + UVic_mon_mean_sat(n,:,:);
  count = count + 1;
  if (count == 13)
    count = 1;
  end
end

%average monthly means
UVic_mon_mean_sat_base(:,:,:) = UVic_mon_mean_sat_base(:,:,:)/num_years;