

%% Get monthly NCEP fields, interpolated to UVic grid
%% Load monthly mean temperature

clear all

%set years over which to obtain a base climate to add anomalies to.
start_year = 1971;
num_years = 30.;

%% 
% Load NCEP data.
% Interpolate to UVic grid.
% Get monthly averages from start_year, for num_years.

[NCEP_mon_mean_sat] = NCEP_UVic_interp(start_year,num_years);

%%
% Load UVic data.
% Assign base climate to calculate anomalies from.

[UVic_mon_mean_sat,UVic_mon_mean_sat_base,kmt] = load_UVic(start_year,num_years);

%%
%Calculate climate anomalies.
%Scan base+anomaly climate for significant periods above 0C each year.
%Accumulate yearly PDDs.
%Assume piecewise linear trends from mid-month to mid-month).
a='profiling'
profile clear
profile on
[bsat,isomap,pdd] = iso_calc(NCEP_mon_mean_sat,UVic_mon_mean_sat,UVic_mon_mean_sat_base);
profile viewer

%%
%Save .mat file with pertinent variables

save isostabilis

%%
%Plot 'er up.
test_param=0.5;  %test value
mask=zeros(size(isomap));
for t=1:size(isomap,1);
    array_slice(:,:)=isomap(t,:,:);
    %more than melt period over land    
    i = find(kmt < 0.5 & array_slice > test_param);
    mask(t,i) = 1;   
    %more than melt period over ocean
    i = find(kmt > 0.5 & array_slice > test_param);
    mask(t,i) = 2;  
    %less than melt period over land
    i = find(kmt < 0.5 & array_slice < test_param);
    mask(t,i) = 3;
    %less than melt period over ocean
    i = find(kmt > 0.5 & array_slice < test_param);
    mask(t,i) = 4;
end

%%
count = 0;
for t=1:5:450
    count = count+1;
    movie_array(count,:,:) = mask(t,:,:);
end    
movie_maker(movie_array,'cum_time');
    