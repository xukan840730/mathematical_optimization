function[year_avg_array] = mon2yearavg(mon_avg_array)

%returns annual averages of an input file of monthly averages (number of
%months must divide into 12, i.e. only full years accepted).

dim = size(mon_avg_array);
for y=1:(dim(1)/12);
  s = (y-1)*12. + 1;
  e = y*12;
  for i=1:dim(2);
    for j=1:dim(3);
      year_avg_array(y,i,j) = sum(mon_avg_array(s:e,i,j))/12.;
    end
  end
end
