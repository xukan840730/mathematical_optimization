function[dtdec] = temp_trend(sat)

%Determine dT/decade for all grid points
%Requires: -arbitrarily sized input array ('sat') 2D annual average surface temperature (over time)

dim = size(sat);
sdat = 1;
slat = 1;
slon = 1;
edat = dim(1);
elat = dim(2);
elon = dim(3);
dtdec = zeros(dim(2),dim(3));
temp = 1:dim(1);

for i = slat:elat
  for j = slon:elon
    %get data for one point over time
    data = sat(sdat:edat,i,j);
    %get coefficients of best-fit line (using least-squares) through
    %temporal data at point    
    P = polyfit(temp,rot90(data),1);
    dtdec(i,j) = P(1).*10;
    %dt/decade, multiply P(1) accordingly.
  end
end
