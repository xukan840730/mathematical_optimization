%% determine melt extent
%requires cum_time_mapper to be previously called

%get yt from tavg file.
ncload('tavg.nc','yt');

%%
%set up area matrix (cubic meters)
dx=400238.893464630;
dy=200119.446732315;
for j=1:100;
    area(j,1:100) = dx*cos(yt(j)*pi/180.)*dy/1000000.;
end
temp(1:100,1:100)=tmsk(2:101,2:101);
meltextent(1:990)=0.;
for i=1:17
  for j=1:100
    if (temp(i,j) < 0.5)      
      for n=1:990
        if (posdd(n,i,j) > 0.)
          meltextent(n)=meltextent(n)+area(i,j);
        end
      end
    end
  end
end

%%
year=[1850:1:2839];
avg=0.
i=0.
for n=138:158;
  avg=avg+meltextent(n);
  i=i+1.;
end
avg=avg/i
for n=1:990
  if (meltextent(n) >= 297500)
    year(n)
    meltextent(n)
    break
  end
end

%%

plot(year,meltextent)



%% plot trend in PDD days over period 1965:2005 to compare against Vaughan (2006)
%note that length of timeseries corresponds to length of observational
%dataset as quoted in Vaughan (2006).
close all
figure
hold on

%RECENT OBSERVED:
%plot(posdd(140:155,15,85), 'r', 'linewidth', 5) %Esperanza
%plot(posdd(100:155,15,83), 'b', 'linewidth', 5) %Faraday
%plot(posdd(115:155,16,84), 'k', 'linewidth', 5) %Bellingshausen
%plot(posdd(130:155,13,81), 'g', 'linewidth', 5) %Rothera
%plot(posdd(130:155,12,100), 'g', 'linewidth', 5) %Fimbul
%plot(posdd(130:155,12,4), 'g', 'linewidth', 5) %Amery
%plot(posdd(130:155,96,87), 'g', 'linewidth', 5) %Ayles

%NEAR FUTURE:
% plot(posdd(156:196,15,85), 'r', 'linewidth', 5) %Esperanza
% plot(posdd(156:196,15,83), 'b', 'linewidth', 5) %Faraday
% plot(posdd(156:196,16,84), 'k', 'linewidth', 5) %Bellingshausen
% plot(posdd(156:196,13,81), 'g', 'linewidth', 5) %Rothera
% plot(posdd(156:196,12,100), 'g', 'linewidth', 5) %Fimbul
%  plot(posdd(156:196,12,4), 'g', 'linewidth', 5) %Amery
plot(posdd(156:196,96,87), 'g', 'linewidth', 5) %Ayles

hold off




















