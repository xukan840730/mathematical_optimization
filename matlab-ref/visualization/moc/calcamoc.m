%plot barotropic stream fucntion for Atlantic
ncload('~/srv_ccrc/data01/z3263455/escm27/pdetop5/100free/tavg.03001.01.01.nc','zw');

%calc dz: zdiff
zw=[0 zw'];
for z = 1:19
    zdiff(z) = zw(z+1)-zw(z);
end

ncload('~/srv_ccrc/data01/z3263455/escm27/pdetop5/100free/tavg.03001.01.01.nc','v');
ncload('~/srv_ccrc/data01/z3263455/escm27/pdetop5/100free/tavg.03001.01.01.nc','time');
ncload('~/srv_ccrc/data01/z3263455/escm27/pdetop5/100free/tavg.03001.01.01.nc','yu');
ncload('~/srv_ccrc/data01/z3263455/escm27/pdetop5/100free/tavg.03001.01.01.nc','mskhr');
ncload('~/srv_ccrc/data01/z3263455/escm27/pdetop5/100free/tavg.03001.01.01.nc','xu');
free100_yu=yu;
free100_xu = xu;

%calc dx = acos(lat)dlon = 6400000m *1.8deg/180*pi
for i=1:100
xdiff(i) = 6400000*cos(yu(i)/180*pi)*(xu(3)-xu(2))/180*pi;
end

tmp = find(v>9.96e36);
v(tmp)=NaN;

My = zeros(19,100);
%calculate Mx
for z=1:19
    for y = 1:100
    	for x = 1:100
            if(mskhr(y,x) == 1 && (isnan(v(z,y,x)) ~= 1))
                My(z,y) = My(z,y) + v(z,y,x)*xdiff(y);
            end
        end
    end
end

amoc = zeros(19,100);
for y = 100:-1:1
    sumpsi = 0;
    for z = 1:19
        %if(mskhr(z,y) == 1 && (isnan(v(1,y,x)) ~= 1))
	    sumpsi = sumpsi+My(z,y)*zdiff(z);
       	    amoc(z,y) = sumpsi;
        %else
	 %   amoc(z,y) = 0;
	 %   sumpsi = 0;		
	%end
    end
end

save 100free.mat amoc -APPEND;
save 100free.mat xu -APPEND;
save 100free.mat yu -APPEND;
save 100free.mat time -APPEND;
