%%%%%%%%%%%%%%%
function o2sat=o2sat(sst,sss)

%sst in degC
%sss in psu
%from ./common/gasbc.F of UVic_ESCM
f1 = log((298.15 - sst)./(273.15 + sst));
f2 = f1.*f1;
f3 = f2.*f1;
f4 = f3.*f1;
f5 = f4.*f1;
o2sat = exp (2.00907 + 3.22014.*f1 + 4.05010.*f2...
     + 4.94457.*f3 - 2.56847E-1.*f4 + 3.88767.*f5...
     + sss.*(-6.24523e-3 - 7.37614e-3.*f1 - 1.03410e-2.*f2...
     - 8.17083E-3.*f3) - 4.88682E-7.*sss.*sss);
%           Convert from ml/l to mol/m^3
o2sat = o2sat/22391.6.*1000.0;
%%%%%%%%%%%%%%%%%%%%%
