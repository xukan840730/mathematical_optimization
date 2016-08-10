function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/10/18 03:25:55 $

minfo.F0 = this.F0;
minfo.Q = this.Q;
minfo.BW = minfo.F0/minfo.Q;

minfo.BWpass = [];
minfo.BWstop = [];
if this.NormalizedFrequency,
    [Flow,Fhigh] = parameqbandedge(minfo.F0*pi,minfo.BW*pi,0);
    minfo.Flow   = Flow/pi;
    minfo.Fhigh  = Fhigh/pi;
else
    Fs = this.Fs;
    [Flow,Fhigh] = parameqbandedge(pi*minfo.F0/(Fs/2),pi*minfo.BW/(Fs/2),0);
    minfo.Flow   = Flow*Fs/(2*pi);
    minfo.Fhigh  = Fhigh*Fs/(2*pi);
end

minfo.GBW    = 10*log10(.5);

if isprop(this,'Apass')
    minfo.Apass = this.Apass;
else
    minfo.Apass = [];
end

if isprop(this,'Astop')
    minfo.Astop = this.Astop;
else
    minfo.Astop = [];
end

% [EOF]
