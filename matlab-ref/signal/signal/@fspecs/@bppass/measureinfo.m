function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/06/16 08:29:06 $

minfo.Fstop1 = [];
minfo.Fcutoff1 = [];
minfo.Fpass1 = this.Fpass1;
minfo.Fpass2 = this.Fpass2;
minfo.Fcutoff2 = [];
minfo.Fstop2 = [];
minfo.Astop1 = [];
minfo.Apass  = this.Apass;
minfo.Astop2 = [];

% [EOF]
