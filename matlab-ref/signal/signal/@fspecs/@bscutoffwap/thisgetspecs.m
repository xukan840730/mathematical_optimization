function specs = thisgetspecs(this)
%THISGETSPECS   

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2005/06/16 08:29:36 $

specs.Fpass1 = this.F3dB1;
specs.Fstop1 = this.F3dB1;
specs.Fstop2 = this.F3dB2;
specs.Fpass2 = this.F3dB2;
specs.Apass1 = this.Apass;
specs.Astop  = NaN;
specs.Apass2 = this.Apass;

% [EOF]
