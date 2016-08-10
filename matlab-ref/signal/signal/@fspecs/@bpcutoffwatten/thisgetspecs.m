function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs. This is used by FVTOOL for drawing the mask. 

%   Copyright 2008 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2008/05/12 21:36:30 $

specs.Fstop1 = this.Fcutoff1;
specs.Fpass1 = this.Fcutoff1;
specs.Fpass2 = this.Fcutoff2;
specs.Fstop2 = this.Fcutoff2;
specs.Astop1 = this.Astop1;
specs.Apass = this.Apass;
specs.Astop2 = this.Astop2;

% [EOF]
