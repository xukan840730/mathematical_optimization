function s = getdesignpanelstate(this)
%GETDESIGNPANELSTATE   Get the designpanelstate.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2005/06/16 08:32:48 $

s = aswofs_getdesignpanelstate(this);

s.Components{1}.Tag   = 'fdadesignpanel.hpfreqpass';
s.Components{1}.Fpass = sprintf('%g', this.Fpass);

s.Components{3}.Tag      = 'fdadesignpanel.hpmag'
s.Components{3}.magUnits = 'dB';
s.Components{3}.Apass    = sprintf('%g', this.Apass);
s.Components{3}.Astop    = sprintf('%g', this.Astop);

% [EOF]
