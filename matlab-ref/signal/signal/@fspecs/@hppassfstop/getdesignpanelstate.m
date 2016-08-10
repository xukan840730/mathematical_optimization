function s = getdesignpanelstate(this)
%GETDESIGNPANELSTATE   Get the designpanelstate.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2005/06/16 08:32:51 $

s = aswofs_getdesignpanelstate(this);

s.Components{1}.Tag = 'fdadesignpanel.hpfreqpassstop';
s.Components{1}.Fstop = sprintf('%g', this.Fstop);
s.Components{1}.Fpass = sprintf('%g', this.Fpass);

s.Components{3}.Tag = 'fdadesignpanel.hpmagpass';
s.Components{3}.magUnits = 'dB';
s.Components{3}.Apass = sprintf('%g', this.Apass);

% [EOF]
