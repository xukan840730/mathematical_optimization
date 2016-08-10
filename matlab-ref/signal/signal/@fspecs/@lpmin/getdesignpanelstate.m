function s = getdesignpanelstate(this)
%GETDESIGNPANELSTATE   Get the designpanelstate.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2005/06/16 08:34:19 $

c = aswfs_getdesignpanelstate(this);
c = c.Components;

c{1}.Tag   = 'fdadesignpanel.lpfreqpassstop';
c{1}.Fpass = sprintf('%g', this.Fpass);
c{1}.Fstop = sprintf('%g', this.Fstop);

c{2}.Tag   = 'siggui.filterorder';
c{2}.order = '10';
c{2}.mode  = 'minimum';

c{3}.Tag      = 'fdadesignpanel.lpmag';
c{3}.magUnits = 'dB';
c{3}.Apass    = sprintf('%d', this.Apass);
c{3}.Astop    = sprintf('%d', this.Astop);

s.Components = c;

% [EOF]
