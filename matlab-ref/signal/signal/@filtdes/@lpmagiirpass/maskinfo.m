function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2002/06/25 13:18:02 $

if isdb(d), apass = get(d, 'Apass');
else,       apass = get(d, 'Epass'); end

cmd{1}.magfcn     = 'cpass';
cmd{1}.amplitude  = apass;
cmd{1}.filtertype = 'lowpass';
cmd{1}.magunits   = get(d, 'magUnits');

cmd{2}            = [];

% [EOF]
