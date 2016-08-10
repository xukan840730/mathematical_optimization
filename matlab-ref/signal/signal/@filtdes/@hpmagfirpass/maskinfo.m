function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:08:37 $

if isdb(d),
    apass = get(d, 'Apass');
else
    apass = get(d, 'Dpass');
end

cmd{1}.magfcn     = 'pass';
cmd{1}.amplitude  = apass;
cmd{1}.filtertype = 'lowpass';
cmd{1}.magunits   = get(d, 'magUnits');

cmd{2} = [];

% [EOF]
