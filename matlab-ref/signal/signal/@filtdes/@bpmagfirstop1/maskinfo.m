function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:02:09 $

if isdb(d),
    astop = get(d, 'Astop1');
else
    astop = get(d, 'Dstop1');
end

cmd{1}.magfcn     = 'stop';
cmd{1}.amplitude  = astop;
cmd{1}.filtertype = 'highpass';

cmd{2} = [];
cmd{3} = [];

% [EOF]
