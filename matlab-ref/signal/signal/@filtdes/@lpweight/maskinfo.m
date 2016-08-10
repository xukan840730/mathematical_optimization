function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2002/06/25 13:18:04 $

cmd{1}.magfcn     = 'wpass';
cmd{1}.filtertype = 'lowpass';
cmd{1}.weight     = get(d, 'Wpass');

cmd{2}.magfcn     = 'wstop';
cmd{2}.filtertype = 'lowpass';
cmd{2}.weight     = get(d, 'Wstop');

% [EOF]
