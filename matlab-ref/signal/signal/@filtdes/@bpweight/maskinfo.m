function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2002/06/25 13:17:51 $

cmd{1}.magfcn     = 'wstop';
cmd{1}.filtertype = 'highpass';
cmd{1}.tag        = 'stop1';
cmd{1}.weight     = get(d, 'Wstop1');

cmd{2}.magfcn     = 'wpass';
cmd{2}.filtertype = 'bandpass';
cmd{2}.weight     = get(d, 'Wpass');

cmd{3}.magfcn     = 'wstop';
cmd{3}.filtertype = 'lowpass';
cmd{3}.tag        = 'stop2';
cmd{3}.weight     = get(d, 'Wstop2');

% [EOF]
