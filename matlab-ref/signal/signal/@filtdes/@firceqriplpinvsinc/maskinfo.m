function cmd = maskinfo(h,d)
%MASKINFO

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2002/07/11 14:56:08 $

cmd = firceqriplp_maskinfo(h,d);

cmd.bands{1}.magfcn     = 'invsinc';
cmd.bands{1}.frequency  = cmd.bands{1}.frequency(2);
cmd.bands{1}.freqfactor = get(d, 'invSincFreqFactor');
cmd.bands{1}.sincpower  = get(d, 'InvSincPower');

% [EOF]
