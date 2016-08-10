function cmd = fs2_maskinfo(hObj, d)
%FS2_MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2002/06/25 13:17:26 $

cmd{1}.frequency = [0 get(d, 'Fstop1')];
cmd{2}.frequency = [get(d, 'Fstop1') get(d, 'Fstop2')];
cmd{3}.frequency = [get(d, 'Fstop2') getnyquist(d)];

% [EOF]
