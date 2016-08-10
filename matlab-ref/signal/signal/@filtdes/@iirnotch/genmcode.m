function b = genmcode(h, d)
%GENMCODE Generate M-code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:09:17 $

% Make sure we have the magnitude units in dB.
mu = get(d, 'MagUnits');
set(d, 'MagUnits', 'dB');

[params, values, descs, str] = getbandwidth(d);

params = {'N', params, 'Apass'};
values = {getmcode(d, 'Order'), values, getmcode(d, 'Apass')};
descs  = {'', descs, 'Bandwidth Attenuation'};

set(d, 'MagUnits', mu);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams(params, values, descs));
b.addcr(str, designdesc(d));
b.addcr('[b, a] = iircomb(N, BW%s, Apass);', getfsstr(d));
b.add('Hd     = dfilt.df2(b, a);');

% [EOF]