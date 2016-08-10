function b = genmcode(h, d)
%GENMCODE Generate M-code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:09:26 $

% Make sure we have the magnitude units in dB.
mu = get(d, 'MagUnits');
set(d, 'MagUnits', 'dB');

[params, values, descs, str] = getbandwidth(d);

apass = getmcode(d, 'Apass');

set(d, 'MagUnits', mu);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', params, 'Apass'}, ...
    {getmcode(d, 'Order'), values, apass}, ...
    {'', descs, 'Bandwidth Attenuation'}));
b.addcr(str, designdesc(d));
b.addcr('[b, a] = iircomb(N, BW%s, Apass, ''peak'');', getfsstr(d));
b.add('Hd     = dfilt.df2(b, a);');

% [EOF]
