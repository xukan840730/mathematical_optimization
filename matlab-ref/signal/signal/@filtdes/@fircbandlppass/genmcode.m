function [params, values, descs, iargs] = genmcode(h, d)
%GENMCODE Generate M-code that designs a lowpass filter using GREMEZ

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:06:07 $

mu = get(d, 'magUnits'); set(d, 'magUnits', 'linear');

params = {'N', 'Fpass', 'Fstop', 'Dpass'};
values = {getmcode(d, 'Order'), getmcode(d, 'Fpass'), getmcode(d, 'Fstop'), ...
        getmcode(d, 'Dpass')};
descs  = {'', '', '', ''};

set(d, 'magUnits', mu);

[fsstr, fs] = getfsstr(d);

iargs = sprintf('[0 Fpass Fstop %s]%s, [1 1 0 0], [Dpass 1]', fs, fsstr);

% [EOF]
