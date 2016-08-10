function b = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2008/05/31 23:26:52 $

[params, values, descs, str] = fir1_genmcode(d);

% Throw out the flag, it isn't used.
b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', params{2:end}, 'L'}, ...
    {getmcode(d, 'order'), values{2:end}, getmcode(d, 'band')}, ...
    {'', descs{2:end}, ''}));
b.cr;
b.addcr(str);
b.cr;
b.addcr('% Calculate the coefficients with the FIRNYQUIST function.');
b.addcr('b  = firnyquist(N, L, win);');
b.add('Hd = dfilt.dffir(b);');

% [EOF]