function b = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:08:23 $

[params, values, descs] = abstract_genmcode(h, d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', params{:}}, {getmcode(d, 'Order'), values{:}}, ...
    {'', descs{:}}));
b.cr;
b.addcr('b  = %s(N, F%s, A, W, ''hilbert'');', designfunction(h, d), getfsstr(d));
b.add('Hd = dfilt.dffir(b);');

% [EOF]