function b = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:08:35 $

[aparams, avalues, adescs]   = abstract_genmcode(h, d);
[params, values, descs, str] = fir1_genmcode(d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({aparams{:}, params{2:end}}, ...
    {avalues{:}, values{2:end}}, ...
    {adescs{:},  descs{2:end}}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = firhalfband(''minorder'', 1-Fpass%s, Dpass, ''kaiser'', ''high'');', getfsstr(d));
b.add('Hd = dfilt.dffir(b);');

% [EOF]
