function b = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2004/12/26 22:12:39 $

[params, values, descs, iargs] = abstract_genmcode(h, d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', params{:}, 'dens'}, ...
    {getmcode(d, 'Order'), values{:}, getmcode(d, 'DensityFactor')}, {'', descs{:}, ''}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = firpm(N, %s, {dens});', iargs);
b.add('Hd = dfilt.dffir(b);');

% [EOF]