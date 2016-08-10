function b = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2004/12/26 22:12:36 $

[params, values, descs] = abstract_genmcode(h, d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({params{:}, 'dens'}, ...
    {values{:}, getmcode(d, 'DensityFactor')}, {descs{:}, ''}));
b.cr;
b.addcr(minorddesc(h, 'firpmord'));
b.addcr('[N, Fo, Ao, W] = firpmord([Fstop, Fpass]%s, [0 1], [Dstop, Dpass]);', getfsstr(d));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = firpm(N, Fo, Ao, W, {dens});');
b.add('Hd = dfilt.dffir(b);');

% [EOF]