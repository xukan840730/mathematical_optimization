function b = genmcode(h, d)
%GENMCODE Returns the MCode necessary to generate the filter.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2004/12/26 22:12:32 $

[params, values, descs] = abstract_genmcode(h, d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({params{:}, 'dens'}, {values{:}, getmcode(d, 'DensityFactor')}, ...
    {descs{:}, ''}));
b.cr;
b.addcr(minorddesc(h, 'firpmord'));
b.addcr('[N, Fo, Ao, W] = firpmord(%s%s, %s);', '[Fstop1 Fpass1 Fpass2 Fstop2]', ...
    getfsstr(d), '[0 1 0], [Dstop1 Dpass Dstop2]');
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = firpm(N, Fo, Ao, W, {dens});');
b.add('Hd = dfilt.dffir(b);');

% [EOF]