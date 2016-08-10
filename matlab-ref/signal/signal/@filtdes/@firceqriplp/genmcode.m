function b = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:06:22 $

[params, values, descs, args] = abstract_genmcode(h, d);
[Dpass, Dstop]                = getdesignspecs(h, d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({params{:}, 'Dstop', 'Dpass'}, ...
    {values{:}, getmcode(d, Dstop), getmcode(d, Dpass)}, ...
    {descs{:}, '', ''}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = firceqrip(%s);', sprintf(args, '[Dpass, Dstop]'));
b.add('Hd = dfilt.dffir(b);');

% [EOF]