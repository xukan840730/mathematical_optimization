function b = genmcode(h, d)
%GENMCODE Generate M-code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:10:47 $

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', 'Fpass'}, {getmcode(d, 'Order'), getmcode(d, 'Fpass')}))
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = firhalfband(N, Fpass%s);', getfsstr(d));
b.add('Hd = dfilt.dffir(b);');

% [EOF]