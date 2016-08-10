function b = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.4.2 $  $Date: 2004/04/13 00:10:48 $

b = abstract_genmcode(h, d);
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = firhalfband(''minorder'', Fpass%s, Dpass);', getfsstr(d));
b.add('Hd = dfilt.dffir(b);');

% [EOF]