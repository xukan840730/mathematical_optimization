function b = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.2.6.2 $  $Date: 2004/04/13 00:09:43 $

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', 'Fc'}, {getmcode(d, 'Order'), getmcode(d, 'Fc')}));
b.cr;
b.addcr(designdesc(d));
b.addcr('[b,a,b1,b2,sos_var,g] = maxflat(N, ''sym'', Fc%s);', getfsstr(d));
b.add('Hd                    = dfilt.df2sos(sos_var, g);');

% [EOF]