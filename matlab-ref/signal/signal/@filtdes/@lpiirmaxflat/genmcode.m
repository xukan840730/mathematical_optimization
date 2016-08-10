function b = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.2.6.2 $  $Date: 2004/04/13 00:09:46 $

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'Nb', 'Na', 'Fc'}, ...
    {getmcode(d, 'numOrder'), getmcode(d, 'denOrder'), getmcode(d, 'Fc')}));
b.cr;
b.addcr(designdesc(d));
b.addcr('[b,a,b1,b2,sos_var,g] = maxflat(Nb, Na, Fc%s);', getfsstr(d));
b.add('Hd                    = dfilt.df2sos(sos_var, g);');

% [EOF]