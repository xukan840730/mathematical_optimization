function b = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.7.32.1 $  $Date: 2009/12/17 13:59:37 $

[Fpass, Apass, Astop] = getdesignspecs(h, d);

b = sigcodegen.mcodebuffer;

p = {'N', 'Fpass', 'Apass', 'Astop'};
v = {getmcode(d, 'Order'), getmcode(d, Fpass), getmcode(d, Apass), getmcode(d, Astop)};

b.addcr(b.formatparams(p,v));
b.cr;
b.addcr(designdesc(d));

b.addcr('h  = fdesign.lowpass(''N,Fp,Ap,Ast'', N, Fpass, Apass, Astop%s);', getfsinput(d));
b.add('Hd = design(h, ''ellip'');');

% [EOF]
