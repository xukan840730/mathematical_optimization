function b = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.7.32.1 $  $Date: 2009/12/17 13:59:13 $

[Fstop1, Fstop2, Astop] = getdesignspecs(h, d);

b = sigcodegen.mcodebuffer;

p = {'N', 'Fstop1', 'Fstop2', 'Astop'};
v = {getmcode(d, 'Order'), getmcode(d, Fstop1), ...
    getmcode(d, Fstop2), getmcode(d, Astop)};

b.addcr(b.formatparams(p, v));
b.cr;
b.addcr(designdesc(d));

b.addcr('h  = fdesign.bandstop(''N,Fst1,Fst2,Ast'', N, Fstop1, Fstop2, Astop%s);', getfsinput(d));
b.add('Hd = design(h, ''cheby2'');');

% [EOF]
