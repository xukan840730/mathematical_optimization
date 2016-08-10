function b = genmcode(h,d)
%GENMCODE Generate M-code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:06:49 $

[fs, fsstr] = getfsstr(d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', 'F', 'A', 'U', 'L'}, ...
    {getmcode(d, 'Order'), getmcode(d, 'FrequencyVector'), getmcode(d, 'MagnitudeVector'), ...
        getmcode(d, 'UpperVector'), getmcode(d, 'LowerVector')}, ...
    {'', '', '', 'Upper Vector', 'Lower Vector'}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = fircls(N, F%s, A, U, L);', fs);
b.add('Hd = dfilt.dffir(b);');

% [EOF]
