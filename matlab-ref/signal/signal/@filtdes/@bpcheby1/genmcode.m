function b = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.7.32.1 $  $Date: 2009/12/17 13:59:04 $

[~, ~, Apass] = getdesignspecs(h, d);

b = sigcodegen.mcodebuffer;

p = {'N', 'Fpass1', 'Fpass2', 'Apass'};
v = {getmcode(d, 'Order'), getmcode(d, 'Fpass1'), getmcode(d, 'Fpass2'), ...
        getmcode(d, Apass)};

b.addcr(b.formatparams(p, v));
b.cr;
b.addcr(designdesc(d));
b.addcr('h  = fdesign.bandpass(''N,Fp1,Fp2,Ap'', N, Fpass1, Fpass2, Apass%s);', getfsinput(d));
b.add('Hd = design(h, ''cheby1'');');

% [EOF]
