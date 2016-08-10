function b = genmcode(~, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.7.10.1 $  $Date: 2009/12/17 13:59:26 $

b = sigcodegen.mcodebuffer;

p = {'N', 'Fc'};
v = {getmcode(d, 'Order'), getmcode(d, 'Fc')};

b.addcr(b.formatparams(p, v));
b.cr;
b.addcr(designdesc(d));

b.addcr('h  = fdesign.highpass(''N,F3dB'', N, Fc%s);', getfsinput(d));
b.add('Hd = design(h, ''butter'');');

% [EOF]