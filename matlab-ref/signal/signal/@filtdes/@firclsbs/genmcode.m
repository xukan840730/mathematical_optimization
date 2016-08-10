function b = genmcode(h, d)
%GENMCODE Generate M-code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2003/03/02 10:17:54 $

[fs, fsstr] = getfsstr(d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', 'Fc1', 'Fc2', 'Dpass1U', 'Dpass1L', 'DstopU', ...
        'DstopL', 'Dpass2U', 'Dpass2L'}, ...
    {getmcode(d, 'Order'), getmcode(d, 'Fc1'), getmcode(d, 'Fc2'), ...
        getmcode(d, 'Dpass1Upper'), getmcode(d, 'Dpass1Lower'), getmcode(d, 'DstopUpper'), ...
        getmcode(d, 'DstopLower'), getmcode(d, 'Dpass2Upper'), getmcode(d, 'Dpass2Lower')}, ...
    {'', '', '', '', '', '', '', '', ''}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = fircls(N, [0 Fc1 Fc2 %s]%s, %s, %s, %s);', fsstr, fs, ...
    '[1 0 1]', '[1+Dpass1U DstopU 1+Dpass2U]', '[1-Dpass1L -DstopL 1-Dpass2L]');
b.add('Hd = dfilt.dffir(b);');

% [EOF]