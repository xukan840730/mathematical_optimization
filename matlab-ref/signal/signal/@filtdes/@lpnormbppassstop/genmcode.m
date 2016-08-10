function [params, values, descs, str, args] = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2003/03/02 10:22:36 $

[fsstr, fs] = getfsstr(d);

[params, values, descs] = abstract_genmcode(h, d);

str = sprintf('\nF = [0 Fstop1 Fpass1 Fpass2 Fstop2 %s]%s;\n', fs, fsstr);

args = 'F, F, [0 0 1 1 0 0], [Wstop1 Wstop1 Wpass Wpass Wstop2 Wstop2]';

% [EOF]
