function [params, values, descs, iargs] = genmcode(h, d)
%GENMCODE Generate M-code that designs a lowpass filter using GREMEZ

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:08:14 $

[params, values, descs] = abstract_genmcode(h, d);

[fs, fsstr] = getfsstr(d);

iargs = sprintf('[0 Fpass Fstop %s]%s, [1 1 0 0], [Dpass Dstop]', fsstr, fs);

% [EOF]
