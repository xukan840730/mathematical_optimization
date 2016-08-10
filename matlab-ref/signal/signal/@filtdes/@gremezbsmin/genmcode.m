function [params, values, descs, iargs] = genmcode(h, d)
%GENMCODE Generate M-code that designs a lowpass filter using GREMEZ

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:07:44 $

[params, values, descs]  = abstract_genmcode(h, d);

[fs, fsstr] = getfsstr(d);

iargs = sprintf('[0 Fpass1 Fstop1 Fstop2 Fpass2 %s]%s, %s, %s', fsstr, fs, ...
    '[1 1 0 0 1 1]', '[Dpass1 Dstop Dpass2]');

% [EOF]
