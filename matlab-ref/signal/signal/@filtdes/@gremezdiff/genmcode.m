function [params, values, descs, iargs] = genmcode(h, d)
%GENMCODE Generate M-code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:07:50 $

[params, values, descs] = abstract_genmcode(h,d);

iargs = sprintf('F%s, A, W, ''differentiator''', getfsstr(d));

% [EOF]
