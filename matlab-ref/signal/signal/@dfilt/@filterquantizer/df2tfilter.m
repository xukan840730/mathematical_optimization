function [y,zf] = df2tfilter(q,b,a,x,zi)
% DF2TFILTER Filter for DFILT.DF2T class in double precision mode

%   Author(s): V.Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2004/04/12 23:57:27 $

x = quantizeinput(q,x);

[y,zf] = df2tfilter(b,a,x,zi);
