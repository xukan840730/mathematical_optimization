function [y,zf] = latticecapcfilter(q,Hd,x,zi)
% LATTICECAPCFILTER Filter for DFILT.CALATTICEPC class in double precision mode

%   Author(s): V.Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2004/04/12 23:57:43 $

x = quantizeinput(q,x);

k1 = Hd.Allpass1;
k2 = Hd.Allpass2;
beta = Hd.Beta;

[y,zf] = latticecapcfilter(k1,k2,beta,x,zi);

