function [y,zf] = secfilter(Hd,x,zi)
%SECFILTER Filter this section.
%   [Y,Zf] = SECFILTER(Hd,X,ZI) filters this section.  This function is only
%   intended to be called from DFILT/FILTER.  The initial conditions have
%   already been padded for the C++ implementation.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.2.4.3 $  $Date: 2004/04/12 23:55:45 $
  
q = Hd.filterquantizer;
b = Hd.privnum;
tapIndexi = Hd.TapIndex;
[y,zf,tapIndexf] = dffirfilter(q,b,x,zi,tapIndexi);
Hd.TapIndex = tapIndexf;
