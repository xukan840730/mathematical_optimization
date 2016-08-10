function [y,zf] = secfilter(Hd,x,zi)
%SECFILTER Filter this section.
%   [Y,Zf] = SECFILTER(Hd,X,ZI) filters this section.  This function is only
%   intended to be called from DFILT/FILTER.  The initial conditions have
%   already been padded for the C++ implementation.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.2.4.3 $  $Date: 2004/04/12 23:58:20 $
  
q = Hd.filterquantizer;
k = Hd.privlattice;
kconj = Hd.privconjlattice;
ladder = Hd.privladder;
[y,zf] = latticearmafilter(q,k,kconj,ladder,x,zi);
