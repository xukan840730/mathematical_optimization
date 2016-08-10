function [y,z] = secfilter(Hd,x,z)
%SECFILTER Filter this section.
%   [Y,Zf] = SECFILTER(Hd,X,ZI) filters this section.  This function is only
%   intended to be called from DFILT/FILTER.  
%
%   See also DFILT.   
  
%   Author: R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.2.4.2 $  $Date: 2004/04/12 23:57:15 $

q = Hd.filterquantizer;
bfft = Hd.fftcoeffs; % Have been expanded to multichannel by ziexpand. Do not use fftcoeffs(Hd).
[y,z] = fftfirfilter(q,Hd,bfft,x,z);
