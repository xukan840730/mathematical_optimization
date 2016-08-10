function f = isfir(Hb)
%ISFIR  True for FIR filter.
%   ISFIR(Hb) returns 1 if filter Hd is FIR, and 0 otherwise.
%
%   See also DFILT.   
  
%   Author: J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2002/07/29 21:41:34 $

f = base_is(Hb, 'thisisfir');

% [EOF]
