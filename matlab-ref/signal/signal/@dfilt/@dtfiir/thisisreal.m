function f = thisisreal(Hd)
%THISISREAL  True for filter with real coefficients.
%   THISISREAL(Hd) returns 1 if filter Hd has real coefficients, and 0
%   otherwise. 
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 2002/07/29 21:42:35 $

% This should be private

f = isreal(Hd.Numerator) & isreal(Hd.Denominator);