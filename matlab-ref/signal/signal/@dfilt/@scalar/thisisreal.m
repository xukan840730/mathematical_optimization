function f = thisisreal(Hd)
%THISISREAL  True for filter with real coefficients.
%   THISISREAL(Hd) returns 1 if filter Hd has real coefficients, and 0
%   otherwise. 
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 2002/07/29 21:42:02 $

% This should be private

f = isreal(Hd.Gain);