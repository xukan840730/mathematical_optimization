function n = thisnstates(Hd)
%NSTATES  Number of states in discrete-time filter.
%   NSTATES(Hd) returns the number of states in the
%   discrete-time filter Hd.  The number of states depends on the filter
%   structure and the coefficients.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/12 23:55:11 $

ncoeffs = Hd.ncoeffs;
n = max(ncoeffs)-1;

% If ncoeffs defauls to [], don't have nstates return a negative value.
if n < 0, n = 0; end
