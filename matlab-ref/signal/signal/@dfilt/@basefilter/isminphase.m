function f = isminphase(Hb,varargin)
%ISMINPHASE True if minimum phase.
%   ISMINPHASE(Hb) returns 1 if filter Hb is minimum phase, and 0 otherwise.
%
%   ISMINPHASE(Hb,TOL) uses tolerance TOL to determine when two numbers are
%   close enough to be considered equal.
%
%   See also DFILT.   
  
%   Author: J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.2.4.1 $  $Date: 2007/12/14 15:07:34 $
  
error(nargchk(1,2,nargin,'struct'));

f = base_is(Hb, 'thisisminphase', varargin{:});

% [EOF]
