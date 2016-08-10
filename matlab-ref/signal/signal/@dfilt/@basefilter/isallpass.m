function f = isallpass(Hb, varargin)
%ISALLPASS  True for allpass filter.
%   ISALLPASS(Hb) returns 1 if filter Hb is all-pass, and 0 otherwise.
%
%   ISALLPASS(Hb,TOL) uses tolerance TOL to determine when two numbers are
%   close enough to be considered equal.
%
%   See also DFILT.   
  
%   Author: J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.2.4.1 $  $Date: 2007/12/14 15:07:30 $

error(nargchk(1,2,nargin,'struct'));

f = base_is(Hb, 'thisisallpass', varargin{:});

% [EOF]
