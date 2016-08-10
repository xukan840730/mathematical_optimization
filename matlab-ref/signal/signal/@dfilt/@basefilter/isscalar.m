function f = isscalar(Hb,varargin)
%ISSCALAR  True if scalar filter.
%   ISSCALAR(Hb) returns 1 if Hb is a scalar filter, and 0 otherwise.
%
%   See also DFILT.   
  
%   Author: J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 2002/08/01 20:17:02 $

f = base_is(Hb, 'thisisscalar', varargin{:});

% [EOF]
