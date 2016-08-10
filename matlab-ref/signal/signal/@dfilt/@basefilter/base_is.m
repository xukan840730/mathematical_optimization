function f = base_is(Hb, fcn, varargin)
%BASE_IS Returns an array for the requested fcn.

%   Author: J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.2.4.2 $  $Date: 2004/04/12 23:52:46 $

% This should be private

f = base_num(Hb, fcn, varargin{:});

% Make sure that we return logicals.
f = logical(f);

% [EOF]
