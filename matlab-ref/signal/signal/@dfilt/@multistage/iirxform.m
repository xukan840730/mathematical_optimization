function [Ht,anum,aden] = iirxform(Ho,fun,varargin)
%XFORM Frequency Transformations.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2005/11/18 14:24:32 $

for n = 1:length(Ho.Stage),   
   Ht(n) = feval(fun, Ho.Stage(n), varargin{:});
end

% Call the transformation on a dummy filter to get the allpass coeffs
[num,den,anum,aden] = feval(fun,1,1,varargin{:});

Ht = feval(str2func(class(Ho)),Ht(:));

% [EOF]
