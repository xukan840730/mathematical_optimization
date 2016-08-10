function [I, T] = computeimpz(this, varargin)
%COMPUTEIMPZ   Calculate the impulse response.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2004/07/14 06:46:25 $

[I, T] = impz(this.Numerator, this.Denominator, varargin{:});

% [EOF]
