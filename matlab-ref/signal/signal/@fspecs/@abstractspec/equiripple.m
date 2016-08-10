function varargout = equiripple(this, varargin)
%EQUIRIPPLE   Design an equiripple filter.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2005/06/16 08:27:13 $

[varargout{1:nargout}] = design(this, 'equiripple', varargin{:});

% [EOF]
