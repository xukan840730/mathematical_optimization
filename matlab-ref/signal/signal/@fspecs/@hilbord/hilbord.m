function this = hilbord(varargin)
%HILBORD   Construct a HILBORD object.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/06/16 08:31:33 $

this = fspecs.hilbord;

this.ResponseType = 'Hilbert Transformer with filter order';

this.FilterOrder = 30;

this.setspecs(varargin{:});

% [EOF]
