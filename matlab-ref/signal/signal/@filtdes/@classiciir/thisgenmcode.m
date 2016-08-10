function b = thisgenmcode(d)
%THISGENMCODE Perform the IIR genmcode.

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.2.4.3.54.1 $  $Date: 2009/12/17 13:59:23 $

% Frequencies Have been prenormalized (0 to 1)

% Call type specific design
b = genmcode(d.responseTypeSpecs, d);

% [EOF]