function [Dstop, Dpass] = getdesignspecs(h, d)
%GETDESIGNSPECS Returns the specs for the design

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2002/10/04 18:12:36 $

magUnits = get(d,'magUnits');
set(d,'magUnits','linear'); % Set temporarily to linear
Dstop = get(d,'Dstop');
Dpass = get(d,'Dpass');
set(d,'magUnits',magUnits); % Reset units

% [EOF]
