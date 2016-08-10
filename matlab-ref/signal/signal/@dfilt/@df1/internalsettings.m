function s = internalsettings(h)
%INTERNALSETTINGS Returns the fixed-point settings viewed by the algorithm.  

%   Author(s): P. Costa
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2004/01/25 23:03:01 $

s = internalsettings(h.filterquantizer);

% [EOF]
